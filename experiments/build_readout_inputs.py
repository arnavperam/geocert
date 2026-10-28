"""Reconstruct crop-level gene counts with PyArrow and verify source-mask labels.

Usage: python experiments/build_readout_inputs.py --source-root SOURCE_BASE
SOURCE_BASE contains data/native/roi_01 and roi_02. NumPy and PyArrow are needed
only for this raw-crop check; the downstream readout audit is standard-library.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse, hashlib, json
import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=ROOT/'data/readout_inputs.json');args=ap.parse_args()
    charts=json.loads((ROOT/'data/charts.json').read_text());rois={};hashes={};checks=0
    for roi in sorted({c['roi'] for c in charts}):
        path=args.source_root/'data/native'/roi
        meta=json.loads((path/'metadata.json').read_text());mask=np.load(path/'masks.npz')['cell_mask']
        rows=pq.read_table(path/'transcripts.parquet').to_pylist()
        scale=F(float(meta['transform'][0][0]));bx,by,*_=meta['bbox_pixels']
        byid={};counts=Counter()
        for r in rows:
            x=F(float(r['x_location']))*scale-bx;y=F(float(r['y_location']))*scale-by
            col=x.numerator//x.denominator;row=y.numerator//y.denominator
            assert int(mask[row,col])==r['mask_label']
            assert r['is_gene'] and r['qv']>=20
            tid=str(r['transcript_id']);assert tid not in byid;byid[tid]=r
            counts[(int(r['mask_label']),str(r['feature_name']))]+=1
        labels=sorted({lab for c in charts if c['roi']==roi for lab in c['cell_labels']})
        edges=set(map(int,np.concatenate([mask[0],mask[-1],mask[:,0],mask[:,-1]])))
        cells={str(lab):{'counts':dict(sorted((gene,n) for (cell,gene),n in counts.items() if cell==lab)),
                        'touches_crop_edge':lab in edges,'mask_pixels':int((mask==lab).sum()),
                        'mask_area_um2':str(F(int((mask==lab).sum()))/(scale*scale)),
                        'cell_id':meta['label_to_cell_id'].get(str(lab))} for lab in labels}
        for c in charts:
            if c['roi']!=roi:continue
            for pt in c['points']:
                r=byid[pt['transcript_id']]
                assert r['feature_name']==pt['gene'] and r['mask_label']==pt['mask_label']==c['cell_labels'][pt['source']]
                gx,gt=(F(float(r['x_location'])),F(float(r['y_location']))) if c['orientation']==0 else (F(float(r['y_location'])),F(float(r['x_location'])))
                ox=(bx if c['orientation']==0 else by)+c['left'];ot=F(by if c['orientation']==0 else bx)+c['row']+F(1,2)
                assert gx-F(ox)/scale==F(pt['x']) and gt-ot/scale==F(pt['t'])
                checks+=1
        rois[roi]={'rna_records':len(rows),'cells':cells}
        for name in ['transcripts.parquet','masks.npz','metadata.json']:
            hashes[f'{roi}/{name}']=hashlib.sha256((path/name).read_bytes()).hexdigest()
    result={'scope':'Q20 gene counts in the supplied crops, using mask_label; outside-corridor assignments remain fixed in each separate query.',
            'parser':'PyArrow '+__import__('pyarrow').__version__,'source_hashes':hashes,
            'source_rna_records':sum(r['rna_records'] for r in rois.values()),'chart_records_checked':checks,'rois':rois}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['parser','source_rna_records','chart_records_checked']}))

if __name__=='__main__':main()
