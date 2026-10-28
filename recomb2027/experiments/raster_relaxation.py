"""Whole supplied crop count bounds under maximum label-transport displacement.

This is a different, explicit uncertainty model from summed corridor area. It
admits any source raster topology; common geometry is relaxed, so no attaining
segmentation is claimed for a marginal endpoint. Outside-crop candidate windows
include every observed label plus an unknown exterior label conservatively.
"""
from pathlib import Path
import sys,json,csv,math,time
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.relaxation import candidate_count_intervals

def main():
    t0=time.perf_counter();obs=json.loads((ROOT/'data/observed_crop_counts.json').read_text());summary=[];outputs=[]
    for roi,rr in obs['rois'].items():
        meta=json.loads((ROOT/'data/source_crops'/roi/'metadata.json').read_text());mask=np.load(ROOT/'data/source_crops'/roi/'masks.npz')['cell_mask']
        transform=F(str(meta['transform'][0][0]));x0,y0=meta['bbox_pixels'][:2]
        coords=[]
        for row in rr['records']:
            x=F(row['x_um'])*transform-x0;y=F(row['y_um'])*transform-y0
            assert math.floor(x)==row['pixel_col'] and math.floor(y)==row['pixel_row']
            coords.append((x,y))
        for radius in [F(0),F(1,4),F(1,2),F(1)]:
            candidates=[];edge_cases=0;ambiguous=0;all_labels=set(map(int,np.unique(mask)))|{-1}
            for row,(x,y) in zip(rr['records'],coords):
                rad=radius*transform;l,h=math.floor(x-rad),math.floor(x+rad);a,b=math.floor(y-rad),math.floor(y+rad)
                if min(l,a)<0 or h>=mask.shape[1] or b>=mask.shape[0]:allowed=all_labels;edge_cases+=1
                else:allowed=set(map(int,np.unique(mask[a:b+1,l:h+1])))
                assert row['cell'] in allowed
                ambiguous+=len(allowed)>1
                candidates.append({'source':row['cell'],'gene':row['gene'],'candidates':sorted(allowed)})
            bounds=candidate_count_intervals(candidates)
            for (cell,gene),r in sorted(bounds.items()):outputs.append(dict(roi=roi,radius_um=str(radius),cell_id=cell,gene=gene,**r))
            summary.append({'roi':roi,'radius_um':str(radius),'RNA':len(candidates),'ambiguous_RNA':ambiguous,
                            'outside_context_windows':edge_cases,'cell_gene_entries':len(bounds),
                            'nontrivial_intervals':sum(r['lower']<r['upper'] for r in bounds.values()),
                            'candidate_edges':sum(len(r['candidates']) for r in candidates)})
    with (ROOT/'results/raster_expression_outer.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(outputs[0]));w.writeheader();w.writerows(outputs)
    out={'model':'supremum displacement of a common source-label transport, with independent candidate-owner relaxation; L-infinity physical norm',
         'count_scope':'Only measured retained RNA in the two supplied crops; outside RNA not modeled.',
         'not_the_area_budget_model':True,'table':summary,'seconds':time.perf_counter()-t0}
    (ROOT/'results/raster_outer_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
