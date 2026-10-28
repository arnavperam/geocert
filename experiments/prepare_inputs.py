"""Build the observation table with PyArrow and verify every raster source owner."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import json,hashlib,math,sys
import numpy as np
import pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parents[1]
def prepare(source=None):
    source=Path(source) if source else ROOT/'data/source_crops'
    result={'scope':'Observed Q20 measured-gene counts in supplied crops. Complete-cell flags are geometric containment checks.',
            'reader':'PyArrow standard Parquet reader; exact binary input coordinates retained','rois':{}}
    for d in sorted(source.iterdir()):
        if not d.is_dir() or not (d/'masks.npz').exists():continue
        meta=json.loads((d/'metadata.json').read_text());mask=np.load(d/'masks.npz')['cell_mask']
        if (d/'transcripts.parquet').exists():df=pq.read_table(d/'transcripts.parquet').to_pandas()
        else:
            import pandas as pd
            df=pd.read_csv(d/'transcripts.csv.gz',float_precision='round_trip')
        if 'qv' in df:df=df[df.qv>=20].copy()
        if 'is_gene' in df:df=df[df.is_gene].copy()
        scale=F(str(meta['transform'][0][0]));x0,y0,*_=meta['bbox_pixels']
        ids={str(v):int(k) for k,v in meta['label_to_cell_id'].items()}
        counts=defaultdict(Counter);records=[];disagreements=[]
        for z in df.itertuples(index=False):
            x,y=F(float(z.x_location)),F(float(z.y_location));col=math.floor(x*scale)-x0;row=math.floor(y*scale)-y0
            if not 0<=row<mask.shape[0] or not 0<=col<mask.shape[1]:raise ValueError('outside crop')
            cell=int(mask[row,col]);supplied=ids.get(str(z.cell_id),0)
            if cell!=supplied:disagreements.append(str(z.transcript_id))
            counts[str(cell)][str(z.feature_name)]+=1
            records.append({'id':str(z.transcript_id),'gene':str(z.feature_name),'cell':cell,
                            'pixel_row':row,'pixel_col':col,'x_um':str(float(x)),'y_um':str(float(y))})
        if disagreements:raise ValueError(f'{d.name}: {len(disagreements)} source-mask mismatches')
        result['rois'][d.name]={'counts':{k:dict(v) for k,v in counts.items()},'records':records,
                               'border_labels':list(map(int,np.unique(np.r_[mask[0],mask[-1],mask[:,0],mask[:,-1]]))),
                               'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in d.iterdir() if p.is_file()}}
    result['checked_source_records']=sum(len(v['records']) for v in result['rois'].values())
    (ROOT/'data/observed_crop_counts.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Prepared',result['checked_source_records'],'source records using PyArrow')
if __name__=='__main__':prepare(sys.argv[1] if len(sys.argv)>1 else None)
