"""Reconstruct observed crop counts, retaining transcript provenance."""
from pathlib import Path
import sys,json,collections,hashlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'validation/vendor'))
from parquet_minimal import read_flat_parquet
import numpy as np
out={};checked=0
for d in sorted((ROOT/'data/source_crops').iterdir()):
    df=read_flat_parquet(d/'transcripts.parquet');mask=np.load(d/'masks.npz')['cell_mask']
    counts=collections.defaultdict(collections.Counter);records=[]
    border=set(map(int,np.unique(np.r_[mask[0],mask[-1],mask[:,0],mask[:,-1]])))
    for row in df.itertuples(index=False):
        r,c=int(row.pixel_row),int(row.pixel_col);owner=int(row.mask_label)
        if int(mask[r,c])!=owner:raise ValueError('record/mask mismatch')
        gene=str(row.feature_name);tid=str(int(row.transcript_id))
        counts[owner][gene]+=1;checked+=1
        records.append({'id':tid,'gene':gene,'cell':owner,'pixel_row':r,'pixel_col':c,
                        'x_um':str(float(row.x_location)),'y_um':str(float(row.y_location))})
    out[d.name]={'counts':{str(c):dict(sorted(v.items())) for c,v in counts.items()},
                 'records':records,'border_labels':sorted(border),
                 'source_hashes':{f:hashlib.sha256((d/f).read_bytes()).hexdigest() for f in ['masks.npz','transcripts.parquet','metadata.json']}}
charts=json.loads((ROOT/'data/charts.json').read_text());n=0
for ch in charts:
    byid={r['id']:r for r in out[ch['roi']]['records']}
    for p in ch['points']:
        row=byid[p['transcript_id']]
        assert row['cell']==ch['cell_labels'][p['source']] and row['gene']==p['gene'];n+=1
payload={'scope':'Observed counts within supplied crops; exterior-to-corridor counts fixed. Border-clipped cells are not claimed complete whole-cell profiles.',
         'reader':'Inherited restricted Parquet decoder; shared ingestion dependency.',
         'checked_source_records':checked,'checked_corridor_records':n,'rois':out}
(ROOT/'data/observed_crop_counts.json').write_text(json.dumps(payload,indent=2)+'\n')
print({'source_records':checked,'corridor_records':n})
