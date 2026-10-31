"""Add surrounding transcripts for the original published lung comparison cores."""
import argparse,json,gzip,csv,shutil,hashlib
from pathlib import Path
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

def main(raw,old,dest):
 genes={t['type']['data']['name'] for t in json.loads((raw/'gene_panel.json').read_text())['payload']['targets'] if t['type']['descriptor']=='gene'}
 metas={r:json.loads((old/'data/source_crops'/r/'metadata.json').read_text()) for r in ['roi_01','roi_02']};parts={r:[] for r in metas}
 for b in pq.ParquetFile(raw/'transcripts.parquet').iter_batches(batch_size=262144):
  x=b.column(b.schema.get_field_index('x_location')).to_numpy().astype(float);y=b.column(b.schema.get_field_index('y_location')).to_numpy().astype(float);q=b.column(b.schema.get_field_index('qv')).to_numpy()
  for roi,m in metas.items():
   x0,y0,x1,y1=m['selection']['bounds_um'];take=np.flatnonzero((x>=x0-300)&(x<x1+300)&(y>=y0-300)&(y<y1+300)&(q>=20))
   if len(take):parts[roi].append(pa.Table.from_batches([b]).take(pa.array(take)))
 for roi,m in metas.items():
  out=dest/roi;out.mkdir(parents=True,exist_ok=True);(out/'data/source_crops').mkdir(parents=True,exist_ok=True);(out/'results').mkdir(exist_ok=True);shutil.copytree(old/'data/source_crops'/roi,out/'data/source_crops'/roi,dirs_exist_ok=True);shutil.copy(old/'data'/f'atlas_{roi}.json',out/'data'/f'atlas_{roi}.json')
  for name in ['cell_spectra.json','patch_spectra.json']:
   d=[v for v in json.loads((old/'results'/name).read_text()) if v['roi']==roi];(out/'results'/name).write_text(json.dumps(d)+'\n')
  table=pa.concat_tables(parts[roi]);rows=table.to_pylist();fields=['transcript_id','cell_id','overlaps_nucleus','feature_name','x_location','y_location','z_location','qv','mask_label','source_cell_id']
  lookup={v:int(k) for k,v in m['label_to_cell_id'].items()};ids=sorted(set(str(r['cell_id']) for r in rows)-set(lookup)-{'UNASSIGNED'});lookup.update({v:1000000+i for i,v in enumerate(ids)});lookup['UNASSIGNED']=0;count=0
  with gzip.open(out/'context_transcripts.csv.gz','wt',newline='') as f:
   w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
   for r in rows:
    g=r['feature_name'];g=g.decode() if isinstance(g,bytes) else str(g)
    if g not in genes:continue
    cid=str(r['cell_id']);new={k:r.get(k,'') for k in fields};new.update(feature_name=g,cell_id=cid,mask_label=lookup[cid],source_cell_id=cid);w.writerow(new);count+=1
  shutil.copy(old/'results/proseg'/roi/'input.csv.gz',out/'data/source_crops'/roi/'transcripts.csv.gz')
  (out/'context_identity_map.json').write_text(json.dumps(lookup)+'\n')
  print(roi,'context_RNA',count,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--old',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();main(a.raw,a.old,a.output)
