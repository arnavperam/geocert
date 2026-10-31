"""Independently reselect core transcripts from original Parquet and check masks."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,json,math,hashlib,csv,gzip
import pyarrow.parquet as pq
import numpy as np

def main(raw,base):
 selection=json.loads((base/raw.name/'selection.json').read_text());p=F(selection['pixel_size_um']);genes={t['type']['data']['name'] for t in json.loads((raw/'gene_panel.json').read_text())['payload']['targets'] if t['type']['descriptor']=='gene'};truth={r['roi']:{} for r in selection['regions']}
 for b in pq.ParquetFile(raw/'transcripts.parquet').iter_batches(batch_size=262144):
  x=b.column(b.schema.get_field_index('x_location')).to_numpy().astype(float);y=b.column(b.schema.get_field_index('y_location')).to_numpy().astype(float);q=b.column(b.schema.get_field_index('qv')).to_numpy()
  for r in selection['regions']:
   lo,top,hi,bot=r['bbox_pixels'];take=np.flatnonzero((x>=float((lo-1)*p))&(x<float((hi+1)*p))&(y>=float((top-1)*p))&(y<float((bot+1)*p))&(q>=20))
   for i in take:
    xx,yy=F(float(x[i])),F(float(y[i]));px,py=math.floor(xx/p),math.floor(yy/p)
    if not(lo<=px<hi and top<=py<bot):continue
    name=b.column(b.schema.get_field_index('feature_name'))[int(i)].as_py();name=name.decode() if isinstance(name,bytes) else name
    if name not in genes:continue
    ident=str(b.column(b.schema.get_field_index('transcript_id'))[int(i)].as_py());owner=b.column(b.schema.get_field_index('cell_id'))[int(i)].as_py();assert ident not in truth[r['roi']];truth[r['roi']][ident]=(name,xx,yy,px,py,owner)
 results=[]
 for r in selection['regions']:
  roi=base/raw.name/r['roi'];d=json.loads((roi/'data/observed_crop_counts.json').read_text())['rois'][r['roi']];records={x['id']:x for x in d['records']};assert records.keys()==truth[r['roi']].keys();z=np.load(roi/'data/source_crops'/r['roi']/'masks.npz');m=z['cell_mask'];meta=json.loads((roi/'data/source_crops'/r['roi']/'metadata.json').read_text());counts=Counter();disagreements=0
  for ident,(gene,x,y,px,py,vendor) in truth[r['roi']].items():
   rec=records[ident];label=int(m[py-r['bbox_pixels'][1],px-r['bbox_pixels'][0]]);assert rec['gene']==gene and F(float(rec['x_um']))==x and F(float(rec['y_um']))==y and rec['cell']==label;counts[(label,gene)]+=1;source=0 if vendor==-1 else vendor;disagreements+=source!=label
  assert dict(counts)=={(int(k),g):n for k,v in d['counts'].items() for g,n in v.items()};assert disagreements==meta['source_assignment_disagreements']
  with gzip.open(roi/'context_transcripts.csv.gz','rt') as f:context={v['transcript_id']:v for v in csv.DictReader(f)}
  assert set(records)<=context.keys();assert len(context)==meta['context_q20_gene_rows']
  result={'roi':r['roi'],'raw_source_rows_reselected':len(records),'source_counts_reconciled':len(counts),'context_rows':len(context),'mask_vendor_assignment_disagreements':disagreements,'failures':0};results.append(result);(roi/'results/source_input_audit.json').write_text(json.dumps(result,indent=2)+'\n')
 print(results,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--analysis',type=Path,required=True);a=ap.parse_args();main(a.raw,a.analysis)
