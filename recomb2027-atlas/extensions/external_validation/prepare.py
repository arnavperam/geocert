"""Extract prespecified regions from native Xenium masks and measured transcripts."""
from pathlib import Path
import argparse,json,hashlib,math,gzip,csv,time
from fractions import Fraction as F
from collections import Counter,defaultdict
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import zarr

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--only-roi');ap.add_argument('--context-only',action='store_true');args=ap.parse_args()
 root=Path(__file__).resolve().parent;protocol=json.loads((root/'protocol.json').read_text());donor=args.raw.name
 dest=args.output/donor;dest.mkdir(parents=True,exist_ok=True)
 cells=pq.read_table(args.raw/'cells.parquet');ids=cells['cell_id'].to_numpy();xx=cells['x_centroid'].to_numpy();yy=cells['y_centroid'].to_numpy()
 store=zarr.storage.ZipStore(str(args.raw/'cells.zarr.zip'),mode='r');group=zarr.open_group(store=store,mode='r');transform=group['masks/homogeneous_transform'][:].tolist();p=1/F(float(transform[0][0]));h,w=group['masks/1'].shape
 assert np.array_equal(group['cell_id'][:],ids);n=math.floor(F(200)/p);nc=math.floor(F(800)/p)
 regions=[]
 for i,q in enumerate([.25,.5,.75],1):
  target=[float(np.quantile(xx,q)),float(np.quantile(yy,.5))];j=int(np.argmin((xx-target[0])**2+(yy-target[1])**2));cx,cy=float(xx[j]),float(yy[j]);x0,y0=math.floor((F(cx)-100)/p),math.floor((F(cy)-100)/p);xc,yc=math.floor((F(cx)-400)/p),math.floor((F(cy)-400)/p)
  assert 0<=x0<x0+n<=w and 0<=y0<y0+n<=h,'analysis crop outside raster'
  context=[max(0,xc),max(0,yc),min(w,xc+nc),min(h,yc+nc)]
  regions.append({'roi':f'{donor.lower()}_{i:02d}','target_um':target,'nearest_cell_id':str(ids[j]),'center_um':[cx,cy],'bbox_pixels':[x0,y0,x0+n,y0+n],'context_bbox_pixels':context,'context_clipped':context!=[xc,yc,xc+nc,yc+nc]})
 selection={'donor':donor,'protocol_sha256':sha(root/'protocol.json'),'pixel_size_um':str(p),'regions':regions,'selected_before_transcript_read':True,'mask_shape':[h,w]};(dest/'selection.json').write_text(json.dumps(selection,indent=2)+'\n');print('SELECTION',json.dumps(selection),flush=True)
 genes={t['type']['data']['name'] for t in json.loads((args.raw/'gene_panel.json').read_text())['payload']['targets'] if t['type']['descriptor']=='gene'}
 chunks=[[] for _ in regions];pf=pq.ParquetFile(args.raw/'transcripts.parquet');print('SCHEMA',str(pf.schema_arrow),flush=True);start=time.time()
 for bi,b in enumerate(pf.iter_batches(batch_size=262144)):
  x=b.column(b.schema.get_field_index('x_location')).to_numpy().astype(np.float64);y=b.column(b.schema.get_field_index('y_location')).to_numpy().astype(np.float64);qv=b.column(b.schema.get_field_index('qv')).to_numpy()
  for i,r in enumerate(regions):
   x0,y0,x1,y1=[float(p*v) for v in r['context_bbox_pixels']];take=np.flatnonzero((x>=x0-float(p))&(x<x1+float(p))&(y>=y0-float(p))&(y<y1+float(p))&(qv>=20))
   if len(take):chunks[i].append(pa.Table.from_batches([b]).take(pa.array(take)))
  if bi%40==0:print(donor,'batch',bi,'seconds',round(time.time()-start,1),flush=True)
 for r,cc in zip(regions,chunks):
  if args.only_roi and r['roi']!=args.only_roi:continue
  roi=r['roi'];out=dest/roi;crop=out/'data/source_crops'/roi;crop.mkdir(parents=True,exist_ok=True);(out/'results').mkdir(exist_ok=True);(out/'docs').mkdir(exist_ok=True)
  table=pa.concat_tables(cc);names=[g.decode() if isinstance(g,bytes) else str(g) for g in table['feature_name'].to_pylist()];take=[i for i,g in enumerate(names) if g in genes];table=table.take(pa.array(take,type=pa.int64()));raw=table.to_pydict();raw['feature_name']=[g.decode() if isinstance(g,bytes) else str(g) for g in raw['feature_name']]
  x0,y0,x1,y1=r['bbox_pixels'];xc,yc,xe,ye=r['context_bbox_pixels'];cm=group['masks/1'][yc:ye,xc:xe];nm=group['masks/0'][yc:ye,xc:xe];cell=cm[y0-yc:y1-yc,x0-xc:x1-xc];nuc=nm[y0-yc:y1-yc,x0-xc:x1-xc]
  if not args.context_only:np.savez_compressed(crop/'masks.npz',cell_mask=cell,nucleus_mask=nuc)
  labelmap={str(int(k)):str(ids[int(k)-1]) for k in np.unique(cm) if k>0};rows=[];records=[];counts=defaultdict(Counter);disagreements=0;contextdisagreements=0;contextrows=0
  fieldnames=['transcript_id','cell_id','overlaps_nucleus','feature_name','x_location','y_location','z_location','qv','mask_label','source_cell_id']
  with gzip.open(out/'context_transcripts.csv.gz','wt',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=fieldnames);writer.writeheader()
   for i in range(len(table)):
    x,y=float(raw['x_location'][i]),float(raw['y_location'][i]);xp,yp=math.floor(F(x)/p),math.floor(F(y)/p)
    if not (xc<=xp<xe and yc<=yp<ye):continue
    label=int(cm[yp-yc,xp-xc]);owner='UNASSIGNED' if label==0 else labelmap[str(label)];vendor=str(raw['cell_id'][i]);vendor='UNASSIGNED' if vendor in ['-1','0','UNASSIGNED'] else vendor
    nuclear=int(nm[yp-yc,xp-xc]);ident=str(raw['transcript_id'][i]);gene=str(raw['feature_name'][i]);contextdisagreements+=owner!=vendor
    row={'transcript_id':ident,'cell_id':owner,'overlaps_nucleus':int(nuclear==label and label>0),'feature_name':gene,'x_location':repr(x),'y_location':repr(y),'z_location':repr(float(raw['z_location'][i])),'qv':repr(float(raw['qv'][i])),'mask_label':label,'source_cell_id':vendor};writer.writerow(row);contextrows+=1
    if x0<=xp<x1 and y0<=yp<y1:
     rows.append(row);disagreements+=owner!=vendor;records.append({'id':ident,'gene':gene,'cell':label,'pixel_row':yp-y0,'pixel_col':xp-x0,'x_um':repr(x),'y_um':repr(y)});counts[label][gene]+=1
  if args.context_only:
   print('CONTEXT_REBUILT',roi,contextrows,flush=True);continue
  assert len({r['id'] for r in records})==len(records)
  with gzip.open(crop/'transcripts.csv.gz','wt',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=fieldnames);writer.writeheader();writer.writerows(rows)
  border=set(map(int,np.unique(np.r_[cell[0],cell[-1],cell[:,0],cell[:,-1]])))
  meta={**r,'pixel_size_um':float(p),'transform':transform,'bounds_um':[float(p*v) for v in r['bbox_pixels']],'label_to_cell_id':labelmap,'complete_source_labels':[int(k) for k in np.unique(cell) if k>0 and k not in border],'q20_gene_rows':len(records),'context_q20_gene_rows':contextrows,'source_assignment_disagreements':disagreements,'context_source_assignment_disagreements':contextdisagreements,'selection_protocol_sha256':selection['protocol_sha256'],'source_archive_sha256':None,'source_archive_hash_note':'Archive members retrieved by HTTP range; full archive not hashed. Verified member SHA-256 values are recorded.','source_provenance':json.loads((args.raw/'provenance.json').read_text()),'source_segmentation':'DAPI nuclei plus vendor 15-um expansion; native raster cell mask','scope':'External biological tissue replication; vendor segmentation is a reference, not truth'}
  (crop/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n');hashes={name:sha(crop/name) for name in ['metadata.json','masks.npz','transcripts.csv.gz']}
  obs={'scope':meta['scope'],'checked_source_records':len(records),'rois':{roi:{'counts':dict(counts),'records':records,'source_hashes':hashes}}};(out/'data/observed_crop_counts.json').write_text(json.dumps(obs,separators=(',',':'))+'\n')
  (out/'docs/ATLAS_PROTOCOL.json').write_text((root/'protocol.json').read_text())
  print('PREPARED',roi,len(records),len(table),len(meta['complete_source_labels']),'vendor_disagreements',disagreements,flush=True)
 store.close()
if __name__=='__main__':main()
