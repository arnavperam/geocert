"""Audit Proseg outputs and compare paired runs on the same retained core IDs."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import argparse,csv,gzip,json,sys,hashlib
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from atlas.model import make_patch

def read(p):
 with gzip.open(p,'rt') as f:return list(csv.DictReader(f))
def audit(roi,dest):
 import pyarrow.dataset as pads
 run=json.loads((dest/'run.json').read_text())
 if run['returncode']!=0:return None
 inp=roi/'context_transcripts.csv.gz' if run['input_side_um']==800 else roi/'data/source_crops'/roi.name/'transcripts.csv.gz'
 if (dest/'input.csv.gz').exists():inp=dest/'input.csv.gz'
 source={r['transcript_id']:r for r in read(inp)};assert hashlib.sha256(inp.read_bytes()).hexdigest()==run['input_sha256']
 identity=json.loads((roi/'context_identity_map.json').read_text()) if (roi/'context_identity_map.json').exists() else None
 cellmap={r['cell']:(identity[r['original_cell_id']] if identity is not None else int(r['original_cell_id'])) for r in read(dest/'cell-metadata.csv.gz')}
 counts=Counter();target={};observedresid=0
 native=dest/'proseg-output.zarr/points/transcripts/points.parquet'
 assert native.exists(),'Complete native transcript output required'
 def output_rows():
  for batch in pads.dataset(native,format='parquet').scanner(batch_size=65536).to_batches():
   for row in batch.to_pylist():yield {k:str(v) for k,v in row.items()}
 for r in output_rows():
  ident=r['transcript_id'];assert ident in source and ident not in target;assert r['gene']==source[ident]['feature_name']
  observedresid=max(observedresid,*[abs(float(r['observed_'+x])-float(source[ident][x+'_location'])) for x in ['x','y']]);assert observedresid<=.002
  owner=0 if r['background'].lower()=='true' else cellmap[r['assignment']];target[ident]=owner
  if owner:counts[(int(r['assignment']),r['gene'])]+=1
 genes=[r['gene'] for r in read(dest/'gene-metadata.csv.gz')];matrix={}
 with gzip.open(dest/'counts.mtx.gz','rt') as f:
  lines=(s for s in f if not s.startswith('%'));shape=list(map(int,next(lines).split()))
  for l in lines:
   i,j,v=l.split();matrix[(int(i)-1,genes[int(j)-1])]=float(v)
 assert dict(counts)=={k:v for k,v in matrix.items() if v!=0},'count matrix/assignment discrepancy'
 csv_complete=True
 try:
  with gzip.open(dest/'transcript-metadata.csv.gz','rb') as f:
   while f.read(1048576):pass
 except EOFError:csv_complete=False
 record={'transcript_output_source':'native spatialdata Parquet, independently reconciled with counts.mtx.gz','optional_csv_export_complete':csv_complete,'input_records':len(source),'output_records':len(target),'excluded_records':len(source.keys()-target.keys()),'count_matrix_nonzeros_verified':len(counts),'max_observed_xy_residual_um':observedresid,'failures':0}
 (dest/'audit.json').write_text(json.dumps(record,indent=2)+'\n');return target

def main(roi,saved_core=False):
 atlas=json.loads(next((roi/'data').glob('atlas_*.json')).read_text());source={r['id']:r for r in atlas['records']};complete=set(atlas['complete_cells']);states={(int(r['cell']),r['gene']):r for r in json.loads((roi/'results/cell_spectra.json').read_text())};targets={}
 for dest in sorted((roi/'proseg').iterdir()):
  if not (dest/'run.json').exists():continue
  run=json.loads((dest/'run.json').read_text())
  if run['returncode']!=0:continue
  if saved_core:
   with gzip.open(dest/'core_assignments.json.gz','rt') as f:t=json.load(f)
   assert set(t)<=source.keys()
  else:
   t=audit(roi,dest)
   with gzip.open(dest/'core_assignments.json.gz','wt') as f:json.dump({i:t[i] for i in source if i in t},f,sort_keys=True,separators=(',',':'))
  targets[dest.name]=t
 common=set(source)
 for t in targets.values():common &= t.keys()
 removed=set(source)-common;removedcounts=Counter((source[i]['cell'],source[i]['gene']) for i in removed)
 excluded=set();patchbyid={}
 for p in atlas['patches']:
  for i,r in enumerate(p['points']):
   patchbyid[r['transcript_id']]=(p,i)
   if r['transcript_id'] in removed:excluded.update((k,r['gene']) for k in p['cell_labels'] if k>0)
 src=Counter((source[i]['cell'],source[i]['gene']) for i in common);summary=[];outrows=[]
 fixed_keys=set(states)|{(k,g) for k,g in src if k in complete}|{(target[i],source[i]['gene']) for target in targets.values() for i in common if target[i] in complete}
 for name,target in targets.items():
  dst=Counter((target[i],source[i]['gene']) for i in common);transfers=Counter();changedfixed=unavailable=badpins=0;cost=F(0)
  for i in common:
   old,new=source[i]['cell'],target[i]
   if old!=new:
    transfers['cell_to_cell' if old and new else 'to_background' if old else 'from_background']+=1
    if i not in patchbyid:changedfixed+=1
    elif new not in patchbyid[i][0]['cell_labels']:unavailable+=1
  for p in atlas['patches']:
   pins={i:p['cell_labels'].index(target[r['transcript_id']]) for i,r in enumerate(p['points']) if r['transcript_id'] in common and target[r['transcript_id']] in p['cell_labels']}
   if not pins:continue
   ans=make_patch(p).clipped(pins)
   if not ans['feasible']:badpins+=1
   else:cost+=F(ans['area'])*F(p['area_scale_um2'])
  keys=fixed_keys;cnt=Counter();bgkeys=set()
  for i in common:
   old,new=source[i]['cell'],target[i]
   if old!=new and (not old or not new):bgkeys.add((old or new,source[i]['gene']))
  for k,g in sorted(keys):
   if (k,g) in excluded:cnt['unavailable_queries_due_to_excluded_mutable_RNA']+=1;continue
   observed,wanted=src[(k,g)],dst[(k,g)];s=states.get((k,g));spec={int(v)-removedcounts[(k,g)]:F(d) for v,d in s['minima'].items()} if s else {observed:F(0)};assert observed in spec
   amin=spec.get(wanted);included=amin is not None and amin<=4;cnt['queries']+=1;changed=observed!=wanted;cnt['changed_counts']+=changed;cnt['changed_counts_reachable_at4']+=changed and included;cnt['changed_counts_involving_background']+=changed and (k,g) in bgkeys
   cnt['invariant_counts']+=len([v for v,d in spec.items() if d<=4])==1;cnt['invariant_counts_disagreeing_with_proseg']+=len([v for v,d in spec.items() if d<=4])==1 and changed
   outrows.append({'roi':roi.name,'run':name,'cell':k,'gene':g,'source_count':observed,'alternative_count':wanted,'changed':changed,'involves_background_transfer':(k,g) in bgkeys,'attainable_at4':included,'minimum_area_um2':str(amin) if amin is not None else 'inf'})
  summary.append({'roi':roi.name,'run':name,'common_retained_core_RNA':len(common),'core_input_RNA':len(source),'common_excluded_core_RNA':len(removed),'background_core_RNA':sum(target[i]==0 for i in common),'changed_core_RNA':sum(transfers.values()),'transfer_categories':dict(transfers),**cnt,'changed_fixed_RNA':changedfixed,'unavailable_destinations':unavailable,'infeasible_pinned_patches':badpins,'full_retained_assignment_representable':not(changedfixed or unavailable or badpins),'feasible_patch_area_um2':str(cost)})
 (roi/'results/proseg_comparison.json').write_text(json.dumps({'scope':'Descriptive alternative-segmentation comparison, not accuracy or biological calibration. All completed variants use the same retained core transcript IDs. Unavailable queries remain counted.','runs':summary},indent=2)+'\n')
 with gzip.open(roi/'results/proseg_count_comparison.csv.gz','wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(outrows[0]));w.writeheader();w.writerows(outrows)
 print(json.dumps(summary),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('roots',nargs='+',type=Path);ap.add_argument('--saved-core',action='store_true',help='Recompute core comparisons from the published assignment extracts; does not repeat the native-output/count-matrix audit.');a=ap.parse_args()
 for r in a.roots:main(r,a.saved_core)
