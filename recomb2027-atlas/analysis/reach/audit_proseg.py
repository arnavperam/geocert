"""Verify the supplied Proseg execution artifacts against the original observations.

This reads existing stochastic runs; it does not execute Proseg or treat its
assignments as ground truth. Background overrides the overlapping cell index.
"""
from pathlib import Path
from collections import Counter
import csv,gzip,json,hashlib,math
R=Path(__file__).resolve().parents[2]
def rows(path):
 with gzip.open(path,'rt') as f:return list(csv.DictReader(f))
def main():
 observed=json.loads((R/'data/observed_crop_counts.json').read_text())['rois'];all_results=[]
 for roi in ('roi_01','roi_02'):
  d=R/'results/proseg'/roi;run=json.loads((d/'run.json').read_text())
  inp=rows(d/'input.csv.gz');tr=rows(d/'transcript-metadata.csv.gz');cells=rows(d/'cell-metadata.csv.gz');genes=rows(d/'gene-metadata.csv.gz')
  if hashlib.sha256((d/'input.csv.gz').read_bytes()).hexdigest()!=run['input_csv_sha256']:raise ValueError('input hash mismatch')
  original={r['id']:r for r in observed[roi]['records']};inputs={r['transcript_id']:r for r in inp}
  assert len(inputs)==len(inp)==len(original) and set(inputs)==set(original)
  cellmap={}
  for ident,r in inputs.items():
   o=original[ident]
   assert r['feature_name']==o['gene'] and int(r['mask_label'])==o['cell']
   for axis in ('x','y'):assert abs(float(r[f'{axis}_location'])-float(o[f'{axis}_um']))<=0.001
   if r['cell_id'] in cellmap:assert cellmap[r['cell_id']]==o['cell']
   cellmap[r['cell_id']]=o['cell']
  output_to_reference={r['cell']:cellmap[r['original_cell_id']] for r in cells}
  target={};gene_counts=Counter();max_resid=0.
  for r in tr:
   ident=r['transcript_id'];assert ident in original and ident not in target
   o=original[ident];assert r['gene']==o['gene']
   for axis in ('x','y'):
    err=abs(float(r[f'observed_{axis}'])-float(o[f'{axis}_um']));max_resid=max(max_resid,err);assert err<=0.001
   background=r['background'].lower()=='true'
   owner=0 if background else output_to_reference[r['assignment']]
   target[ident]=owner
   if not background:gene_counts[(int(r['assignment']),r['gene'])]+=1
  excluded=sorted(set(original)-set(target));assert excluded==sorted(run['comparison']['excluded_source_ids'])
  # Check the saved sparse count matrix against transcript-level final assignments.
  gene_names=[g['gene'] for g in genes];matrix={}
  with gzip.open(d/'counts.mtx.gz','rt') as f:
   lines=(s for s in f if not s.startswith('%'));shape=list(map(int,next(lines).split()))
   for line in lines:
    i,j,v=line.split();matrix[(int(i)-1,gene_names[int(j)-1])]=float(v)
  discrepancies=[(k,v,gene_counts[k]) for k,v in matrix.items() if v!=gene_counts[k]]
  discrepancies += [(k,0,v) for k,v in gene_counts.items() if k not in matrix]
  assert not discrepancies,discrepancies[:3]
  changes=sum(target[i]!=original[i]['cell'] for i in target)
  assert changes==run['comparison']['changed_owner']
  result={'roi':roi,'source_records':len(original),'retained_output_records':len(target),'excluded_records':len(excluded),'changed_owner':changes,'background_records':sum(v==0 for v in target.values()),'sparse_count_entries_verified':len(matrix),'maximum_observed_xy_rounding_residual_um':max_resid,'failures':0,'source_to_output_cell_map':output_to_reference,'target_owners':target,'excluded_ids':excluded,'file_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(d.iterdir()) if p.is_file()}}
  (R/'results/reach'/f'proseg_join_{roi}.json').write_text(json.dumps(result,indent=2)+'\n')
  all_results.append({k:v for k,v in result.items() if k not in ['target_owners','excluded_ids','source_to_output_cell_map','file_hashes']})
 (R/'results/reach/proseg_audit.json').write_text(json.dumps({'runs':all_results,'scope':'Independent ingestion and arithmetic verification of supplied runs, not a fresh Proseg execution, five-run replication or accuracy validation.'},indent=2)+'\n')
 print(all_results)
if __name__=='__main__':main()
