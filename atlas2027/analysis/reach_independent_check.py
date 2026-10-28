"""Independent read-only consistency check of saved reach tables and witnesses.
Does not rerun optimization or atlas extraction; it checks identities, scoring,
recorded feasible witnesses and exact score-membership logic.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import csv,gzip,json,sys
import numpy as np
from scipy.ndimage import distance_transform_edt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from atlas.model import verify_witness

def readcsv(p):
 with gzip.open(p,'rt') if p.suffix=='.gz' else p.open() as f:return list(csv.DictReader(f))
def refs(a):
 records={r['id']:r for r in a['records']};d=ROOT/'results/proseg'/a['roi']
 ins={r['transcript_id']:r for r in readcsv(d/'input.csv.gz')};assert len(ins)==len(records)
 labels={r['cell_id']:int(r['mask_label']) for r in ins.values()};cmap={r['cell']:labels[r['original_cell_id']] for r in readcsv(d/'cell-metadata.csv.gz')}
 proseg={}
 for r in readcsv(d/'transcript-metadata.csv.gz'):
  tid=r['transcript_id'];assert tid in records and tid not in proseg;assert r['gene']==records[tid]['gene']
  assert int(ins[tid]['mask_label'])==records[tid]['cell']
  for axis in ['x','y']:assert abs(float(records[tid][axis+'_um'])-float(r['observed_'+axis]))<=.002
  proseg[tid]=0 if r['assignment']=='' else cmap[r['assignment']]
 z=np.load(ROOT/'data/source_crops'/a['roi']/'masks.npz');p=float(F(a['pixel_size_um']));mask=z['cell_mask'];nuc=np.where(z['nucleus_mask']>0,mask,0);H,W=mask.shape
 dist,ix=distance_transform_edt(nuc==0, sampling=p,return_indices=True)
 nuclear={}
 for tid,r in records.items():
  row,col=r['pixel_row'],r['pixel_col']
  if min(row+.5,col+.5,H-row-.5,W-col-.5)*p>=5:nuclear[tid]=0 if dist[row,col]>5 else int(nuc[ix[0,row,col],ix[1,row,col]])
 return {'proseg':proseg,'nuclear':nuclear}
summary=json.loads((ROOT/'results/revision_reach/summary.json').read_text());checks=[]
for cap in [4,12,28]:
 d=ROOT/'results/revision_reach'/f'cap_{cap}'
 for af in sorted((d/'data').glob('atlas_roi_*.json')):
  a=json.loads(af.read_text());records={r['id']:r for r in a['records']};patches={p['id']:p for p in a['patches']}
  for name,ref in refs(a).items():
   source=Counter((records[tid]['cell'],records[tid]['gene']) for tid in ref);target=Counter((owner,records[tid]['gene']) for tid,owner in ref.items())
   queries=json.loads((d/'results'/f'{a["roi"]}_{name}_queries.json').read_text());qmap={(q['cell'],q['gene']):q for q in queries}
   changed={k for k in source.keys()|target.keys() if k[0] in a['complete_cells'] and source[k]!=target[k]}
   assert changed=={(q['cell'],q['gene']) for q in queries if q['changed']}
   with gzip.open(d/'results'/f'{a["roi"]}_{name}_witnesses.json.gz','rt') as f:local={(v['patch'],v['gene']):v for v in json.load(f)}
   physical_checks=0
   for (pid,g),sp in local.items():
    patch=patches[pid]
    for value,w in sp['witnesses'].items():
     assert F(w['area_um2'])==F(sp['minima'][value])
     verify_witness(patch,w['owners'],w['boundaries'],F(w['area_um2']));physical_checks+=1
     count=sum(r['gene']==g and r['transcript_id'] in ref and w['owners'][i]==0 for i,r in enumerate(patch['points']))
     assert count==int(value)
   globalchecks=0
   for q in queries:
    k,g=q['cell'],q['gene'];assert q['source_count']==source[k,g] and q['reference_count']==target[k,g]
    assert q['attainable']==(str(target[k,g]) in q['minima'])
    if q['attainable']:assert F(q['minimum_area_um2'])==F(q['minima'][str(target[k,g])])
    for value,cost in q['minima'].items():
     trace=q['witness_traces'][value];total=F(0);delta=0
     for pid,change in trace.items():
      patch=patches[pid];side=patch['cell_labels'].index(k);original=sum(r['gene']==g and r['transcript_id'] in ref and r['source']==0 for r in patch['points']);localcount=original+int(change)*(1 if side==0 else -1)
      total+=F(local[pid,g]['minima'][str(localcount)]);delta+=int(change)
     assert source[k,g]+delta==int(value) and total==F(cost);globalchecks+=1
   rows=readcsv(d/'results'/f'{a["roi"]}_{name}_intervals.csv');hole_cases=0
   for r in rows:
    q=qmap[int(r['cell']),r['gene']];budget=F(r['budget_um2']);scores=[int(v) for v,c in q['minima'].items() if F(c)<=budget]
    assert int(r['lower'])==min(scores) and int(r['upper'])==max(scores)
    exact=q['reference_count'] in scores
    assert (r['reference_count_attainable']=='True')==exact
    if min(scores)<=q['reference_count']<=max(scores) and not exact:hole_cases+=1
   agg=next(x for x in summary['aggregate'] if x['cap']==cap and x['reference']==name)
   checks.append({'cap':cap,'roi':a['roi'],'reference':name,'shared_transcripts':len(ref),'changed_queries':len(changed),'attainable':sum(qmap[k]['attainable'] for k in changed),'local_witness_checks':physical_checks,'global_trace_checks':globalchecks,'csv_rows':len(rows),'interval_contains_but_unattainable':hole_cases})
for cap in [4,12,28]:
 for name in ['nuclear','proseg']:
  cc=[x for x in checks if x['cap']==cap and x['reference']==name];aa=next(x for x in summary['aggregate'] if x['cap']==cap and x['reference']==name)
  assert sum(x['changed_queries'] for x in cc)==aa['changed_queries'] and sum(x['attainable'] for x in cc)==aa['attainable_unlimited']
result={'scope':'Saved-witness and mapping consistency; no optimizer or extraction rerun','failures':0,'checks':checks}
(ROOT/'results/revision_reach/independent_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
