"""Project unchanged baseline spectra onto the reach study's common budget grid."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import csv,json,os,sys,time
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from atlas.model import make_patch,verify_witness
BUDGETS=list(map(F,json.loads((R/'results/reach/plan_identity.json').read_text())['budgets_um2']))
def main():
 for ds in ['lung','breast']:
  base=R if ds=='lung' else R/'external';singles=defaultdict(list);checked=0;t=time.time()
  for f in sorted((base/'data').glob('atlas_*.json')):
   a=json.loads(f.read_text())
   for p in a['patches']:
    if not p['points']:continue
    e=make_patch(p);scale=F(p['area_scale_um2'])
    for i,r in enumerate(p['points']):
     q=e.clipped({i:1-r['source']});checked+=1
     if not q['feasible']:continue
     cost=F(q['area'])*scale;verify_witness(p,q['owners'],q['boundaries'],cost)
     for k in p['cell_labels']:
      if k>0:singles[(a['roi'],k,r['gene'])].append((p['cell_labels'][r['source']],cost))
  rows=[]
  for s in json.loads((base/'results/cell_spectra.json').read_text()):
   n=s['source_count'];key=s['roi'],s['cell'],s['gene'];values={int(v):F(d) for v,d in s['minima'].items()}
   for b in BUDGETS:
    feasible=[v for v,d in values.items() if d<=b];lo,hi=min(feasible),max(feasible)
    low=n-sum(c==s['cell'] and d<=b for c,d in singles[key]);high=n+sum(c!=s['cell'] and d<=b for c,d in singles[key]);assert low<=lo<=hi<=high
    rows.append({'atlas_id':s['atlas_id'],'roi':s['roi'],'cell':s['cell'],'gene':s['gene'],'budget_um2':str(b),'source_count':n,'lower':lo,'upper':hi,'status':'exact','mutable_patches':len(s['patches']),'independent_lower':low,'independent_upper':high,'lower_cost_um2':str(values[lo]),'upper_cost_um2':str(values[hi])})
  old=list(csv.DictReader((base/'results/cell_expression_intervals.csv').open()))
  index={(z['roi'],str(z['cell']),z['gene'],z['budget_um2']):z for z in rows}
  for z in old:
   q=index[(z['roi'],z['cell'],z['gene'],z['budget_um2'])]
   for k,v in z.items():assert str(q[k])==v,(ds,k,q[k],v)
  out=R/'results/reach'/f'baseline_{ds}_fine_intervals.csv'
  with out.open('w',newline='') as f:
   w=csv.DictWriter(f,rows[0]);w.writeheader();w.writerows(rows)
  result={'dataset':ds,'rows':len(rows),'old_rows_unchanged':len(old),'singleton_queries':checked,'seconds':time.time()-t,'changes':'Only budget grid; original geometry, spectra and original CSV unchanged.'}
  (out.with_suffix('.json')).write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)
if __name__=='__main__':main()
