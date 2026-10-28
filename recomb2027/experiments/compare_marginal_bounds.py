"""Post-census, same-query attribution: sound marginal envelope vs joint budget.

An additive sum of singleton radii is also tested, but is explicitly a feasible
INNER approximation with common feasibility, never a sound count enclosure.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine
from mixed import make_problem
from data_api import variant

def main():
 t0=time.perf_counter();rows=[];examples=[];patterns=0
 for ch in json.loads((ROOT/'data/charts.json').read_text()):
  p0,pins,n=make_problem(ch)
  for setting,p in [('baseline',p0),('maximal',variant(p0))]:
   data=json.loads((ROOT/'results/native'/f'{ch["id"]}_{setting}.json').read_text())
   en=ExpressionEngine(p,pins);singleton={}
   for g in range(n):
    singleton[g]={}
    for o in range(p.m+1):
     a=en.clipped({g:o})
     if a['feasible']:singleton[g][o]=F(a['area'])
   for q in data['queries']:
    gs=q['local_indices'];cell=q['cell_index'];outside=q['fixed_outside_count'];added={}
    for oo in product(*(singleton[g] for g in gs)):
     patterns+=1
     if not en.source_engine.closure(dict(zip(gs,oo)))['feasible']:continue
     val=outside+sum(o==cell for o in oo);cost=sum((singleton[g][o] for g,o in zip(gs,oo)),F(0))
     added[val]=min(cost,added.get(val,cost))
    for r in data['rows']:
     if r['gene']!=q['gene'] or r['cell_index']!=cell:continue
     d=F(r['budget_um2']);domains=[[o for o,v in singleton[g].items() if v<=d] for g in gs]
     lo=outside+sum(dom==[cell] for dom in domains);hi=outside+sum(cell in dom for dom in domains)
     vals=[v for v,c in added.items() if c<=d];alo,ahi=min(vals),max(vals)
     assert lo<=r['lower']<=alo<=ahi<=r['upper']<=hi
     rec={k:r[k] for k in ['corridor','roi','setting','cell_id','cell_index','gene','budget_um2','source_count','lower','upper']}
     rec.update(marginal_lower=lo,marginal_upper=hi,additive_inner_lower=alo,additive_inner_upper=ahi)
     rows.append(rec)
     if (lo,hi)!=(r['lower'],r['upper']) or (alo,ahi)!=(r['lower'],r['upper']):examples.append(rec)
 table=[]
 for s in ['baseline','maximal']:
  rr=[r for r in rows if r['setting']==s]
  table.append({'setting':s,'questions':len(rr),
                'joint_tighter_than_sound_marginals':sum((r['lower'],r['upper'])!=(r['marginal_lower'],r['marginal_upper']) for r in rr),
                'additive_inner_misses_attainable_endpoint':sum((r['lower'],r['upper'])!=(r['additive_inner_lower'],r['additive_inner_upper']) for r in rr)})
 out={'protocol':'post-census attribution on every existing count question; no geometry/budget selection',
      'scope':'marginal model is a sound outer envelope; additive-feasible model is an inner approximation and cannot certify exclusion',
      'table':table,'target_patterns':patterns,'changed_questions':examples,'seconds':time.perf_counter()-t0,'failures':0}
 (ROOT/'results/marginal_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='changed_questions'},indent=2))
if __name__=='__main__':main()
