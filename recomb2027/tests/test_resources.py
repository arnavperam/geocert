from pathlib import Path
from fractions import Fraction as F
from itertools import product
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'tests')]
from test_budget import model
from geocert_expression.budget import ExpressionEngine
from geocert_expression.resources import resource_spectrum,resource_bound,changed_weight
from layered import realize
checks=0;bounds=0;patterns=0
for seed in range(1010,1030):
 p,r=model(seed);engine=ExpressionEngine(p);w=[1+(g%2) for g in range(p.n)]
 score=[[int(g<3 and c==0) for c in range(p.m+1)] for g in range(p.n)]
 truth=[]
 for oo in product(range(p.m+1),repeat=p.n):
  patterns+=1;a=realize(p,list(oo))
  if a['feasible']:truth.append((sum(row[o] for row,o in zip(score,oo)),a['area'],changed_weight(engine,oo,w)))
 for k in [0,1,3,6]:
  sp=resource_spectrum(engine,score,w,k)
  exact={}
  for v,d,rr in truth:
   if rr<=k:exact[v]=min(exact.get(v,d),d)
  assert sp.minima==exact;checks+=1
  for delta in [F(0),F(1,4),F(2)]:
   vals=[v for v,d,rr in truth if d<=delta and rr<=k]
   for sign in [1,-1]:
    a=resource_bound(engine,score,w,delta,k,F(1,2),F(2,3),sign=sign)
    assert a['value']>=max(vals) if sign==1 else a['value']<=min(vals);bounds+=1
out={'models':20,'full_owner_patterns':patterns,'exact_two_budget_spectra':checks,'rational_two_resource_bounds':bounds,'failures':0}
(ROOT/'results/resource_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
