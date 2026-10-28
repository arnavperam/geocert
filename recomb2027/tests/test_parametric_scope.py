from pathlib import Path
import sys,json
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'tests'))
from test_budget import model
from geocert_expression.budget import ExpressionEngine
from pinned_compiler import solve_pinned
from layered import Layers,Box
checks=0
for seed in range(1000,1020):
    p,_=model(seed);e=ExpressionEngine(p);last=set()
    for reward in [F(0),F(1,10),F(1,2),F(1),F(2),F(10)]:
        costs=[[-reward if owner==0 and g%2==0 else 0 for owner in range(p.m+1)] for g in range(p.n)]
        ans=solve_pinned(e.compiled,costs,area_weight=1,reconstruct=False)
        side=set(ans['certificate']['source_side']);assert last<=side;last=side;checks+=1
p=Layers([[(0,0),(10,0)],[(0,2),(10,2)]],[Box(-F(2,5),-F(2,5),5,5,0),Box(F(12,5),F(12,5),5,5,1)],
         [0,2],1,[F(1,2)],[F(7,2)],clearance=F(1,20),bands=[(-1,1),(1,3)],mode='nominal')
e=ExpressionEngine(p);answers=[]
for t in [0,1]:
 a=solve_pinned(e.compiled,[[0,-t,0] for _ in range(2)],area_weight=1,reconstruct=False)
 answers.append({'reward':t,'owners':a['owners'],'bits':p.encode(a['owners'])})
a,b=[{i for i,v in enumerate(r['bits']) if v} for r in answers]
assert not a<=b and not b<=a
out={'leftmost_cell_nested_cut_checks':checks,'middle_cell_counterexample':answers,'failures':0,
     'scope':'The standard monotone-terminal conditions apply to extremal-cell nonnegative count rewards after scaling. Arbitrary middle-cell scores are not automatically covered; no optimized parametric-flow implementation is claimed.'}
(ROOT/'results/parametric_scope_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
