from pathlib import Path
import sys,subprocess,json,random,itertools,time
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.cpp_backend import flow,solve
from geocert_expression.budget import ExpressionEngine
from data_api import reassignment_costs,variant
from mixed import make_problem
def load_cases():
 charts=json.loads((ROOT/'data/charts.json').read_text());expected=json.loads((ROOT/'data/expected_optional.json').read_text())
 for i,ch in enumerate(charts):
  p,pins,n=make_problem(ch)
  for j,(tag,q) in enumerate([('baseline',p),('maximal',variant(p))]):
   yield i,ch,tag,q,pins,n,expected[2*i+j]
(ROOT/'run').mkdir(exist_ok=True)
binary=ROOT/'run/exact_flow';subprocess.run(['g++','-std=c++17','-O2',str(ROOT/'cpp/exact_flow.cpp'),'-o',str(binary)],check=True)
r=random.Random(92302);truthchecks=0;t0=time.perf_counter()
for case in range(100):
 n=r.randrange(2,9);edges=[]
 for u in range(n):
  for v in range(n):
   if r.random()<.2:edges.append((u,v,r.randrange(0,20)*(10**60 if case%10==0 else 1)))
 graph={'n':n,'s':0,'t':n-1,'edges':edges};cert=flow(graph,binary)
 best=min(sum(c for u,v,c in edges if u in side and v not in side) for bits in itertools.product([0,1],repeat=n-2) for side in [{0}|{i+1 for i,b in enumerate(bits) if b}])
 assert cert['value']==best;truthchecks+=1
native=0
for i,ch,tag,p,pins,n,expected in load_cases():
 e=ExpressionEngine(p,pins);costs=reassignment_costs(p,n)
 for lam in [F(1,4),F(1),F(4)]:
  a=solve(e.compiled,costs,lam,binary)
  assert a['feasible'] and a['objective']==F(expected['objectives'][str(lam)]);native+=1
out={'exhaustive_graph_optima':truthchecks,'huge_capacity_graphs':10,'native_optional_objectives':native,'independent_flow_certificate_check':True,'failures':0,'seconds_excluding_compilation':time.perf_counter()-t0,
     'scope':'standard arbitrary-precision Dinic, fresh flows; no dynamic residual reuse or state-of-the-art maxflow speed claim'}
(ROOT/'results/cpp_flow_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
