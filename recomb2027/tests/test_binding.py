from pathlib import Path
import sys,json,copy
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'tests')]
from test_budget import model
from geocert_expression.budget import ExpressionEngine
from geocert_expression.verification import verify_penalized
from layered import Layers,Box
p,_=model(1000);e=ExpressionEngine(p);scores=[[int(g==0 and c==0) for c in range(2)] for g in range(p.n)]
a=e.penalized(scores,1);assert verify_penalized(e,scores,a)
mutations=[]
for key,value in [('score',999),('objective','999'),('area','999'),('multiplier','999'),('query_key','wrong')]:
 b=copy.deepcopy(a);b[key]=value;mutations.append(b)
b=copy.deepcopy(a);b['owners'][0]=1-b['owners'][0];mutations.append(b)
b=copy.deepcopy(a);b['certificate']['arcs'][0]=list(b['certificate']['arcs'][0]);b['certificate']['arcs'][0][2]+=1;mutations.append(b)
rejected=0
for b in mutations:
 try:verify_penalized(e,scores,b)
 except (AssertionError,ValueError):rejected+=1
assert rejected==len(mutations)
other=[[int(g==1 and c==0) for c in range(2)] for g in range(p.n)]
try:e.enclosure(other,1,queries=[a])
except ValueError:rejected+=1
else:raise AssertionError('wrong-score proof accepted')
# A source-inconsistent region can require positive deformation even at source
# owners. The exact root bound remains valid; source-clipping refinement is refused.
p=Layers([[(0,0),(2,0)]],[Box(-F(1,10),F(1,10),1,1,0)],[0],1,[],clearance=F(1,20),bands=[(-1,1)],mode='universal')
e=ExpressionEngine(p);assert e.source_engine is None
assert e.enclosure([[1,0]],0)['status']=='empty_uncertainty_set'
b=e.enclosure([[1,0]],1,refine=True);assert (b['lower'],b['upper'])==(0,1)
out={'deliberate_changes_rejected':rejected,'positive_minimum_deformation_handled':True,'failures':0}
(ROOT/'results/binding_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
