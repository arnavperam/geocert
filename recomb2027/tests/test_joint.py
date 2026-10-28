from pathlib import Path
import sys,json,itertools
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'tests'))
from geocert_expression.budget import ExpressionEngine
from geocert_expression.joint import outcome_spectrum,predicate_status,all_optimal_assignments
from test_budget import model
from layered import Layers,Box,realize
from pinned_compiler import solve_pinned
checks=domain_checks=0
for seed in range(1000,1020):
    p,r=model(seed);e=ExpressionEngine(p)
    channels=[]
    for gene in [0,1]:channels.append([[int(g%2==gene and o==0) for o in range(p.m+1)] for g in range(p.n)])
    got=outcome_spectrum(e,channels);truth={}
    for own in itertools.product(range(p.m+1),repeat=p.n):
        a=realize(p,list(own))
        if a['feasible']:
            v=tuple(sum(row[o] for row,o in zip(s,own)) for s in channels)
            truth[v]=min(a['area'],truth.get(v,a['area']))
    assert {tuple(r['counts']):F(r['area']) for r in got['outcomes']}==truth;checks+=len(truth)
    costs=[[-int(g%2==0 and o==0) for o in range(p.m+1)] for g in range(p.n)]
    sol=solve_pinned(e.compiled,costs,area_weight=1);domains=all_optimal_assignments(e.compiled,sol)
    opt=sol['objective'];minimizers=[]
    for own in itertools.product(range(p.m+1),repeat=p.n):
        a=realize(p,list(own))
        if a['feasible'] and a['area']+sum(row[o] for row,o in zip(costs,own))==opt:minimizers.append(own)
    assert domains['marginal_owner_domains']==[sorted({o[g] for o in minimizers}) for g in range(p.n)]
    domain_checks+=p.n
# Opposing assignments of two distinct genes: individually possible, jointly
# excluded by minimum thickness. No biological signaling assertion is made.
p=Layers([[(0,0),(10,0)],[(0,1),(10,1)]],[Box(F(2,5),F(2,5),5,5,0),Box(F(3,5),F(3,5),5,5,1)],
         [1,1],1,[F(3,5)],[None],clearance=F(1,20),bands=[(-1,1),(0,2)],mode='nominal')
e=ExpressionEngine(p);channels=[[[1,0,0],[0,0,0]],[[0,0,0],[0,0,1]]]
sp=outcome_spectrum(e,channels);st=predicate_status(sp,1,lambda x:x[0]>=1 and x[1]>=1)
assert st['status']=='never' and any(r['counts'][0]==1 for r in sp['outcomes']) and any(r['counts'][1]==1 for r in sp['outcomes'])
(ROOT/'results/joint_presence_example.json').write_text(json.dumps({'spectrum':sp,'status':st,'interpretation':'two individually possible gene transfers cannot coexist'},default=str,indent=2)+'\n')
res={'models':20,'joint_outcome_minima':checks,'all_optimal_owner_domains':domain_checks,'marginal_joint_counterexample':True,'failures':0}
(ROOT/'results/joint_tests.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
