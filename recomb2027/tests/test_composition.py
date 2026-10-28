"""Disjoint-module composition and the succinct weighted-count knapsack embedding."""
from pathlib import Path
import sys,json,itertools,random,time
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine,convolve_score_tables
from layered import Layers,Box
rng=random.Random(73101);cache={};checks=0;owners=0;t0=time.perf_counter()
for case in range(120):
    n=rng.randrange(2,9);weights=[rng.randrange(1,10) for _ in range(n)];profits=[rng.randrange(1,13) for _ in range(n)]
    parts=[]
    for a,v in zip(weights,profits):
        if (a,v) not in cache:
            p=Layers([[(-a,0),(a,0)]],[Box(F(9,10),F(9,10),0,0,0)],[1],F(1,a),[],clearance=F(1,10),bands=[(-2,2)],mode='nominal')
            sp=ExpressionEngine(p).spectrum([[v,0]])
            assert sp.minima=={0:F(0),v:F(a)};cache[(a,v)]=sp
        parts.append(cache[(a,v)])
    joint=convolve_score_tables(parts);truth={}
    for z in itertools.product([0,1],repeat=n):
        cost=sum(a*x for a,x in zip(weights,z));value=sum(v*x for v,x in zip(profits,z));truth[value]=min(truth.get(value,F(cost)),F(cost));owners+=1
    assert joint==truth
    for budget in range(sum(weights)+1):
        assert max(v for v,d in joint.items() if d<=budget)==max(v for v,d in truth.items() if d<=budget);checks+=1
out={'instances':120,'independent_owner_combinations':owners,'budget_checks':checks,'local_geometric_embeddings':len(cache),'failures':0,'seconds':time.perf_counter()-t0,
     'scope':'disjoint modules with certified compatible cores; not stitching the overlapping native corridors',
     'hardness_scope':'succinct binary-encoded integer transcript multiplicities; not a hardness proof for fixed-label planar unit-transcript segmentation'}
(ROOT/'results/composition_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
