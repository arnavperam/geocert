from pathlib import Path
import sys,json,itertools,random
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.raster import label_candidates
rng=random.Random(8321);checks=0
for d,shape in [(2,(6,7)),(3,(4,5,6))]:
    labels=np.array([rng.randrange(5) for _ in range(int(np.prod(shape)))],dtype=np.int64).reshape(shape)
    pts=[tuple(F(rng.randrange(n)*2+1,2) for n in shape) for _ in range(20)]
    for radius in [F(0),F(1,4),F(5,4)]:
        result=label_candidates(labels,pts,['G']*len(pts),spacing=[1]*d,radius=radius)
        for p,row in zip(pts,result['records']):
            for idx in itertools.product(*(range(n) for n in shape)):
                q=tuple(F(2*j+1,2) for j in idx)
                if max(abs(a-b) for a,b in zip(p,q))<=radius:
                    assert int(labels[idx]) in row['candidates'];checks+=1
        if radius==0:assert all(row['candidates']==[row['source']] for row in result['records'])
out={'dimensions_tested':[2,3],'transport_membership_checks':checks,'zero_radius_identity':True,'failures':0}
(ROOT/'results/raster_adapter_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
