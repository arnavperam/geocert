from pathlib import Path
import sys,json,itertools,random
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.relaxation import candidate_count_intervals,expression_presence
rng=random.Random(4091);checks=0;assignments=0
for seed in range(120):
    records=[]
    for i in range(rng.randrange(1,7)):
        src=rng.randrange(3);allowed={src}|{j for j in range(3) if rng.random()<.3}
        records.append({'source':src,'gene':rng.randrange(2),'candidates':allowed})
    for k in range(len(records)+1):
        bs=candidate_count_intervals(records,k);truth={key:[] for key in bs}
        for owners in itertools.product(*(r['candidates'] for r in records)):
            if sum(c!=r['source'] for r,c in zip(records,owners))>k:continue
            assignments+=1
            for key in truth:truth[key].append(sum((c,r['gene'])==key for r,c in zip(records,owners)))
        for key,v in truth.items():
            assert (bs[key]['lower'],bs[key]['upper'])==(min(v),max(v));checks+=1
# One RNA cannot simultaneously be in both cells, despite both upper bounds=1.
b=candidate_count_intervals([{'source':0,'gene':'G','candidates':[0,1]}])
assert expression_presence(b,[(0,'G'),(1,'G')])=='joint_presence_unresolved'
try:candidate_count_intervals([{'source':0,'gene':1,'candidates':[1]}])
except ValueError:pass
else:raise AssertionError('incompatible source accepted')
r={'instances':120,'budget_key_checks':checks,'feasible_assignments':assignments,'failures':0}
(ROOT/'results/relaxation_tests.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
