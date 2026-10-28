from common_test import *
from hard_scc import *
from scc_solver import prepare,solve_compiled
from itertools import product
import random,importlib.util
rng=random.Random(720261027);out={'graphs':0,'cuts_checked':0,'models':0,'optional_objectives':0,'failures':0}
for _ in range(200):
    n=rng.randint(3,9);hard=[(u,v) for u in range(n) for v in range(n) if u!=v and rng.random()<.12]
    finite=[(u,v,rng.randint(1,7)) for u in range(n) for v in range(n) if u!=v and rng.random()<.10]
    big=sum(w for u,v,w in finite)+1
    g={'n':n,'s':0,'t':n-1,'edges':finite+[(u,v,big) for u,v in hard],'big':big}
    plan=build_plan(g);q=contract(g,plan)
    old_vals=[]
    for bits in product([0,1],repeat=n-2):
        side={0}|{i+1 for i,b in enumerate(bits) if b}
        if any(u in side and v not in side for u,v in hard):continue
        old_vals.append(sum(w for u,v,w in g['edges'] if u in side and v not in side))
        assert all((i in side)==(j in side) for i in range(n) for j in range(n) if plan.mapping[i]==plan.mapping[j])
        reduced={plan.mapping[i] for i in side}
        assert sum(w for u,v,w in q['edges'] if u in reduced and v not in reduced)==old_vals[-1]
        out['cuts_checked']+=1
    new_vals=[]
    other=[i for i in range(q['n']) if i not in (q['s'],q['t'])]
    for bits in product([0,1],repeat=len(other)):
        ss={q['s']}|{i for i,b in zip(other,bits) if b};side=lift_side(ss,plan)
        if any(u in side and v not in side for u,v in hard):continue
        new_vals.append(sum(w for u,v,w in q['edges'] if u in ss and v not in ss))
    assert sorted(old_vals)==sorted(new_vals)
    out['graphs']+=1
import model_generator as mod
for seed in range(24):
    p,rng=mod.make(270100+seed,n=1+seed%3,m=1+seed%2);p.mode=['nominal','universal','possible_outer'][seed%3]
    pins={0:p.source[0]} if seed%7==0 else {};compiled={}
    for r in ['direct','height']:
        c,plan,_=prepare(p,pins,r,True);compiled[r]=(c,plan)
    for j,lam in enumerate([F(0),F(1,4),F(3)]):
        costs=[[F(rng.randint(-5,5),3) for _ in range(p.m+1)] for _ in range(p.n)]
        answers=[]
        for r,(c,plan) in compiled.items():
            for use in [False,True]:
                a=solve_compiled(c,r,costs,lam,plan if use else None)
                answers.append((a['feasible'],a.get('objective')))
        assert len(set(answers))==1,(seed,j,answers)
        out['optional_objectives']+=1
    out['models']+=1
save7('scc_validation.json',out);print(out)
