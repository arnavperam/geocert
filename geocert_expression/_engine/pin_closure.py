"""Exact source-clipped pinned deformation queries; no integrated-cost compiler.

The full geometry implementation is reused for feasibility and fixed-owner
reconstruction, not for analytic activation layers or weighted cut optimization.
The theorem requires a zero-deformation source satisfying ALL ownership obligations.
"""
from __future__ import annotations
from collections import deque
from dataclasses import dataclass
from fractions import Fraction as F
from time import perf_counter
from pinned_compiler import feasibility
from layered import realize
from mixed import verify_mixed

class SourceNotFeasible(ValueError):
    """Source is not a valid zero-deformation ownership realization."""

@dataclass
class PinClosure:
    problem: object
    structural_pins: dict
    def __post_init__(self):
        t0=perf_counter();p=self.problem
        self.family=feasibility(p,self.structural_pins)
        self.inconsistent=self.family['inconsistent']
        self.out=[[] for _ in range(p.n*p.m)]
        self.rev=[[] for _ in range(p.n*p.m)]
        for a,b in self.family['edges']:
            self.out[a].append(b);self.rev[b].append(a)
        self.bits=p.source_bits()
        if not self.inconsistent:
            if any(self.bits[v]!=z for v,z in self.family['fixed'].items()) or any(self.bits[a]>self.bits[b] for a,b in self.family['edges']):
                raise SourceNotFeasible('source ownership violates geometric implications or fixed pins')
            source_segments=[[(a,b,u,v) for (a,u),(b,v) in zip(base,base[1:])] for base in p.baselines]
            try:verify_mixed(p,p.source,source_segments,F(0))
            except (AssertionError, ValueError) as e:
                raise SourceNotFeasible('source does not satisfy every retained ownership/region constraint at zero deformation') from e
        self.setup_seconds=perf_counter()-t0

    def closure(self,pins):
        p=self.problem
        if self.inconsistent:return {'status':'model_infeasible','feasible':False}
        initial_one={v for v,z in self.family['fixed'].items() if z}
        initial_zero={v for v,z in self.family['fixed'].items() if not z}
        for g,o in pins.items():
            if type(g) is not int or type(o) is not int or not 0<=g<p.n or not 0<=o<=p.m:
                raise ValueError('invalid group/owner pin')
            for a in range(p.m):
                (initial_one if o<=a else initial_zero).add(p.bit(g,a))
        # One forward proof path into a forced-zero bit suffices for a conflict.
        pred={v:None for v in initial_one};q=deque(sorted(initial_one))
        while q:
            u=q.popleft()
            if u in initial_zero:
                chain=[]
                while u is not None:chain.append(u);u=pred[u]
                chain.reverse()
                return {'status':'pin_infeasible','feasible':False,'contradiction':chain}
            for v in self.out[u]:
                if v not in pred:pred[v]=u;q.append(v)
        T=set(pred);Fset=set(initial_zero);q=deque(sorted(initial_zero))
        while q:
            u=q.popleft()
            for v in self.rev[u]:
                if v not in Fset:Fset.add(v);q.append(v)
        assert T.isdisjoint(Fset)
        bits=list(self.bits)
        for i in T:bits[i]=1
        for i in Fset:bits[i]=0
        owners=p.decode(bits)
        assert all(bits[a]<=bits[b] for a,b in self.family['edges'])
        assert all(owners[g]==o for g,o in pins.items())
        return {'status':'optimal','feasible':True,'owners':owners,
                'forced_ones':sorted(T),'forced_zeros':sorted(Fset)}

    def query(self,pins,*,check=True):
        t0=perf_counter();a=self.closure(pins);t1=perf_counter()
        a['closure_seconds']=t1-t0
        if not a['feasible']:
            a.update(reconstruct_seconds=0.,check_seconds=0.,seconds=perf_counter()-t0)
            return a
        w=realize(self.problem,a['owners']);t2=perf_counter()
        if not w['feasible']:raise AssertionError('closure/reconstruction contradiction')
        a.update(area=w['area'],objective=w['area'],boundaries=w['boundaries'])
        if check:verify_mixed(self.problem,a['owners'],w['boundaries'],w['area'])
        a.update(reconstruct_seconds=t2-t1,check_seconds=perf_counter()-t2,seconds=perf_counter()-t0)
        return a


def verify_closure_answer(engine,pins,answer):
    """Set-based independent implication replay, not the query's BFS routine."""
    fam=engine.family;p=engine.problem
    one={i for i,v in fam['fixed'].items() if v}
    zero={i for i,v in fam['fixed'].items() if not v}
    for g,o in pins.items():
        for a in range(p.m):(one if o<=a else zero).add(p.bit(g,a))
    changed=True
    while changed:
        old=(len(one),len(zero))
        for u,v in fam['edges']:
            if u in one:one.add(v)
            if v in zero:zero.add(u)
        changed=old!=(len(one),len(zero))
    feasible=not fam['inconsistent'] and not (one&zero)
    if feasible!=answer['feasible']:raise AssertionError('closure feasibility mismatch')
    if feasible:
        expected=p.source_bits()
        for i in one:expected[i]=1
        for i in zero:expected[i]=0
        if p.encode(answer['owners'])!=expected:raise AssertionError('wrong source-clipped owner vector')
        verify_mixed(p,answer['owners'],answer['boundaries'],answer['area'])
    elif answer.get('status')=='pin_infeasible':
        chain=answer['contradiction'];edges=set(fam['edges'])
        initial_one={i for i,v in fam['fixed'].items() if v};initial_zero={i for i,v in fam['fixed'].items() if not v}
        for g,o in pins.items():
            for a in range(p.m):(initial_one if o<=a else initial_zero).add(p.bit(g,a))
        if not chain or chain[0] not in initial_one or chain[-1] not in initial_zero or any((u,v) not in edges for u,v in zip(chain,chain[1:])):
            raise AssertionError('invalid contradiction path')
    return True
