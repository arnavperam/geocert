"""Target-restricted exact spectra via the published source-clipping routine."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from dataclasses import dataclass
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'core'))
from pin_closure import PinClosure,verify_closure_answer

@dataclass
class Spectrum:
    minima:dict
    witnesses:dict
    complete:bool
    examined:int
    total_patterns:int

class LocalEngine:
    def __init__(self,p,pins):
        self.problem=p;self.source_engine=PinClosure(p,pins);self.cache={}
    def clipped(self,pins):
        key=tuple(sorted(pins.items()))
        if key not in self.cache:
            a=self.source_engine.query(pins)
            verify_closure_answer(self.source_engine,pins,a)
            self.cache[key]=a
        return self.cache[key]
    def spectrum(self,scores,constant=0,*,pattern_limit=65536):
        p=self.problem
        if len(scores)!=p.n or any(len(r)!=p.m+1 for r in scores):raise ValueError('score shape')
        if any(type(v)!=int for r in scores for v in r):raise TypeError('integer score required')
        relevant=[i for i,r in enumerate(scores) if len(set(r))>1]
        base=constant+sum(r[0] for r in scores if len(set(r))==1)
        domains=[[o for o in range(p.m+1) if self.source_engine.closure({i:o})['feasible']] for i in relevant]
        total=1
        for d in domains:total*=len(d)
        costs={};wits={};examined=0
        for owners in product(*domains):
            if examined>=pattern_limit:break
            examined+=1;pins=dict(zip(relevant,owners));a=self.clipped(pins)
            if not a['feasible']:continue
            count=base+sum(scores[i][o] for i,o in pins.items());area=F(a['area'])
            if count not in costs or area<costs[count]:
                costs[count]=area;wits[count]={'pins':pins,'owners':a['owners'],'area':str(area),'boundaries':a['boundaries']}
        return Spectrum(costs,wits,examined==total,examined,total)
