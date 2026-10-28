"""Pinned structural regions keep universal semantics in every RNA mode.

This adapter is not a new min-cut theorem. Each obligation is still an interval
lower/upper constraint already covered by the predecessor proof. Geometry pins
must be supplied to solve; the helper does so explicitly.
"""
from pathlib import Path
import sys
from dataclasses import dataclass
from fractions import Fraction as F
from layered import Box, Layers, compile_layers, solve, optimal_owner_domains, value
from verify import verify
@dataclass(frozen=True)
class TaggedBox(Box):
    fixed_region:bool=False
class MixedLayers(Layers):
    def obligation(self,z,k=0):
        if getattr(z,'fixed_region',False):return z.tlo,z.thi,z.xhi+self.clearance,z.xlo-self.clearance
        return super().obligation(z,k)
def verify_mixed(p,owners,bounds,reported_area=None):
    # Geometry and RNA check, then separately check all structural rectangles
    # universally through the independent, original-coordinate checker.
    import copy
    z=copy.copy(p);z.boxes=[b for b in p.boxes if not getattr(b,'fixed_region',False)]
    ans=verify(z,owners,bounds,reported_area)
    z.boxes=[b for b in p.boxes if getattr(b,'fixed_region',False)];z.mode='universal'
    verify(z,owners,bounds,reported_area)
    return ans

def make_problem(record,mode='nominal',radius=F(0)):
    n=len(record['points']);src=[int(z['source']) for z in record['points']];bs=[]
    lo=F(record['baselines'][0][0][0]);hi=F(record['baselines'][0][-1][0])
    for g,z in enumerate(record['points']):
        x,t=F(z['x']),F(z['t']);dy=min(radius,t-lo,hi-t)
        bs.append(TaggedBox(x-radius,x+radius,t-dy,t+dy,g,False))
    pins={}
    for k,rects in enumerate(record.get('nuclear_boxes',[])):
        if not rects:continue
        g=len(src);src.append(k);pins[g]=k
        for x0,x1,t0,t1 in merge_rectangles(rects):bs.append(TaggedBox(F(x0),F(x1),F(t0),F(t1),g,True))
    p=MixedLayers(record['baselines'],bs,src,record['slopes'],record['minimum_gaps'],record['maximum_gaps'],clearance=F(1,1000),bands=record['bands'],mode=mode)
    return p,pins,n


def merge_rectangles(rects):
    byx={}
    for x0,x1,t0,t1 in rects:byx.setdefault((F(x0),F(x1)),[]).append((F(t0),F(t1)))
    out=[]
    for (x0,x1),rr in sorted(byx.items()):
        merged=[]
        for lo,hi in sorted(rr):
            if merged and lo<=merged[-1][1]:merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
            else:merged.append((lo,hi))
        out.extend((x0,x1,lo,hi) for lo,hi in merged)
    return out
