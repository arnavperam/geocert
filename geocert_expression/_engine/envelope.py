"""Exact continuous L1-minimum graph-interface realization.

Classical Lipschitz envelopes are used explicitly, not claimed as new mathematics.
All arithmetic in the geometric oracle is rational. The result is exact ONLY
for the stated single-valued, slope-bounded physical interface family.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as Q
from typing import Sequence

@dataclass(frozen=True)
class Atom:
    x: Q
    y: Q
    owner: int                 # 1: left/A, 0: right/B
    def __post_init__(self):
        object.__setattr__(self,'x',Q(self.x));object.__setattr__(self,'y',Q(self.y))
        if self.owner not in (0,1): raise ValueError('owner must be 0 or 1')

@dataclass
class Realization:
    feasible: bool
    area: Q|None
    segments: list[tuple[Q,Q,Q,Q]] # y0,y1,f(y0),f(y1)
    reason: str=''
    breakpoints: int=0

def _line(a,b):
    x0,y0=a;x1,y1=b
    m=(y1-y0)/(x1-x0)
    return m,y0-m*x0

def _at(line,y): return line[0]*y+line[1]

def evaluate(segments,y):
    y=Q(y)
    for a,b,u,v in segments:
        if a<=y<=b:return u+(v-u)*(y-a)/(b-a)
    raise ValueError('outside interval')

def realize(atoms: Sequence[Atom], *, slope, radius, baseline, clearance=0)->Realization:
    """Minimize integral |f-f0| for |f|<=radius and Lip(f)<=slope.

    Baseline is ordered (y,f0(y)) with slope <= slope; both endpoints are fixed.
    owner A demands f(y_i)>=x_i+clearance, B demands f(y_i)<=x_i-clearance.
    Positive clearance avoids ambiguity of atoms lying exactly on the boundary.
    A closed, zero-clearance relaxation is allowed but labeled as such by caller.
    """
    s=Q(slope);B=Q(radius);rho=Q(clearance)
    if s<=0 or B<=0 or rho<0:raise ValueError('invalid physical parameters')
    base=[(Q(y),Q(f)) for y,f in baseline]
    if len(base)<2 or any(base[j][0]>=base[j+1][0] for j in range(len(base)-1)):
        raise ValueError('baseline knot order')
    if any(abs(f)>B for _,f in base) or any(abs((base[j+1][1]-base[j][1])/(base[j+1][0]-base[j][0]))>s for j in range(len(base)-1)):
        raise ValueError('baseline exceeds family')
    lo,hi=base[0][0],base[-1][0]
    lower={lo:base[0][1],hi:base[-1][1]};upper=dict(lower)
    for z in atoms:
        if not lo<=z.y<=hi:raise ValueError('atom outside longitudinal interval')
        if z.owner:
            lower[z.y]=max(lower.get(z.y,-B),z.x+rho)
        else:upper[z.y]=min(upper.get(z.y,B),z.x-rho)
    if any(v>B for v in lower.values()) or any(v< -B for v in upper.values()):
        return Realization(False,None,[],'outside displacement band')
    knots=sorted(set(lower)|set(upper)|{y for y,_ in base})
    # Prefix/suffix extrema: after the sort, the envelope has <=3 lines per side.
    def scan(values,plus,op,reverse=False):
        out=[None]*len(knots);v=None
        rng=range(len(knots)-1,-1,-1) if reverse else range(len(knots))
        for j in rng:
            t=knots[j]
            if t in values:
                a=values[t]+plus*s*t;v=a if v is None else op(v,a)
            out[j]=v
        return out
    lp=scan(lower,1,max);ls=scan(lower,-1,max,True)
    up=scan(upper,-1,min);us=scan(upper,1,min,True)
    kbase=0;area=Q(0);segs=[]
    for j,(a,b) in enumerate(zip(knots,knots[1:])):
        while kbase+1<len(base)-1 and base[kbase+1][0]<=a:kbase+=1
        bl=_line(base[kbase],base[kbase+1]);low=[(Q(0),-B)];high=[(Q(0),B)]
        if lp[j] is not None:low.append((-s,lp[j]))
        if ls[j+1] is not None:low.append((s,ls[j+1]))
        if up[j] is not None:high.append((s,up[j]))
        if us[j+1] is not None:high.append((-s,us[j+1]))
        lines=low+high+[bl];pts={a,b}
        for u in range(len(lines)):
            for v in range(u):
                if lines[u][0]==lines[v][0]:continue
                z=(lines[v][1]-lines[u][1])/(lines[u][0]-lines[v][0])
                if a<z<b:pts.add(z)
        pts=sorted(pts)
        for c,d in zip(pts,pts[1:]):
            mid=(c+d)/2;ll=max(low,key=lambda z:_at(z,mid));uu=min(high,key=lambda z:_at(z,mid))
            if any(_at(ll,y)>_at(uu,y) for y in (c,mid,d)):
                return Realization(False,None,[],'incompatible Lipschitz ownership constraints')
            f=ll if _at(bl,mid)<_at(ll,mid) else (uu if _at(bl,mid)>_at(uu,mid) else bl)
            fc,fd=_at(f,c),_at(f,d)
            area+=(abs(fc-_at(bl,c))+abs(fd-_at(bl,d)))*(d-c)/2
            segs.append((c,d,fc,fd))
    return Realization(True,area,segs,breakpoints=len(segs)+1)

def verify(result,atoms,*,slope,radius,baseline,clearance=0):
    """Direct rational witness checks, independent of envelope construction."""
    if not result.feasible:return False
    s=Q(slope);B=Q(radius);rho=Q(clearance);base=[(Q(y),Q(f)) for y,f in baseline]
    assert result.segments and result.segments[0][0]==base[0][0] and result.segments[-1][1]==base[-1][0]
    area=Q(0)
    for j,(a,b,u,v) in enumerate(result.segments):
        assert a<b and abs(v-u)<=s*(b-a) and max(abs(u),abs(v))<=B
        if j:assert result.segments[j-1][1:2]==(a,) and result.segments[j-1][3]==u
        f0=next(_line(base[t],base[t+1]) for t in range(len(base)-1) if base[t][0]<=a and b<=base[t+1][0])
        aa=u-_at(f0,a);bb=v-_at(f0,b)
        # Oracle splits at every intersection with f0. Validate this property.
        assert aa*bb>=0
        area+=(abs(aa)+abs(bb))*(b-a)/2
    assert area==result.area
    assert result.segments[0][2]==base[0][1] and result.segments[-1][3]==base[-1][1]
    for z in atoms:
        val=evaluate(result.segments,z.y)
        assert val>=z.x+rho if z.owner else val<=z.x-rho
    return True
