"""Exact rational finite-cut representation of continuous Lipschitz ownership.

The compiler, not maximum-flow, is the project-specific derivation.
Scores are separable unary terms. Hard budgets/nonlinear rankings are NOT
solved by the polynomial-time core. Input units and slope must be physical.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from bisect import bisect_right
from collections import deque
from math import lcm
from time import perf_counter
import sys
from pathlib import Path
from envelope import Atom,realize,verify as verify_boundary

@dataclass(frozen=True)
class Site:
    x:F
    y:F
    group:int
    def __post_init__(self):
        object.__setattr__(self,'x',F(self.x));object.__setattr__(self,'y',F(self.y))

@dataclass
class Interface:
    sites:list[Site]
    source:list[int]
    baseline:list[tuple[F,F]]
    slope:F
    radius:F
    clearance:F=F(1,100)
    def __post_init__(self):
        self.baseline=[(F(y),F(x)) for y,x in self.baseline]
        self.slope=F(self.slope);self.radius=F(self.radius);self.clearance=F(self.clearance)
        if any(x not in (0,1) for x in self.source):raise ValueError('binary source required')
        if self.clearance<=0:raise ValueError('strictly positive ownership clearance required')
        if any(type(z.group) is not int or not 0<=z.group<len(self.source) for z in self.sites):raise ValueError('invalid group')
        if set(z.group for z in self.sites)!=set(range(len(self.source))):raise ValueError('every group must contain at least one RNA site')
        realize([],**self.params())
        for z in self.sites:
            if not self.baseline[0][0]<=z.y<=self.baseline[-1][0]:raise ValueError('site outside interval')
            b=self.at(z.y)
            if (self.source[z.group] and b<z.x+self.clearance) or (not self.source[z.group] and b>z.x-self.clearance):
                raise ValueError('source fails declared ownership clearance')
    def at(self,y):
        j=max(0,min(len(self.baseline)-2,bisect_right([a for a,_ in self.baseline],y)-1))
        a,u=self.baseline[j];b,v=self.baseline[j+1]
        return u+(v-u)*(y-a)/(b-a)
    def params(self):return dict(baseline=self.baseline,slope=self.slope,radius=self.radius,clearance=self.clearance)
    def atoms(self,owners):return [Atom(z.x,z.y,owners[z.group]) for z in self.sites]
    def weights(self):return [sum(z.group==g for z in self.sites) for g in range(len(self.source))]

@dataclass
class AreaTerm:
    mask:int
    parent:int
    added:int
    area:F=F(0)

@dataclass
class Compiled:
    problem:Interface
    implies:tuple[tuple[int,int],...]
    forced:dict[int,int]
    positive:dict[int,AreaTerm]
    negative:dict[int,AreaTerm]
    stats:dict
    def feasible(self,x):
        return all(x[a]<=x[b] for a,b in self.implies) and all(x[g]==v for g,v in self.forced.items())
    def area(self,x):
        yes=sum((1<<i) for i,v in enumerate(x) if v)
        allmask=(1<<len(x))-1
        return sum((z.area for z in self.positive.values() if z.mask&yes),F(0))+sum((z.area for z in self.negative.values() if z.mask&(allmask^yes)),F(0))


def _base_line(p,y):
    bs=p.baseline;j=max(0,min(len(bs)-2,bisect_right([a for a,_ in bs],y)-1))
    a,u=bs[j];b,v=bs[j+1];m=(v-u)/(b-a)
    return m,u-m*a

def _affine_value(line,t):return line[0]*t+line[1]

def _area_terms(p,source):
    """Integrate max cones using exact arrangement slabs; merge equal coverage sets.

    source=0: positive (B->A) excursion; source=1: negative (A->B).
    Each OR-set is constructed from an earlier proper subset and one input.
    This is classical coverage-to-cut reduction applied to exact analytic layers.
    """
    sites=[z for z in p.sites if p.source[z.group]==source]
    if not sites:return {},0,0
    lo,hi=p.baseline[0][0],p.baseline[-1][0];s=p.slope;rho=p.clearance
    # h_i is tent height before subtracting/adding the baseline.
    heights=[(z.x if source==0 else -z.x)+rho for z in sites]
    points={y for y,_ in p.baseline}|{z.y for z in sites}
    basepoints=sorted(points)
    def lines_at(mid):
        bm,bc=_base_line(p,mid);sign=-1 if source==0 else 1
        return [((-s if z.y<mid else s)+sign*bm,
                 h+(s*z.y if z.y<mid else -s*z.y)+sign*bc)
                for z,h in zip(sites,heights)]
    # Pair order changes at at most one ordinate between their centers.
    for i in range(len(sites)):
        for j in range(i):
            a,b=(j,i) if sites[j].y<=sites[i].y else (i,j)
            y=(heights[a]-heights[b]+s*(sites[a].y+sites[b].y))/(2*s)
            if sites[a].y<y<sites[b].y and lo<y<hi:
                height=heights[a]-s*abs(y-sites[a].y)+(-p.at(y) if source==0 else p.at(y))
                if height>0:points.add(y)
    # Add exact zeros where a cone crosses the source boundary.
    for a,b in zip(basepoints,basepoints[1:]):
        for m,c in lines_at((a+b)/2):
            if m:
                z=-c/m
                if a<z<b:points.add(z)
    points=sorted(points);terms={};raw_layers=0
    for a,b in zip(points,points[1:]):
        mid=(a+b)/2;lines=lines_at(mid)
        order=sorted((i for i,l in enumerate(lines) if _affine_value(l,mid)>0),key=lambda i:(-_affine_value(lines[i],mid),i))
        mask=0
        for k,i in enumerate(order):
            nxt=mask|(1<<sites[i].group)
            if nxt not in terms:terms[nxt]=AreaTerm(nxt,mask,sites[i].group)
            mask=nxt
            low=lines[order[k+1]] if k+1<len(order) else (F(0),F(0))
            top=lines[i]
            w=((_affine_value(top,a)-_affine_value(low,a))+(_affine_value(top,b)-_affine_value(low,b)))*(b-a)/2
            if w<0:raise ArithmeticError('incorrect analytic cone arrangement')
            terms[mask].area+=w;raw_layers+=1
    return terms,len(points)-1,raw_layers


def compile_interface(p:Interface)->Compiled:
    t0=perf_counter();edges=set();fixed={};s=p.slope;rho=p.clearance
    for i,z in enumerate(p.sites):
        lower=z.x+rho;upper=z.x-rho
        canA=lower<=p.radius and all(lower-f<=s*abs(z.y-y) for y,f in (p.baseline[0],p.baseline[-1]))
        canB=upper>=-p.radius and all(f-upper<=s*abs(z.y-y) for y,f in (p.baseline[0],p.baseline[-1]))
        for allowed,label in ((canA,0),(canB,1)):
            if not allowed:
                if z.group in fixed and fixed[z.group]!=label:raise ValueError('inconsistent source')
                fixed[z.group]=label
        for u in p.sites:
            if z.group!=u.group and lower-(u.x-rho)>s*abs(z.y-u.y):edges.add((z.group,u.group))
    pos,sp,lp=_area_terms(p,0);neg,sn,ln=_area_terms(p,1)
    out=Compiled(p,tuple(sorted(edges)),fixed,pos,neg,dict(compile_seconds=perf_counter()-t0,sites=len(p.sites),groups=len(p.source),implications=len(edges),positive_terms=len(pos),negative_terms=len(neg),arrangement_slabs=sp+sn,raw_layers=lp+ln))
    assert out.feasible(p.source) and out.area(p.source)==0
    return out

class Dinic:
    """Ordinary exact integer max flow; not claimed as new."""
    def __init__(self,n):self.g=[[] for _ in range(n)];self.records=[]
    def add(self,u,v,c):
        if c<0:raise ValueError('negative capacity')
        if u==v:return
        e=[v,int(c),len(self.g[v])];r=[u,0,len(self.g[u])]
        self.g[u].append(e);self.g[v].append(r);self.records.append((u,e,c))
    def run(self,s,t,deadline):
        n=len(self.g);total=0
        sys.setrecursionlimit(max(sys.getrecursionlimit(),2*n+100))
        while True:
            if perf_counter()>deadline:raise TimeoutError('max-flow resource limit: no optimum returned')
            level=[-1]*n;level[s]=0;q=deque([s])
            while q:
                u=q.popleft()
                for v,c,_ in self.g[u]:
                    if c and level[v]<0:level[v]=level[u]+1;q.append(v)
            if level[t]<0:break
            it=[0]*n
            def push(u,f):
                if u==t:return f
                while it[u]<len(self.g[u]):
                    e=self.g[u][it[u]];v,c,rev=e
                    if c and level[v]==level[u]+1:
                        h=push(v,min(f,c))
                        if h:e[1]-=h;self.g[v][rev][1]+=h;return h
                    it[u]+=1
                return 0
            inf=sum(e[1] for e in self.g[s])+1
            while True:
                h=push(s,inf)
                if not h:break
                total+=h
                if perf_counter()>deadline:raise TimeoutError('max-flow resource limit: no optimum returned')
        seen={s};q=deque([s])
        while q:
            u=q.popleft()
            for v,c,_ in self.g[u]:
                if c and v not in seen:seen.add(v);q.append(v)
        records=[(u,e[0],int(c),int(c-e[1])) for u,e,c in self.records]
        return total,seen,records


def graph_for(c:Compiled,cost0,cost1,area_weight=1,pins=None):
    n=len(c.problem.source);lam=F(area_weight)
    if lam<0:raise ValueError('area multiplier must be nonnegative')
    if len(cost0)!=n or len(cost1)!=n:raise ValueError('one unary pair per group')
    S=n;T=n+1;nextid=n+2;finite=[];hard=[];constant=F(0)
    for i,(u,v) in enumerate(zip(cost0,cost1)):
        u,v=F(u),F(v);m=min(u,v);constant+=m
        if u>m:finite.append((S,i,u-m))
        if v>m:finite.append((i,T,v-m))
    for direction,terms in ((0,c.positive),(1,c.negative)):
        ids={}
        for mask in sorted(terms,key=lambda x:(x.bit_count(),x)):
            z=terms[mask];v=nextid;nextid+=1;ids[mask]=v
            if lam*z.area:finite.append((v,T,lam*z.area) if direction==0 else (S,v,lam*z.area))
            deps=[z.added]+([ids[z.parent]] if z.parent else [])
            hard.extend((d,v) if direction==0 else (v,d) for d in deps)
    hard.extend(c.implies)
    for i,v in list(c.forced.items())+list((pins or {}).items()):
        if not 0<=i<n or v not in (0,1):raise ValueError('invalid pin')
        hard.append((S,i) if v else (i,T))
    scale=1
    for _,_,v in finite:scale=lcm(scale,v.denominator)
    fin=[(u,v,int(w*scale)) for u,v,w in finite]
    big=sum(w for _,_,w in fin)+1
    all_edges=fin+[(u,v,big) for u,v in hard if u!=v]
    return dict(n=nextid,s=S,t=T,edges=all_edges,scale=scale,constant=constant,big=big,hard=hard)


def solve(c:Compiled,cost0,cost1,*,area_weight=1,pins=None,seconds=60,certificate=True,fast_verify=False):
    t0=perf_counter();data=graph_for(c,cost0,cost1,area_weight,pins);d=Dinic(data['n'])
    for u,v,w in data['edges']:d.add(u,v,w)
    value,side,records=d.run(data['s'],data['t'],t0+seconds)
    feasible=value<data['big']
    out=dict(feasible=feasible,seconds=perf_counter()-t0,nodes=data['n'],arcs=len(data['edges']))
    if certificate:
        out['certificate']=dict(n=data['n'],s=data['s'],t=data['t'],arcs=records,source_side=sorted(side),value=value,scale=data['scale'],constant=str(data['constant']))
    if not feasible:return out
    x=[int(i in side) for i in range(len(c.problem.source))]
    area=c.area(x);obj=F(area_weight)*area+sum((F(cost1[i]) if x[i] else F(cost0[i]) for i in range(len(x))),F(0))
    assert obj==F(value,data['scale'])+data['constant']
    assert c.feasible(x) and all(x[i]==v for i,v in (pins or {}).items())
    geo=realize(c.problem.atoms(x),**c.problem.params())
    assert geo.feasible and geo.area==area
    if fast_verify:
        from check_boundary import verify_boundary_indexed
        verify_boundary_indexed(geo,c.problem.atoms(x),**c.problem.params())
    else:
        verify_boundary(geo,c.problem.atoms(x),**c.problem.params())
    out.update(owners=x,area=area,objective=obj,segments=geo.segments)
    return out
