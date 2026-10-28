"""Exact pin-aware compiler; archived implementation is retained unchanged.

This is constraint propagation and partial evaluation, not a new min-cut theorem.
Active fixed obligations may require positive deformation; their cost must not be
silently deleted. See docs/OWNERWISE_THEORY.md for the preservation argument.
"""
from __future__ import annotations
from pathlib import Path
import sys
from fractions import Fraction as F
from dataclasses import dataclass
from collections import deque
from copy import deepcopy
from itertools import combinations
from time import perf_counter

from layered import (Layers,Tent,AreaTerm,CompiledLayers,_terms,line_at,value,
    dist_intervals,graph_for as original_graph_for,realize,Dinic,check_certificate)


def feasibility(p, pins=None):
    """Exact implication family, with proof causes and propagated fixed bits."""
    m=p.m; edges=set(); fixed={}; causes={}; inconsistent=False
    forced_log=[]; edge_causes={}
    def force(i,v,reason):
        nonlocal inconsistent
        forced_log.append((i,v,reason))
        if i in fixed and fixed[i]!=v: inconsistent=True
        if i not in fixed: causes[i]=reason
        fixed[i]=v
    def edge(u,v,reason):
        if u!=v:
            edges.add((u,v));edge_causes.setdefault((u,v),reason)
    for g in range(p.n):
        for a in range(m-1):edge(p.bit(g,a),p.bit(g,a+1),{'kind':'nesting','group':g})
    # Precompute obligations; do not repeatedly construct identical fractions.
    ob=[[(p.bit(z.group,a),*p.obligation(z,a)) for a in range(m)] for z in p.boxes]
    for zi,z in enumerate(p.boxes):
        for a in range(m):
            vid,zl,zh,L,U=ob[zi][a]
            for b in range(m):
                for d,s in p.kernels[a][b]:
                    h=L+d
                    if h>p.bands[b][1]: force(vid,0,{'kind':'propagated_upper_band','site':zi,'origin':a,'boundary':b,'margin':h-p.bands[b][1]})
                    for t,x in (p.baselines[b][0],p.baselines[b][-1]):
                        margin=h-x-dist_intervals(zl,zh,t,t)*s
                        if margin>0:force(vid,0,{'kind':'propagated_upper_endpoint','site':zi,'origin':a,'boundary':b,'t':t,'margin':margin})
                for d,s in p.kernels[b][a]:
                    h=U-d
                    if h<p.bands[b][0]:force(vid,1,{'kind':'propagated_lower_band','site':zi,'origin':a,'boundary':b,'margin':p.bands[b][0]-h})
                    for t,x in (p.baselines[b][0],p.baselines[b][-1]):
                        margin=x-h-dist_intervals(zl,zh,t,t)*s
                        if margin>0:force(vid,1,{'kind':'propagated_lower_endpoint','site':zi,'origin':a,'boundary':b,'t':t,'margin':margin})
    for i,zi in enumerate(p.boxes):
        for j,zj in enumerate(p.boxes):
            for a in range(m):
                u,il,ih,L,_=ob[i][a]
                for b in range(m):
                    v,jl,jh,_,U=ob[j][b]
                    if u==v:continue
                    dd=dist_intervals(il,ih,jl,jh)
                    for d,s in p.kernels[a][b]:
                        margin=L+d-U-s*dd
                        if margin>0:
                            edge(u,v,{'kind':'ownership_implication','lower_site':i,'upper_site':j,'origin':a,'boundary':b,'margin':margin})
                            break
    for g,owner in (pins or {}).items():
        if not 0<=g<p.n or not 0<=owner<=m:raise ValueError('invalid owner pin')
        for a in range(m):force(p.bit(g,a),int(owner<=a),{'kind':'pin','group':g,'owner':owner})
    base_fixed=dict(fixed)
    out=[[] for _ in range(m*p.n)];rev=[[] for _ in range(m*p.n)]
    for u,v in sorted(edges):out[u].append(v);rev[v].append(u)
    queue=deque(sorted(fixed));processed=set()
    while queue:
        u=queue.popleft()
        if u in processed:continue
        processed.add(u); val=fixed[u]
        for v in (out[u] if val else rev[u]):
            if v in fixed:
                if fixed[v]!=val:inconsistent=True
            else:
                fixed[v]=val;causes[v]={'kind':'implication_propagation','from':u,'value':val,
                                     'edge':[u,v] if val else [v,u]};queue.append(v)
    return {'edges':tuple(sorted(edges)), 'fixed':fixed,'base_fixed':base_fixed,
            'inconsistent':inconsistent,'causes':causes,'edge_causes':edge_causes,'initial_forcing':forced_log}


def mandatory_envelope(base,tents,s,sign):
    """Return shifted reference and exact integrated unavoidable excursion.

    sign=-1: new reference=max(base, mandatory lower cones).
    sign=+1: new reference=min(base, negative mandatory upper cones).
    """
    if not tents:return list(base),F(0)
    pts=sorted({t for t,_ in base}|{t for z in tents for t in (z.lo,z.hi)})
    segments=[];area=F(0)
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; bm,bc=line_at(base,mid)
        # Mandatory nonnegative excursion = max(0, z+sign*base).
        by_slope={F(0):F(0)}
        for z in tents:
            mm,cc=z.line(mid,s);mm+=sign*bm;cc+=sign*bc
            by_slope[mm]=max(by_slope.get(mm,cc),cc)
        lines=list(by_slope.items());cuts={a,b}
        for (mi,ci),(mj,cj) in combinations(lines,2):
            if mi!=mj:
                t=(cj-ci)/(mi-mj)
                if a<t<b:cuts.add(t)
        cuts=sorted(cuts)
        for aa,bb in zip(cuts,cuts[1:]):
            t=(aa+bb)/2;mm,cc=max(lines,key=lambda z:z[0]*t+z[1])
            ea,eb=mm*aa+cc,mm*bb+cc
            assert ea>=0 and eb>=0
            area+=(ea+eb)*(bb-aa)/2
            ma=bm-sign*mm;ca=bc-sign*cc
            ua,ub=ma*aa+ca,ma*bb+ca
            if segments and segments[-1][1]==aa and segments[-1][3]==ua and (segments[-1][3]-segments[-1][2])/(segments[-1][1]-segments[-1][0])==ma:
                segments[-1]=(segments[-1][0],bb,segments[-1][2],ub)
            else:segments.append((aa,bb,ua,ub))
    out=[(segments[0][0],segments[0][2])]+[(b,v) for a,b,u,v in segments]
    return out,area

@dataclass
class PinnedCompiled(CompiledLayers):
    constant_area:F=F(0)
    pins:dict|None=None
    closure_record:dict|None=None
    def area(self,owners):return self.constant_area+super().area(owners)


def compile_pinned(problem:Layers,pins=None,*,reduce_area=True):
    t0=perf_counter();p=deepcopy(problem);m=p.m
    ff=feasibility(p,pins); t1=perf_counter()
    positive={};negative={};slabs=raw=cones=0;inactive=active=0;constant=F(0)
    def merge(dst,src):
        for mask,t in src.items():
            if mask not in dst:dst[mask]=AreaTerm(mask,t.parent,t.added)
            dst[mask].area+=t.area
    if not ff['inconsistent']:
        for k in range(m):
            pos=[];neg=[]
            for z in p.boxes:
                for a in range(m):
                    v=p.bit(z.group,a);zl,zh,L,U=p.obligation(z,a)
                    pos.extend(Tent(zl,zh,L+d,v,s) for d,s in p.kernels[a][k])
                    neg.extend(Tent(zl,zh,-U+d,v,s) for d,s in p.kernels[k][a])
            for tents,sign,dst,on in ((pos,-1,positive,1),(neg,1,negative,0)):
                free=[];mandatory=[]
                for z in tents:
                    if not reduce_area or z.variable not in ff['fixed']:free.append(z)
                    elif ff['fixed'][z.variable]==on:mandatory.append(z);active+=1
                    else:inactive+=1
                base,basearea=mandatory_envelope(p.baselines[k],mandatory,p.slopes[k],sign)
                constant+=basearea
                terms,ss,rr=_terms(free,base,p.slopes[k],sign)
                merge(dst,terms);slabs+=ss;raw+=rr;cones+=len(free)
    stats={'groups':p.n,'sites':len(p.boxes),'source_knots':sum(len(b) for b in p.baselines),
           'boundaries':m,'owner_bits':p.n*m,'implications':len(ff['edges']),
           'forced_bits':len(ff['fixed']),'initial_forced_bits':len(ff['base_fixed']),
           'positive_terms':len(positive),'negative_terms':len(negative),
           'propagated_free_cones':cones,'fixed_active_cones':active,'fixed_inactive_cones':inactive,
           'slabs':slabs,'raw_layers':raw,'constant_area':str(constant),'area_substitution':bool(reduce_area),
           'feasibility_seconds':t1-t0,'area_seconds':perf_counter()-t1,'compile_seconds':perf_counter()-t0}
    return PinnedCompiled(p,ff['edges'],ff['fixed'],positive,negative,stats,ff['inconsistent'],constant,dict(pins or {}),ff)


def graph_for(c,costs,area_weight=1,pins=None):
    data=original_graph_for(c,costs,area_weight,pins)
    data['constant']+=F(area_weight)*c.constant_area
    return data


def solve_pinned(c,costs,*,area_weight=1,pins=None,seconds=60,backend='dinic',certificate=True,reconstruct=True):
    t0=perf_counter();data=graph_for(c,costs,area_weight,pins);S,T=data['s'],data['t']
    if backend=='dinic':
        d=Dinic(data['n'])
        for u,v,w in data['edges']:d.add(u,v,w)
        vv,side,records=d.run(S,T,t0+seconds)
    elif backend=='networkx':
        import networkx as nx
        G=nx.DiGraph();G.add_nodes_from(range(data['n']))
        for u,v,w in data['edges']:
            if G.has_edge(u,v):G[u][v]['capacity']+=w
            else:G.add_edge(u,v,capacity=w)
        R=nx.algorithms.flow.preflow_push(G,S,T);vv=R.graph['flow_value'];side={S};queue=deque([S])
        while queue:
            u=queue.popleft()
            for v,attr in R[u].items():
                if attr['capacity']-attr['flow']>0 and v not in side:side.add(v);queue.append(v)
        records=[(u,v,attr['capacity'],max(0,R[u][v]['flow'])) for u,v,attr in G.edges(data=True)]
    else:raise ValueError('unknown backend')
    cert={'n':data['n'],'s':S,'t':T,'arcs':records,'source_side':sorted(side),'value':int(vv),'scale':data['scale'],'constant':str(data['constant'])}
    check_certificate(cert)
    ans={'feasible':vv<data['big'],'nodes':data['n'],'arcs':len(data['edges']),'seconds':perf_counter()-t0}
    if certificate:ans['certificate']=cert
    if not ans['feasible']:return ans
    bits=[int(i in side) for i in range(c.problem.n*c.problem.m)];owners=c.problem.decode(bits)
    area=c.area(owners);obj=F(area_weight)*area+sum((F(row[q]) for row,q in zip(costs,owners)),F(0))
    assert obj==F(vv,data['scale'])+data['constant']
    ans.update(owners=owners,area=area,objective=obj)
    if reconstruct:
        ww=realize(c.problem,owners)
        assert ww['feasible'] and ww['area']==area
        ans['boundaries']=ww['boundaries']
    ans['seconds']=perf_counter()-t0
    return ans
