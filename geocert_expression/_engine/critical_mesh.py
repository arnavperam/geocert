"""Exact critical-mesh/ordered-height reference construction.

This is an implemented alternative derived from the reviewer's canonical-mesh
argument, NOT a reimplementation or benchmark of a published imaging package.
It shares primitive model semantics, propagated kernels, feasibility/pin closure,
and the exact flow backend with the direct representation; it does not call the
analytic ownership-cost compiler to form its graph. Every node cost is the
trapezoidal absolute-displacement weight on the ownership-independent mesh.

Safe pruning removes crossings and candidate values outside the envelopes forced
by geometry and structural pins. This gives the reference the same pin-closure
information as the direct compiler. No molecular cost influences compilation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from bisect import bisect_left,bisect_right
from time import perf_counter
from math import lcm
from collections import deque
from pinned_compiler import feasibility
from continuum_cut import Dinic
from check_flow import check_certificate
from mixed import verify_mixed
from layered import Tent, line_at, value, realize

@dataclass
class MeshCompiled:
    problem: object
    mesh: list
    candidates: list
    ids: list
    area_coefs: dict
    area_constant: F
    hard: tuple
    nvertices: int
    source: int
    sink: int
    pins: dict
    stats: dict
    inconsistent: bool=False
    fixed_answer: dict|None=None


def compile_mesh(p,pins=None,*,seconds=60,node_cap=300000,prune=True):
    start=perf_counter();deadline=start+seconds;pins=dict(pins or {})
    def guard():
        if perf_counter()>deadline:raise TimeoutError('critical-mesh compilation limit; no exact result returned')
    ff=feasibility(p,pins);guard();S=p.m*p.n;T=S+1
    if ff['inconsistent']:
        return MeshCompiled(p,[],[],[],{},F(0),((S,T),),S+2,S,T,pins,{'inconsistent':True},True)
    # Do not handicap the reference on models with no optional ownership left.
    # One fixed-owner envelope reconstruction replaces any height expansion.
    if len(ff['fixed'])==p.m*p.n:
        owners=p.decode([ff['fixed'][i] for i in range(p.m*p.n)])
        w=realize(p,owners)
        if not w['feasible']:raise AssertionError('fixed closure/geometry mismatch')
        w['owners']=owners
        hard=set(ff['edges'])
        for i,v in ff['fixed'].items():hard.add((S,i) if v else (i,T))
        stats={'fixed_owner_fast_path':True,'mesh_nodes':0,'height_bits':0,'owner_bits':p.n*p.m,
               'graph_nodes':S+2,'hard_arcs':len(hard),'compile_seconds':perf_counter()-start}
        return MeshCompiled(p,[],[],[],{},w['area'],tuple(sorted(hard)),S+2,S,T,pins,stats,False,w)
    # Every possible canonical piece, plus forced lower/upper pieces.
    # Tuples encode (cone, output sign); constants have cone=None.
    fun=[];fixed_lo=[];fixed_hi=[];initial={p.lo,p.hi}
    for b in p.baselines:initial.update(t for t,_ in b)
    for k in range(p.m):
        fs=[];ll=[];uu=[]
        def add(z,sgn,fixed=None):
            q=(z,sgn);fs.append(q)
            if fixed=='low':ll.append(q)
            if fixed=='high':uu.append(q)
            if isinstance(z,Tent):initial.update((z.lo,z.hi))
        for a in range(p.m):
            if p.closure[a][k] is not None:add(p.bands[a][0]+p.closure[a][k],1,'low')
            if p.closure[k][a] is not None:add(p.bands[a][1]-p.closure[k][a],1,'high')
            for d,s in p.kernels[a][k]:
                for t,x in (p.baselines[a][0],p.baselines[a][-1]):add(Tent(t,t,x+d,-1,s),1,'low')
            for d,s in p.kernels[k][a]:
                for t,x in (p.baselines[a][0],p.baselines[a][-1]):add(Tent(t,t,-x+d,-1,s),-1,'high')
        for z in p.boxes:
            for a in range(p.m):
                v=p.bit(z.group,a);lo,hi,L,U=p.obligation(z,a);f=ff['fixed'].get(v)
                if f!=0:
                    for d,s in p.kernels[a][k]:add(Tent(lo,hi,L+d,-1,s),1,'low' if f==1 else None)
                if f!=1:
                    for d,s in p.kernels[k][a]:add(Tent(lo,hi,-U+d,-1,s),-1,'high' if f==0 else None)
        fun.append(list(dict.fromkeys(fs)));fixed_lo.append(list(dict.fromkeys(ll)));fixed_hi.append(list(dict.fromkeys(uu)))
    initial=sorted(initial)
    def eval_fun(f,t,k):
        z,sgn=f
        return sgn*(z.at(t,p.slopes[k]) if isinstance(z,Tent) else z)
    def linear(f,t,k):
        z,sgn=f
        if isinstance(z,Tent):a,b=z.line(t,p.slopes[k]);return sgn*a,sgn*b
        return F(0),sgn*z
    critical=set(initial);pair_checks=0;retained=0
    for k in range(p.m):
        for l,r in zip(initial,initial[1:]):
            guard();mid=(l+r)/2
            ls=set(linear(f,mid,k) for f in fun[k]);ls.add(line_at(p.baselines[k],mid))
            low_lines=[linear(f,mid,k) for f in fixed_lo[k]];high_lines=[linear(f,mid,k) for f in fixed_hi[k]]
            for (a,b),(c,d) in combinations(ls,2):
                pair_checks+=1
                if a==c:continue
                t=(d-b)/(a-c)
                if not l<t<r:continue
                h=a*t+b
                if prune and (h<max(u*t+v for u,v in low_lines) or h>min(u*t+v for u,v in high_lines)):continue
                critical.add(t);retained+=1
            guard()
    mesh=sorted(critical);mesh_time=perf_counter()-start
    vals=[];ids=[];nextid=T+1;area_coefs={};constant=F(0);hard=set(ff['edges'])
    for v,z in ff['fixed'].items():hard.add((S,v) if z else (v,T))
    for k in range(p.m):
        row=[];ii=[]
        for j,t in enumerate(mesh):
            if j%32==0:guard()
            lo=max(eval_fun(f,t,k) for f in fixed_lo[k]);hi=min(eval_fun(f,t,k) for f in fixed_hi[k])
            choices=sorted({eval_fun(f,t,k) for f in fun[k]}|{value(p.baselines[k],t)})
            choices=[x for x in choices if lo<=x<=hi]
            if not choices:raise AssertionError('source-independent canonical candidate set empty')
            row.append(choices);xx=[S]
            for h in range(1,len(choices)):
                xx.append(nextid);nextid+=1
                if nextid>node_cap:raise MemoryError('critical-mesh node cap; no exact result returned')
                if h>1:hard.add((xx[h],xx[h-1]))
            ii.append(xx)
            weight=((t-mesh[j-1]) if j else F(0))+((mesh[j+1]-t) if j+1<len(mesh) else F(0));weight/=2
            y=value(p.baselines[k],t);prev=abs(choices[0]-y)*weight;constant+=prev
            for h in range(1,len(choices)):
                cost=abs(choices[h]-y)*weight;coef=cost-prev
                if coef:area_coefs[xx[h]]=coef
                prev=cost
        vals.append(row);ids.append(ii)
    def lower_token(k,j,h):
        q=bisect_left(vals[k][j],h)
        return T if q==len(vals[k][j]) else ids[k][j][q]
    def strictly_above(k,j,h):
        q=bisect_right(vals[k][j],h)
        return T if q==len(vals[k][j]) else ids[k][j][q]
    def add_arc(u,v):
        if u!=T and v!=S and u!=v:hard.add((u,v))
    def difference(a,i,b,j,d):
        # H_b(j) >= H_a(i)+d, including the unconditional minimum level.
        for v,h in zip(ids[a][i],vals[a][i]):add_arc(v,lower_token(b,j,h+d))
    for k in range(p.m):
        for j in range(len(mesh)-1):
            dd=-p.slopes[k]*(mesh[j+1]-mesh[j]);difference(k,j,k,j+1,dd);difference(k,j+1,k,j,dd)
    for k in range(p.m-1):
        for j in range(len(mesh)):
            difference(k,j,k+1,j,p.minimum_gaps[k])
            if p.maximum_gaps[k] is not None:difference(k+1,j,k,j,-p.maximum_gaps[k])
    for z in p.boxes:
        for k in range(p.m):
            lo,hi,L,U=p.obligation(z,k);bit=p.bit(z.group,k)
            for j in range(bisect_left(mesh,lo),bisect_right(mesh,hi)):
                add_arc(bit,lower_token(k,j,L));add_arc(strictly_above(k,j,U),bit)
    guard()
    stats={'mesh_nodes':len(mesh),'initial_mesh_nodes':len(initial),'height_bits':nextid-(T+1),
           'owner_bits':p.n*p.m,'graph_nodes':nextid,'hard_arcs':len(hard),
           'area_unaries':len(area_coefs),'pair_checks':pair_checks,'retained_crossing_occurrences':retained,
           'total_candidates':sum(len(z) for row in vals for z in row),
           'fixed_range_pruning':prune,'mesh_seconds':mesh_time,'compile_seconds':perf_counter()-start}
    return MeshCompiled(p,mesh,vals,ids,area_coefs,constant,tuple(sorted(hard)),nextid,S,T,pins,stats)


def graph_mesh(c,costs,area_weight=1,pins=None):
    p=c.problem;lam=F(area_weight)
    if lam<0:raise ValueError('nonnegative area penalty required')
    if len(costs)!=p.n or any(len(r)!=p.m+1 for r in costs):raise ValueError('bad unary costs')
    S,T=c.source,c.sink;finite=[];hard=set(c.hard);constant=lam*c.area_constant
    def unary(v,a):
        nonlocal constant
        if a>=0:
            if a:finite.append((v,T,a))
        else:constant+=a;finite.append((S,v,-a))
    for v,a in c.area_coefs.items():unary(v,lam*a)
    for g,r in enumerate(costs):
        r=list(map(F,r));constant+=r[-1]
        for a in range(p.m):unary(p.bit(g,a),r[a]-r[a+1])
    for g,o in (pins or {}).items():
        if not 0<=g<p.n or not 0<=o<=p.m:raise ValueError('bad owner pin')
        for a in range(p.m):hard.add((S,p.bit(g,a)) if o<=a else (p.bit(g,a),T))
    scale=1
    for u,v,a in finite:scale=lcm(scale,a.denominator)
    edges=[(u,v,int(a*scale)) for u,v,a in finite];big=sum(a for u,v,a in edges)+1
    edges += [(u,v,big) for u,v in sorted(hard) if u!=v]
    return dict(n=c.nvertices,s=S,t=T,edges=edges,scale=scale,constant=constant,big=big)


def solve_mesh(c,costs,*,area_weight=1,pins=None,seconds=60,certificate=False,backend='dinic'):
    start=perf_counter();g=graph_mesh(c,costs,area_weight,pins);t1=perf_counter()
    if backend=='dinic':
        f=Dinic(g['n'])
        for u,v,a in g['edges']:f.add(u,v,a)
        opt,side,records=f.run(g['s'],g['t'],start+seconds)
    elif backend=='networkx':
        import networkx as nx
        G=nx.DiGraph();G.add_nodes_from(range(g['n']))
        for u,v,a in g['edges']:
            if G.has_edge(u,v):G[u][v]['capacity']+=a
            else:G.add_edge(u,v,capacity=a)
        residual=nx.algorithms.flow.preflow_push(G,g['s'],g['t']);opt=residual.graph['flow_value'];side={g['s']};q=deque(side)
        while q:
            u=q.popleft()
            for v,d in residual[u].items():
                if d['capacity']>d['flow'] and v not in side:side.add(v);q.append(v)
        records=[(u,v,d['capacity'],max(0,residual[u][v]['flow'])) for u,v,d in G.edges(data=True)]
    else:raise ValueError('unknown backend')
    t2=perf_counter();proof={'n':g['n'],'s':g['s'],'t':g['t'],'arcs':records,'source_side':sorted(side),'value':int(opt),'scale':g['scale'],'constant':str(g['constant'])}
    check_certificate(proof);t3=perf_counter()
    result={'feasible':opt<g['big'],'nodes':g['n'],'arcs':len(g['edges']),'network_seconds':t1-start,'flow_seconds':t2-t1,'flow_check_seconds':t3-t2}
    if certificate:result['certificate']=proof
    if not result['feasible']:result['seconds']=perf_counter()-start;return result
    p=c.problem;owners=p.decode([int(v in side) for v in range(p.n*p.m)])
    if c.fixed_answer is not None:
        w=c.fixed_answer;area=w['area'];boundaries=w['boundaries']
        objective=F(area_weight)*area+sum(F(row[o]) for row,o in zip(costs,owners))
        if owners!=w['owners'] or objective!=F(opt,g['scale'])+g['constant']:raise AssertionError('fixed outcome mismatch')
        t4=perf_counter();verify_mixed(p,owners,boundaries,area);t5=perf_counter()
        result.update(owners=owners,area=area,objective=objective,trapezoid_area=area,boundaries=boundaries,
                      reconstruct_seconds=t4-t3,boundary_check_seconds=t5-t4,seconds=t5-start)
        return result
    heights=[[max(x for x,v in zip(c.candidates[k][j],c.ids[k][j]) if v==g['s'] or v in side) for j in range(len(c.mesh))] for k in range(p.m)]
    boundaries=[[(a,b,x,y) for a,b,x,y in zip(c.mesh,c.mesh[1:],h,h[1:])] for h in heights]
    trap=F(0);actual=F(0)
    for k,h in enumerate(heights):
        for a,b,x,y in zip(c.mesh,c.mesh[1:],h,h[1:]):
            u=x-value(p.baselines[k],a);v=y-value(p.baselines[k],b)
            trap+=(b-a)*(abs(u)+abs(v))/2
            actual+=(b-a)*(u*u+v*v)/(2*(abs(u)+abs(v))) if u*v<0 else (b-a)*(abs(u)+abs(v))/2
    objective=F(area_weight)*actual+sum(F(row[o]) for row,o in zip(costs,owners))
    node_objective=F(opt,g['scale'])+g['constant']
    if F(area_weight)>0 and objective!=node_objective:raise AssertionError(('mesh/continuum cost mismatch',objective,node_objective))
    if actual>trap:raise AssertionError('invalid trapezoid upper bound')
    t4=perf_counter();verify_mixed(p,owners,boundaries,actual);t5=perf_counter()
    result.update(owners=owners,area=actual,objective=objective,trapezoid_area=trap,boundaries=boundaries,
                  reconstruct_seconds=t4-t3,boundary_check_seconds=t5-t4,seconds=t5-start)
    return result
