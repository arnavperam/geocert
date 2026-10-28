"""Exact event-driven integration of maximum-envelope activation costs.

The geometry, feasibility rules and fixed-owner reconstruction are unchanged.
This module replaces only the slab-by-slab sorting/integration step. It uses
standard arrangement events and adjacent-order updates; no new max-flow method.
"""
from __future__ import annotations
from bisect import bisect_right
from collections import defaultdict, deque
from fractions import Fraction as F
from itertools import combinations
from copy import deepcopy
from time import perf_counter
from layered import Tent, AreaTerm
from pinned_compiler import feasibility, mandatory_envelope, PinnedCompiled


def _line_arrays(knots):
    ts=[t for t,_ in knots];lines=[]
    for (a,u),(b,v) in zip(knots,knots[1:]):
        s=(v-u)/(b-a);lines.append((s,u-s*a))
    return ts,lines


def event_terms(tents,base,default_rate,sign,*,debug=False):
    """Return integrated positive OR terms and exact update-count diagnostics.

    All tent plateaus and the reference are continuous. Stable total orders at
    event endpoints use right slopes, then input index. Group-label duplication
    is handled by prefix bitsets rather than independent molecular decisions.
    """
    lo,hi=base[0][0],base[-1][0];bt,bl=_line_arrays(base)
    def base_line(t):return bl[max(0,min(len(bl)-1,bisect_right(bt,t)-1))]
    # Safe inactive-cone pruning: its affine difference can attain a maximum
    # only at a source knot or one of its own two plateau endpoints.
    kept=[]
    for z in tents:
        candidates=set(bt)|{max(lo,min(hi,z.lo)),max(lo,min(hi,z.hi))}
        if any(z.at(t,default_rate)+sign*(base_line(t)[0]*t+base_line(t)[1])>0 for t in candidates):kept.append(z)
    tents=kept;n=len(tents);ref=n
    if not n:return {},{'cones':0,'events':0,'swaps':0,'edge_integrations':0,'prefixes':0}
    events=defaultdict(set)
    events[lo].update(range(n+1));events[hi]
    for i,z in enumerate(tents):
        for t in (z.lo,z.hi):
            if lo<=t<=hi:events[t].add(i)
    for t in bt:events[t].add(ref)
    for i,j in combinations(range(n),2):
        zi,zj=tents[i],tents[j]
        cuts=sorted({lo,hi}|{t for t in (zi.lo,zi.hi,zj.lo,zj.hi) if lo<t<hi})
        for a,b in zip(cuts,cuts[1:]):
            mi,ci=zi.line((a+b)/2,default_rate);mj,cj=zj.line((a+b)/2,default_rate)
            if mi!=mj:
                t=(cj-ci)/(mi-mj)
                if a<=t<=b:events[t].update((i,j))
    # Crossings with the common reference govern positivity, not cone ordering.
    for i,z in enumerate(tents):
        cuts=sorted(set(bt)|{t for t in (z.lo,z.hi) if lo<t<hi})
        for a,b in zip(cuts,cuts[1:]):
            mid=(a+b)/2;mi,ci=z.line(mid,default_rate);bm,bc=base_line(mid)
            mm,cc=mi+sign*bm,ci+sign*bc
            if mm:
                t=-cc/mm
                if a<=t<=b:events[t].update((i,ref))
    times=sorted(events);mid=(times[0]+times[1])/2
    lines=[z.line(mid,default_rate) for z in tents]
    bm,bc=base_line(mid);lines.append((-sign*bm,-sign*bc))
    def key(i,t):
        m,c=lines[i];return (-(m*t+c),-m,i)
    order=sorted(range(n+1),key=lambda i:key(i,lo));pos={v:i for i,v in enumerate(order)}
    bits=[1<<z.variable for z in tents]+[0];prefix=[];terms={};mask=0
    def register(new,parent,node):
        if new and new not in terms:
            assert bits[node] and new==parent|bits[node] and new!=parent
            terms[new]=AreaTerm(new,parent,tents[node].variable)
    for node in order:
        parent=mask;mask|=bits[node];register(mask,parent,node);prefix.append(mask)
    # State for each adjacent-order gap. Its integral is postponed until that
    # gap's affine expression, activation prefix, or positivity status changes.
    state=[None]*n;integrations=0;swaps=0
    def install(k,t):
        if 0<=k<n:
            if k<pos[ref]:
                i,j=order[k],order[k+1];mi,ci=lines[i];mj,cj=lines[j]
                state[k]=(prefix[k],mi-mj,ci-cj,t)
            else:state[k]=None
    def flush(k,t):
        nonlocal integrations
        old=state[k]
        if old is None:return
        mask,m,c,a=old
        if t>a:
            va,vb=m*a+c,m*t+c
            if va<0 or vb<0:raise ArithmeticError('negative ordered gap')
            terms[mask].area+=(va+vb)*(t-a)/2;integrations+=1
    for k in range(n):install(k,lo)
    for idx,t in enumerate(times[1:-1],1):
        affected=events[t];mid=(t+times[idx+1])/2;dirty=set();queue=deque();queued=set()
        def enqueue(k):
            if 0<=k<n and k not in queued:queue.append(k);queued.add(k)
        for node in affected:
            r=pos[node];dirty.update(k for k in (r-1,r) if 0<=k<n);enqueue(r-1);enqueue(r)
            if node==ref:
                bm,bc=base_line(mid);lines[node]=(-sign*bm,-sign*bc)
            else:lines[node]=tents[node].line(mid,default_rate)
        while queue:
            k=queue.popleft();queued.remove(k)
            i,j=order[k],order[k+1]
            if key(i,t)>key(j,t):
                order[k],order[k+1]=j,i;pos[i],pos[j]=k+1,k;swaps+=1
                parent=prefix[k-1] if k else 0;new=parent|bits[j]
                register(new,parent,j);prefix[k]=new
                for q in (k-1,k,k+1):
                    if 0<=q<n:dirty.add(q);enqueue(q)
        if debug:
            assert order==sorted(range(n+1),key=lambda i:key(i,t))
            mm=0
            for k,i in enumerate(order):mm|=bits[i];assert prefix[k]==mm
        for k in dirty:flush(k,t);install(k,t)
    for k in range(n):flush(k,hi)
    # Retain positive terms and the zero-weight predecessors they need.
    needed={mask for mask,z in terms.items() if z.area}
    todo=list(needed)
    while todo:
        par=terms[todo.pop()].parent
        if par and par not in needed:needed.add(par);todo.append(par)
    result={mask:terms[mask] for mask in needed}
    return result,{'cones':n,'events':len(times),'swaps':swaps,'edge_integrations':integrations,'prefixes':len(terms),'retained_terms':len(result)}


def compile_events(problem,pins=None,*,debug=False):
    start=perf_counter();p=deepcopy(problem);ff=feasibility(p,pins);after_feas=perf_counter()
    positive={};negative={};constant=F(0);pieces=[]
    def merge(dst,src):
        for mask,z in src.items():
            if mask not in dst:dst[mask]=AreaTerm(mask,z.parent,z.added)
            dst[mask].area+=z.area
    if not ff['inconsistent']:
        for k in range(p.m):
            pos=[];neg=[]
            for z in p.boxes:
                for a in range(p.m):
                    v=p.bit(z.group,a);zl,zh,L,U=p.obligation(z,a)
                    pos.extend(Tent(zl,zh,L+d,v,s) for d,s in p.kernels[a][k])
                    neg.extend(Tent(zl,zh,-U+d,v,s) for d,s in p.kernels[k][a])
            for tt,sign,dst,on in [(pos,-1,positive,1),(neg,1,negative,0)]:
                free=[];mandatory=[]
                for z in tt:
                    if z.variable not in ff['fixed']:free.append(z)
                    elif ff['fixed'][z.variable]==on:mandatory.append(z)
                base,area=mandatory_envelope(p.baselines[k],mandatory,p.slopes[k],sign);constant+=area
                new,stats=event_terms(free,base,p.slopes[k],sign,debug=debug);merge(dst,new);pieces.append(stats)
    stats={'compiler':'event-driven','groups':p.n,'sites':len(p.boxes),'boundaries':p.m,'source_knots':sum(len(b) for b in p.baselines),
           'positive_terms':len(positive),'negative_terms':len(negative),'events':sum(z['events'] for z in pieces),'swaps':sum(z['swaps'] for z in pieces),
           'edge_integrations':sum(z['edge_integrations'] for z in pieces),'direction_details':pieces,'feasibility_seconds':after_feas-start,
           'area_seconds':perf_counter()-after_feas,'compile_seconds':perf_counter()-start}
    return PinnedCompiled(p,ff['edges'],ff['fixed'],positive,negative,stats,ff['inconsistent'],constant,dict(pins or {}),ff)
