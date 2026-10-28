"""Exact count spectra by a resource-augmented critical-height chain.

For one boundary, point-valued free RNA groups and fixed structural regions,
count is a local nodal statistic. A count-indexed shortest-path dynamic program
therefore obtains every minimum area without enumerating ownership vectors.
The canonical-envelope theorem supplies an exact mesh, not an approximation.

This module does not use min-cut, the activation-layer compiler or pin-clipping.
It shares the physical model and direct final boundary checker.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from bisect import bisect_right, bisect_left
from collections import defaultdict, deque
from itertools import combinations
from time import perf_counter
from geocert_expression import ROOT as _ENGINE_BOOTSTRAP
from mixed import verify_mixed
from layered import value


@dataclass
class Function:
    pieces: list[tuple[F,F,F,F]]
    def at(self,t):
        starts=[a for a,b,m,c in self.pieces]
        j=max(0,min(len(self.pieces)-1,bisect_right(starts,t)-1))
        a,b,m,c=self.pieces[j]
        return m*t+c


def cone(a,b,height,slope,lo,hi,upper=False):
    a,b,height,slope,lo,hi=map(F,(a,b,height,slope,lo,hi));sgn=-1 if upper else 1
    knots=sorted({lo,hi}|{t for t in (a,b) if lo<t<hi});pieces=[]
    for u,v in zip(knots,knots[1:]):
        t=(u+v)/2
        if t<a:m,c=sgn*slope,height-sgn*slope*a
        elif t>b:m,c=-sgn*slope,height+sgn*slope*b
        else:m,c=F(0),height
        pieces.append((u,v,m,c))
    return Function(pieces)


@dataclass
class ChainMesh:
    problem: object
    pins: dict
    times: list[F]
    heights: list[list[F]]
    owners_at: list[list[dict[int,int]]]
    point_groups: list[list[int]]
    node_costs: list[list[F]]
    stats: dict
    feasible: bool


def compile_chain(problem,pins=None,*,seconds=120,candidate_cap=1000000):
    start=perf_counter();deadline=start+seconds;p=problem;pins=dict(pins or {})
    if p.m!=1:raise ValueError('count-chain requires exactly one boundary')
    if any(g<0 or g>=p.n or o not in (0,1) for g,o in pins.items()):raise ValueError('invalid owner pin')
    bygroup=defaultdict(list)
    for z in p.boxes:bygroup[z.group].append(z)
    free_points={}
    for g in range(p.n):
        if g in pins:continue
        if len(bygroup[g])!=1:raise ValueError('each unpinned group must be one point observation')
        z=bygroup[g][0]
        if z.xlo!=z.xhi or z.tlo!=z.thi:
            raise ValueError('count-chain requires explicitly point-valued free RNA; location boxes are unsupported')
        a,b,L,U=p.obligation(z,0)
        if a!=b or L<=U:
            raise ValueError('free point owner thresholds must be strictly separated at one node')
        free_points[g]=(a,L,U)
    lo,hi,s=p.lo,p.hi,p.slopes[0];functions=[];lower=[];upper=[]
    def fixed(fn,is_lower):
        functions.append(fn);(lower if is_lower else upper).append(fn)
    lb,ub=p.bands[0]
    fixed(Function([(lo,hi,F(0),lb)]),True);fixed(Function([(lo,hi,F(0),ub)]),False)
    for t,h in (p.baselines[0][0],p.baselines[0][-1]):
        fixed(cone(t,t,h,s,lo,hi),True);fixed(cone(t,t,h,s,lo,hi,True),False)
    for z in p.boxes:
        a,b,L,U=p.obligation(z,0);g=z.group
        if g not in pins or pins[g]==0:
            fn=cone(a,b,L,s,lo,hi);functions.append(fn)
            if g in pins:lower.append(fn)
        if g not in pins or pins[g]==1:
            fn=cone(a,b,U,s,lo,hi,True);functions.append(fn)
            if g in pins:upper.append(fn)
    basepieces=[]
    for (a,x),(b,y) in zip(p.baselines[0],p.baselines[0][1:]):
        m=(y-x)/(b-a);basepieces.append((a,b,m,x-m*a))
    functions.append(Function(basepieces))
    pieces=sorted(set(q for fn in functions for q in fn.pieces))
    events={lo,hi}|{t for q in pieces for t in q[:2]}
    events.update(t for t,L,U in free_points.values())
    pairchecks=0;retained=0
    for i,(a,b,m,c) in enumerate(pieces):
        if i%8==0 and perf_counter()>deadline:raise TimeoutError('critical-chain mesh limit')
        for aa,bb,mm,cc in pieces[:i]:
            left,right=max(a,aa),min(b,bb)
            if left>=right or m==mm:continue
            pairchecks+=1;t=(cc-c)/(m-mm)
            if not left<t<right:continue
            h=m*t+c
            if not lb<=h<=ub:continue
            if h<max(fn.at(t) for fn in lower) or h>min(fn.at(t) for fn in upper):continue
            events.add(t);retained+=1
    times=sorted(events);at=defaultdict(list)
    for g,(t,L,U) in free_points.items():at[t].append(g)
    heights=[];owners=[];costs=[];total=0
    for j,t in enumerate(times):
        if j%32==0 and perf_counter()>deadline:raise TimeoutError('critical-chain candidate limit')
        L=max(fn.at(t) for fn in lower);U=min(fn.at(t) for fn in upper)
        vals=sorted({fn.at(t) for fn in functions});vv=[];oo=[]
        for h in vals:
            if not L<=h<=U:continue
            assignment={};valid=True
            for g in at[t]:
                _,reqL,reqU=free_points[g]
                if h>=reqL:assignment[g]=0
                elif h<=reqU:assignment[g]=1
                else:valid=False;break
            if valid:vv.append(h);oo.append(assignment)
        heights.append(vv);owners.append(oo);total+=len(vv)
        if total>candidate_cap:raise MemoryError('critical-chain candidate cap')
        weight=((t-times[j-1]) if j else F(0))+((times[j+1]-t) if j+1<len(times) else F(0));weight/=2
        costs.append([weight*abs(h-value(p.baselines[0],t)) for h in vv])
    stats={'mesh_nodes':len(times),'candidate_heights':total,'max_heights_per_node':max(map(len,heights),default=0),
           'affine_pieces':len(pieces),'intersection_tests':pairchecks,'retained_intersection_occurrences':retained,
           'compile_seconds':perf_counter()-start}
    return ChainMesh(p,pins,times,heights,owners,[at[t] for t in times],costs,stats,all(heights))


def spectrum(mesh,scores,constant=0,*,seconds=120,state_cap=2000000):
    start=perf_counter();p=mesh.problem;deadline=start+seconds
    if len(scores)!=p.n or any(len(row)!=2 or any(type(v) is not int or v<0 for v in row) for row in scores):
        raise ValueError('nonnegative integer score for each group and owner required')
    if type(constant) is not int:raise TypeError('integer fixed count required')
    if not mesh.feasible:return {'complete':True,'minima':{},'witnesses':{},'status':'infeasible_model'}
    offset=constant+sum(scores[g][o] for g,o in mesh.pins.items())
    contrib=[[sum(scores[g][o] for g,o in assignment.items()) for assignment in oo] for oo in mesh.owners_at]
    state={(h,offset+contrib[0][h]):cost for h,cost in enumerate(mesh.node_costs[0])}
    parents=[{}];visited=len(state);transitions=0
    for j in range(1,len(mesh.times)):
        if perf_counter()>deadline:raise TimeoutError('count-chain DP limit; no completed spectrum')
        allowed=p.slopes[0]*(mesh.times[j]-mesh.times[j-1]);oldh=mesh.heights[j-1];newh=mesh.heights[j]
        nxt={};parent={};qs=sorted({q for h,q in state})
        for q in qs:
            end=0;queue=deque()
            for h,x in enumerate(newh):
                while end<len(oldh) and oldh[end]<=x+allowed:
                    c=state.get((end,q))
                    if c is not None:
                        while queue and queue[-1][1]>c:queue.pop()
                        queue.append((end,c))
                    end+=1
                while queue and oldh[queue[0][0]]<x-allowed:queue.popleft()
                if not queue:continue
                prev,c=queue[0];v=q+contrib[j][h];key=(h,v);cost=c+mesh.node_costs[j][h];transitions+=1
                if key not in nxt or cost<nxt[key]:nxt[key]=cost;parent[key]=(prev,q)
        state=nxt;parents.append(parent);visited+=len(nxt)
        if visited>state_cap:raise MemoryError('count-chain state cap')
        if not state:break
    minima={};endkeys={}
    for key,cost in state.items():
        h,q=key
        if q not in minima or cost<minima[q]:minima[q]=cost;endkeys[q]=key
    witnesses={}
    for q,cost in sorted(minima.items()):
        key=endkeys[q];indices=[]
        for j in range(len(mesh.times)-1,-1,-1):
            indices.append(key[0])
            if j:key=parents[j][key]
        indices.reverse();hs=[row[i] for row,i in zip(mesh.heights,indices)];owners=[None]*p.n
        for g,o in mesh.pins.items():owners[g]=o
        for j,i in enumerate(indices):
            for g,o in mesh.owners_at[j][i].items():owners[g]=o
        if any(o is None for o in owners):raise AssertionError('group counted zero or multiple times')
        bounds=[[(a,b,u,v) for a,b,u,v in zip(mesh.times,mesh.times[1:],hs,hs[1:])]]
        # The direct geometric checker integrates absolute displacement exactly.
        verify_mixed(p,owners,bounds,cost)
        if constant+sum(row[o] for row,o in zip(scores,owners))!=q:raise AssertionError('count mismatch')
        if any(owners[g]!=o for g,o in mesh.pins.items()):raise AssertionError('fixed owner violated')
        witnesses[q]={'owners':owners,'boundaries':bounds,'area':cost}
    return {'complete':True,'status':'exact' if minima else 'infeasible_model','minima':minima,'witnesses':witnesses,
            'DP_states':visited,'feasible_transition_evaluations':transitions,'seconds':perf_counter()-start}
