"""Exact continuous ordered-cell ownership with robust rectangular RNA sites.

Heterogeneous physical slope bounds; rational piecewise-linear source surfaces;
constant lower/upper separation bounds. Costs are additive in each group's cell
label plus a nonnegative weighted sum of swept interface areas. This is NOT an
arbitrary-cell tessellation solver, probabilistic uncertainty calibration, or
hard-budget optimizer. Coverage/cut and max-flow machinery are established.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from fractions import Fraction as F
from bisect import bisect_right
from itertools import combinations
from math import lcm
from time import perf_counter
from collections import deque
from copy import deepcopy
from pathlib import Path
import sys
from continuum_cut import Dinic, AreaTerm
from check_flow import check_certificate


def line_at(knots, t):
    j=max(0,min(len(knots)-2,bisect_right([x for x,_ in knots],t)-1))
    a,u=knots[j];b,v=knots[j+1];m=(v-u)/(b-a)
    return m,u-m*a

def value(knots,t):
    m,c=line_at(knots,t);return m*t+c

def dist_intervals(a,b,c,d):
    return max(F(0),c-b,a-d)

@dataclass(frozen=True)
class Box:
    xlo:F
    xhi:F
    tlo:F
    thi:F
    group:int
    def __post_init__(self):
        for a in ('xlo','xhi','tlo','thi'):object.__setattr__(self,a,F(getattr(self,a)))
        if self.xlo>self.xhi or self.tlo>self.thi:raise ValueError('reversed location region')

@dataclass
class Layers:
    baselines:list
    boxes:list[Box]
    source:list[int]
    slope:F|list|tuple
    minimum_gaps:list
    maximum_gaps:list|None=None
    clearance:F=F(1,100)
    bands:list|None=None
    mode:str="universal"
    closure:list=field(init=False,repr=False)
    slopes:list=field(init=False,repr=False)
    kernels:list=field(init=False,repr=False)
    def __post_init__(self):
        if self.mode not in ("universal","nominal","possible_outer"):raise ValueError("invalid uncertainty quantifier")
        self.baselines=[[(F(t),F(x)) for t,x in b] for b in self.baselines]
        self.slopes=[F(v) for v in self.slope] if isinstance(self.slope,(list,tuple)) else [F(self.slope)]*self.m
        if len(self.slopes)!=self.m or any(v<=0 for v in self.slopes):raise ValueError("positive slope per boundary required")
        self.slope=tuple(self.slopes);self.clearance=F(self.clearance)
        self.minimum_gaps=list(map(F,self.minimum_gaps));m=self.m
        if m<1 or self.clearance<=0:raise ValueError('positive slope/clearance and nonempty layers required')
        if len(self.minimum_gaps)!=m-1 or any(x<=0 for x in self.minimum_gaps):raise ValueError('strictly positive adjacent gaps required')
        if self.maximum_gaps is None:self.maximum_gaps=[None]*(m-1)
        if len(self.maximum_gaps)!=m-1:raise ValueError('wrong maximum gap count')
        self.maximum_gaps=[None if d is None else F(d) for d in self.maximum_gaps]
        if any(hi is not None and hi<lo for lo,hi in zip(self.minimum_gaps,self.maximum_gaps)):raise ValueError('empty separation interval')
        if self.bands is None:
            low=min(x for b in self.baselines for _,x in b)-10
            high=max(x for b in self.baselines for _,x in b)+10
            self.bands=[(low,high)]*m
        self.bands=[(F(a),F(b)) for a,b in self.bands]
        if len(self.bands)!=m:raise ValueError('wrong band count')
        for b,(lo,hi),slope in zip(self.baselines,self.bands,self.slopes):
            if len(b)<2 or any(t>=u for (t,_),(u,_) in zip(b,b[1:])):raise ValueError('baseline knots not increasing')
            if b[0][0]!=self.lo or b[-1][0]!=self.hi:raise ValueError('different longitudinal domains')
            if lo>=hi or any(not lo<=x<=hi for _,x in b):raise ValueError('baseline outside band')
            if any(abs(v-u)>slope*(bb-a) for (a,u),(bb,v) in zip(b,b[1:])):raise ValueError('source exceeds its declared slope')
        if self.lo>=self.hi:raise ValueError('empty domain')
        D=[[None]*m for _ in range(m)]
        for k in range(m):D[k][k]=F(0)
        for k,d in enumerate(self.minimum_gaps):
            D[k][k+1]=d
            if self.maximum_gaps[k] is not None:D[k+1][k]=-self.maximum_gaps[k]
        for k in range(m):
            for a in range(m):
                if D[a][k] is None:continue
                for b in range(m):
                    if D[k][b] is None:continue
                    w=D[a][k]+D[k][b]
                    if D[a][b] is None or w>D[a][b]:D[a][b]=w
        if any(D[k][k]>0 for k in range(m)):raise ValueError('positive difference-constraint cycle')
        self.closure=D
        # A path may travel longitudinally on a smoother neighboring surface.
        # Its exact best propagation is max_c D[a,c]+D[c,b]-s[c]*distance.
        kernels=[[[] for _ in range(m)] for _ in range(m)]
        for a in range(m):
            for b in range(m):
                by_rate={}
                for c in range(m):
                    if D[a][c] is None or D[c][b] is None:continue
                    h=D[a][c]+D[c][b];rate=self.slopes[c]
                    by_rate[rate]=max(by_rate.get(rate,h),h)
                candidates=[(h,rate) for rate,h in by_rate.items()]
                kernels[a][b]=[(h,rate) for h,rate in candidates if not any(hh>=h and ss<=rate and (hh,ss)!=(h,rate) for hh,ss in candidates)]
        self.kernels=kernels
        knots=sorted({t for b in self.baselines for t,_ in b})
        for t in knots:
            vals=[value(b,t) for b in self.baselines]
            if any(D[a][b] is not None and vals[b]-vals[a]<D[a][b] for a in range(m) for b in range(m)):
                raise ValueError('source violates separation constraints')
        if any(type(c) is not int or not 0<=c<=m for c in self.source):raise ValueError('invalid source cell')
        if any(type(z.group) is not int or not 0<=z.group<len(self.source) for z in self.boxes):raise ValueError('invalid group')
        if {b.group for b in self.boxes}!=set(range(len(self.source))):raise ValueError('group must contain a site')
        for z in self.boxes:
            if not self.lo<=z.tlo<=z.thi<=self.hi:raise ValueError('box outside domain')
            # The source must match the observed center, NOT every uncertain
            # position. Regions may cross source boundaries; universal feasibility
            # may then require deformation or be impossible for every owner.
            center_x=(z.xlo+z.xhi)/2;center_t=(z.tlo+z.thi)/2
            for k,b in enumerate(self.baselines):
                f=value(b,center_t)
                if self.source[z.group]<=k:
                    if f<center_x:raise ValueError('source label disagrees with nominal site')
                elif f>center_x:raise ValueError('source label disagrees with nominal site')
    @property
    def m(self):return len(self.baselines)
    @property
    def n(self):return len(self.source)
    @property
    def lo(self):return self.baselines[0][0][0]
    @property
    def hi(self):return self.baselines[0][-1][0]
    def bit(self,g,k):return g*self.m+k
    def obligation(self,z,k=0):
        if self.mode=="universal":return z.tlo,z.thi,z.xhi+self.clearance,z.xlo-self.clearance
        t=(z.tlo+z.thi)/2
        if self.mode=="nominal":
            x=(z.xlo+z.xhi)/2
            return t,t,x+self.clearance,x-self.clearance
        d=(z.thi-z.tlo)/2*self.slopes[k]
        return t,t,z.xlo+self.clearance-d,z.xhi-self.clearance+d
    def source_bits(self):return [int(c<=k) for c in self.source for k in range(self.m)]
    def encode(self,owners):
        if len(owners)!=self.n or any(type(c) is not int or not 0<=c<=self.m for c in owners):raise ValueError('invalid target owners')
        return [int(c<=k) for c in owners for k in range(self.m)]
    def decode(self,bits):
        ans=[]
        for g in range(self.n):
            bs=bits[g*self.m:(g+1)*self.m]
            if any(a>b for a,b in zip(bs,bs[1:])):raise ValueError('nonmonotone owner encoding')
            ans.append(sum(1-b for b in bs))
        return ans

@dataclass(frozen=True)
class Tent:
    lo:F
    hi:F
    height:F
    variable:int
    rate:F|None=None
    def line(self,t,s):
        if self.rate is not None:s=self.rate
        if t<self.lo:return s,self.height-s*self.lo
        if t>self.hi:return -s,self.height+s*self.hi
        return F(0),self.height
    def at(self,t,s):
        a,b=self.line(t,s);return a*t+b


def _terms(tents,base,s,sign):
    """Integrate arbitrary source-relative flat-top cones exactly.

    sign=-1 gives lower-cone minus source, sign=1 negated-upper plus source.
    This is a reference exact arrangement, not the previous point-support bound.
    Boxes can make support disconnected when the source oscillates inside them.
    """
    if not tents:return {},0,0
    lo,hi=base[0][0],base[-1][0]
    # Remove curves that are never positive; retain all crossings for the others.
    kept=[]
    for z in tents:
        ts=sorted({lo,hi,z.lo,z.hi}|{t for t,_ in base})
        if any(z.at(t,s)+sign*value(base,t)>0 for t in ts):kept.append(z)
    tents=kept
    if not tents:return {},0,0
    points={lo,hi}|{t for t,_ in base}
    for z in tents:
        points.update((z.lo,z.hi))
        ts=sorted({lo,hi,z.lo,z.hi}|{t for t,_ in base})
        for a,b in zip(ts,ts[1:]):
            mid=(a+b)/2;tm,tc=z.line(mid,s);bm,bc=line_at(base,mid)
            mm=tm+sign*bm;cc=tc+sign*bc
            if mm:
                t=-cc/mm
                if a<t<b:points.add(t)
    for zi,zj in combinations(tents,2):
        ts=sorted({lo,hi,zi.lo,zi.hi,zj.lo,zj.hi})
        for a,b in zip(ts,ts[1:]):
            mid=(a+b)/2;mi,ci=zi.line(mid,s);mj,cj=zj.line(mid,s)
            if mi!=mj:
                t=(cj-ci)/(mi-mj)
                if a<t<b and zi.at(t,s)+sign*value(base,t)>0:points.add(t)
    points=sorted(points);out={};raw=0
    for a,b in zip(points,points[1:]):
        mid=(a+b)/2;bm,bc=line_at(base,mid)
        lines=[(mm+sign*bm,cc+sign*bc) for mm,cc in (z.line(mid,s) for z in tents)]
        active=sorted([i for i,(mm,cc) in enumerate(lines) if mm*mid+cc>0],key=lambda i:(-(lines[i][0]*mid+lines[i][1]),i))
        mask=0
        for j,i in enumerate(active):
            nxt=mask|1<<tents[i].variable
            if nxt not in out:out[nxt]=AreaTerm(nxt,mask,tents[i].variable)
            mask=nxt
            mm,cc=lines[i]
            ml,cl=lines[active[j+1]] if j+1<len(active) else (F(0),F(0))
            weight=((mm-ml)*(a+b)+2*(cc-cl))*(b-a)/2
            if weight<0:raise ArithmeticError('negative arrangement layer')
            out[mask].area+=weight;raw+=1
    return out,len(points)-1,raw

@dataclass
class CompiledLayers:
    problem:Layers
    implies:tuple
    forced:dict
    positive:dict
    negative:dict
    stats:dict
    inconsistent:bool=False
    def feasible(self,owners):
        x=self.problem.encode(owners)
        return not self.inconsistent and all(x[a]<=x[b] for a,b in self.implies) and all(x[i]==v for i,v in self.forced.items())
    def area(self,owners):
        x=self.problem.encode(owners);mask=sum(1<<i for i,v in enumerate(x) if v);allmask=(1<<len(x))-1
        return sum((z.area for z in self.positive.values() if z.mask&mask),F(0))+sum((z.area for z in self.negative.values() if z.mask&(mask^allmask)),F(0))


def compile_layers(p:Layers):
    start=perf_counter();p=deepcopy(p);m=p.m;D=p.closure;edges=set();fixed={};inconsistent=False
    def force(i,v):
        nonlocal inconsistent
        if i in fixed and fixed[i]!=v:inconsistent=True
        fixed[i]=v
    for g in range(p.n):
        for k in range(m-1):edges.add((p.bit(g,k),p.bit(g,k+1)))
    # Exact cross-surface path kernels, including heterogeneous slope effects.
    for z in p.boxes:
        for a in range(m):
            vid=p.bit(z.group,a);zl,zh,lower,upper=p.obligation(z,a)
            for b in range(m):
                for offset,rate in p.kernels[a][b]:
                    h=lower+offset
                    if h>p.bands[b][1] or any(h-end>dist_intervals(zl,zh,t,t)*rate for t,end in (p.baselines[b][0],p.baselines[b][-1])):force(vid,0)
                for offset,rate in p.kernels[b][a]:
                    h=upper-offset
                    if h<p.bands[b][0] or any(end-h>dist_intervals(zl,zh,t,t)*rate for t,end in (p.baselines[b][0],p.baselines[b][-1])):force(vid,1)
    for zi in p.boxes:
        for zj in p.boxes:
            for a in range(m):
                il,ih,lower,_=p.obligation(zi,a)
                for b in range(m):
                    jl,jh,_,upper=p.obligation(zj,b)
                    distance=dist_intervals(il,ih,jl,jh)
                    if any(lower+offset-upper>rate*distance for offset,rate in p.kernels[a][b]):
                        u,v=p.bit(zi.group,a),p.bit(zj.group,b)
                        if u!=v:edges.add((u,v))
    pos={};neg={};slabs=raw=0;cone_count=0
    def merge(dst,src):
        for mask,t in src.items():
            if mask not in dst:dst[mask]=AreaTerm(mask,t.parent,t.added)
            dst[mask].area+=t.area
    for k in range(m):
        positives=[];negatives=[]
        for z in p.boxes:
            for a in range(m):
                v=p.bit(z.group,a);zl,zh,lower,upper=p.obligation(z,a)
                positives.extend(Tent(zl,zh,lower+offset,v,rate) for offset,rate in p.kernels[a][k])
                negatives.extend(Tent(zl,zh,-upper+offset,v,rate) for offset,rate in p.kernels[k][a])
        for tents,sign,dst in ((positives,-1,pos),(negatives,1,neg)):
            tt,ss,rr=_terms(tents,p.baselines[k],p.slopes[k],sign)
            merge(dst,tt);slabs+=ss;raw+=rr;cone_count+=len(tents)
    out=CompiledLayers(p,tuple(sorted(edges)),fixed,pos,neg,{'groups':p.n,'sites':len(p.boxes),'boundaries':m,'owner_bits':m*p.n,'implications':len(edges),'forced_bits':len(fixed),'positive_terms':len(pos),'negative_terms':len(neg),'propagated_cones':cone_count,'slabs':slabs,'raw_layers':raw,'compile_seconds':perf_counter()-start},inconsistent)
    return out


def graph_for(c:CompiledLayers,costs,area_weight=1,pins=None):
    p=c.problem;n=p.n*p.m;lam=F(area_weight)
    if lam<0 or len(costs)!=p.n or any(len(row)!=p.m+1 for row in costs):raise ValueError('invalid owner costs or area weight')
    S=n;T=n+1;nextid=n+2;finite=[];hard=list(c.implies);constant=F(0)
    if c.inconsistent:hard.append((S,T))
    for g,row in enumerate(costs):
        row=list(map(F,row));constant+=row[-1]
        for k in range(p.m):
            coef=row[k]-row[k+1];v=p.bit(g,k)
            if coef>=0:finite.append((v,T,coef))
            else:constant+=coef;finite.append((S,v,-coef))
    for polarity,terms in ((0,c.positive),(1,c.negative)):
        ids={}
        for mask in sorted(terms,key=lambda x:(x.bit_count(),x)):
            z=terms[mask];v=nextid;nextid+=1;ids[mask]=v
            if lam*z.area:finite.append((v,T,lam*z.area) if polarity==0 else (S,v,lam*z.area))
            deps=[z.added]+([ids[z.parent]] if z.parent else [])
            hard.extend((u,v) if polarity==0 else (v,u) for u in deps)
    for i,v in c.forced.items():hard.append((S,i) if v else (i,T))
    for g,owner in (pins or {}).items():
        if not 0<=g<p.n or not 0<=owner<=p.m:raise ValueError('invalid owner pin')
        for k in range(p.m):
            i=p.bit(g,k);hard.append((S,i) if owner<=k else (i,T))
    scale=1
    for _,_,a in finite:scale=lcm(scale,a.denominator)
    finite=[(u,v,int(a*scale)) for u,v,a in finite if a]
    big=sum(a for _,_,a in finite)+1
    edges=finite+[(u,v,big) for u,v in hard if u!=v]
    return {'n':nextid,'s':S,'t':T,'scale':scale,'constant':constant,'big':big,'edges':edges,'hard':hard}


def solve(c,costs,*,area_weight=1,pins=None,seconds=60,backend='dinic',certificate=True,reconstruct=True):
    start=perf_counter();data=graph_for(c,costs,area_weight,pins);S,T=data['s'],data['t']
    if backend=='dinic':
        d=Dinic(data['n'])
        for u,v,w in data['edges']:d.add(u,v,w)
        vv,side,records=d.run(S,T,start+seconds)
    elif backend=='networkx':
        import networkx as nx
        G=nx.DiGraph();G.add_nodes_from(range(data['n']))
        for u,v,w in data['edges']:
            if G.has_edge(u,v):G[u][v]['capacity']+=w
            else:G.add_edge(u,v,capacity=w)
        R=nx.algorithms.flow.preflow_push(G,S,T)
        vv=R.graph['flow_value'];side={S};queue=deque([S])
        while queue:
            u=queue.popleft()
            for v,attr in R[u].items():
                if attr['capacity']-attr['flow']>0 and v not in side:side.add(v);queue.append(v)
        records=[(u,v,attr['capacity'],max(0,R[u][v]['flow'])) for u,v,attr in G.edges(data=True)]
    else:raise ValueError('unknown backend')
    cert={'n':data['n'],'s':S,'t':T,'arcs':records,'source_side':sorted(side),'value':int(vv),'scale':data['scale'],'constant':str(data['constant'])}
    check_certificate(cert)
    result={'feasible':vv<data['big'],'nodes':data['n'],'arcs':len(data['edges']),'seconds':perf_counter()-start}
    if certificate:result['certificate']=cert
    if not result['feasible']:return result
    bits=[int(i in side) for i in range(c.problem.m*c.problem.n)];owners=c.problem.decode(bits)
    area=c.area(owners);obj=F(area_weight)*area+sum((F(row[q]) for row,q in zip(costs,owners)),F(0))
    if obj!=F(vv,data['scale'])+data['constant']:raise AssertionError('cut/objective disagreement')
    result.update(owners=owners,area=area,objective=obj)
    if reconstruct:
        witness=realize(c.problem,owners)
        if not witness['feasible'] or witness['area']!=area:raise AssertionError('cut/continuum disagreement')
        result['boundaries']=witness['boundaries']
    result['seconds']=perf_counter()-start
    return result


def realize(p:Layers,owners):
    """Direct fixed-owner all-cone oracle, separately from finite graph layers."""
    p.encode(owners);m=p.m;D=p.closure;allbounds=[];area=F(0)
    for k in range(m):
        lowers=[];uppers=[];lowconst=[];upconst=[]
        for a in range(m):
            if D[a][k] is not None:lowconst.append(p.bands[a][0]+D[a][k])
            if D[k][a] is not None:upconst.append(p.bands[a][1]-D[k][a])
            for offset,rate in p.kernels[a][k]:
                lowers.extend(Tent(t,t,v+offset,-1,rate) for t,v in (p.baselines[a][0],p.baselines[a][-1]))
            for offset,rate in p.kernels[k][a]:
                uppers.extend(Tent(t,t,-v+offset,-1,rate) for t,v in (p.baselines[a][0],p.baselines[a][-1]))
        for z in p.boxes:
            for a in range(m):
                zl,zh,lower,upper=p.obligation(z,a)
                if owners[z.group]<=a:lowers.extend(Tent(zl,zh,lower+offset,-1,rate) for offset,rate in p.kernels[a][k])
                else:uppers.extend(Tent(zl,zh,-upper+offset,-1,rate) for offset,rate in p.kernels[k][a])
        base=p.baselines[k];pts={t for t,_ in base}
        for z in lowers+uppers:pts.update((z.lo,z.hi))
        pts=sorted(pts);out=[]
        for a,b in zip(pts,pts[1:]):
            mid=(a+b)/2
            ll=[z.line(mid,p.slopes[k]) for z in lowers]+[(F(0),max(lowconst))]
            ul=[(-mm,-cc) for mm,cc in (z.line(mid,p.slopes[k]) for z in uppers)]+[(F(0),min(upconst))]
            # Same-slope lines are globally dominated on this entire slab.
            ldict={};udict={}
            for mm,cc in ll:ldict[mm]=max(ldict.get(mm,cc),cc)
            for mm,cc in ul:udict[mm]=min(udict.get(mm,cc),cc)
            ll=list(ldict.items());ul=list(udict.items())
            bl=line_at(base,mid);lines=ll+ul+[bl];cuts={a,b}
            for (mi,ci),(mj,cj) in combinations(lines,2):
                if mi!=mj:
                    t=(cj-ci)/(mi-mj)
                    if a<t<b:cuts.add(t)
            cuts=sorted(cuts)
            for aa,bb in zip(cuts,cuts[1:]):
                cc=(aa+bb)/2
                l=max(ll,key=lambda z:z[0]*cc+z[1]);u=min(ul,key=lambda z:z[0]*cc+z[1])
                if l[0]*cc+l[1]>u[0]*cc+u[1]:return {'feasible':False}
                fs=l if bl[0]*cc+bl[1]<l[0]*cc+l[1] else u if bl[0]*cc+bl[1]>u[0]*cc+u[1] else bl
                va=fs[0]*aa+fs[1];vb=fs[0]*bb+fs[1]
                delta_a=va-value(base,aa);delta_b=vb-value(base,bb)
                area+=(abs(delta_a)+abs(delta_b))*(bb-aa)/2
                if out and (out[-1][3]-out[-1][2])/(out[-1][1]-out[-1][0])==fs[0] and out[-1][1]==aa and out[-1][3]==va:
                    out[-1]=(out[-1][0],bb,out[-1][2],vb)
                else:out.append((aa,bb,va,vb))
        allbounds.append(out)
    return {'feasible':True,'area':area,'boundaries':allbounds}


def optimal_owner_domains(c,solution):
    """Exact coordinate-wise label domains across ALL minimum cuts.

    Uses the classical residual-closure characterization. Marginal domains do not
    mean every combination of allowed labels is jointly optimal or feasible.
    """
    if not solution.get('feasible') or 'certificate' not in solution:raise ValueError('feasible proof required')
    cert=solution['certificate'];check_certificate(cert);adj=[set() for _ in range(cert['n'])]
    for u,v,cap,f in cert['arcs']:
        if cap-f:adj[u].add(v)
        if f:adj[v].add(u)
    p=c.problem;domains=[]
    for g in range(p.n):
        allowed=[]
        for cell in range(p.m+1):
            ones={cert['s']}|{p.bit(g,k) for k in range(cell,p.m)}
            zeros={cert['t']}|{p.bit(g,k) for k in range(cell)}
            seen=set(ones);queue=deque(ones)
            while queue:
                u=queue.popleft()
                for v in adj[u]:
                    if v not in seen:seen.add(v);queue.append(v)
            if not seen&zeros:allowed.append(cell)
        if not allowed:raise AssertionError('optimum has no owner')
        domains.append(allowed)
    return domains
