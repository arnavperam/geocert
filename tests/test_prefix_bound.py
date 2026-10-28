"""Independent rational test of the supplied prefix-family size refinement.
Truth integration is written here and does not use the production area oracle.
Tests cover simultaneous crossings, coincident functions, repeated activation
labels and many source knots. They do not establish priority or native anatomy.
"""
from common_test import *
from layered import Tent,_terms
from itertools import combinations
import random

def funline(z,t):
    if t<z.lo:return z.rate,z.height-z.rate*z.lo
    if t>z.hi:return -z.rate,z.height+z.rate*z.hi
    return F(0),z.height

def baseln(base,t):
    for (a,u),(b,v) in zip(base,base[1:]):
        if a<=t<=b:
            slope=(v-u)/(b-a);return slope,u-slope*a
    raise ValueError('outside source')

def partitions(tents,base=None,sign=1):
    grid={F(0),F(10)}|{x for z in tents for x in (z.lo,z.hi)}
    if base is not None:grid|={t for t,v in base}
    knots=sorted(grid);more=set(grid)
    for a,b in zip(knots,knots[1:]):
        t=(a+b)/2;lines=[funline(z,t) for z in tents]
        if base is not None:
            bm,bc=baseln(base,t);lines=[(m+sign*bm,c+sign*bc) for m,c in lines]+[(F(0),F(0))]
        for (m,c),(u,v) in combinations(lines,2):
            if m!=u:
                x=(v-c)/(m-u)
                if a<x<b:more.add(x)
    return sorted(more)

def direct_integral(tents,base,sign,mask):
    selected=[z for z in tents if mask>>z.variable&1];grid=partitions(selected,base,sign);area=F(0)
    for a,b in zip(grid,grid[1:]):
        t=(a+b)/2;bm,bc=baseln(base,t)
        lines=[(F(0),F(0))]+[(m+sign*bm,c+sign*bc) for m,c in [funline(z,t) for z in selected]]
        m,c=max(lines,key=lambda q:q[0]*t+q[1]);area+=(b-a)*(m*(a+b)/2+c)
    return area

rng=random.Random(2710287);out={'arrangements':0,'assignment_checks':0,'family_bound_checks':0,'helper_checks':0,'failures':0,'records':[]}
for seed in range(84):
    M=1+seed%8;G=1+seed%4;tents=[]
    for j in range(M):
        lo=F(rng.randrange(0,17),2);hi=min(F(10),lo+F(rng.randrange(0,4),2))
        tents.append(Tent(lo,hi,F(rng.randrange(-4,17),4),rng.randrange(G),F(rng.randrange(1,5),2)))
    if M>=3 and seed%3==0:tents[-1]=Tent(tents[0].lo,tents[0].hi,tents[0].height,(tents[0].variable+1)%G,tents[0].rate)
    if M>=4 and seed%4==0:tents[-2]=Tent(F(4),F(6),F(1),0,F(1))
    J=[2,5,21,41][seed%4];base=[(F(10*j,J-1),F(rng.randrange(-8,9),4)) for j in range(J)]
    grid=partitions(tents);perms=[];prefix=set();swaps=0;last=None
    for a,b in zip(grid,grid[1:]):
        t=(a+b)/2;perm=sorted(range(M),key=lambda i:(-(funline(tents[i],t)[0]*t+funline(tents[i],t)[1]),i))
        if last is not None:
            pos={v:i for i,v in enumerate(last)};swaps+=sum(pos[perm[i]]>pos[perm[j]] for i in range(M) for j in range(i+1,M))
        last=perm;mask=0
        for i in perm:mask|=1<<tents[i].variable;prefix.add(mask)
    assert len(prefix)<=M+swaps
    assert swaps<=9*M*(M-1)//2
    out['family_bound_checks']+=1
    for sign in [-1,1]:
        terms,slabs,raw=_terms(tents,base,F(1),sign)
        assert set(terms)<=prefix
        for mask,z in terms.items():
            assert z.parent==0 or z.parent in terms
            assert not z.parent>>z.added&1 and mask==(z.parent|(1<<z.added))
            assert z.area>=0;out['helper_checks']+=1
        for mask in range(1<<G):
            result=sum((z.area for S,z in terms.items() if S&mask),F(0))
            assert result==direct_integral(tents,base,sign,mask),(seed,sign,mask)
            out['assignment_checks']+=1
        out['arrangements']+=1
        out['records'].append({'seed':seed,'cones':M,'source_knots':J,'polarity':sign,'slabs':slabs,'raw_layers':raw,'helpers':len(terms),'all_prefix_sets':len(prefix),'adjacent_swaps':swaps})
save7('prefix_validation.json',out);print({k:v for k,v in out.items() if k!='records'})
