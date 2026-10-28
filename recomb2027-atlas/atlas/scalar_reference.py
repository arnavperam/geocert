"""Independent scalar reference for source-flat, two-region atlas patches.

No production geometry, implication, envelope, integration or cut modules are
imported. Used for validation, not as the primary query implementation.
"""
from fractions import Fraction as F
from itertools import product


def distance(a,b,c,d):return max(F(0),c-b,a-d)


def prepare(p):
    n=len(p['points']);rho=F(p['clearance_pixels']);T=F(p['domain_pixels'][1]);s=F(p['slope']);band=F(p['band_halfwidth_pixels'])
    low=[(F(0),F(0),F(0)),(T,T,F(0)),(F(0),T,-band)]
    up=[(F(0),F(0),F(0)),(T,T,F(0)),(F(0),T,band)]
    for side,v,a,b in p['nuclear_constraints']:(low if side==0 else up).append((F(a),F(b),F(v)))
    for a,b,v in low:
        for c,d,w in up:
            if v-s*distance(a,b,c,d)>w:raise ValueError('inconsistent structural model')
    pts=[(F(r['t']),F(r['x'])) for r in p['points']]
    ones=set();zeros=set();edges=[]
    for i,(t,x) in enumerate(pts):
        if any(x+rho-s*distance(t,t,a,b)>v for a,b,v in up):zeros.add(i)
        if any(v-s*distance(a,b,t,t)>x-rho for a,b,v in low):ones.add(i)
        for j,(tt,xx) in enumerate(pts):
            if i!=j and x+rho-s*abs(t-tt)>xx-rho:edges.append((i,j))
    return n,rho,T,s,pts,ones,zeros,edges


def tents_area(terms,T,s):
    if not terms:return F(0)
    ts=sorted({F(0),T}|{t for t,h in terms if 0<t<T})
    area=F(0)
    for a,b in zip(ts,ts[1:]):
        mid=(a+b)/2;lines=[(F(0),F(0))]
        for t,h in terms:
            lines.append((s,h-s*t) if mid<t else (-s,h+s*t))
        events={a,b}
        for i,(m,c) in enumerate(lines):
            for mm,cc in lines[:i]:
                if m!=mm:
                    u=(cc-c)/(m-mm)
                    if a<u<b:events.add(u)
        es=sorted(events)
        for u,v in zip(es,es[1:]):
            m,c=max(lines,key=lambda mc:mc[0]*(u+v)/2+mc[1])
            area+=(v-u)*(m*(u+v)/2+c)
    return area


def query(p,pins,prepared=None):
    n,rho,T,s,pts,initial1,initial0,edges=prepared or prepare(p)
    one=set(initial1);zero=set(initial0)
    for i,owner in pins.items():(one if owner==0 else zero).add(i)
    change=True
    while change:
        old=(len(one),len(zero))
        for u,v in edges:
            if u in one:one.add(v)
            if v in zero:zero.add(u)
        change=(len(one),len(zero))!=old
    if one&zero:return None
    own=[r['source'] for r in p['points']]
    for i in one:own[i]=0
    for i in zero:own[i]=1
    plus=[(t,x+rho) for (t,x),o in zip(pts,own) if o==0 and x+rho>0]
    minus=[(t,rho-x) for (t,x),o in zip(pts,own) if o==1 and rho-x>0]
    area=(tents_area(plus,T,s)+tents_area(minus,T,s))*F(p['area_scale_um2'])
    return own,area


def spectrum(p,gene):
    pre=prepare(p);ids=[i for i,r in enumerate(p['points']) if r['gene']==gene]
    out={};checked=0
    for oo in product([0,1],repeat=len(ids)):
        ans=query(p,dict(zip(ids,oo)),pre);checked+=1
        if ans is None:continue
        owners,area=ans;value=sum(owners[i]==0 for i in ids)
        out[value]=min(out.get(value,area),area)
    return out,checked
