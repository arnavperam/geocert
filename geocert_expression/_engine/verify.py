"""Original-coordinate verifier; does not import the compiler or flow solver."""
from fractions import Fraction as F
from bisect import bisect_right

def at_segments(segs,t):
    starts=[a for a,_,_,_ in segs]
    j=max(0,min(len(segs)-1,bisect_right(starts,t)-1))
    a,b,u,v=segs[j];return u+(v-u)*(t-a)/(b-a)

def baseline_at(base,t):
    for (a,u),(b,v) in zip(base,base[1:]):
        if a<=t<=b:return u+(v-u)*(t-a)/(b-a)
    raise ValueError('outside')

def verify(p,owners,bounds,reported_area=None):
    if len(bounds)!=len(p.baselines) or len(owners)!=len(p.source):raise ValueError('incorrect dimensions')
    if any(type(c) is not int or not 0<=c<=len(p.baselines) for c in owners):raise ValueError('invalid owner label')
    allknots={t for base in p.baselines for t,_ in base};area=F(0)
    for base,segs,band,slope in zip(p.baselines,bounds,p.bands,p.slopes):
        if not segs or segs[0][0]!=base[0][0] or segs[-1][1]!=base[-1][0]:raise ValueError('incorrect domain')
        if segs[0][2]!=base[0][1] or segs[-1][3]!=base[-1][1]:raise ValueError('endpoint change')
        for j,(a,b,u,v) in enumerate(segs):
            if a>=b or abs(v-u)>slope*(b-a):raise ValueError('slope/domain error')
            if min(u,v)<band[0] or max(u,v)>band[1]:raise ValueError('band violation')
            if j and (segs[j-1][1]!=a or segs[j-1][3]!=u):raise ValueError('discontinuity')
            allknots.update((a,b))
        ts=sorted({t for t,_ in base}|{x for a,b,_,_ in segs for x in (a,b)})
        for a,b in zip(ts,ts[1:]):
            d1=at_segments(segs,a)-baseline_at(base,a);d2=at_segments(segs,b)-baseline_at(base,b)
            if d1*d2>=0:area+=(abs(d1)+abs(d2))*(b-a)/2
            else:area+=(d1*d1+d2*d2)*(b-a)/(2*(abs(d1)+abs(d2)))
    for t in allknots:
        vals=[at_segments(segs,t) for segs in bounds]
        for k,d in enumerate(p.minimum_gaps):
            if vals[k+1]-vals[k]<d:raise ValueError('minimum thickness violation')
            if p.maximum_gaps[k] is not None and vals[k+1]-vals[k]>p.maximum_gaps[k]:raise ValueError('maximum thickness violation')
    for box in p.boxes:
        if p.mode=='universal':lo,hi,lower,upper=box.tlo,box.thi,box.xhi+p.clearance,box.xlo-p.clearance
        else:
            lo=hi=(box.tlo+box.thi)/2
            if p.mode=='nominal':lower=(box.xlo+box.xhi)/2+p.clearance;upper=(box.xlo+box.xhi)/2-p.clearance
            else:
                dy=(box.thi-box.tlo)/2
                lower=box.xlo+p.clearance-p.slopes[0]*dy;upper=box.xhi-p.clearance+p.slopes[0]*dy
        ts=sorted({lo,hi}|{t for t in allknots if lo<t<hi})
        for k,segs in enumerate(bounds):
            if p.mode=="possible_outer":
                dy=(box.thi-box.tlo)/2
                lower=box.xlo+p.clearance-p.slopes[k]*dy;upper=box.xhi-p.clearance+p.slopes[k]*dy
            vals=[at_segments(segs,t) for t in ts]
            if owners[box.group]<=k:
                if min(vals)<lower:raise ValueError('RNA lower ownership violation')
            elif max(vals)>upper:raise ValueError('RNA upper ownership violation')
    if reported_area is not None and area!=reported_area:raise ValueError('incorrect physical cost')
    return area


def changed_tissue_area(p,bounds):
    """Actual area whose ordered cell index changes, counted once.

    Unlike summed swept interface area this does not count multi-interface
    crossings repeatedly. Exact analytic integration; not used as cut objective.
    """
    old=[[(a,b,u,v) for (a,u),(b,v) in zip(base,base[1:])] for base in p.baselines]
    fs=old+bounds
    ts=sorted({x for segs in fs for a,b,_,_ in segs for x in (a,b)})
    out=F(0)
    for a,b in zip(ts,ts[1:]):
        lines=[]
        for segs in fs:
            u=at_segments(segs,a);v=at_segments(segs,b);m=(v-u)/(b-a);lines.append((m,u-m*a))
        extra={a,b}
        for i,(m,c) in enumerate(lines):
            for mm,cc in lines[:i]:
                if m!=mm:
                    t=(cc-c)/(m-mm)
                    if a<t<b:extra.add(t)
        extra=sorted(extra)
        for aa,bb in zip(extra,extra[1:]):
            t=(aa+bb)/2
            ids=sorted(range(len(lines)),key=lambda i:lines[i][0]*t+lines[i][1])
            for i,j in zip(ids,ids[1:]):
                x=((lines[i][0]+lines[j][0])*t+lines[i][1]+lines[j][1])/2
                oldcell=sum(lines[k][0]*t+lines[k][1]<x for k in range(p.m))
                newcell=sum(lines[k+p.m][0]*t+lines[k+p.m][1]<x for k in range(p.m))
                if oldcell!=newcell:
                    m=lines[j][0]-lines[i][0];c=lines[j][1]-lines[i][1]
                    out+=(m*(aa+bb)+2*c)*(bb-aa)/2
    return out
