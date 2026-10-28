"""Independent rational scalar-envelope audit of the baseline native showcase.

Only primitive input extraction uses the supplied model ingestion and obligation
semantics. Feasibility, envelope integration, full-owner enumeration, singleton
costs, and objective/penalty comparisons below use Python's standard library.
No production feasibility kernels, realization, area compiler, flow, or saved
singleton exclusions determine the independently enumerated candidate set.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
from bisect import bisect_right
import sys,json,hashlib

source=Path(sys.argv[1]).resolve();output=Path(sys.argv[2])
sys.path.insert(0,str(source/'core'))
from data_api import load_cases
ix,ch,tag,p,pins,n,prior=next(z for z in load_cases() if z[1]['id']=='4037abe9fea38e' and z[2]=='baseline')
assert p.m==1
primitive={'chart':ch['id'],'source_charts_sha256':hashlib.sha256((source/'data/charts.json').read_bytes()).hexdigest(),
           'slope':p.slopes[0],'band':p.bands[0],'base':p.baselines[0],
           'source_owners':p.source,'structural_pins':pins,'n_rna':n,'n_groups':p.n,
           'obligations':[{'group':z.group,'interval_lower_upper':p.obligation(z,0)} for z in p.boxes],
           'genes':[z['gene'] for z in ch['points']]}

def distance(t,lo,hi):return max(lo-t,F(0),t-hi)
base=primitive['base'];slope=primitive['slope'];band=primitive['band']
tmin,tmax=base[0][0],base[-1][0]
base_times=[t for t,x in base]
knots=sorted({t for t,x in base}|{v for o in primitive['obligations'] for v in o['interval_lower_upper'][:2]})
def source_line(t):
    j=min(len(base)-2,max(0,bisect_right(base_times,t)-1))
    (a,x),(b,y)=base[j:j+2];m=(y-x)/(b-a)
    return m,x-m*a
def line(cone,a,b):
    lo,hi,height,polarity=cone
    va=height+polarity*slope*distance(a,lo,hi);vb=height+polarity*slope*distance(b,lo,hi)
    m=(vb-va)/(b-a);return m,va-m*a
def val(line,t):return line[0]*t+line[1]

def cones_for(owners):
    lower=[(t,t,x,-1) for t,x in (base[0],base[-1])]
    upper=[(t,t,x,1) for t,x in (base[0],base[-1])]
    for o in primitive['obligations']:
        if o['group'] not in owners:continue
        lo,hi,L,U=o['interval_lower_upper']
        if owners[o['group']]==0:lower.append((lo,hi,L,-1))
        else:upper.append((lo,hi,U,1))
    return lower,upper

def evaluate(owners):
    lower,upper=cones_for(owners);segments=[];area=F(0)
    for a,b in zip(knots,knots[1:]):
        ls=[(F(0),band[0])]+[line(c,a,b) for c in lower]
        us=[(F(0),band[1])]+[line(c,a,b) for c in upper]
        # Same-slope domination is exact throughout this affine interval.
        lm={};um={}
        for m,z in ls:lm[m]=max(lm.get(m,z),z)
        for m,z in us:um[m]=min(um.get(m,z),z)
        ls=list(lm.items());us=list(um.items());src=source_line((a+b)/2)
        lines=ls+us+[src];cuts={a,b}
        for (m,z),(mm,zz) in combinations(lines,2):
            if m!=mm:
                t=(zz-z)/(m-mm)
                if a<t<b:cuts.add(t)
        cuts=sorted(cuts)
        for l,r in zip(cuts,cuts[1:]):
            mid=(l+r)/2;L=max(ls,key=lambda z:val(z,mid));U=min(us,key=lambda z:val(z,mid))
            if val(L,mid)>val(U,mid) or val(L,l)>val(U,l) or val(L,r)>val(U,r):return None
            f=L if val(src,mid)<val(L,mid) else U if val(src,mid)>val(U,mid) else src
            dl=val(f,l)-val(src,l);dr=val(f,r)-val(src,r)
            assert dl*dr>=0
            area+=(abs(dl)+abs(dr))*(r-l)/2
            if segments and segments[-1][1]==l and segments[-1][2:]==f:
                segments[-1]=(segments[-1][0],r,*f)
            else:segments.append((l,r,*f))
    return {'D':area,'segments':segments}

# Obtain safe candidate domains independently: retain a label unless it conflicts
# even with the relaxed model containing just the structural and that RNA's pins.
domains=[];excluded=[]
for g in range(primitive['n_groups']):
    if g in pins:domains.append([pins[g]]);continue
    choices=[]
    for owner in [0,1]:
        if evaluate({**pins,g:owner}) is not None:choices.append(owner)
        else:excluded.append([g,owner])
    domains.append(choices)
rows=[];total=0
for owners in product(*domains):
    total+=1;answer=evaluate(dict(enumerate(owners)))
    if answer is None:continue
    R=sum(owners[g]!=primitive['source_owners'][g] for g in range(n))
    rows.append({'owners':list(owners),'R':R,**answer})
assert rows
single={}
for g in range(n):
    for o in [0,1]:
        costs=[r['D'] for r in rows if r['owners'][g]==o]
        single[g,o]=min(costs) if costs else None
for r in rows:
    r['D_add']=sum(single[g,r['owners'][g]] for g in range(n))
    assert r['D']<=r['D_add']

saved=json.loads((source/'data/showcase_scalar_audit.json').read_text())
saved_map={tuple(r['owners']):(F(r['D']),r['R'],F(r['D_add'])) for r in saved['rows']}
fresh_map={tuple(r['owners']):(r['D'],r['R'],r['D_add']) for r in rows}
assert saved_map==fresh_map

def optimality_interval(winner,cost):
    lo=F(0);hi=None
    for r in rows:
        delta=winner[cost]-r[cost];gain=winner['R']-r['R']
        if delta>0:
            bound=gain/delta;hi=bound if hi is None else min(hi,bound)
        elif delta<0:lo=max(lo,gain/delta)
        else:assert gain>=0
    assert hi is None or lo<=hi
    return lo,hi

queries=[]
for lam in [F(1,4),F(1),F(4)]:
    best=min(lam*r['D']-r['R'] for r in rows)
    exact=[r for r in rows if lam*r['D']-r['R']==best]
    apr=min(lam*r['D_add']-r['R'] for r in rows)
    ar=[r for r in rows if lam*r['D_add']-r['R']==apr]
    additive=min(ar,key=lambda r:lam*r['D']-r['R'])
    q={'lambda':lam,'exact_objective':best,'additive_primary':apr,
       'best_additive_true_objective':lam*additive['D']-additive['R'],
       'exact_ties':len(exact),'additive_primary_ties':len(ar),
       'exact_owners':exact[0]['owners'],'additive_owners':additive['owners']}
    if lam==1:
        ew=exact[0];aw=additive
        er=optimality_interval(ew,'D');ar=optimality_interval(aw,'D_add')
        common_lo=max(er[0],ar[0]);common_hi=min(x for x in [er[1],ar[1]] if x is not None)
        q.update(exact_optimality_interval=er,additive_optimality_interval=ar,
                 same_two_winners_interval=[common_lo,common_hi],
                 regret_affine_intercept=ew['R']-aw['R'],regret_affine_slope=aw['D']-ew['D'],
                 interior_min_regret=min(aw['D']*t-aw['R']-(ew['D']*t-ew['R']) for t in [common_lo,common_hi]))
    queries.append(q)

# Build a COMMON, ownership-independent event mesh from every possible cone and
# source segment, not from a selected witness. At these events the canonical
# optimum's trapezoidal cost is exact for every feasible owner vector.
allcones=[(t,t,x,pol) for t,x in (base[0],base[-1]) for pol in [-1,1]]
for o in primitive['obligations']:
    lo,hi,L,U=o['interval_lower_upper'];allcones.extend([(lo,hi,L,-1),(lo,hi,U,1)])
critical=set(knots)
for a,b in zip(knots,knots[1:]):
    lines=list(set([line(c,a,b) for c in allcones]+[(F(0),band[0]),(F(0),band[1]),source_line((a+b)/2)]))
    for (m,z),(mm,zz) in combinations(lines,2):
        if m!=mm:
            t=(zz-z)/(m-mm)
            if a<t<b:critical.add(t)
critical=sorted(critical)
def trapezoidal(row,mesh):
    segs=row['segments'];starts=[s[0] for s in segs]
    def displacement(t):
        j=min(len(segs)-1,max(0,bisect_right(starts,t)-1));l,r,m,b=segs[j]
        return abs(m*t+b-val(source_line(t),t))
    v=[displacement(t) for t in mesh]
    return sum((u+w)*(b-a)/2 for a,b,u,w in zip(mesh,mesh[1:],v,v[1:]))
coarse_gaps=[]
for r in rows:
    refined=trapezoidal(r,critical);coarse=trapezoidal(r,knots)
    assert refined==r['D'] and coarse>=r['D']
    coarse_gaps.append(coarse-r['D'])
result={'scope':'Independent scalar-envelope arithmetic and all candidate owner enumeration; shared primitive input ingestion only. One two-cell native baseline example, not a general multi-boundary implementation.',
        'primitive':primitive,'candidate_vectors':total,'feasible_vectors':len(rows),'excluded_by_relaxed_single_pin_geometry':excluded,
        'all_saved_owner_costs_match':True,'queries':queries,
        'critical_mesh':{'coarse_nodes':len(knots),'event_nodes':len(critical),'all_26_costs_exact_on_common_event_mesh':True,
                         'coarse_mesh_overestimates':sum(v>0 for v in coarse_gaps),'max_coarse_excess':max(coarse_gaps)},
        'rows':rows}
output.write_text(json.dumps(result,indent=2,default=lambda v:str(v) if isinstance(v,F) else None))
print(json.dumps({k:v for k,v in result.items() if k not in ['rows','primitive']},indent=2,default=str))
