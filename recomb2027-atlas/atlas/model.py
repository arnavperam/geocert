"""Prepared exact patch models and rational physical-witness checks."""
from fractions import Fraction as F
from copy import copy
from bisect import bisect_right
import hashlib,json
from geocert_expression.budget import ExpressionEngine
from layered import Box, Layers
from mixed import TaggedBox, MixedLayers
from verify import at_segments
from .geometry import evaluate_path


def make_patch(patch):
    """Encode chart-displacement constraints without changing the exact engine.

    The transverse coordinate is displacement from the *exact* source raster arc.
    For a pinned left nuclear bound h>=v, the synthetic box is at v-rho; for a
    pinned right bound h<=v it is at v+rho. Adding the engine's RNA clearance
    cancels this shift. Thus structural interiors have zero extra clearance.
    These boxes are constraint encodings, not transformed RNA observations.
    """
    rho=F(patch['clearance_pixels']);T=F(patch['domain_pixels'][1]);band=F(patch['band_halfwidth_pixels'])
    boxes=[];source=[];pins={}
    for i,r in enumerate(patch['points']):
        x,t=F(r['x']),F(r['t']);boxes.append(TaggedBox(x,x,t,t,i,False));source.append(r['source'])
    groups={}
    for side,value,a,b in patch.get('nuclear_constraints',[]):
        if side not in groups:
            i=len(source);groups[side]=i;source.append(side);pins[i]=side
        i=groups[side];x=F(value)+(-rho if side==0 else rho)
        boxes.append(TaggedBox(x,x,F(a),F(b),i,True))
    p=MixedLayers([[(0,0),(T,0)]],boxes,source,[F(patch['slope'])],[],[],
                  clearance=rho,bands=[(-band,band)],mode='nominal')
    engine=ExpressionEngine(p,pins)
    if engine.source_engine is None or engine.source_engine.inconsistent:
        raise ValueError('reference patch does not satisfy all source and structural owners')
    return engine


def count_scores(patch,engine,cell,gene):
    out=[[0,0] for _ in range(engine.problem.n)]
    for i,r in enumerate(patch['points']):
        if r['gene']==gene:
            out[i]=[int(k==cell) for k in patch['cell_labels']]
    return out


def serialized(obj):
    return json.loads(json.dumps(obj,default=lambda v:str(v) if isinstance(v,F) else (_ for _ in ()).throw(TypeError(type(v).__name__))))


def verify_witness(patch, owners, bounds, area_um2):
    """Checks chart geometry, physical area, RNA and true structural constraints.

    No compiler, cut optimizer, implication routine, or envelope constructor.
    Topological validity additionally relies on the separately verified atlas.
    """
    seg=[tuple(map(F,z)) for z in bounds[0]]
    T=F(patch['domain_pixels'][1]);width=F(patch['halfwidth_pixels']);band=F(patch['band_halfwidth_pixels'])
    slope=F(patch['slope']);rho=F(patch['clearance_pixels'])
    if not seg or seg[0][0]!=0 or seg[-1][1]!=T or seg[0][2]!=0 or seg[-1][3]!=0:raise ValueError('trace mismatch')
    area=F(0)
    for j,(a,b,u,v) in enumerate(seg):
        if a>=b or abs(v-u)>slope*(b-a) or min(u,v)<-band or max(u,v)>band:raise ValueError('invalid boundary geometry')
        if not -width<u<width or not -width<v<width:raise ValueError('homeomorphism loses strict rails')
        if j and seg[j-1][1]!=a:raise ValueError('gap in segments')
        if j and seg[j-1][3]!=u:raise ValueError('discontinuous segments')
        if u*v>=0:area+=(b-a)*(abs(u)+abs(v))/2
        else:area+=(b-a)*(u*u+v*v)/(2*(abs(u)+abs(v)))
    for i,r in enumerate(patch['points']):
        x,t=F(r['x']),F(r['t']);f=at_segments(seg,t)
        if owners[i] not in [0,1]:raise ValueError('bad RNA owner')
        if owners[i]==0 and f<x+rho:raise ValueError('RNA outside requested left cell')
        if owners[i]==1 and f>x-rho:raise ValueError('RNA outside requested right cell')
    for side,bound,a,b in patch.get('nuclear_constraints',[]):
        a,b,bound=F(a),F(b),F(bound)
        ts={a,b}|{t for aa,bb,_,_ in seg for t in [aa,bb] if a<t<b}
        vals=[at_segments(seg,t) for t in ts]
        if side==0 and min(vals)<bound:raise ValueError('nuclear left region displaced out of cell')
        if side==1 and max(vals)>bound:raise ValueError('nuclear right region displaced out of cell')
    physical=area*F(patch['area_scale_um2'])
    if physical!=F(area_um2):raise ValueError('physical area/Jacobian mismatch')
    return physical


def compose(parts, constant=0):
    """Exact scalar min-plus composition, with attaining patch-state indices."""
    state={constant:(F(0),{})}
    for pid,minima in parts:
        nxt={}
        for old,(cost,trace) in state.items():
            for delta,localcost in minima.items():
                value=old+int(delta);candidate=cost+F(localcost)
                if value not in nxt or candidate<nxt[value][0]:
                    nxt[value]=(candidate,{**trace,pid:int(delta)})
        state=nxt
    return state


def interval(state,budget):
    feasible=[v for v,(d,_) in state.items() if d<=F(budget)]
    if not feasible:return {'status':'empty_uncertainty_set','budget':str(budget)}
    lo,hi=min(feasible),max(feasible)
    return {'status':'exact','lower':lo,'upper':hi,'budget':str(budget),
            'lower_cost':str(state[lo][0]),'upper_cost':str(state[hi][0]),
            'lower_trace':state[lo][1],'upper_trace':state[hi][1]}
