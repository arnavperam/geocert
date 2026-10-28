"""Joint count outcomes and residual descriptions of penalized optima.

Marginal interval endpoints and marginal optimal owners cannot be multiplied into
an attainable joint outcome. This module retains the shared assignment explicitly.
"""
from itertools import product
from fractions import Fraction as F
from .budget import check_scores
from layered import optimal_owner_domains
from check_flow import check_certificate

def outcome_spectrum(engine, channels, constants=None, pattern_limit=100000):
    channels=[check_scores(engine.problem,s) for s in channels]
    constants=list(constants or [0]*len(channels))
    if len(constants)!=len(channels) or any(type(c) is not int for c in constants):raise ValueError('integer offset per channel required')
    targets=sorted({g for s in channels for g,row in enumerate(s) if len(set(row))>1})
    if engine.source_engine is None:raise ValueError('zero-cost source required for joint target projection')
    domains=[[o for o in range(engine.problem.m+1) if engine.source_engine.closure({g:o})['feasible']] for g in targets]
    total=1
    for d in domains:total*=len(d)
    values={};examined=0
    for owners in product(*domains):
        if examined>=pattern_limit:break
        examined+=1;pins=dict(zip(targets,owners));a=engine.clipped(pins)
        if not a['feasible']:continue
        key=tuple(c+sum(row[o] for row,o in zip(s,a['owners'])) for c,s in zip(constants,channels))
        if key not in values or F(a['area'])<F(values[key]['area']):
            values[key]={'area':str(a['area']),'owners':a['owners'],'boundaries':a['boundaries'],'pins':pins}
    return {'complete':examined==total,'examined':examined,'patterns':total,
            'outcomes':[{'counts':list(k),**v} for k,v in sorted(values.items())]}

def predicate_status(spectrum,budget,predicate):
    feasible=[r for r in spectrum['outcomes'] if F(r['area'])<=F(budget)]
    truth=[bool(predicate(r['counts'])) for r in feasible]
    if not spectrum['complete']:
        return {'status':'unresolved','true_witness':next((r for r,t in zip(feasible,truth) if t),None),
                'false_witness':next((r for r,t in zip(feasible,truth) if not t),None)}
    if not feasible:return {'status':'empty_uncertainty_set'}
    return {'status':'always' if all(truth) else 'never' if not any(truth) else 'variable',
            'true_witness':next((r for r,t in zip(feasible,truth) if t),None),
            'false_witness':next((r for r,t in zip(feasible,truth) if not t),None)}

def all_optimal_assignments(compiled,solution):
    """Residual SCC DAG represents all mincuts; owner domains are only marginals."""
    cert=solution.get('certificate')
    if not solution.get('feasible') or cert is None:raise ValueError('feasible flow certificate required')
    check_certificate(cert);n=cert['n'];adj=[set() for _ in range(n)];rev=[set() for _ in range(n)]
    for u,v,cap,flow in cert['arcs']:
        if cap>flow:adj[u].add(v);rev[v].add(u)
        if flow:adj[v].add(u);rev[u].add(v)
    seen=set();order=[]
    for v in range(n):
        if v in seen:continue
        stack=[(v,False)]
        while stack:
            u,done=stack.pop()
            if done:order.append(u);continue
            if u in seen:continue
            seen.add(u);stack.append((u,True));stack.extend((x,False) for x in sorted(adj[u],reverse=True) if x not in seen)
    component=[-1]*n;groups=[]
    for v in reversed(order):
        if component[v]>=0:continue
        cid=len(groups);stack=[v];gg=[];component[v]=cid
        while stack:
            u=stack.pop();gg.append(u)
            for x in rev[u]:
                if component[x]<0:component[x]=cid;stack.append(x)
        groups.append(sorted(gg))
    arcs=sorted({(component[u],component[v]) for u in range(n) for v in adj[u] if component[u]!=component[v]})
    return {'scope':'all optima of this penalized objective, not all hard-budget feasible assignments',
            'components':groups,'implications':arcs,'source_component':component[cert['s']],
            'sink_component':component[cert['t']],'owner_bit_components':component[:compiled.problem.n*compiled.problem.m],
            'marginal_owner_domains':optimal_owner_domains(compiled,solution),
            'semantics':'choose a forward-closed component set containing source and excluding sink; threshold bits decode owners'}
