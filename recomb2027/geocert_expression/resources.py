"""Optional simultaneous physical-deformation and all-RNA change budgets."""
from fractions import Fraction as F
from itertools import product
from time import perf_counter
from .budget import check_scores,ScoreSpectrum,floor,ceil
from .scope import digest
from pinned_compiler import solve_pinned
from mixed import verify_mixed

def _weights(engine,weights):
    if len(weights)!=engine.problem.n or any(type(w) is not int or w<0 for w in weights):
        raise ValueError('nonnegative integer reassignment weight for every group required')
    return tuple(weights)

def changed_weight(engine,owners,weights):
    return sum(w for w,o,s in zip(weights,owners,engine.problem.source) if o!=s)

def resource_spectrum(engine,scores,weights,rna_budget,constant=0,pattern_limit=100000):
    """Exact target projection at a second hard RNA-change budget.

    Source clipping minimizes D and source-owner change weight simultaneously
    within each target pin assignment. Across target assignments it still keeps
    the minimum D separately at the supplied RNA budget.
    """
    rows=check_scores(engine.problem,scores);weights=_weights(engine,weights)
    if engine.source_engine is None:raise ValueError('zero-cost source required')
    if type(rna_budget) is not int or rna_budget<0 or type(constant) is not int or type(pattern_limit) is not int or pattern_limit<1:raise ValueError('invalid RNA budget, integer offset or pattern limit')
    t0=perf_counter();targets=[g for g,row in enumerate(rows) if len(set(row))>1]
    domains=[[o for o in range(engine.problem.m+1) if engine.source_engine.closure({g:o})['feasible']] for g in targets]
    total=1
    for d in domains:total*=len(d)
    costs={};witnesses={};examined=0
    for oo in product(*domains):
        if examined>=pattern_limit:break
        examined+=1;pins=dict(zip(targets,oo));a=engine.clipped(pins)
        if not a['feasible']:continue
        r=changed_weight(engine,a['owners'],weights)
        if r>rna_budget:continue
        v=constant+sum(row[o] for row,o in zip(rows,a['owners']));d=F(a['area'])
        if v not in costs or d<costs[v]:
            costs[v]=d;witnesses[v]={'pins':pins,'owners':a['owners'],'area':str(d),'rna_changes':r,'boundaries':a['boundaries']}
    return ScoreSpectrum(costs,witnesses,examined==total,examined,total,perf_counter()-t0,
                         {'model_key': engine.model_key, 'score_key': engine.score_key(rows),
                          'resource_key': digest([weights, rna_budget]), 'constant': constant,
                          'scope': 'supplied_geometric_model_and_resource_budget_only'})

def resource_bound(engine,scores,weights,delta,rna_budget,mu,nu,constant=0,sign=1):
    """One exact rational Lagrange endpoint bound for two resource constraints."""
    rows=check_scores(engine.problem,scores);weights=_weights(engine,weights)
    delta,mu,nu=F(delta),F(mu),F(nu)
    if min(delta,mu,nu)<0 or type(rna_budget) is not int or rna_budget<0 or type(constant) is not int or sign not in (-1,1):raise ValueError('invalid resources')
    costs=[[nu*w*(o!=s)-sign*v for o,v in enumerate(row)] for row,w,s in zip(rows,weights,engine.problem.source)]
    binding={'model_key':engine.model_key,'score_key':engine.score_key(rows),
             'resource_key':digest([weights,rna_budget]),'constant':constant,
             'scope':'supplied_geometric_model_and_resource_budget_only'}
    a=solve_pinned(engine.compiled,costs,area_weight=mu,certificate=True,reconstruct=True)
    if not a['feasible']:return {'status':'model_infeasible','certificate':a.get('certificate'),'binding':binding}
    verify_mixed(engine.problem,a['owners'],a['boundaries'],a['area'])
    raw=mu*delta+nu*rna_budget-F(a['objective'])
    bound=constant+floor(raw) if sign==1 else constant+ceil(-raw)
    return {'status':'bound','binding':binding,'direction':'upper' if sign==1 else 'lower','value':bound,
            'mu':str(mu),'nu':str(nu),'rna_budget':rna_budget,'budget_um2':str(delta),
            'witness_score':constant+sum(row[o] for row,o in zip(rows,a['owners'])),
            'witness_within_budgets':a['area']<=delta and changed_weight(engine,a['owners'],weights)<=rna_budget,
            'solution':a}
