"""One common solving/verification wrapper for the SCC ablation."""
from core_support import *
from critical_mesh import compile_mesh,graph_mesh
from pinned_compiler import graph_for as graph_direct
from layered import value
from hard_scc import build_plan,contract,lift_side
from time import perf_counter


def prepare(p,pins,representation,use_scc=False,seconds=30):
    start=perf_counter()
    c=compile_mesh(p,pins,seconds=seconds,node_cap=150000) if representation=='height' else compile_pinned(p,pins)
    compiled=perf_counter();plan=None
    if use_scc:
        graph=(graph_mesh if representation=='height' else graph_direct)(c,[[0]*(p.m+1) for _ in range(p.n)],1)
        plan=build_plan(graph)
    return c,plan,{'compile_seconds':compiled-start,'scc_prepare_seconds':perf_counter()-compiled,'prepare_seconds':perf_counter()-start}


def solve_compiled(c,representation,costs,lam,plan=None,seconds=30,certificate=False):
    start=perf_counter();original=(graph_mesh if representation=='height' else graph_direct)(c,costs,lam)
    g=contract(original,plan) if plan is not None else original;t1=perf_counter()
    d=Dinic(g['n'])
    for u,v,w in g['edges']:d.add(u,v,w)
    opt,side,records=d.run(g['s'],g['t'],start+seconds);t2=perf_counter()
    cert={'n':g['n'],'s':g['s'],'t':g['t'],'arcs':records,'source_side':sorted(side),'value':int(opt),'scale':g['scale'],'constant':str(g['constant'])}
    check_certificate(cert);t3=perf_counter()
    lifted=lift_side(side,plan) if plan is not None else side
    # Establish exact cut-cost preservation, not only objective agreement.
    old_cut=sum(w for u,v,w in original['edges'] if u in lifted and v not in lifted)
    if old_cut!=opt:raise AssertionError('lifted cut value differs')
    ans={'feasible':opt<g['big'],'nodes':g['n'],'arcs':len(g['edges']),
         'original_nodes':original['n'],'original_arcs':len(original['edges']),
         'network_seconds':t1-start,'flow_seconds':t2-t1,'flow_check_seconds':t3-t2}
    if certificate:ans['certificate']=cert;ans['mapping']=None if plan is None else list(plan.mapping)
    if not ans['feasible']:ans['seconds']=perf_counter()-start;return ans
    p=c.problem;owners=p.decode([int(i in lifted) for i in range(p.n*p.m)])
    if representation=='direct':
        area=c.area(owners);w=realize(p,owners)
        if not w['feasible'] or w['area']!=area:raise AssertionError('direct realization disagreement')
        boundaries=w['boundaries']
    elif c.fixed_answer is not None:
        if owners!=c.fixed_answer['owners']:raise AssertionError('fixed owner disagreement')
        area=c.fixed_answer['area'];boundaries=c.fixed_answer['boundaries']
    else:
        heights=[[max(x for x,v in zip(c.candidates[k][j],c.ids[k][j]) if v==original['s'] or v in lifted) for j in range(len(c.mesh))] for k in range(p.m)]
        boundaries=[[(a,b,x,y) for a,b,x,y in zip(c.mesh,c.mesh[1:],h,h[1:])] for h in heights]
        area=F(0);trap=F(0)
        for k,h in enumerate(heights):
            for a,b,x,y in zip(c.mesh,c.mesh[1:],h,h[1:]):
                u=x-value(p.baselines[k],a);v=y-value(p.baselines[k],b)
                trap+=(b-a)*(abs(u)+abs(v))/2
                area+=(b-a)*(u*u+v*v)/(2*(abs(u)+abs(v))) if u*v<0 else (b-a)*(abs(u)+abs(v))/2
        if area>trap:raise AssertionError('invalid quadrature bound')
    objective=F(lam)*area+sum(F(row[o]) for row,o in zip(costs,owners))
    if objective!=F(opt,g['scale'])+g['constant']:raise AssertionError('physical objective disagreement')
    t4=perf_counter();verify_mixed(p,owners,boundaries,area);t5=perf_counter()
    ans.update(owners=owners,area=area,objective=objective,boundaries=boundaries,
               reconstruct_seconds=t4-t3,boundary_check_seconds=t5-t4,seconds=t5-start)
    return ans
