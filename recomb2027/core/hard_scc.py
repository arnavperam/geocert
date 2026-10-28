"""Cost-preserving contraction of hard implication SCCs.

An independent, standard graph simplification, applied identically to both
representations. Hard arcs are explicitly supplied by their big-capacity value;
this value exceeds the sum of all finite capacities in the original graph.
No finite-cost residual arcs enter the SCC calculation. No transitive reduction,
forced-bit removal, or query-dependent pin relaxation is performed here.
"""
from dataclasses import dataclass

@dataclass(frozen=True)
class SCCPlan:
    old_n: int
    mapping: tuple[int,...]
    new_n: int
    hard_edges: tuple[tuple[int,int],...]
    terminal_conflict: bool=False


def hard_edges(g):
    # Every finite capacity is strictly smaller than big by its construction.
    return tuple(sorted(set((u,v) for u,v,c in g['edges'] if c==g['big'] and u!=v)))


def build_plan(g):
    n=g['n'];hard=hard_edges(g);adj=[[] for _ in range(n)];rev=[[] for _ in range(n)]
    for u,v in hard:adj[u].append(v);rev[v].append(u)
    seen=[False]*n;order=[]
    for start in range(n):
        if seen[start]:continue
        seen[start]=True;stack=[(start,0)]
        while stack:
            u,k=stack[-1]
            if k<len(adj[u]):
                v=adj[u][k];stack[-1]=(u,k+1)
                if not seen[v]:seen[v]=True;stack.append((v,0))
            else:order.append(u);stack.pop()
    mapping=[-1]*n;cid=0
    for start in reversed(order):
        if mapping[start]>=0:continue
        mapping[start]=cid;stack=[start]
        while stack:
            u=stack.pop()
            for v in rev[u]:
                if mapping[v]<0:mapping[v]=cid;stack.append(v)
        cid+=1
    conflict=mapping[g['s']]==mapping[g['t']]
    # Retain the original infeasible network when its terminals would coincide,
    # so the ordinary flow proof still reports infeasibility in that case.
    if conflict:return SCCPlan(n,tuple(range(n)),n,hard,True)
    return SCCPlan(n,tuple(mapping),cid,hard)


def contract(g,plan):
    if g['n']!=plan.old_n or hard_edges(g)!=plan.hard_edges:
        raise ValueError('SCC plan does not match hard constraints; recompile after new pins')
    edge={};p=plan.mapping
    for u,v,w in g['edges']:
        a,b=p[u],p[v]
        if a!=b and w:edge[a,b]=edge.get((a,b),0)+w
    q=dict(g);q.update(n=plan.new_n,s=p[g['s']],t=p[g['t']],edges=[(u,v,w) for (u,v),w in sorted(edge.items())])
    return q


def lift_side(side,plan):
    return {i for i,c in enumerate(plan.mapping) if c in side}
