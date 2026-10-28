"""Picard--Queyranne residual-closure representation of all penalized optima.

This is established minimum-cut theory, not a new GeoCert theorem. Projection
onto ownership bits can identify several cuts with the same assignment.
"""
from itertools import product
from check_flow import check_certificate
from hard_scc import build_plan


def reachable(starts, adjacency):
    seen = set(starts)
    stack = list(seen)
    while stack:
        for v in adjacency[stack.pop()]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def all_min_cuts(certificate, owner_bits, original_to_flow=None):
    check_certificate(certificate)
    n, s, t = (certificate[k] for k in ('n', 's', 't'))
    residual = set()
    for u, v, cap, flow in certificate['arcs']:
        if cap > flow:
            residual.add((u, v))
        if flow > 0:
            residual.add((v, u))
    # All residual edges are closure implications; no cost arc is discarded.
    graph = {'n': n, 's': s, 't': t, 'big': 1,
             'edges': [(u, v, 1) for u, v in sorted(residual)]}
    plan = build_plan(graph)
    if plan.terminal_conflict:
        raise AssertionError("A verified maximum flow cannot have an augmenting residual path")
    component = plan.mapping
    edges = sorted({(component[u], component[v]) for u, v in residual if component[u] != component[v]})
    adj = [[] for _ in range(plan.new_n)]
    rev = [[] for _ in range(plan.new_n)]
    for u, v in edges:
        adj[u].append(v)
        rev[v].append(u)
    forced_one = reachable([component[s]], adj)
    forced_zero = reachable([component[t]], rev)
    assert forced_one.isdisjoint(forced_zero)
    mapping = list(range(n)) if original_to_flow is None else original_to_flow
    return {'components': plan.new_n, 'edges': edges,
            'forced_source': sorted(forced_one), 'forced_sink': sorted(forced_zero),
            'flow_vertex_component': list(component),
            'owner_bit_component': [component[mapping[i]] for i in range(owner_bits)],
            'interpretation': 'Successor-closed component sets containing forced_source and excluding forced_sink represent exactly all minimum cuts.'}


def owner_options(problem, representation):
    """Possible owners at each group; choices across groups remain dependent."""
    r = representation
    adj = [[] for _ in range(r['components'])]
    for u, v in r['edges']:
        adj[u].append(v)
    result = []
    for g in range(problem.n):
        allowed = []
        for owner in range(problem.m + 1):
            one, zero = set(r['forced_source']), set(r['forced_sink'])
            for a in range(problem.m):
                v = r['owner_bit_component'][problem.bit(g, a)]
                (one if owner <= a else zero).add(v)
            if not (reachable(one, adj) & zero):
                allowed.append(owner)
        result.append(allowed)
    return result


def enumerate_optimal_owners(problem, representation, max_free_components=20):
    r = representation
    one, zero = set(r['forced_source']), set(r['forced_sink'])
    free = sorted(set(range(r['components'])) - one - zero)
    if len(free) > max_free_components:
        raise ValueError('All-optima enumeration cap exceeded; compact DAG is still complete')
    owners = set()
    for values in product([0, 1], repeat=len(free)):
        side = one | {v for v, z in zip(free, values) if z}
        if any(u in side and v not in side for u, v in r['edges']):
            continue
        bits = [int(c in side) for c in r['owner_bit_component']]
        owners.add(tuple(problem.decode(bits)))
    return sorted(owners)
