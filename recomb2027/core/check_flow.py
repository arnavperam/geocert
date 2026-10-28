"""Independent integer primal-dual max-flow/min-cut certificate checker.

Does not import compiler or flow solver. This certifies the supplied finite graph;
the continuum equivalence is a separate mathematical/compiler obligation.
"""
def check_certificate(cert):
    n=cert['n'];s=cert['s'];t=cert['t'];S=set(cert['source_side'])
    if not(0<=s<n and 0<=t<n and s!=t and s in S and t not in S):raise ValueError('bad terminals/cut')
    bal=[0]*n;cut=0
    for u,v,cap,flow in cert['arcs']:
        if not all(isinstance(q,int) for q in (u,v,cap,flow)) or not(0<=u<n and 0<=v<n and 0<=flow<=cap):raise ValueError('bad arc capacity or flow')
        bal[u]-=flow;bal[v]+=flow
        if u in S and v not in S:cut+=cap
    if any(bal[v] for v in range(n) if v not in (s,t)):raise ValueError('flow conservation')
    if not -bal[s]==bal[t]==cut==cert['value']:raise ValueError('primal-dual gap')
    return True
