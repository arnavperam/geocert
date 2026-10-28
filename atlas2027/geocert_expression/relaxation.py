"""Sparse, topology-independent outer intervals from possible cell labels.

A candidate set must contain every possible true owner under the caller's model.
The routine then drops spatial coherence. Its endpoints are exact for independent
unit-RNA assignments, but only outer bounds for a coherent segmentation family.
"""
from collections import defaultdict

def candidate_count_intervals(records, reassignment_budget=None):
    counts=defaultdict(lambda:[0,0,0])
    n=0
    for rec in records:
        source=rec['source'];gene=rec['gene'];allowed=set(rec['candidates'])
        if source not in allowed:raise ValueError('source must be a candidate')
        n+=1;counts[(source,gene)][0]+=1
        if len(allowed)>1:counts[(source,gene)][1]+=1
        for cell in allowed-{source}:counts[(cell,gene)][2]+=1
    k=n if reassignment_budget is None else reassignment_budget
    if type(k) is not int or k<0:raise ValueError('nonnegative integer RNA budget required')
    return {(cell,gene):{'source_count':v[0],'lower':v[0]-min(k,v[1]),'upper':v[0]+min(k,v[2]),
                        'removable':v[1],'potential_incoming':v[2],
                        'status':'outer_bound_for_geometry'} for (cell,gene),v in counts.items()}


def robust_linear_contrast(intervals, weights, offset=0):
    """Sound interval for a fixed linear contrast, without statistical DE claims."""
    lower=upper=offset
    for key,w in weights.items():
        lo,hi=intervals[key]['lower'],intervals[key]['upper']
        lower+=w*(lo if w>=0 else hi);upper+=w*(hi if w>=0 else lo)
    return {'lower':lower,'upper':upper,'sign':'positive' if lower>0 else 'negative' if upper<0 else 'unresolved'}


def expression_presence(intervals, keys):
    """Safe conjunction status. Separate positive upper bounds do not prove coexistence."""
    if all(intervals[k]['lower']>=1 for k in keys):return 'present_under_every_allowed_assignment'
    if any(intervals[k]['upper']==0 for k in keys):return 'joint_presence_impossible'
    return 'joint_presence_unresolved'
