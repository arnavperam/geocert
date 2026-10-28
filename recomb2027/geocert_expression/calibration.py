"""Split-conformal calibration at the independent region/slide unit.

This supplies a coverage theorem only under exchangeability and a family fitted
without calibration/test labels. Scores may be infinite: model failures must not
be discarded. Calibration on simulation is not calibration on biological tissue.
"""
from fractions import Fraction as F
from math import ceil, sqrt

def threshold(scores, alpha=F(1,10)):
    alpha=F(alpha)
    if not 0<alpha<1 or not scores:raise ValueError('nonempty calibration and 0<alpha<1 required')
    n=len(scores);rank=ceil((n+1)*(1-alpha))
    ordered=sorted((None if x is None else F(x) for x in scores),key=lambda x:(x is None,x or F(0)))
    if any(x is not None and x<0 for x in ordered):raise ValueError('area nonconformity scores cannot be negative')
    value=None if rank>n else ordered[rank-1]
    return {'calibration_units':n,'rank':rank,'alpha':str(alpha),'threshold':None if value is None else str(value),
            'infinite_scores':sum(x is None for x in ordered),
            'status':'unbounded' if value is None else 'finite',
            'guarantee':'marginal simultaneous count coverage at an exchangeable new region, conditional on the fitted family and valid optimization bounds; not per-donor conditional coverage'}

def wilson(successes,total,z=1.959963984540054):
    if not 0<=successes<=total or total<1:raise ValueError('invalid binomial counts')
    p=successes/total;den=1+z*z/total
    mid=(p+z*z/(2*total))/den
    half=z*sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return [mid-half,mid+half]
