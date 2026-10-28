"""A split-conformal budget utility, not a fitted GeoCert uncertainty model.

Scores must be constructed independently on exchangeable calibration units after
all geometry-fitting choices are frozen. None denotes an unbounded score. A
nonfinite calibrated threshold returns abstention, never a false finite guarantee.
"""
from fractions import Fraction as F
from math import ceil


def split_budget(scores, alpha=F(1, 20), calibration_ids=None, evaluation_ids=None):
    alpha = F(alpha)
    if not 0 < alpha < 1 or not scores:
        raise ValueError('Nonempty calibration set and 0 < alpha < 1 required')
    if calibration_ids is not None:
        if len(calibration_ids) != len(scores) or len(set(calibration_ids)) != len(scores):
            raise ValueError('Each independent calibration unit needs one unique identifier')
        if set(calibration_ids) & set(evaluation_ids or []):
            raise ValueError('Calibration and evaluation units overlap')
    values = [None if s is None else F(s) for s in scores]
    if any(s is not None and s < 0 for s in values):
        raise ValueError('Nonnegative nonconformity scores required')
    n = len(values)
    rank = ceil((n + 1) * (1 - alpha))
    ordered = sorted((s for s in values if s is not None))
    if rank > len(ordered):
        return {'status': 'abstain', 'rank': rank, 'units': n, 'alpha': alpha, 'budget': None}
    return {'status': 'finite', 'rank': rank, 'units': n, 'alpha': alpha, 'budget': ordered[rank - 1],
            'scope': 'Marginal new-unit coverage requires exchangeability and a valid nested-family score; not cellwise conditional coverage.'}
