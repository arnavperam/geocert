"""Hard-budget count bounds and exact target-projected cost spectra.

All arithmetic in certificates is rational/integer. A penalized optimum is only
an outer bound for a hard budget. Optional exact refinement enumerates scored
owners, not unscored collateral owners, and requires a zero-cost admissible source.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from time import perf_counter
from typing import Sequence
from copy import deepcopy
import hashlib, json
from . import ROOT
from pin_closure import PinClosure, SourceNotFeasible, verify_closure_answer
from event_layers import compile_events
from pinned_compiler import solve_pinned
from mixed import verify_mixed
from layered import realize


def floor(x: F) -> int:
    return x.numerator // x.denominator


def ceil(x: F) -> int:
    return -((-x.numerator) // x.denominator)


def check_scores(problem, scores: Sequence[Sequence[int]]):
    if len(scores) != problem.n or any(len(row) != problem.m + 1 for row in scores):
        raise ValueError('one score per group and possible cell is required')
    if any(type(v) is not int for row in scores for v in row):
        raise TypeError('integer scores required; scale rational contrasts explicitly')
    return tuple(tuple(row) for row in scores)


@dataclass
class ScoreSpectrum:
    minima: dict[int, F]
    witnesses: dict[int, dict]
    complete: bool
    examined: int
    total_patterns: int
    elapsed: float

    def interval(self, budget) -> dict:
        budget = F(budget)
        if budget < 0:
            raise ValueError('negative budget')
        if not self.complete:
            raise ValueError('an incomplete spectrum cannot certify exact extrema')
        values = sorted(v for v, d in self.minima.items() if d <= budget)
        if not values:
            return {'status': 'empty_uncertainty_set', 'budget': str(budget)}
        return {'status': 'exact', 'budget': str(budget), 'lower': values[0],
                'upper': values[-1], 'minimum_witness': self.witnesses[values[0]],
                'maximum_witness': self.witnesses[values[-1]]}

    def serializable(self):
        return {'minima': {str(k):str(v) for k,v in sorted(self.minima.items())},
                'witnesses': self.witnesses, 'complete': self.complete,
                'examined': self.examined, 'total_patterns': self.total_patterns,
                'elapsed': self.elapsed}


class ExpressionEngine:
    def __init__(self, problem, structural_pins=None, *, cpp_binary=None):
        self.problem = deepcopy(problem)
        problem = self.problem
        self.pins = dict(structural_pins or {})
        self.cpp_binary = cpp_binary
        self._compiled = None
        self._minimum_model = None
        physical = (problem.baselines, problem.boxes, problem.source, problem.slopes, problem.minimum_gaps, problem.maximum_gaps, problem.bands, problem.clearance, problem.mode, sorted(self.pins.items()))
        self.model_key = hashlib.sha256(repr(physical).encode()).hexdigest()
        try:
            self.source_engine = PinClosure(problem, self.pins)
        except SourceNotFeasible:
            self.source_engine = None
        self._geometric_cache = {}
        self._boundary_cache = {}

    @property
    def compiled(self):
        if self._compiled is None:
            self._compiled = compile_events(self.problem, self.pins)
        return self._compiled

    def clipped(self, pins):
        if self.source_engine is None:
            raise SourceNotFeasible('target projection requires a source satisfying all ownership and structural constraints')
        key = tuple(sorted(pins.items()))
        if key not in self._geometric_cache:
            ans = self.source_engine.query(pins)
            verify_closure_answer(self.source_engine, pins, ans)
            self._geometric_cache[key] = ans
            if ans['feasible']:
                self._boundary_cache[tuple(ans['owners'])] = {'area':ans['area'],'boundaries':ans['boundaries']}
        return self._geometric_cache[key]

    def spectrum(self, scores, constant=0, *, pattern_limit=100000):
        """Exact least deformation for each score, including unsupported outcomes.

        Unscored owners are free collateral variables. For every scored assignment,
        source-clipped completion minimizes cost among ALL of its completions.
        """
        scores = check_scores(self.problem, scores)
        if type(constant) is not int or pattern_limit < 1:
            raise ValueError('invalid fixed contribution or pattern limit')
        if self.source_engine is None:
            raise SourceNotFeasible('source-feasible target projection unavailable')
        t0 = perf_counter()
        variable = [g for g,row in enumerate(scores) if len(set(row)) > 1]
        base = constant + sum(row[0] for row in scores if len(set(row)) == 1)
        domains = [[o for o in range(self.problem.m+1)
                    if self.source_engine.closure({g:o})['feasible']] for g in variable]
        total = 1
        for dom in domains: total *= len(dom)
        minima = {}; witnesses = {}; examined = 0
        for owners in product(*domains):
            if examined >= pattern_limit: break
            examined += 1
            pins = dict(zip(variable, owners)); answer = self.clipped(pins)
            if not answer['feasible']: continue
            value = base + sum(scores[g][o] for g,o in pins.items()); cost = F(answer['area'])
            if value not in minima or cost < minima[value]:
                minima[value] = cost
                witnesses[value] = {'pins': pins, 'owners': answer['owners'],
                                    'area': str(cost), 'boundaries': answer['boundaries']}
        return ScoreSpectrum(minima, witnesses, examined == total, examined, total, perf_counter()-t0)

    def score_key(self, scores):
        return hashlib.sha256(json.dumps([self.model_key, scores], separators=(",",":")).encode()).hexdigest()

    def penalized(self, scores, multiplier, *, sign=1, certificate=True):
        scores = check_scores(self.problem, scores); mu = F(multiplier)
        if mu < 0 or sign not in (-1,1): raise ValueError('invalid multiplier or direction')
        costs = [[-sign*v for v in row] for row in scores]
        if self.cpp_binary is None:
            ans = solve_pinned(self.compiled, costs, area_weight=mu,
                               certificate=certificate, reconstruct=False)
        else:
            from .cpp_backend import solve as solve_cpp
            ans = solve_cpp(self.compiled, costs, mu, self.cpp_binary, reconstruct=False)
            if not certificate:
                ans.pop('certificate', None)
        if ans['feasible']:
            key = tuple(ans['owners'])
            if key not in self._boundary_cache:
                witness = realize(self.problem, ans['owners'])
                if not witness['feasible'] or witness['area'] != ans['area']:
                    raise AssertionError('compiled cost differs from reconstructed geometry')
                verify_mixed(self.problem, ans['owners'], witness['boundaries'], witness['area'])
                self._boundary_cache[key] = witness
            witness = self._boundary_cache[key]
            if witness['area'] != ans['area']:
                raise AssertionError('cached physical area differs')
            ans['boundaries'] = witness['boundaries']
            ans['score'] = sum(row[o] for row,o in zip(scores, ans['owners']))
        ans.update(multiplier=str(mu), sign=sign, query_key=self.score_key(scores))
        return ans

    def dual_queries(self, scores, multipliers, *, certificate=True):
        return [self.penalized(scores, mu, sign=s, certificate=certificate)
                for s in (1,-1) for mu in sorted(set(map(F,multipliers)))]

    def enclosure(self, scores, budget, constant=0, *, queries=None,
                  multipliers=(0,F(1,4),1,4,16), refine=False, pattern_limit=100000):
        """Return sound outer endpoints and attained inner bounds.

        L_bound <= min S <= L_witness; U_witness <= max S <= U_bound.
        Numerical tolerances never enter endpoint rounding. Empty models and
        unresolved feasibility are distinct from preservation certificates.
        """
        scores = check_scores(self.problem, scores); delta = F(budget)
        if delta < 0 or type(constant) is not int: raise ValueError('invalid budget or offset')
        if queries is None: queries = self.dual_queries(scores, multipliers)
        low = constant + sum(min(row) for row in scores)
        high = constant + sum(max(row) for row in scores)
        attained = []
        if self.source_engine is not None:
            ans = self.clipped({})
            if ans['feasible'] and F(ans['area']) <= delta:
                attained.append({'value':constant+sum(row[o] for row,o in zip(scores,ans['owners'])),
                                 'area':str(ans['area']),'owners':ans['owners'],'boundaries':ans['boundaries']})
        if self.source_engine is None:
            if self._minimum_model is None:
                self._minimum_model = solve_pinned(self.compiled, [[0]*(self.problem.m+1) for _ in range(self.problem.n)], area_weight=1)
            a = self._minimum_model
            if not a['feasible']:
                return {'status':'model_infeasible','budget':str(delta)}
            if F(a['area'])>delta:
                return {'status':'empty_uncertainty_set','budget':str(delta),'minimum_feasible_deformation':str(a['area'])}
            verify_mixed(self.problem, a['owners'], a['boundaries'], a['area'])
            attained.append({'value':constant+sum(row[o] for row,o in zip(scores,a['owners'])), 'area':str(a['area']),'owners':a['owners'],'boundaries':a['boundaries']})
        for ans in queries:
            if ans.get('query_key') != self.score_key(scores):
                raise ValueError('penalized answer belongs to another physical model or molecular score')
            if not ans['feasible']: return {'status':'model_infeasible','budget':str(delta)}
            mu = F(ans['multiplier']); m = F(ans['objective'])
            if ans['sign'] == 1: high = min(high, constant + floor(mu*delta-m))
            else: low = max(low, constant + ceil(m-mu*delta))
            if F(ans['area']) <= delta:
                attained.append({'value':constant+ans['score'],'area':str(ans['area']),
                                 'owners':ans['owners'],'boundaries':ans['boundaries']})
        if low > high: return {'status':'empty_uncertainty_set','budget':str(delta)}
        if refine and self.source_engine is not None:
            sp = self.spectrum(scores, constant, pattern_limit=pattern_limit)
            if sp.complete:
                result = sp.interval(delta)
                result.update(dual_outer=[low,high],target_patterns=sp.examined)
                return result
            for value,cost in sp.minima.items():
                if cost <= delta:
                    attained.append(dict(sp.witnesses[value],value=value))
        if not attained:
            return {'status':'feasibility_unresolved','budget':str(delta),'lower':low,'upper':high}
        lower_w = min(attained,key=lambda z:z['value']); upper_w = max(attained,key=lambda z:z['value'])
        if not low <= lower_w['value'] <= upper_w['value'] <= high:
            raise AssertionError('invalid dual enclosure or primal witness')
        return {'status':'exact' if (low,high)==(lower_w['value'],upper_w['value']) else 'bounded',
                'budget':str(delta),'lower':low,'upper':high,
                'minimum_witness_value':lower_w['value'],'maximum_witness_value':upper_w['value'],
                'lower_gap':lower_w['value']-low,'upper_gap':high-upper_w['value'],
                'minimum_witness':lower_w,'maximum_witness':upper_w}


def compose_spectra(parts, constant=0):
    """Min-plus composition, conditional on independent gluing/disjoint RNA.

    This function does not certify that overlapping native corridors can be glued.
    Its geometric precondition must be established separately.
    """
    state = {constant:F(0)}
    for part in parts:
        if not part.complete: raise ValueError('complete local spectra required')
        new = {}
        for old,area in state.items():
            for val,cost in part.minima.items():
                v=old+val; d=area+cost
                if v not in new or d < new[v]: new[v]=d
        state=new
    return state
