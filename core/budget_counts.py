"""Conditional hard-budget count intervals, exact witnesses, and dual bounds.

The uncertainty family is the supplied ordered-corridor model. Counts outside
its retained sites are held fixed. All arithmetic used in guarantees is rational.
This module does not supply anatomical calibration or statistical coverage.
"""
from dataclasses import dataclass, asdict
from fractions import Fraction as F
from itertools import product
from math import ceil, floor, lcm, prod
import hashlib
import json

from pin_closure import PinClosure, verify_closure_answer
from scc_solver import prepare, solve_compiled


@dataclass(frozen=True)
class CountTarget:
    """Integer molecular weights at selected groups, a target owner, and offset."""
    groups: tuple
    owner: int
    outside: int = 0
    weights: tuple = ()

    def checked_weights(self, problem):
        weights = self.weights or (1,) * len(self.groups)
        if len(set(self.groups)) != len(self.groups):
            raise ValueError("A group may occur only once")
        if not 0 <= self.owner <= problem.m or type(self.outside) is not int or self.outside < 0:
            raise ValueError("Invalid owner or fixed count")
        if len(weights) != len(self.groups) or any(type(w) is not int or w <= 0 for w in weights):
            raise ValueError("Positive integer weights are required")
        if any(type(g) is not int or not 0 <= g < problem.n for g in self.groups):
            raise ValueError("Invalid group")
        return weights

    def count(self, owners):
        weights = self.weights or (1,) * len(self.groups)
        return self.outside + sum(w for g, w in zip(self.groups, weights) if owners[g] == self.owner)

    def unary_scores(self, problem):
        scores = [[0] * (problem.m + 1) for _ in range(problem.n)]
        for g, w in zip(self.groups, self.checked_weights(problem)):
            scores[g][self.owner] = w
        return scores


@dataclass(frozen=True)
class IntegerContrast:
    """Signed integer linear count contrast; no statistical-DE interpretation."""
    groups: tuple
    scores: tuple
    outside: int = 0

    def checked_weights(self, problem):
        if len(set(self.groups)) != len(self.groups) or len(self.groups) != len(self.scores):
            raise ValueError('Distinct groups and one score row per group required')
        if type(self.outside) is not int or any(type(g) is not int or not 0 <= g < problem.n for g in self.groups):
            raise ValueError('Integer offset and valid groups required')
        if any(len(row) != problem.m + 1 or any(type(v) is not int for v in row) for row in self.scores):
            raise ValueError('One integer score per owner required')
        return tuple(max(row) - min(row) for row in self.scores)

    def count(self, owners):
        return self.outside + sum(row[owners[g]] for g, row in zip(self.groups, self.scores))

    def unary_scores(self, problem):
        self.checked_weights(problem)
        scores = [[0] * (problem.m + 1) for _ in range(problem.n)]
        for g, row in zip(self.groups, self.scores):
            scores[g] = list(row)
        return scores


def witness(answer):
    return {k: answer[k] for k in ("owners", "area", "boundaries")}


def enumerate_gene_states(engine, groups, max_assignments=100000):
    """Exact minimum-cost completion for each feasible choice of target owners.

    Requires the source-admissibility check in PinClosure. The exponential factor
    involves only these groups. No partial enumeration is reported as complete.
    """
    p = engine.problem
    groups = tuple(groups)
    CountTarget(groups, 0).checked_weights(p)
    choices = [[o for o in range(p.m + 1) if engine.closure({g: o})['feasible']]
               for g in groups]
    combinations = prod(map(len, choices))
    if combinations > max_assignments:
        raise ValueError(f"Target enumeration needs {combinations} combinations; cap is {max_assignments}")
    states = []
    for selected in product(*choices):
        pins = dict(zip(groups, selected))
        answer = engine.query(pins)
        verify_closure_answer(engine, pins, answer)
        if answer['feasible']:
            states.append({"selected": list(selected), **witness(answer)})
    return {"groups": list(groups), "choices": choices, "combinations_checked": combinations,
            "complete": True, "source_admissible": not engine.inconsistent, "states": states}


def count_frontier(problem, target, enumeration):
    """Minimum deformation for each attainable count, retaining witness indices."""
    target.checked_weights(problem)
    if not enumeration['complete'] or tuple(enumeration['groups']) != target.groups:
        raise ValueError("Complete matching target enumeration required")
    best = {}
    for i, state in enumerate(enumeration['states']):
        q = target.count(state['owners'])
        if q not in best or F(state['area']) < best[q]['minimum_deformation']:
            best[q] = {"count": q, "minimum_deformation": F(state['area']), "state_index": i}
    return [best[q] for q in sorted(best)]


def exact_interval(frontier, budget):
    budget = F(budget)
    if budget < 0:
        raise ValueError("Budget must be nonnegative")
    feasible = [r for r in frontier if F(r['minimum_deformation']) <= budget]
    if not feasible:
        return {"status": "empty_budget_family", "budget": budget}
    lo = min(feasible, key=lambda r: r['count'])
    hi = max(feasible, key=lambda r: r['count'])
    return {"status": "exact", "budget": budget, "lower": lo['count'], "upper": hi['count'],
            "lower_state_index": lo['state_index'], "upper_state_index": hi['state_index']}


class LagrangianCounts:
    """Repeated exact cuts give sound outer intervals, including duality gaps.

    supported_path uses dichotomic oracle calls. It is not a GGT implementation
    and does not assume nested ownership sets. A query limit returns sound bounds
    but never silently labels an incomplete path complete.
    """
    def __init__(self, problem, pins):
        self.problem = problem
        self.pins = dict(pins)
        self.model_signature = self._model_signature()
        self.compiled, self.plan, _ = prepare(problem, pins, 'direct', True)
        self.minimum = solve_compiled(self.compiled, 'direct',
                                      [[0] * (problem.m + 1) for _ in range(problem.n)],
                                      F(1), self.plan, certificate=True)

    def _model_signature(self):
        p = self.problem
        payload = {k: getattr(p, k) for k in ('baselines', 'source', 'slopes',
                   'minimum_gaps', 'maximum_gaps', 'clearance', 'bands', 'mode')}
        payload.update(boxes=[asdict(b) for b in p.boxes],
                       model_class=type(p).__module__ + '.' + type(p).__qualname__,
                       pins=sorted(self.pins.items()))
        return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()

    def _binding(self, target):
        if self._model_signature() != self.model_signature:
            raise ValueError('Model changed after compilation; construct a new oracle')
        payload = {'model': self.model_signature, 'outside': target.outside,
                   'scores': target.unary_scores(self.problem)}
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def supported_path(self, target, direction=1, max_queries=256):
        p = self.problem
        weights = target.checked_weights(p)
        binding = self._binding(target)
        if direction not in (-1, 1):
            raise ValueError("Direction must be +1 (maximum) or -1 (minimum)")
        if not self.minimum['feasible']:
            return {"status": "model_infeasible", "direction": direction,
                    "binding": binding, "complete": True, "cuts": []}
        if max_queries < 3:
            raise ValueError("At least three oracle calls must be allowed")
        costs = [[-direction * value for value in row] for row in target.unary_scores(p)]
        cache = {}

        def solve(mu):
            mu = F(mu)
            if mu in cache:
                return cache[mu]
            if len(cache) >= max_queries:
                raise OverflowError("Supported-path oracle limit")
            a = solve_compiled(self.compiled, 'direct', costs, mu, self.plan, certificate=True)
            if not a['feasible']:
                raise AssertionError("Fixed-family feasibility changed with objective")
            q = target.count(a['owners'])
            r = direction * (q - target.outside)
            assert a['objective'] == mu * a['area'] - r
            row = {"mu": mu, "objective": a['objective'], "count": q, "reward": r,
                   **witness(a), "certificate": a['certificate'], "mapping": a['mapping']}
            cache[mu] = row
            return row

        # Exact lexicographic endpoints using explicit range and lattice bounds.
        area_bound = (p.hi - p.lo) * sum(hi - lo for lo, hi in p.bands)
        denominator = lcm(F(self.compiled.constant_area).denominator,
                          *(z.area.denominator for z in self.compiled.positive.values()),
                          *(z.area.denominator for z in self.compiled.negative.values()))
        solve(F(0))
        right = solve(F(1) / (area_bound + 1))  # largest integer reward, then least area
        left = solve(F((sum(weights) + 1) * denominator))  # least area, then largest reward
        stack = [(left, right)]
        complete = True
        try:
            while stack:
                a, b = stack.pop()
                if a['area'] == b['area']:
                    assert a['reward'] == b['reward']
                    continue
                assert a['area'] < b['area'] and a['reward'] < b['reward']
                mu = F(b['reward'] - a['reward']) / (b['area'] - a['area'])
                z = solve(mu)
                chord = mu * a['area'] - a['reward']
                if z['objective'] < chord:
                    assert a['area'] < z['area'] < b['area']
                    stack.extend([(a, z), (z, b)])
                else:
                    assert z['objective'] == chord
        except OverflowError:
            complete = False
        return {"status": "bounded", "direction": direction, "complete": complete,
                "binding": binding,
                "cuts": [cache[k] for k in sorted(cache)], "oracle_calls": len(cache),
                "count_lattice": 1, "area_lattice_denominator": denominator,
                "minimum_area": self.minimum['area'], "minimum_witness": witness(self.minimum)}

    def enclosure(self, target, budget, maximum_path, minimum_path):
        target.checked_weights(self.problem)
        budget = F(budget)
        if budget < 0:
            raise ValueError("Budget must be nonnegative")
        binding = self._binding(target)
        if any(path.get('binding') != binding for path in (maximum_path, minimum_path)):
            raise ValueError('Paths must use the same model, pins, and target as this query')
        if maximum_path['direction'] != 1 or minimum_path['direction'] != -1:
            raise ValueError("Mismatched path directions")
        if not self.minimum['feasible'] or self.minimum['area'] > budget:
            return {"status": "empty_budget_family", "budget": budget,
                    "minimum_area": self.minimum.get('area'), "certificate": self.minimum['certificate']}
        uppers = [(floor(target.outside + c['mu'] * budget - c['objective']), i)
                  for i, c in enumerate(maximum_path['cuts'])]
        lowers = [(ceil(target.outside + c['objective'] - c['mu'] * budget), i)
                  for i, c in enumerate(minimum_path['cuts'])]
        upper, ui = min(uppers)
        lower, li = max(lowers)
        feasible = [self.minimum] + [c for path in [maximum_path, minimum_path]
                                    for c in path['cuts'] if c['area'] <= budget]
        lo_w = min(feasible, key=lambda c: target.count(c['owners']))
        hi_w = max(feasible, key=lambda c: target.count(c['owners']))
        a, b = target.count(lo_w['owners']), target.count(hi_w['owners'])
        assert lower <= a <= b <= upper
        return {"status": "exact" if lower == a and b == upper else "outer_bound",
                "budget": budget, "lower": lower, "upper": upper,
                "attained_lower": a, "attained_upper": b,
                "lower_gap": a - lower, "upper_gap": upper - b,
                "lower_cut_index": li, "upper_cut_index": ui,
                "complete_dual_paths": maximum_path['complete'] and minimum_path['complete'],
                "lower_witness": witness(lo_w), "upper_witness": witness(hi_w)}
