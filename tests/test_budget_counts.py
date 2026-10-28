"""Independent small-state checks of budget bounds and all-optima projection."""
from pathlib import Path
import sys, json
from fractions import Fraction as F
from itertools import product

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'core'), str(ROOT / 'tests')]
from layered import Layers, Box, realize
from pin_closure import PinClosure, SourceNotFeasible
from budget_counts import (CountTarget, IntegerContrast, LagrangianCounts, enumerate_gene_states,
                           count_frontier, exact_interval)
from optimal_cuts import all_min_cuts, owner_options, enumerate_optimal_owners
from model_generator import make
from data_api import load_cases, encode
from calibration import split_budget


def toy():
    p = Layers([[(F(0), F(0)), (F(3), F(0))]],
               [Box(F(1, 5), F(1, 5), t, t, i) for i, t in enumerate([F(1), F(5, 4)])],
               [1, 1], F(1), [], [], clearance=F(1, 20), bands=[(F(-1, 2), F(1, 2))])
    return p


def brute(p, pins=None):
    rows = []
    for owners in product(range(p.m + 1), repeat=p.n):
        if any(owners[g] != o for g, o in (pins or {}).items()):
            continue
        w = realize(p, owners)
        if w['feasible']:
            rows.append({'owners': owners, 'area': w['area']})
    return rows


def check_problem(p, target, rows):
    engine = PinClosure(p, {})
    states = enumerate_gene_states(engine, target.groups)
    frontier = count_frontier(p, target, states)
    expected = {}
    for r in rows:
        q = target.count(r['owners'])
        expected[q] = min(expected.get(q, r['area']), r['area'])
    assert {r['count']: r['minimum_deformation'] for r in frontier} == expected
    oracle = LagrangianCounts(p, {})
    mx, mn = oracle.supported_path(target), oracle.supported_path(target, -1)
    assert mx['complete'] and mn['complete']
    for direction, path in [(1, mx), (-1, mn)]:
        for r in path['cuts']:
            opt = min(r['mu'] * w['area'] - direction * (target.count(w['owners']) - target.outside) for w in rows)
            assert opt == r['objective']
    levels = sorted({r['area'] for r in rows})
    budgets = sorted(set(levels + [(a + b) / 2 for a, b in zip(levels, levels[1:])] + [F(0), F(10)]))
    for delta in budgets:
        e = exact_interval(frontier, delta)
        valid = [target.count(r['owners']) for r in rows if r['area'] <= delta]
        assert (e['lower'], e['upper']) == (min(valid), max(valid))
        outer = oracle.enclosure(target, delta, mx, mn)
        assert outer['lower'] <= e['lower'] <= e['upper'] <= outer['upper']
        if outer['status'] == 'exact':
            assert (outer['lower'], outer['upper']) == (e['lower'], e['upper'])
    # Verify all-optima projection, not just possible labels independently.
    a = mx['cuts'][0]
    dag = all_min_cuts(a['certificate'], p.n * p.m, a['mapping'])
    got = set(enumerate_optimal_owners(p, dag, max_free_components=20))
    best = min(a['mu'] * r['area'] - (target.count(r['owners']) - target.outside) for r in rows)
    wanted = {tuple(r['owners']) for r in rows if a['mu'] * r['area'] - (target.count(r['owners']) - target.outside) == best}
    assert got == wanted
    opts = owner_options(p, dag)
    assert opts == [sorted({c[g] for c in wanted}) for g in range(p.n)]
    return len(budgets), len(rows), len(got)


def main():
    totals = {'budget_comparisons': 0, 'complete_owner_states': 0, 'optimal_vectors_checked': 0}
    for seed in range(16):
        p, _ = make(seed, n=3, m=1 + seed % 2)
        target = CountTarget((0, 2), seed % (p.m + 1), outside=7, weights=(1, 2))
        nb, ns, no = check_problem(p, target, brute(p))
        for k, v in zip(totals, [nb, ns, no]): totals[k] += v
    for seed in range(4):
        p, _ = make(seed + 20, n=3, m=2)
        target = IntegerContrast((0, 2), ((1, -1, 0), (2, 0, -2)), outside=-3)
        nb, ns, no = check_problem(p, target, brute(p))
        for k, v in zip(totals, [nb, ns, no]): totals[k] += v
    assert split_budget(list(range(19)))['budget'] == 18
    assert split_budget(list(range(5)))['status'] == 'abstain'
    assert split_budget([None] * 19)['status'] == 'abstain'
    try:
        split_budget([1, 2], calibration_ids=['a', 'b'], evaluation_ids=['b'])
        raise AssertionError('Calibration/evaluation leakage was accepted')
    except ValueError:
        pass

    p = toy(); target = CountTarget((0, 1), 0)
    states = enumerate_gene_states(PinClosure(p, {}), target.groups)
    frontier = count_frontier(p, target, states)
    assert {r['count']: r['minimum_deformation'] for r in frontier} == {0: F(0), 1: F(1, 16), 2: F(7, 64)}
    oracle = LagrangianCounts(p, {}); mx = oracle.supported_path(target); mn = oracle.supported_path(target, -1)
    delta = F(3, 50)
    exact = exact_interval(frontier, delta)
    outer = oracle.enclosure(target, delta, mx, mn)
    for other_oracle, other_target in [(oracle, CountTarget((0,), 0)),
                                        (LagrangianCounts(p, {0: 1}), target)]:
        try:
            other_oracle.enclosure(other_target, delta, mx, mn)
            raise AssertionError('Foreign objective/model paths were accepted')
        except ValueError:
            pass
    assert exact['upper'] == 0 and outer['upper'] == 1 and outer['complete_dual_paths']
    # Count one is a hard-budget optimum at 1/16 but is never penalty optimal.
    assert exact_interval(frontier, F(1, 16))['upper'] == 1
    assert all(c['count'] != 1 for c in mx['cuts'])
    tie = next(c for c in mx['cuts'] if c['mu'] == F(128, 7))
    dag = all_min_cuts(tie['certificate'], p.n * p.m, tie['mapping'])
    assert enumerate_optimal_owners(p, dag) == [(0, 0), (1, 1)]
    try:
        enumerate_gene_states(PinClosure(p, {}), target.groups, max_assignments=3)
        raise AssertionError('Enumeration cap was ignored')
    except ValueError:
        pass
    assert not oracle.supported_path(target, max_queries=3)['complete']
    bad = LagrangianCounts(p, {0: 0})
    assert bad.minimum['area'] == F(1, 16)
    a, b = bad.supported_path(target), bad.supported_path(target, -1)
    assert bad.enclosure(target, F(0), a, b)['status'] == 'empty_budget_family'
    try:
        PinClosure(p, {0: 0})
        raise AssertionError('Source-inconsistent shortcut was accepted')
    except SourceNotFeasible:
        pass

    # Native frontier against the previously independent scalar implementation.
    _, ch, _, p, pins, n, _ = next(z for z in load_cases() if z[1]['id'] == '4037abe9fea38e' and z[2] == 'baseline')
    rows = json.loads((ROOT / 'data/showcase_scalar_audit.json').read_text())['rows']
    native_checks = 0
    for gene in sorted({r['gene'] for r in ch['points']}):
        groups = tuple(i for i, r in enumerate(ch['points']) if r['gene'] == gene)
        enumeration = enumerate_gene_states(PinClosure(p, pins), groups)
        for owner in range(p.m + 1):
            t = CountTarget(groups, owner)
            f = count_frontier(p, t, enumeration)
            expected = {}
            for r in rows:
                q = t.count(r['owners']); area = F(r['D'])
                expected[q] = min(expected.get(q, area), area)
            assert {r['count']: r['minimum_deformation'] for r in f} == expected
            native_checks += 1
    sys.path.insert(0, str(ROOT))
    from geocert_counts import calculate
    try:
        calculate({'corridor': '4037abe9fea38e', 'setting': 'baseline', 'cell': 153124,
                   'gene': '__UNMEASURED_GENE__', 'budget': '1', 'method': 'exact',
                   'max_assignments': 100000})
        raise AssertionError('Unknown assay gene was certified as zero')
    except ValueError as e:
        assert 'vocabulary' in str(e)
    result = {'passed': True, 'random_count_models': 16, 'random_contrast_models': 4,
              **totals, 'native_gene_cell_frontiers': native_checks, 'calibration_utility_checks': 4,
              'foreign_path_rejections': 2, 'unknown_gene_rejections': 1,
              'dual_gap_example': {'budget': delta, 'exact_interval': [exact['lower'], exact['upper']],
                                  'best_integer_dual_interval': [outer['lower'], outer['upper']],
                                  'complete_dual_path': True, 'tie_multiplier': F(128, 7)},
              'unsupported_budget_optimum': {'budget': F(1, 16), 'maximum': 1},
              'all_optima_at_tie': [(0, 0), (1, 1)]}
    out = ROOT / 'results/budgets/tests.json'; out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, default=encode) + '\n')
    print(json.dumps(result, indent=2, default=encode))


if __name__ == '__main__':
    main()
