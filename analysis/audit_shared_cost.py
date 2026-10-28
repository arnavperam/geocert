"""Reproduce the joint-versus-additive energy ablation on all published models.

This is a representation ablation, not a segmentation-accuracy benchmark.
The additive model receives the same feasible assignments and its best possible
joint-model tie choice. Run with Python 3.12+ and assertions enabled.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import argparse
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'core'))
from data_api import load_cases, reassignment_costs, encode
from pin_closure import PinClosure
from scc_solver import prepare, solve_compiled


def run(output: Path) -> None:
    if not __debug__:
        raise RuntimeError('Run without -O: verification assertions are required.')
    rows = []
    start = time.perf_counter()
    for _, chart, setting, problem, pins, n_rna, expected in load_cases():
        closure = PinClosure(problem, pins)
        singletons = [[F(0)] * (problem.m + 1) for _ in range(problem.n)]
        for group in range(n_rna):
            for owner in range(problem.m + 1):
                if owner == problem.source[group]:
                    continue
                answer = closure.query({group: owner})
                # Infeasible labels remain excluded by the unchanged hard graph.
                if answer['feasible']:
                    singletons[group][owner] = answer['area']
        compiled, plan, _ = prepare(problem, pins, 'direct', True)
        rewards = reassignment_costs(problem, n_rna)
        bound = (problem.hi - problem.lo) * sum(hi - lo for lo, hi in problem.bands)
        for penalty in (F(1, 4), F(1), F(4)):
            joint = solve_compiled(compiled, 'direct', rewards, penalty, plan)
            additive_costs = [[penalty * d + r for d, r in zip(ds, rs)]
                              for ds, rs in zip(singletons, rewards)]
            additive = solve_compiled(compiled, 'direct', additive_costs, F(0), plan)
            scale = lcm(*(v.denominator for row in additive_costs for v in row))
            dominance = int(penalty * bound + n_rna) + 2
            tie_costs = [[dominance * scale * a + r for a, r in zip(aa, rr)]
                         for aa, rr in zip(additive_costs, rewards)]
            best_tie = solve_compiled(compiled, 'direct', tie_costs, penalty, plan)
            assert joint['feasible'] and additive['feasible'] and best_tie['feasible']
            owners = best_tie['owners']
            primary = sum(additive_costs[g][o] for g, o in enumerate(owners))
            assert primary == additive['objective']
            gain = sum(owners[g] != problem.source[g] for g in range(n_rna))
            evaluated = penalty * best_tie['area'] - gain
            regret = evaluated - joint['objective']
            assert regret >= 0
            rows.append({
                'chart': chart['id'], 'setting': setting, 'penalty': str(penalty),
                'joint_objective': str(joint['objective']), 'joint_area': str(joint['area']),
                'joint_owners': joint['owners'][:n_rna],
                'joint_changed': sum(joint['owners'][g] != problem.source[g] for g in range(n_rna)),
                'additive_primary': str(primary), 'additive_joint_objective': str(evaluated),
                'additive_area': str(best_tie['area']), 'additive_owners': owners[:n_rna],
                'additive_changed': gain, 'regret': str(regret),
            })
        print(chart['id'], setting, len(rows), flush=True)
    summary = {}
    for setting in ('baseline', 'maximal'):
        subset = [r for r in rows if r['setting'] == setting]
        affected = [r for r in subset if F(r['regret']) > 0]
        summary[setting] = {'queries': len(subset), 'positive_regret': len(affected),
                            'affected_corridors': len({r['chart'] for r in affected})}
    assert summary['baseline']['positive_regret'] == 1
    assert summary['maximal']['positive_regret'] == 16
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({'summary': summary, 'rows': rows,
        'elapsed_seconds': time.perf_counter() - start}, indent=2, default=encode) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'run/shared_cost.json')
    run(parser.parse_args().output)
