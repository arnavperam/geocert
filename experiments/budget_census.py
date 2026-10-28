"""Complete local cell/gene interval functions on all 84 unchanged native models.

Budgets 0, 1/4, 1, 4, 16 are an illustrative grid, not fitted error magnitudes.
Exact minimum cost at every count additionally represents the entire budget axis.
"""
from pathlib import Path
import argparse, csv, gzip, json, sys
from fractions import Fraction as F
from collections import Counter
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'core')]
from data_api import load_cases, encode
from pin_closure import PinClosure
from budget_counts import CountTarget, LagrangianCounts, enumerate_gene_states, count_frontier, exact_interval
from geocert_counts import target_for


def run(out):
    inputs = json.loads((ROOT / 'data/readout_inputs.json').read_text())
    models, rows, stats = [], [], Counter()
    start = perf_counter()
    for _, ch, setting, p, pins, n, _ in load_cases():
        engine = PinClosure(p, pins)
        assert not engine.inconsistent
        oracle = None
        singleton_costs = {}
        for g in range(n):
            singleton_costs[g] = {}
            for owner in range(p.m + 1):
                a = engine.query({g: owner})
                stats['singleton_queries'] += 1
                if a['feasible']:
                    singleton_costs[g][owner] = a['area']
        cells = inputs['rois'][ch['roi']]['cells']
        genes = sorted({pt['gene'] for pt in ch['points']} |
                       {gene for cell in ch['cell_labels'] for gene in cells[str(cell)]['counts']})
        model = {'id': ch['id'], 'roi': ch['roi'], 'setting': setting, 'cell_labels': ch['cell_labels'],
                 'source_admissible': True, 'genes': [], 'default_other_gene_interval': [0, 0]}
        for gene in genes:
            groups = tuple(i for i, pt in enumerate(ch['points']) if pt['gene'] == gene)
            enumeration = enumerate_gene_states(engine, groups)
            stats['gene_models'] += 1
            stats['combinations_checked'] += enumeration['combinations_checked']
            stats['max_combinations_per_gene_model'] = max(stats['max_combinations_per_gene_model'], enumeration['combinations_checked'])
            gdata = {'gene': gene, 'enumeration': enumeration, 'cells': []}
            for owner, cell in enumerate(ch['cell_labels']):
                target, source = target_for(ch, cell, gene, inputs)
                frontier = count_frontier(p, target, enumeration)
                stats['cell_gene_functions'] += 1
                # Sound zero-flow shortcut: every group's membership in this cell
                # is fixed by individual feasibility, before examining any costs.
                lo = target.outside + sum(options == [owner] for options in enumeration['choices'])
                hi = target.outside + sum(owner in options for options in enumeration['choices'])
                fixed = lo == hi
                mx = mn = None
                if not fixed:
                    stats['variable_cell_gene_functions'] += 1
                    if oracle is None:
                        oracle = LagrangianCounts(p, pins)
                    mx, mn = oracle.supported_path(target), oracle.supported_path(target, -1)
                    assert mx['complete'] and mn['complete']
                    stats['cut_calls'] += mx['oracle_calls'] + mn['oracle_calls']
                cdata = {'cell': cell, 'reference_count': source, 'outside_count': target.outside,
                         'frontier': frontier, 'fixed_by_feasibility': fixed}
                if not fixed:
                    cdata['dual_paths'] = {name: {'complete': path['complete'],
                            'cuts': [{k: cut[k] for k in ['mu', 'objective', 'count', 'area']} for cut in path['cuts']]}
                            for name, path in [('maximum', mx), ('minimum', mn)]}
                areas = sorted({f['minimum_deformation'] for f in frontier})
                grid = [F(0), F(1, 4), F(1), F(4), F(16)]
                probes = sorted(set(grid + areas + [(a + b) / 2 for a, b in zip(areas, areas[1:])]))
                for delta in probes:
                    exact = exact_interval(frontier, delta)
                    if fixed:
                        bound = {'lower': lo, 'upper': hi, 'attained_lower': lo, 'attained_upper': hi}
                    else:
                        bound = oracle.enclosure(target, delta, mx, mn)
                    assert bound['lower'] <= exact['lower'] <= exact['upper'] <= bound['upper']
                    stats['dual_soundness_checks'] += 1
                    endpoint_gap = (exact['lower'] - bound['lower']) + (bound['upper'] - exact['upper'])
                    stats['probes_with_strict_dual_gap'] += bool(endpoint_gap)
                    if delta not in grid:
                        continue
                    stats['grid_intervals'] += 1
                    stats['grid_nontrivial_intervals'] += exact['lower'] < exact['upper']
                    stats['grid_exact_dual_endpoints'] += endpoint_gap == 0
                    allowed = [[o for o, area in singleton_costs[g].items() if area <= delta] for g in groups]
                    single_lower = target.outside + sum(options == [owner] for options in allowed)
                    single_upper = target.outside + sum(owner in options for options in allowed)
                    assert single_lower <= exact['lower'] <= exact['upper'] <= single_upper
                    stats['grid_joint_tighter_than_singleton'] += (single_lower < exact['lower'] or exact['upper'] < single_upper)
                    rows.append({'corridor': ch['id'], 'setting': setting, 'cell': cell, 'gene': gene,
                                 'budget': str(delta), 'reference_count': source, 'lower': exact['lower'],
                                 'upper': exact['upper'], 'dual_lower': bound['lower'], 'dual_upper': bound['upper'],
                                 'singleton_lower': single_lower, 'singleton_upper': single_upper,
                                 'penalty_witness_min': bound['attained_lower'], 'penalty_witness_max': bound['attained_upper'],
                                 'detection': 'always_detected' if exact['lower'] > 0 else 'always_zero' if exact['upper'] == 0 else 'variable'})
                gdata['cells'].append(cdata)
            model['genes'].append(gdata)
        models.append(model)
        print(ch['id'], setting, 'complete', flush=True)
    out.mkdir(parents=True, exist_ok=True)
    artifact = {'schema': 'geocert.native_budget_functions.v1',
                'scope': 'Separate local corridor families, Q20 counts, all other crop assignments fixed. No cross-corridor stitching.',
                'budget_unit': 'square micrometres of summed interface displacement',
                'gene_vocabulary': sorted({gene for roi in inputs['rois'].values() for cell in roi['cells'].values() for gene in cell['counts']}),
                'default_interval_scope': 'Defaults apply only to this observed-gene vocabulary; assay membership of other genes is unknown.',
                'calibrated': False, 'illustrative_budgets': ['0', '1/4', '1', '4', '16'],
                'complete_interval_functions_for_every_nonnegative_budget': True,
                'statistics': dict(stats), 'models': models}
    (out / 'native_frontiers.json').write_text(json.dumps(artifact, indent=2, default=encode) + '\n')
    with (out / 'native_frontiers.json.gz').open('wb') as raw:
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as compressed:
            compressed.write(json.dumps(artifact, separators=(',', ':'), default=encode).encode())
    with (out / 'native_intervals.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    summary = {'models': len(models), 'statistics': dict(stats), 'seconds': perf_counter() - start}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, default=ROOT / 'run/budgets')
    run(ap.parse_args().output_dir)
