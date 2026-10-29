"""Independently reconcile saved local witnesses, cell traces and CSV endpoints.

This checks the algebra and completeness of exported results. Geometric witness
validity is checked separately by the original rational validator and producer.
It does not call the optimizer, dynamic program or composition implementation.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from fractions import Fraction as F
import json
from pathlib import Path
import time

BUDGETS = {F(s) for s in ['0', '1/20', '1/10', '1/4', '1/2', '1', '2', '4']}


def main():
    if not __debug__:
        raise RuntimeError('Run without Python -O; all scientific assertions are required')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prepared-root', type=Path, required=True)
    args = ap.parse_args()
    root = args.prepared_root
    started = time.perf_counter()
    atlases = [json.loads(p.read_text()) for p in (root / 'data').glob('atlas_*.json')]
    atlas_index = {a['roi']: a for a in atlases}
    assert len(atlas_index) == len(atlases)
    patches = {(a['roi'], p['id']): p for a in atlases for p in a['patches']}
    expected_local = {(roi, pid, gene) for (roi, pid), p in patches.items()
                      for gene in {point['gene'] for point in p['points']}}
    expected_cells = {(a['roi'], int(cell), gene) for a in atlases
                      for cell, genes in a['source_counts'].items() if int(cell) in set(a['complete_cells'])
                      for gene, number in genes.items() if number > 0}
    expected_parts = defaultdict(set)
    for roi, pid, gene in expected_local:
        for cell in patches[roi, pid]['cell_labels']:
            if cell > 0:
                expected_cells.add((roi, cell, gene))
                expected_parts[roi, cell, gene].add(pid)
    local_rows = json.loads((root / 'results/patch_spectra.json').read_text())
    local = {}
    witness_count = 0
    for row in local_rows:
        key = row['roi'], row['patch'], row['gene']
        assert key not in local
        patch = patches[key[:2]]
        points = patch['points']
        assert row['source_left'] == sum(p['gene'] == row['gene'] and p['source'] == 0 for p in points)
        assert F(row['minima'][str(row['source_left'])]) == 0
        assert all(int(v) >= 0 and F(c) >= 0 for v,c in row['minima'].items())
        assert set(row['minima']) == set(row['witnesses'])
        for value, witness in row['witnesses'].items():
            # Structural variables follow RNA, in first-encountered side order.
            sides = list(dict.fromkeys(c[0] for c in patch.get('nuclear_constraints', [])))
            assert len(witness['owners']) == len(points) + len(sides)
            assert witness['owners'][len(points):] == sides
            assert all(owner in (0, 1) for owner in witness['owners'])
            assert sum(p['gene'] == row['gene'] and owner == 0 for p, owner in zip(points, witness['owners'])) == int(value)
            assert F(witness['physical_area_um2']) == F(row['minima'][value])
            witness_count += 1
        local[key] = row
    assert set(local) == expected_local, 'Missing or extra local gene spectra'
    states = json.loads((root / 'results/cell_spectra.json').read_text())
    state_index = {}
    traces_checked = 0
    for state in states:
        key = state['roi'], state['cell'], state['gene']
        assert key not in state_index
        atlas = atlas_index[state['roi']]
        assert state['atlas_id'] == atlas['atlas_id']
        assert state['source_count'] == atlas['source_counts'].get(str(state['cell']), {}).get(state['gene'], 0)
        assert F(state['minima'][str(state['source_count'])]) == 0
        assert all(int(v) >= 0 and F(c) >= 0 for v,c in state['minima'].items())
        assert len(state['patches']) == len(set(state['patches']))
        assert set(state['patches']) == expected_parts[key]
        assert set(state['minima']) == set(state['witness_traces'])
        # Recompute scalar min-plus costs without calling the producer's compose.
        costs = {state['source_count']: F(0)}
        for pid in sorted(expected_parts[key]):
            entry = local[state['roi'], pid, state['gene']]
            direction = 1 if patches[state['roi'], pid]['cell_labels'][0] == state['cell'] else -1
            next_costs = {}
            for old_value, old_cost in costs.items():
                for left, local_cost in entry['minima'].items():
                    value = old_value + direction * (int(left) - entry['source_left'])
                    candidate = old_cost + F(local_cost)
                    if value not in next_costs or candidate < next_costs[value]:
                        next_costs[value] = candidate
            costs = next_costs
        assert costs == {int(v): F(c) for v,c in state['minima'].items()}
        for value, trace in state['witness_traces'].items():
            assert set(trace) == set(state['patches'])
            assert state['source_count'] + sum(trace.values()) == int(value)
            total_cost = F(0)
            for pid, delta in trace.items():
                patch = patches[state['roi'], pid]
                entry = local[state['roi'], pid, state['gene']]
                side = patch['cell_labels'].index(state['cell'])
                left_count = entry['source_left'] + delta * (1 if side == 0 else -1)
                total_cost += F(entry['minima'][str(left_count)])
            assert total_cost == F(state['minima'][value])
            traces_checked += 1
        state_index[key] = (state, {int(k): F(v) for k,v in state['minima'].items()})
    assert set(state_index) == expected_cells, 'Missing or extra cell gene spectra'
    seen = set()
    row_count = 0
    with (root / 'results/cell_expression_intervals.csv').open() as handle:
        for row in csv.DictReader(handle):
            key = row['roi'], int(row['cell']), row['gene']
            budget = F(row['budget_um2'])
            full_key = (*key, budget)
            assert full_key not in seen and budget in BUDGETS
            seen.add(full_key)
            state, minima = state_index[key]
            possible = [v for v,cost in minima.items() if cost <= budget]
            low, high = min(possible), max(possible)
            assert row['atlas_id'] == state['atlas_id']
            assert int(row['source_count']) == state['source_count']
            assert row['status'] == 'exact'
            assert (int(row['lower']), int(row['upper'])) == (low, high)
            assert (F(row['lower_cost_um2']), F(row['upper_cost_um2'])) == (minima[low], minima[high])
            row_count += 1
    assert row_count == len(state_index) * len(BUDGETS), 'Missing exported cell–gene–budget rows'
    summary = {'status': 'passed', 'local_gene_spectra': len(local), 'local_witnesses_checked': witness_count,
               'cell_gene_spectra': len(state_index), 'cell_state_traces_checked': traces_checked,
               'interval_rows_checked': row_count, 'budgets_per_question': len(BUDGETS),
               'seconds': time.perf_counter()-started,
               'scope': 'Complete local/cell census, structural owner bits, witness counts, independent scalar min-plus composition, trace costs and endpoints; local optimality is supplied by the exact local solver and geometric witnesses by its verifier'}
    (root / 'results/output_reconciliation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
