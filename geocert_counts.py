#!/usr/bin/env python3
"""Exact local expression-count intervals under an explicit area budget.

All other crop assignments remain fixed. Budgets are sensitivity parameters,
not estimated segmentation errors or statistical confidence levels.
"""
from pathlib import Path
import argparse, hashlib, json, sys
from fractions import Fraction as F

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'core'))
if not __debug__:
    raise SystemExit('Keep assertions enabled; do not run with -O.')
from data_api import load_cases, encode
from pin_closure import PinClosure
from budget_counts import CountTarget, LagrangianCounts, enumerate_gene_states, count_frontier, exact_interval


def target_for(ch, cell, gene, inputs):
    owner = ch['cell_labels'].index(cell)
    groups = tuple(i for i, pt in enumerate(ch['points']) if pt['gene'] == gene)
    source = inputs['rois'][ch['roi']]['cells'][str(cell)]['counts'].get(gene, 0)
    inside = sum(ch['points'][g]['source'] == owner for g in groups)
    if source < inside:
        raise ValueError('Source count is inconsistent with retained sites')
    return CountTarget(groups, owner, source - inside), source


def calculate(query):
    matches = [z for z in load_cases() if z[1]['id'] == query['corridor'] and z[2] == query['setting']]
    if len(matches) != 1:
        raise ValueError('Unknown corridor/setting')
    _, ch, _, p, pins, n, _ = matches[0]
    inputs = json.loads((ROOT / 'data/readout_inputs.json').read_text())
    observed_vocabulary = {gene for roi in inputs['rois'].values()
                           for cell in roi['cells'].values() for gene in cell['counts']}
    if query['gene'] not in observed_vocabulary:
        raise ValueError('Gene is outside the released observed-gene vocabulary; assay membership is not established.')
    target, reference = target_for(ch, query['cell'], query['gene'], inputs)
    record = {'schema': 'geocert.count_interval.v1', 'query': query,
              'scope': 'Conditional local corridor; all other supplied-crop assignments fixed. Q20 counts, not latent expression.',
              'reference_count': reference, 'outside_count': target.outside,
              'gene_vocabulary_scope': 'Genes observed in the released relevant source-cell records; unmeasured genes are not certified as absent.',
              'input_hashes': {name: hashlib.sha256((ROOT / 'data' / name).read_bytes()).hexdigest()
                               for name in ['charts.json', 'readout_inputs.json']}}
    if query['method'] == 'exact':
        states = enumerate_gene_states(PinClosure(p, pins), target.groups, query['max_assignments'])
        frontier = count_frontier(p, target, states)
        record.update(enumeration=states, frontier=frontier,
                      interval=exact_interval(frontier, F(query['budget'])))
    else:
        oracle = LagrangianCounts(p, pins)
        mx, mn = oracle.supported_path(target), oracle.supported_path(target, -1)
        record.update(maximum_path=mx, minimum_path=mn,
                      interval=oracle.enclosure(target, F(query['budget']), mx, mn))
    return json.loads(json.dumps(record, default=encode))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    p = sub.add_parser('interval')
    p.add_argument('--corridor', required=True)
    p.add_argument('--setting', choices=['baseline', 'maximal'], default='baseline')
    p.add_argument('--cell', type=int, required=True)
    p.add_argument('--gene', required=True)
    p.add_argument('--budget', required=True, help='Nonnegative exact rational area in square micrometres')
    p.add_argument('--method', choices=['exact', 'dual'], default='exact')
    p.add_argument('--max-assignments', type=int, default=100000)
    p.add_argument('--output', default='run/count_interval.json')
    v = sub.add_parser('verify')
    v.add_argument('input', type=Path)
    a = ap.parse_args()
    if a.command == 'verify':
        original = json.loads(a.input.read_text())
        fresh = calculate(original['query'])
        if original != fresh:
            raise ValueError('Certificate replay differs from the supplied record')
        print(json.dumps({'verified': True, 'scope': 'Full deterministic replay; cut proofs and boundaries checked during execution.'}))
        return
    q = {k: getattr(a, k) for k in ['corridor', 'setting', 'cell', 'gene', 'budget', 'method', 'max_assignments']}
    q['budget'] = str(F(q['budget']))
    result = calculate(q)
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + '\n')
    display = {k: v for k, v in result['interval'].items() if not k.endswith('_witness')}
    print(json.dumps({'query': q, 'reference_count': result['reference_count'], 'interval': display, 'record': str(out)}, indent=2))


if __name__ == '__main__':
    main()
