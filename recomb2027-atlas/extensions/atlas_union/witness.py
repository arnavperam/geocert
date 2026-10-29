"""Materialize and check an attaining boundary for one finite-union count."""
import argparse
import json
import hashlib
from fractions import Fraction
from pathlib import Path

import union
from atlas.model import verify_witness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--results-dir', type=Path, required=True)
    parser.add_argument('--dataset', choices=union.DATASETS, required=True)
    parser.add_argument('--roi', required=True)
    parser.add_argument('--cell', type=int, required=True)
    parser.add_argument('--gene', required=True)
    parser.add_argument('--count', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data_root = args.data_root.resolve()
    union.require(not args.output.resolve().is_relative_to(data_root), 'witness output must be outside immutable inputs')
    union.require(args.output.resolve() != (args.results_dir / 'union_spectra.json.gz').resolve(), 'witness output would replace a union input')
    union.check_pinned_inputs(data_root)
    key = args.roi, args.cell, args.gene
    result_rows = union.read_json(args.results_dir / 'union_spectra.json.gz')
    rows = [row for row in result_rows if row['dataset'] == args.dataset and union.key_of(row) == key]
    union.require(len(rows) == 1, 'requested cell-gene key is absent or duplicated')
    result = rows[0]
    value = str(args.count)
    union.require(value in result['minima'], 'requested count is unattainable in this union family')
    claimed = result['witnesses'][value]
    costs = []
    members = []
    member_ids = {}
    selected_row = None
    selected_atlas = None
    for position, setting in enumerate(union.SETTINGS):
        root = union.member_root(data_root, args.dataset, setting)
        atlases = {item['roi']: item for item in (union.read_json(p) for p in sorted((root / 'data').glob('atlas_*.json')))}
        atlas = atlases[args.roi]
        members.append({'atlases': atlases})
        member_ids[setting] = {roi: item['atlas_id'] for roi, item in atlases.items()}
        candidates = [row for row in union.read_json(root / 'results/cell_spectra.json') if union.key_of(row) == key]
        union.require(len(candidates) <= 1, 'duplicate component key')
        if candidates:
            cost = candidates[0]['minima'].get(value)
        else:
            # The pinned complete census permits implicit fixed source rows.
            source = union.source_count(atlas, args.cell, args.gene)
            cost = '0' if source == args.count else None
        if cost is not None:
            costs.append((Fraction(cost), position, setting))
        if setting == claimed['setting']:
            selected_row = candidates[0] if candidates else None
            selected_atlas = atlas
    reference_digest = union.check_common_reference(members)
    model_id = hashlib.sha256(union.canonical({'schema': union.SCHEMA, 'dataset': args.dataset,
                    'reference_digest': reference_digest, 'members': member_ids})).hexdigest()
    union.require(result['model_id'] == model_id, 'union model identity mismatch')
    union.require(bool(costs), 'count has no component witness')
    best = min(costs)
    union.require(best[0] == Fraction(result['minima'][value]) and best[2] == claimed['setting'],
                  'selected witness is not the deterministic exact union minimum')
    union.require(selected_atlas['atlas_id'] == claimed['atlas_id'], 'witness atlas identity mismatch')
    trace = selected_row['witness_traces'][value] if selected_row else {}
    union.require(trace == claimed['trace'], 'union witness differs from its component trace')
    root = union.member_root(data_root, args.dataset, claimed['setting'])
    union.require(claimed['source_spectrum_path'] == (root / 'results/cell_spectra.json').relative_to(data_root).as_posix(), 'source spectrum path mismatch')
    patches = {patch['id']: patch for patch in selected_atlas['patches']}
    local = {(row['roi'], row['patch'], row['gene']): row for row in union.read_json(root / 'results/patch_spectra.json')}
    output_patches = []
    area = Fraction(0)
    for pid, delta in trace.items():
        patch = patches[pid]
        state = local[(args.roi, pid, args.gene)]
        sign = 1 if patch['cell_labels'][0] == args.cell else -1
        left = state['source_left'] + sign * delta
        witness = state['witnesses'][str(left)]
        physical = verify_witness(patch, witness['owners'], witness['boundaries'], witness['physical_area_um2'])
        area += physical
        output_patches.append({'patch_id': pid, 'cell_count_change': delta, 'left_gene_count': left,
                               'physical_area_um2': str(physical), 'owners': witness['owners'],
                               'boundaries': witness['boundaries'], 'patch': patch})
    source = union.source_count(selected_atlas, args.cell, args.gene)
    union.require(source + sum(trace.values()) == args.count and area == best[0],
                  'materialized boundary does not attain the requested count and area')
    output = {'schema': union.SCHEMA, 'dataset': args.dataset, 'roi': args.roi, 'cell': args.cell,
              'gene': args.gene, 'source_count': source, 'attained_count': args.count,
              'minimum_area_um2': str(area), 'setting': claimed['setting'],
              'atlas_id': claimed['atlas_id'], 'model_id': result['model_id'],
              'member_atlas_path': (root / 'data' / ('atlas_' + args.roi + '.json')).relative_to(data_root).as_posix(),
              'source_spectrum_path': claimed['source_spectrum_path'],
              'scope': 'One complete member atlas. All patches absent from this trace stay at their zero-displacement source.',
              'patch_witnesses': output_patches, 'continuous_witness_check': 'passed'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    union.write_json(args.output, output)
    print(json.dumps({'status': 'passed', 'atlas_id': claimed['atlas_id'], 'setting': claimed['setting'],
                      'attained_count': args.count, 'minimum_area_um2': str(area), 'patches': len(output_patches)}))


if __name__ == '__main__':
    main()
