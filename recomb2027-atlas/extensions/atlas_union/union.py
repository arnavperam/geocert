"""Exact count spectra over a finite union of the five published atlas families.

A witness selects one complete member atlas. Local patches from different
members are never composed. The CLI consumes the checksum-pinned 1.1.0
scientific records, with the primary and four fixed-width members unchanged.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[2]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))
from atlas.binding import verify_prepared_atlas

SETTINGS = ('baseline', '0.5', '1', '2', '3')
DATASETS = ('lung', 'breast')
SCHEMA = 'geocert-atlas-union-1'


class ContractError(ValueError):
    """A required identity, complete spectrum or witness contract is violated."""


def require(condition, message):
    if not condition:
        raise ContractError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as stream:
        return json.load(stream)


def write_json(path, value):
    path = Path(path)
    if path.suffix == '.gz':
        with path.open('wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as packed:
                with io.TextIOWrapper(packed, encoding='utf-8') as stream:
                    json.dump(value, stream, separators=(',', ':'), sort_keys=True)
                    stream.write('\n')
    else:
        path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def write_csv(path, fields, rows):
    path = Path(path)
    if path.suffix == '.gz':
        raw = path.open('wb')
        packed = gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0)
        stream = io.TextIOWrapper(packed, encoding='utf-8', newline='')
    else:
        raw = packed = None
        stream = path.open('w', newline='', encoding='utf-8')
    try:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    finally:
        stream.close()
        if packed is not None:
            packed.close()
        if raw is not None:
            raw.close()


def member_root(data_root, dataset, setting):
    if setting == 'baseline':
        return data_root if dataset == 'lung' else data_root / 'external'
    return data_root / 'reach' / dataset / ('r' + setting.replace('.', 'p'))


def key_of(record):
    return record['roi'], int(record['cell']), record['gene']


def parse_minima(values, context):
    require(isinstance(values, dict) and values, f'{context}: empty or missing spectrum')
    result = {}
    for value, cost in values.items():
        require(isinstance(value, str) and str(int(value)) == value and int(value) >= 0,
                f'{context}: count must be a canonical nonnegative integer')
        c = Fraction(cost)
        require(c >= 0, f'{context}: negative area')
        result[int(value)] = c
    return result


def source_count(atlas, cell, gene):
    return atlas['source_counts'].get(str(cell), {}).get(gene, 0)


def expected_questions(atlases):
    """Every source-positive or incident editable gene must have an exact row."""
    expected = set()
    parts = defaultdict(set)
    for roi, atlas in atlases.items():
        eligible = set(atlas['complete_cells'])
        for cell, genes in atlas['source_counts'].items():
            if int(cell) in eligible:
                expected.update((roi, int(cell), gene) for gene, count in genes.items() if count > 0)
        for patch in atlas['patches']:
            for gene in {point['gene'] for point in patch['points']}:
                for cell in patch['cell_labels']:
                    if cell in eligible:
                        key = roi, cell, gene
                        expected.add(key)
                        parts[key].add(patch['id'])
    return expected, parts


def union_minima(member_spectra):
    """Exact minimum per count; deterministic member order breaks area ties."""
    values, choices = {}, {}
    for setting, spectrum in member_spectra:
        for value, cost in spectrum.items():
            if value not in values or cost < values[value]:
                values[value], choices[value] = cost, setting
    require(bool(values), 'union has no certified count state')
    return values, choices


def check_pinned_inputs(data_root):
    pins = read_json(Path(__file__).with_name('frozen_inputs.json'))
    require(pins.get('schema') == SCHEMA, 'unsupported frozen input manifest')
    files = pins['files']
    for relative, expected in files.items():
        path = (data_root / relative).resolve()
        require(path.is_relative_to(data_root), 'source path leaves data root: ' + relative)
        require(path.is_file(), 'missing source record; restore the scientific archive: ' + relative)
        require(digest(path) == expected, 'source checksum mismatch: ' + relative)
    return pins


def check_common_reference(members):
    common = None
    for member in members:
        signature = {}
        for roi, atlas in member['atlases'].items():
            signature[roi] = {key: atlas[key] for key in
                              ('source_archive_sha256', 'pixel_size_um', 'bbox_pixels',
                               'complete_cells', 'source_counts', 'records', 'source_hashes')}
        encoded = canonical(signature)
        if common is None:
            common = encoded
        require(encoded == common, 'members do not share identical reference/count observations')
    return hashlib.sha256(common).hexdigest()


def load_member(data_root, dataset, setting):
    root = member_root(data_root, dataset, setting)
    base = member_root(data_root, dataset, 'baseline')
    atlas_paths = sorted((root / 'data').glob('atlas_*.json'))
    atlases = {}
    geometry_audit = read_json(root / 'results/atlas_geometric_validation.json')
    binding = []
    for path in atlas_paths:
        atlas = read_json(path)
        roi = atlas['roi']
        require(roi not in atlases, 'duplicate atlas ROI')
        binding.append(verify_prepared_atlas(atlas, base / 'data/source_crops' / roi, geometry_audit))
        atlases[roi] = atlas
    require(bool(atlases), f'{dataset}/{setting}: no atlases')
    expected, expected_parts = expected_questions(atlases)
    patch_index = {(roi, patch['id']): patch for roi, atlas in atlases.items() for patch in atlas['patches']}
    local_rows = read_json(root / 'results/patch_spectra.json')
    local = {}
    for row in local_rows:
        key = row['roi'], row['patch'], row['gene']
        require(key not in local, 'duplicate component patch spectrum')
        require(key[:2] in patch_index, 'patch spectrum refers to absent patch')
        patch = patch_index[key[:2]]
        require(row['source_left'] == sum(point['gene'] == row['gene'] and point['source'] == 0
                                          for point in patch['points']), 'local source count mismatch')
        minima = parse_minima(row['minima'], str(key))
        require(minima.get(row['source_left']) == 0, 'local source count is not available at zero area')
        require(set(row['witnesses']) == set(row['minima']), 'incomplete local witness map')
        for value, cost in minima.items():
            witness = row['witnesses'][str(value)]
            require(Fraction(witness['physical_area_um2']) == cost, 'local witness area mismatch')
            structural_sides = list(dict.fromkeys(item[0] for item in patch.get('nuclear_constraints', [])))
            require(len(witness['owners']) == len(patch['points']) + len(structural_sides), 'local owner vector has incorrect length')
            require(witness['owners'][len(patch['points']):] == structural_sides, 'local witness moves a fixed structural group')
            require(all(owner in (0, 1) for owner in witness['owners']), 'invalid local owner')
            attained = sum(point['gene'] == row['gene'] and owner == 0
                           for point, owner in zip(patch['points'], witness['owners']))
            require(attained == value, 'local witness count mismatch')
        local[key] = {'minima': minima, 'source_left': row['source_left']}
    del local_rows
    expected_local = {(roi, patch['id'], gene) for (roi, _), patch in patch_index.items()
                      for gene in {point['gene'] for point in patch['points']}}
    require(set(local) == expected_local, 'component local spectra are incomplete or extraneous')
    spectra = {}
    traces_checked = 0
    for row in read_json(root / 'results/cell_spectra.json'):
        key = key_of(row)
        require(key not in spectra, 'duplicate cell spectrum')
        require(key in expected, 'unexpected component question')
        roi, cell, gene = key
        atlas = atlases[roi]
        require(row['atlas_id'] == atlas['atlas_id'], 'cell spectrum names another atlas')
        n = source_count(atlas, cell, gene)
        require(row['source_count'] == n, 'global source count mismatch')
        require(len(row['patches']) == len(set(row['patches'])) and
                set(row['patches']) == expected_parts[key], 'incident patch list is incomplete or duplicated')
        minima = parse_minima(row['minima'], str(key))
        require(minima.get(n) == 0, 'reference count is not available at zero area')
        # Recompose complete component spectra before taking a union. This also
        # catches an accidentally missing or extra count in a saved global row.
        composed = {n: Fraction(0)}
        for pid in row['patches']:
            patch = patch_index[(roi, pid)]
            loc = local[(roi, pid, gene)]
            sign = 1 if patch['cell_labels'][0] == cell else -1
            require(cell in patch['cell_labels'], 'patch is not incident to the question cell')
            nxt = {}
            for value, cost in composed.items():
                for left, local_cost in loc['minima'].items():
                    target = value + sign * (left - loc['source_left'])
                    candidate = cost + local_cost
                    if target not in nxt or candidate < nxt[target]:
                        nxt[target] = candidate
            composed = nxt
        require(composed == minima, 'global component spectrum does not equal complete patch composition')
        require(set(row['witness_traces']) == set(row['minima']), 'incomplete global witness map')
        for value, cost in minima.items():
            trace = row['witness_traces'][str(value)]
            require(set(trace) == set(row['patches']), 'witness omits or adds incident patch')
            require(n + sum(trace.values()) == value, 'global witness count mismatch')
            summed = Fraction(0)
            for pid, change in trace.items():
                require(type(change) is int, 'patch count change must be an integer')
                patch = patch_index[(roi, pid)]
                loc = local[(roi, pid, gene)]
                sign = 1 if patch['cell_labels'][0] == cell else -1
                left = loc['source_left'] + sign * change
                require(left in loc['minima'], 'global witness refers to absent local state')
                summed += loc['minima'][left]
            require(summed == cost, 'global witness area mismatch')
            traces_checked += 1
        row['_minima'] = minima
        spectra[key] = row
    require(set(spectra) == expected, 'component spectrum census is incomplete')
    relative = str(root.relative_to(data_root))
    spectrum_path = (root / 'results/cell_spectra.json').relative_to(data_root).as_posix()
    interval_path = (data_root / 'results/reach' / f'baseline_{dataset}_fine_intervals.csv'
                     if setting == 'baseline' else root / 'results/cell_expression_intervals.csv')
    return {'setting': setting, 'root': relative, 'atlases': atlases, 'spectra': spectra,
            'source_spectrum_path': spectrum_path, 'interval_path': interval_path,
            'binding': binding, 'global_traces_checked': traces_checked}


def spectrum_or_fixed(member, key):
    if key in member['spectra']:
        return member['spectra'][key]['_minima']
    roi, cell, gene = key
    atlas = member['atlases'][roi]
    require(cell in atlas['complete_cells'], 'implicit source count requested for an ineligible cell')
    # load_member proved all relevant editable and source-positive keys present.
    return {source_count(atlas, cell, gene): Fraction(0)}


def accumulate_marginals(member, keys, budgets, bounds):
    """Union of member-wise marginal relaxations, preserving complete families."""
    lookup = {value: i for i, value in enumerate(budgets)}
    seen = {key: 0 for key in member['spectra']}
    with member['interval_path'].open(newline='') as stream:
        for row in csv.DictReader(stream):
            key = key_of(row)
            require(key in seen, 'interval names absent component spectrum')
            budget = Fraction(row['budget_um2'])
            require(budget in lookup, 'unexpected component budget')
            index = lookup[budget]
            require(not seen[key] & (1 << index), 'duplicate component interval row')
            seen[key] |= 1 << index
            spectrum = member['spectra'][key]
            require(row['atlas_id'] == spectrum['atlas_id'], 'interval atlas identity mismatch')
            require(int(row['source_count']) == spectrum['source_count'], 'interval source count mismatch')
            feasible = [value for value, cost in spectrum['_minima'].items() if cost <= budget]
            lo, hi = min(feasible), max(feasible)
            require(row['status'] == 'exact' and (int(row['lower']), int(row['upper'])) == (lo, hi),
                    'saved component interval does not match its spectrum')
            outer_lo, outer_hi = int(row['independent_lower']), int(row['independent_upper'])
            require(0 <= outer_lo <= lo <= hi <= outer_hi, 'component marginal bound is invalid')
            require(Fraction(row['lower_cost_um2']) == spectrum['_minima'][lo] and
                    Fraction(row['upper_cost_um2']) == spectrum['_minima'][hi], 'component endpoint cost mismatch')
            prev_lo, prev_hi = bounds[key][index]
            bounds[key][index] = min(prev_lo, outer_lo), max(prev_hi, outer_hi)
    require(all(mask == (1 << len(budgets)) - 1 for mask in seen.values()),
            'component interval census is incomplete')


def build_dataset(data_root, dataset, budgets):
    members = [load_member(data_root, dataset, setting) for setting in SETTINGS]
    common_digest = check_common_reference(members)
    by_setting = {member['setting']: member for member in members}
    baseline = by_setting['baseline']
    keys = sorted(set().union(*(set(member['spectra']) for member in members)))
    member_ids = {setting: {roi: atlas['atlas_id'] for roi, atlas in by_setting[setting]['atlases'].items()}
                  for setting in SETTINGS}
    model_id = hashlib.sha256(canonical({'schema': SCHEMA, 'dataset': dataset,
                                        'reference_digest': common_digest, 'members': member_ids})).hexdigest()
    bounds = {}
    results = []
    for key in keys:
        roi, cell, gene = key
        n = source_count(baseline['atlases'][roi], cell, gene)
        bounds[key] = [(n, n) for _ in budgets]
        minima, chosen = union_minima((setting, spectrum_or_fixed(by_setting[setting], key)) for setting in SETTINGS)
        witnesses = {}
        for value, setting in sorted(chosen.items()):
            member = by_setting[setting]
            row = member['spectra'].get(key)
            witnesses[str(value)] = {
                'setting': setting, 'atlas_id': member['atlases'][roi]['atlas_id'],
                'source_spectrum_path': member['source_spectrum_path'],
                'trace': row['witness_traces'][str(value)] if row else {},
                'implicit_fixed': row is None,
            }
        results.append({'dataset': dataset, 'roi': roi, 'cell': cell, 'gene': gene,
                        'source_count': n, 'baseline_question': key in baseline['spectra'],
                        'model_id': model_id,
                        'minima': {str(value): str(cost) for value, cost in sorted(minima.items())},
                        'witnesses': witnesses})
    for member in members:
        accumulate_marginals(member, keys, budgets, bounds)
    summaries = []
    intervals = []
    previous = {}
    ever_tight = {population: set() for population in ('baseline_questions', 'union_questions')}
    ever_eligible = {population: set() for population in ever_tight}
    for bi, budget in enumerate(budgets):
        current = {}
        budget_records = []
        for result in results:
            key = key_of(result)
            minima = {int(value): Fraction(cost) for value, cost in result['minima'].items()}
            feasible = [value for value, cost in minima.items() if cost <= budget]
            lo, hi = min(feasible), max(feasible)
            outer_lo, outer_hi = bounds[key][bi]
            require(outer_lo <= lo <= hi <= outer_hi, 'union interval escapes its union marginal bound')
            baseline_values = [value for value, cost in spectrum_or_fixed(baseline, key).items() if cost <= budget]
            blo, bhi = min(baseline_values), max(baseline_values)
            require(lo <= blo <= bhi <= hi, 'union lost a baseline count endpoint')
            lw, uw = result['witnesses'][str(lo)], result['witnesses'][str(hi)]
            row = {'dataset': dataset, 'model_id': model_id, 'roi': result['roi'], 'cell': result['cell'],
                   'gene': result['gene'], 'baseline_question': result['baseline_question'],
                   'budget_um2': str(budget), 'source_count': result['source_count'],
                   'lower': lo, 'upper': hi, 'status': 'exact',
                   'lower_cost_um2': str(minima[lo]), 'upper_cost_um2': str(minima[hi]),
                   'lower_setting': lw['setting'], 'upper_setting': uw['setting'],
                   'lower_atlas_id': lw['atlas_id'], 'upper_atlas_id': uw['atlas_id'],
                   'independent_lower': outer_lo, 'independent_upper': outer_hi,
                   'baseline_lower': blo, 'baseline_upper': bhi}
            intervals.append(row)
            budget_records.append(row)
            current[key] = lo, hi
        for population in ever_tight:
            rows = [row for row in budget_records if population == 'union_questions' or row['baseline_question']]
            eligible = {key_of(row) for row in rows if row['independent_upper'] - row['independent_lower'] >= 2}
            tight = {key_of(row) for row in rows if (row['lower'], row['upper']) != (row['independent_lower'], row['independent_upper'])}
            ever_tight[population].update(tight)
            ever_eligible[population].update(eligible)
            summaries.append({'dataset': dataset, 'population': population, 'budget_um2': str(budget),
                              'questions': len(rows), 'invariant': sum(row['lower'] == row['upper'] for row in rows),
                              'non_singleton': sum(row['lower'] != row['upper'] for row in rows),
                              'invariant_fraction': sum(row['lower'] == row['upper'] for row in rows) / len(rows),
                              'baseline_invariant': sum(row['baseline_lower'] == row['baseline_upper'] for row in rows),
                              'widened_vs_baseline': sum((row['lower'], row['upper']) != (row['baseline_lower'], row['baseline_upper']) for row in rows),
                              'summed_width': sum(row['upper'] - row['lower'] for row in rows),
                              'maximum_width': max(row['upper'] - row['lower'] for row in rows),
                              'eligible_marginal_width_ge2': len(eligible), 'tightened_questions': len(tight),
                              'endpoint_changes_from_previous_budget': sum(current[key_of(row)] != previous[key_of(row)] for row in rows) if previous else 0})
        previous = current
    summary = {'dataset': dataset, 'model_id': model_id, 'member_atlas_ids': member_ids,
               'common_reference_digest': common_digest,
               'baseline_questions': len(baseline['spectra']), 'union_questions': len(keys),
               'new_questions': len(keys) - len(baseline['spectra']),
               'interval_rows': len(intervals), 'count_states': sum(len(row['minima']) for row in results),
               'component_global_spectra_checked': sum(len(member['spectra']) for member in members),
               'component_global_traces_checked': sum(member['global_traces_checked'] for member in members),
               'ever_tightened_questions': {pop: len(vals) for pop, vals in ever_tight.items()},
               'ever_eligible_marginal_width_ge2': {pop: len(vals) for pop, vals in ever_eligible.items()},
               'binding': {member['setting']: member['binding'] for member in members}}
    return results, intervals, summaries, summary, by_setting


def reference_reach(data_root, dataset, by_setting, union_records, budgets):
    union = {key_of(row): row for row in union_records}
    references = defaultdict(dict)
    with gzip.open(data_root / 'results/reach/reference_counts.csv.gz', 'rt', newline='') as stream:
        for row in csv.DictReader(stream):
            if row['dataset'] != dataset:
                continue
            key = row['reference'], *key_of(row)
            require(row['setting'] in SETTINGS, 'unexpected reference member')
            require(row['setting'] not in references[key], 'duplicate reference key/member')
            references[key][row['setting']] = row
    rows, summary = [], []
    for (reference, roi, cell, gene), records in sorted(references.items()):
        require(set(records) == set(SETTINGS), 'reference comparison does not use all five members')
        identities = {(r['source_count'], r['reference_count'], r['changed']) for r in records.values()}
        require(len(identities) == 1, 'reference comparison keys/counts differ across members')
        source, target, changed = next(iter(identities))
        source, target = int(source), int(target)
        require(changed == str(source != target), 'reference changed flag is inconsistent')
        key = roi, cell, gene
        n = source_count(by_setting['baseline']['atlases'][roi], cell, gene)
        require(source == n, 'reference exclusion guard does not match full retained source count')
        costs = []
        for setting in SETTINGS:
            cost = spectrum_or_fixed(by_setting[setting], key).get(target)
            recorded = records[setting]['minimum_area_um2']
            require((cost is None and recorded == 'inf') or (cost is not None and recorded != 'inf' and cost == Fraction(recorded)),
                    'reference minimum differs from its complete component spectrum')
            if cost is not None:
                costs.append((cost, SETTINGS.index(setting), setting))
        best = min(costs) if costs else None
        union_cost = Fraction(union[key]['minima'][str(target)]) if key in union and str(target) in union[key]['minima'] else (Fraction(0) if key not in union and target == n else None)
        require(union_cost == (best[0] if best else None), 'union reference reach is not exact count membership')
        baseline = records['baseline']['minimum_area_um2']
        rows.append({'dataset': dataset, 'reference': reference, 'roi': roi, 'cell': cell, 'gene': gene,
                     'source_count': source, 'reference_count': target, 'changed': source != target,
                     'minimum_area_um2': str(best[0]) if best else 'inf',
                     'setting': best[2] if best else '',
                     'atlas_id': by_setting[best[2]]['atlases'][roi]['atlas_id'] if best else '',
                     'baseline_minimum_area_um2': baseline})
    for reference in sorted({row['reference'] for row in rows}):
        changed = [row for row in rows if row['reference'] == reference and row['changed']]
        for budget in [*budgets, None]:
            attainable = lambda value: value != 'inf' and (budget is None or Fraction(value) <= budget)
            baseline_n = sum(attainable(row['baseline_minimum_area_um2']) for row in changed)
            union_n = sum(attainable(row['minimum_area_um2']) for row in changed)
            require(union_n >= baseline_n, 'union reduced reference count reach')
            summary.append({'dataset': dataset, 'reference': reference,
                            'budget_um2': str(budget) if budget is not None else 'unbounded',
                            'changed_count_keys': len(changed), 'baseline_attainable': baseline_n,
                            'union_attainable': union_n, 'additional_attainable': union_n - baseline_n,
                            'baseline_reach': baseline_n / len(changed), 'union_reach': union_n / len(changed)})
    return rows, summary


def check_output_location(data_root, destination):
    require(not destination.is_relative_to(data_root) and not data_root.is_relative_to(destination),
            'output directory must be separate from the immutable source tree')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True,
                        help='Expanded and restored immutable GeoCert 1.1.0 scientific package')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    data_root, destination = args.data_root.resolve(), args.output_dir.resolve()
    require(data_root.is_dir(), 'data root is missing')
    check_output_location(data_root, destination)
    start = time.perf_counter()
    pins = check_pinned_inputs(data_root)
    budgets = [Fraction(value) for value in read_json(data_root / 'results/reach/plan_identity.json')['budgets_um2']]
    require(len(budgets) == len(set(budgets)) and budgets == sorted(budgets) and budgets[0] == 0,
            'budgets must be unique, sorted, nonnegative and include zero')
    all_spectra, all_intervals, budget_rows, datasets, reference_rows, reference_summary = [], [], [], [], [], []
    for dataset in DATASETS:
        result, intervals, summaries, summary, members = build_dataset(data_root, dataset, budgets)
        rr, rs = reference_reach(data_root, dataset, members, result, budgets)
        all_spectra.extend(result); all_intervals.extend(intervals); budget_rows.extend(summaries)
        datasets.append(summary); reference_rows.extend(rr); reference_summary.extend(rs)
        print(json.dumps({'dataset': dataset, 'questions': summary['union_questions'],
                          'interval_rows': summary['interval_rows'], 'elapsed_seconds': time.perf_counter() - start}), flush=True)
        del members
    destination.mkdir(parents=True, exist_ok=True)
    write_json(destination / 'union_spectra.json.gz', all_spectra)
    write_csv(destination / 'intervals.csv.gz', list(all_intervals[0]), all_intervals)
    write_csv(destination / 'budget_summary.csv', list(budget_rows[0]), budget_rows)
    write_json(destination / 'budget_summary.json', budget_rows)
    write_csv(destination / 'reference_counts.csv.gz', list(reference_rows[0]), reference_rows)
    write_json(destination / 'reference_summary.json', reference_summary)
    write_csv(destination / 'reference_summary.csv', list(reference_summary[0]), reference_summary)
    manifest = {'schema': SCHEMA, 'scientific_source_commit': pins['scientific_source_commit'],
                'files': pins['files'], 'budgets_um2': list(map(str, budgets)), 'settings': list(SETTINGS),
                'model_ids': {row['dataset']: row['model_id'] for row in datasets}}
    write_json(destination / 'input_manifest.json', manifest)
    summary = {'schema': SCHEMA, 'settings': list(SETTINGS), 'budgets_um2': list(map(str, budgets)),
               'scope': 'Finite union of five complete atlas families on the same reference and retained observations; every witness selects one whole member.',
               'independent_baseline': 'Envelope of the five member-wise independent-molecule outer intervals at the same budget.',
               'population_definition': 'Union census contains every source-positive or editable-in-any-member cell-gene key; baseline_questions retains the original primary census.',
               'datasets': datasets, 'interval_rows': len(all_intervals), 'count_states': sum(len(row['minima']) for row in all_spectra),
               'seconds': time.perf_counter() - start,
               'outputs': {name: digest(destination / name) for name in ('union_spectra.json.gz', 'intervals.csv.gz',
                           'budget_summary.csv', 'budget_summary.json', 'reference_counts.csv.gz',
                           'reference_summary.json', 'reference_summary.csv', 'input_manifest.json')}}
    write_json(destination / 'summary.json', summary)
    print(json.dumps({'status': 'passed', 'questions': len(all_spectra), 'interval_rows': len(all_intervals),
                      'seconds': summary['seconds']}), flush=True)


if __name__ == '__main__':
    main()
