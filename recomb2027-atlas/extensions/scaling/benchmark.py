"""Measure native-crop workload scaling from masks to checked count intervals.

The four existing 200-micrometre crops are processed in fixed increasing
prefixes. Each repeat uses a fresh directory and three fresh processes:
source-mask extraction, independent rational geometry validation, and exact
count spectra with all witnesses. This is an existing-data workload benchmark,
not a larger contiguous specimen or a whole-slide measurement.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

ORDER = [('lung', 'roi_01'), ('lung', 'roi_02'),
         ('breast', 'breast_01'), ('breast', 'breast_02')]
STAGES = [('extraction', 'analysis/build_atlas.py'),
          ('geometry_validation', 'analysis/verify_atlas.py'),
          ('count_spectra', 'analysis/certify_cells.py')]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def base_for(root: Path, dataset: str) -> Path:
    return root if dataset == 'lung' else root / 'external'


def prepare(root: Path, output: Path, prefix: int) -> list[dict]:
    (output / 'data/source_crops').mkdir(parents=True)
    (output / 'results').mkdir()
    (output / 'docs').mkdir()
    (output / 'docs/ATLAS_PROTOCOL.json').symlink_to((root / 'docs/ATLAS_PROTOCOL.json').resolve())
    combined = {'scope': 'Fixed prefix of the four existing primary crops', 'rois': {}}
    members = []
    for dataset, roi in ORDER[:prefix]:
        base = base_for(root, dataset)
        data = json.loads((base / 'data/observed_crop_counts.json').read_text())
        combined['rois'][roi] = data['rois'][roi]
        (output / 'data/source_crops' / roi).symlink_to((base / 'data/source_crops' / roi).resolve(), target_is_directory=True)
        atlas = json.loads((base / 'data' / f'atlas_{roi}.json').read_text())
        members.append({'dataset': dataset, 'roi': roi, 'expected_atlas_id': atlas['atlas_id'],
                        'retained_transcripts': len(atlas['records']),
                        'complete_cells': len(atlas['complete_cells']),
                        'ribbons': len(atlas['patches']),
                        'source_hashes': atlas['source_hashes']})
    (output / 'data/observed_crop_counts.json').write_text(json.dumps(combined, separators=(',', ':')) + '\n')
    return members


def check_run(root: Path, run: Path, members: list[dict]) -> dict:
    rois = {member['roi'] for member in members}
    for member in members:
        atlas = json.loads((run / 'data' / f"atlas_{member['roi']}.json").read_text())
        if atlas['atlas_id'] != member['expected_atlas_id']:
            raise ValueError('Re-extracted atlas differs from frozen input: ' + member['roi'])
        original = base_for(root, member['dataset']) / 'data' / f"atlas_{member['roi']}.json"
        if atlas != json.loads(original.read_text()):
            raise ValueError('Re-extracted atlas body differs: ' + member['roi'])
    expected_spectra, expected_intervals = {}, {}
    for dataset in sorted({member['dataset'] for member in members}):
        base = base_for(root, dataset)
        for row in json.loads((base / 'results/cell_spectra.json').read_text()):
            if row['roi'] in rois:
                expected_spectra[(row['roi'], row['cell'], row['gene'])] = row['minima']
        with (base / 'results/cell_expression_intervals.csv').open() as handle:
            for row in csv.DictReader(handle):
                if row['roi'] in rois:
                    expected_intervals[(row['roi'], row['cell'], row['gene'], row['budget_um2'])] = row
    actual = json.loads((run / 'results/cell_spectra.json').read_text())
    actual_spectra = {(row['roi'], row['cell'], row['gene']): row['minima'] for row in actual}
    if actual_spectra != expected_spectra:
        raise ValueError('Recomputed spectra differ from frozen exact minima')
    with (run / 'results/cell_expression_intervals.csv').open() as handle:
        actual_intervals = {(row['roi'], row['cell'], row['gene'], row['budget_um2']): row for row in csv.DictReader(handle)}
    if actual_intervals != expected_intervals:
        raise ValueError('Recomputed interval rows differ from frozen rows')
    audit = json.loads((run / 'results/atlas_geometric_validation.json').read_text())
    if audit['failures']:
        raise ValueError('Rational geometry validation failed')
    result = json.loads((run / 'results/cell_expression_summary.json').read_text())
    if any(row['count_method'] != 'chain' for row in result['rois']):
        raise ValueError('Count benchmark did not use the requested production method')
    coverage_rows = []
    for dataset in sorted({member['dataset'] for member in members}):
        coverage_rows.extend(row for row in json.loads((base_for(root, dataset) / 'results/coverage_summary.json').read_text())
                             if row['roi'] in rois)
    return {'atlas_identity_matches': len(members), 'count_spectra_match': len(actual_spectra),
            'interval_rows_match': len(actual_intervals),
            'geometry_failures': audit['failures'],
            'resource_limited_queries': 0,
            'resource_limit_semantics': 'Incomplete local spectra raise and fail the run before exact output is accepted',
            'editable_transcripts': audit['editable_RNA'],
            'editable_edges': sum(row['atlas_editable_unique_edges'] for row in coverage_rows),
            'incident_edges': sum(row['all_incident_unique_edges'] for row in coverage_rows),
            'chain_queries': sum(row['chain_queries'] for row in result['rois']),
            'chain_dp_states': sum(row['chain_DP_states'] for row in result['rois']),
            'singleton_questions': sum(row['singleton_questions'] for row in result['rois']),
            'local_gene_spectra': sum(row['local_gene_spectra'] for row in result['rois']),
            'output_sha256': {name: sha(run / 'results' / name) for name in
                              ['cell_spectra.json', 'patch_spectra.json', 'cell_expression_intervals.csv',
                               'atlas_geometric_validation.json']}}


def summarize(output: Path, records: list[dict], expected_repeats: int) -> dict:
    grouped = defaultdict(list)
    for record in records:
        grouped[record['prefix_crops']].append(record)
    rows = []
    for prefix, runs in sorted(grouped.items()):
        passed = [r for r in runs if r['status'] == 'passed']
        row = {'prefix_crops': prefix, 'completed_repeats': len(passed),
               'requested_repeats': expected_repeats, 'failed_repeats': len(runs) - len(passed)}
        if passed:
            first = passed[0]
            for name in ['complete_cells', 'retained_transcripts', 'ribbons']:
                row[name] = first[name]
            for name in ['count_spectra_match', 'interval_rows_match', 'editable_transcripts', 'editable_edges',
                         'incident_edges', 'local_gene_spectra', 'chain_queries', 'chain_dp_states']:
                row[name] = first['verification'][name]
            row['edge_coverage'] = row['editable_edges'] / row['incident_edges']
            for stage, _ in STAGES:
                row[stage + '_median_seconds'] = statistics.median(r['stages'][stage]['elapsed_seconds'] for r in passed)
            totals = [sum(r['stages'][stage]['elapsed_seconds'] for stage, _ in STAGES) for r in passed]
            row['pipeline_median_seconds'] = statistics.median(totals)
            row['pipeline_min_seconds'], row['pipeline_max_seconds'] = min(totals), max(totals)
            row['largest_stage_peak_rss_mib'] = max(r['stages'][stage]['peak_process_rss_kib'] / 1024
                                                   for r in passed for stage, _ in STAGES)
            row['resource_limited_queries'] = sum(r['verification']['resource_limited_queries'] for r in passed)
        rows.append(row)
    result = {'scope': 'Sequential workload scaling over existing disjoint native crops; no added tissue or whole-slide claim',
              'workloads': rows, 'failed_runs': sum(r['status'] != 'passed' for r in records),
              'timing_semantics': 'Fresh-process wall time, including interpreter startup/imports and stage output serialization; staging and final frozen-result comparison excluded',
              'memory_semantics': 'Largest individual stage process RSS, including imports; not total concurrent machine memory',
              'geometry': 'Every repeated extraction must equal the original atlas, including its digest',
              'execution_environment': 'Shared container; other research tasks may run concurrently. Stages and repeats here run sequentially.'}
    (output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    fields = sorted({key for row in rows for key in row})
    with (output / 'summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--stage-timeout-seconds', type=int, default=900)
    parser.add_argument('--prefixes', type=int, nargs='+', default=[1, 2, 3, 4])
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Run without Python -O; source verification uses assertions')
    if args.repeats < 1 or any(p not in range(1, 5) for p in args.prefixes):
        raise ValueError('Positive repetitions and prefixes1–4 required')
    root, output = args.data_root.resolve(), args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Choose an empty output directory')
    output.mkdir(parents=True, exist_ok=True)
    cpu_max = Path('/sys/fs/cgroup/cpu.max')
    protocol = {'created_unix': time.time(), 'crop_order': ORDER, 'prefixes': args.prefixes,
                'repeats': args.repeats, 'stage_timeout_seconds': args.stage_timeout_seconds,
                'scope': 'Existing real disjoint200µmcrops; all biological inputs unchanged',
                'python': sys.version, 'platform': platform.platform(),
                'machine_cpu_count': os.cpu_count(), 'cgroup_cpu_max': cpu_max.read_text().strip() if cpu_max.exists() else None,
                'source_scripts_sha256': {name: sha(root / path) for name, path in STAGES},
                'benchmark_sha256': sha(Path(__file__)),
                'selection': 'Fixed ordered prefixes; no outcome-based patch selection or duplicated/translated crops'}
    (output / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
    records = []
    for prefix in args.prefixes:
        for repeat in range(1, args.repeats + 1):
            run = output / f'crops{prefix}_repeat{repeat}'
            members = prepare(root, run, prefix)
            env = os.environ.copy()
            env.pop('GEOCERT_ROI', None)
            env.update(GEOCERT_DATA_ROOT=str(run), GEOCERT_RESULT_DIR=str(run / 'results'), GEOCERT_COUNT_METHOD='chain',
                       PYTHONHASHSEED='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
            record = {'prefix_crops': prefix, 'repeat': repeat, 'members': members,
                      'status': 'running', 'stages': {}, 'complete_cells': sum(m['complete_cells'] for m in members),
                      'retained_transcripts': sum(m['retained_transcripts'] for m in members),
                      'ribbons': sum(m['ribbons'] for m in members)}
            try:
                for name, script in STAGES:
                    destination = run / f'{name}_resource.json'
                    begin = time.perf_counter()
                    with (run / f'{name}.log').open('w') as log:
                        process = subprocess.run([sys.executable, str(Path(__file__).with_name('stage_runner.py')),
                                                  '--script', str(root / script), '--record', str(destination)],
                                                 env=env, cwd=root, stdout=log, stderr=subprocess.STDOUT,
                                                 timeout=args.stage_timeout_seconds, check=False)
                    phase = json.loads(destination.read_text())
                    phase['elapsed_seconds'] = time.perf_counter() - begin
                    record['stages'][name] = phase
                    if process.returncode or phase['status'] != 'passed':
                        raise RuntimeError(name + ' failed; inspect the preserved log')
                    print(json.dumps({'prefix_crops': prefix, 'repeat': repeat, 'stage': name,
                                      'seconds': phase['elapsed_seconds'], 'peak_rss_mib': phase['peak_process_rss_kib'] / 1024}), flush=True)
                record['verification'] = check_run(root, run, members)
                record['status'] = 'passed'
            except BaseException as exc:
                record['status'] = 'failed'
                record['error'] = f'{type(exc).__name__}: {exc}'
                (run / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
                records.append(record)
                summarize(output, records, args.repeats)
                raise
            record['orchestrator_peak_rss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            (run / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
            records.append(record)
            summary = summarize(output, records, args.repeats)
            print(json.dumps({'completed': [prefix, repeat], 'verification': record['verification']}), flush=True)
    (output / 'runs.json').write_text(json.dumps(records, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
