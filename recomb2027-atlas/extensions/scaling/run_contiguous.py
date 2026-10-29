"""Measure and validate the unchanged pipeline on a prepared contiguous crop."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from benchmark import STAGES, sha


def main():
    if not __debug__:
        raise RuntimeError('Run without Python -O; all result checks are required')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frozen-code-root', type=Path, required=True)
    parser.add_argument('--prepared-root', type=Path, required=True)
    parser.add_argument('--total-timeout-seconds', type=int, default=2400)
    parser.add_argument('--count-method', choices=['chain', 'enumeration'], default='chain')
    parser.add_argument('--enumeration-copy-limit', type=int, default=16)
    parser.add_argument('--reuse-validated-input', type=Path,
                        help='Fresh count-only attempt after a saved whole-run timeout; original inputs and validation remain unchanged')
    args = parser.parse_args()
    source, run = args.frozen_code_root.resolve(), args.prepared_root.resolve()
    previous = None
    if args.reuse_validated_input:
        original = args.reuse_validated_input.resolve()
        previous = json.loads((original / 'measurement.json').read_text())
        if previous['status'] != 'failed' or 'TimeoutExpired' not in previous.get('error', ''):
            raise ValueError('Reuse requires a preserved failed whole-run timeout')
        if previous.get('count_method') != args.count_method:
            raise ValueError('A count-only retry must preserve the exact backend')
        if any(previous['source_scripts_sha256'][name] != sha(source / path) for name, path in STAGES):
            raise ValueError('A count-only retry must preserve all scientific source scripts')
        if run.exists() and any(run.iterdir()):
            raise ValueError('A count-only retry requires a fresh destination')
        run.mkdir(parents=True, exist_ok=True)
        (run / 'data').symlink_to(original / 'data', target_is_directory=True)
        (run / 'docs').symlink_to(original / 'docs', target_is_directory=True)
        (run / 'results').mkdir()
        (run / 'previous_attempt').mkdir()
        for name in ['atlas_geometric_validation.json', 'atlas_extraction.json']:
            shutil.copyfile(original / 'results' / name, run / 'results' / name)
        shutil.copyfile(original / 'preparation_summary.json', run / 'preparation_summary.json')
        for path in original.iterdir():
            if path.is_file() and path.suffix in {'.json', '.log'}:
                shutil.copyfile(path, run / 'previous_attempt' / path.name)
        for name in ['extraction', 'geometry_validation']:
            if previous['stages'][name]['status'] != 'passed':
                raise ValueError('Both original geometric stages must have passed')
        amendment = {'recorded_unix': time.time(), 'reason': 'One additional count-only attempt after the original combined wall-time cap',
                     'maximum_additional_seconds': args.total_timeout_seconds,
                     'input_and_geometry': 'Unchanged source files and atlas, bound to the copied original zero-failure rational audit',
                     'backend': args.count_method, 'original_attempt_preserved': 'previous_attempt/measurement.json',
                     'reporting': 'Preserve the failed attempt and report cumulative measured attempt wall time; completed stage sum alone excludes the interrupted attempt'}
        (run / 'retry_amendment.json').write_text(json.dumps(amendment, indent=2)+'\n')
    if (run / 'measurement.json').exists():
        raise ValueError('Measurement already exists; use a new prepared directory')
    env = os.environ.copy()
    env.pop('GEOCERT_ROI', None)
    env.update(GEOCERT_DATA_ROOT=str(run), GEOCERT_RESULT_DIR=str(run / 'results'), GEOCERT_COUNT_METHOD=args.count_method,
               PYTHONHASHSEED='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    started = time.perf_counter()
    record = {'status': 'running', 'stages': {}, 'scope': 'One genuine central contiguous breast format-test crop',
              'execution_environment': 'Shared Linux container; concurrent research work; each stage is a fresh sequential process',
              'timing_semantics': 'Stage wall time includes startup, imports, computation and output serialization; source preparation and final completeness verification are separate',
              'memory_semantics': 'Maximum stage process RSS in MiB, not aggregate shared-machine memory',
              'source_scripts_sha256': {name: sha(source / path) for name, path in STAGES},
              'count_method': args.count_method,
              'runner_sha256': sha(Path(__file__)), 'stage_wrapper_sha256': sha(Path(__file__).with_name('stage_runner.py')),
              'total_timeout_seconds': args.total_timeout_seconds}
    if previous:
        record['stages'] = {name: previous['stages'][name] for name in ['extraction', 'geometry_validation']}
        record['previous_attempt'] = previous
        record['global_deadline_interruptions'] = 1
        record['pipeline_seconds_semantics'] = 'Sum of the completed extraction, geometry and final count stages; interrupted attempt is additional and included in cumulative_measured_attempt_seconds'
    try:
        stages_to_run = STAGES[2:] if previous else STAGES
        for name, script in stages_to_run + [('coverage', 'analysis/coverage.py')]:
            if name == 'count_spectra' and args.count_method == 'enumeration':
                inputs = [json.loads(p.read_text()) for p in (run / 'data').glob('atlas_*.json')]
                maximum = max((number for atlas in inputs for patch in atlas['patches']
                               for number in Counter(point['gene'] for point in patch['points']).values()), default=0)
                decision = {'selected_unix': time.time(), 'method': 'enumeration',
                            'maximum_local_gene_copies': maximum, 'permitted_maximum': args.enumeration_copy_limit,
                            'pattern_limit': 65536,
                            'rule': 'Select the existing exact enumeration backend from local input multiplicity after geometry validation and before count outputs; require at most 16 copies per local gene',
                            'atlas_ids': [a['atlas_id'] for a in inputs]}
                (run / 'backend_decision.json').write_text(json.dumps(decision, indent=2)+'\n')
                record['backend_decision'] = decision
                if not 0 <= maximum <= min(args.enumeration_copy_limit, 16):
                    raise RuntimeError('Enumeration input exceeds the predeclared copy limit; count stage unrun')
                del inputs
            remaining = args.total_timeout_seconds - (time.perf_counter() - started)
            if remaining <= 0:
                raise TimeoutError('Combined pipeline time budget exhausted')
            stamp = time.perf_counter()
            destination = run / f'{name}_resource.json'
            with (run / f'{name}.log').open('w') as log:
                proc = subprocess.run([sys.executable, str(Path(__file__).with_name('stage_runner.py')),
                                       '--script', str(source / script), '--record', str(destination)],
                                      env=env, cwd=source, stdout=log, stderr=subprocess.STDOUT,
                                      timeout=remaining, check=False)
            phase = json.loads(destination.read_text())
            phase['elapsed_seconds'] = time.perf_counter()-stamp
            record['stages'][name] = phase
            (run / 'measurement_progress.json').write_text(json.dumps(record, indent=2)+'\n')
            print(json.dumps({'stage': name, 'seconds': phase['elapsed_seconds'],
                              'rss_mib': phase['peak_process_rss_kib']/1024}), flush=True)
            if proc.returncode or phase['status'] != 'passed':
                raise RuntimeError(name + ' failed')
        atlases = [json.loads(p.read_text()) for p in (run / 'data').glob('atlas_*.json')]
        assert len(atlases) == 1
        atlas = atlases[0]
        audit = json.loads((run / 'results/atlas_geometric_validation.json').read_text())
        summary = json.loads((run / 'results/cell_expression_summary.json').read_text())
        assert all(r['count_method'] == args.count_method for r in summary['rois'])
        spectra = json.loads((run / 'results/cell_spectra.json').read_text())
        with (run / 'results/cell_expression_intervals.csv').open() as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == sum(r['interval_rows'] for r in summary['rois']), 'Incomplete CSV serialization'
        assert len(spectra) == sum(r['positive_or_editable_cell_gene_keys'] for r in summary['rois'])
        assert len({(r['roi'], r['cell'], r['gene'], r['budget_um2']) for r in rows}) == len(rows)
        assert audit['failures'] == 0 and all(r['status'] == 'exact' for r in rows)
        coverage = json.loads((run / 'results/coverage_summary.json').read_text())[0]
        record.update(status='passed', roi=atlas['roi'], atlas_id=atlas['atlas_id'],
                      complete_cells=len(atlas['complete_cells']), retained_transcripts=len(atlas['records']),
                      ribbons=len(atlas['patches']), count_spectra=len(spectra), interval_rows=len(rows),
                      interval_statuses=dict(Counter(r['status'] for r in rows)), geometry_failures=0,
                      resource_limited_queries=0, editable_transcripts=audit['editable_RNA'],
                      coverage=coverage, count_summary=summary['rois'][0],
                      pipeline_seconds=sum(record['stages'][name]['elapsed_seconds'] for name, _ in STAGES),
                      largest_stage_rss_mib=max(r['peak_process_rss_kib']/1024 for r in record['stages'].values()),
                      output_sha256={name: sha(run / 'results' / name) for name in
                                     ['cell_spectra.json', 'patch_spectra.json', 'cell_expression_intervals.csv',
                                      'atlas_geometric_validation.json', 'coverage_summary.json']})
    except BaseException as exc:
        record['status'] = 'failed'
        record['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        record['outer_elapsed_seconds'] = time.perf_counter()-started
        if previous:
            record['cumulative_measured_attempt_seconds'] = previous['outer_elapsed_seconds'] + record['outer_elapsed_seconds']
        (run / 'measurement.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ['stages', 'coverage', 'count_summary']}), flush=True)


if __name__ == '__main__':
    main()
