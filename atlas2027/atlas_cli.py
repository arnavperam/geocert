#!/usr/bin/env python3
"""Reproduce the exact fixed-seam atlas analysis from prepared source inputs."""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def manifest() -> None:
    manifest_path = ROOT / 'SHA256SUMS.txt'
    if not manifest_path.is_file():
        raise FileNotFoundError('SHA256SUMS.txt is missing from the release')
    count = 0
    for line in manifest_path.read_text().splitlines():
        expected, relative = line.split('  ', 1)
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError('manifest path leaves the release')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'checksum mismatch: {relative}')
        count += 1
    print(json.dumps({'manifest_entries': count, 'status': 'passed'}))


def showcase(output: Path) -> None:
    from atlas.model import verify_witness
    base = ROOT / 'external'
    a = json.loads((base / 'data/atlas_breast_01.json').read_text())
    patches = {p['id']: p for p in a['patches']}
    local = {(r['patch'], r['gene']): r for r in json.loads((base / 'results/patch_spectra.json').read_text()) if r['roi'] == 'breast_01'}
    s = next(r for r in json.loads((base / 'results/cell_spectra.json').read_text()) if r['roi'] == 'breast_01' and r['cell'] == 81 and r['gene'] == 'LDHB')
    minimum = Fraction(s['minima']['17'])
    area = Fraction(0)
    moved: dict[str, dict] = {}
    records = []
    for pid, delta in s['witness_traces']['17'].items():
        p = patches[pid]; rec = local[(pid, 'LDHB')]
        sign = 1 if p['cell_labels'][0] == 81 else -1
        left_count = rec['source_left'] + sign * delta
        w = rec['witnesses'][str(left_count)]
        area += verify_witness(p, w['owners'], w['boundaries'], w['physical_area_um2'])
        for point, owner in zip(p['points'], w['owners']):
            if owner != point['source']:
                ident = point['transcript_id']
                if ident in moved:
                    raise ValueError('one RNA appears in two patch changes')
                moved[ident] = {'gene': point['gene'], 'from_cell': p['cell_labels'][point['source']], 'to_cell': p['cell_labels'][owner]}
        records.append({'patch': pid, 'count_change': delta, 'witness': w})
    source_lookup = {r['id']: r for r in a['records']}
    count = s['source_count']
    for ident, r in moved.items():
        original = source_lookup[ident]
        if original['gene'] == 'LDHB':
            count += int(r['to_cell'] == 81) - int(original['cell'] == 81)
    if count != 17 or area != minimum:
        raise ValueError('showcase count or minimum area mismatch')
    intervals = []
    for budget in [Fraction(1, 4), Fraction(1, 2)]:
        levels = [int(v) for v, d in s['minima'].items() if Fraction(d) <= budget]
        intervals.append({'budget_um2': str(budget), 'lower': min(levels), 'upper': max(levels)})
    result = {'atlas_id': a['atlas_id'], 'roi': 'breast_01', 'cell': 81, 'gene': 'LDHB',
              'source_count': s['source_count'], 'attained_count': count,
              'minimum_area_um2': str(minimum), 'intervals': intervals,
              'changed_molecules': moved, 'patch_witnesses': records,
              'scope': 'Recorded exact spectrum and direct physical witness check; atlas geometry and full optimality have separate validation commands. External format-test data, not biological truth.'}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'count': count, 'area_um2': str(minimum), 'changed_RNA': len(moved), 'output': str(output)}))


def run(script: str, env: dict | None = None) -> None:
    subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT, env=env, check=True)


def reproduce(dataset: str, output: Path) -> None:
    base = ROOT if dataset == 'lung' else ROOT / 'external'
    output = output.resolve()
    if output == (base / 'results').resolve():
        raise ValueError('choose a fresh output directory, not the stored results')
    env = os.environ.copy()
    env.update(GEOCERT_DATA_ROOT=str(base), GEOCERT_RESULT_DIR=str(output), GEOCERT_COUNT_METHOD='chain')
    run('analysis/certify_cells.py', env)
    expected = list(csv.DictReader((base / 'results/cell_expression_intervals.csv').open()))
    actual = list(csv.DictReader((output / 'cell_expression_intervals.csv').open()))
    if expected != actual:
        raise ValueError('regenerated interval table differs from the stored exact table')
    a = json.loads((base / 'results/cell_spectra.json').read_text())
    b = json.loads((output / 'cell_spectra.json').read_text())
    keys = lambda x: (x['roi'], x['cell'], x['gene'])
    expected_spectra = {keys(x): x['minima'] for x in a}
    actual_spectra = {keys(x): x['minima'] for x in b}
    if expected_spectra != actual_spectra:
        raise ValueError('regenerated exact minimum-area spectra differ')
    print(json.dumps({'dataset': dataset, 'spectra': len(a), 'intervals': len(actual), 'status': 'exact match'}))


def main() -> None:
    if not __debug__:
        raise RuntimeError('Scientific checks require assertions: do not run Python with -O')
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    sub.add_parser('verify-manifest')
    show = sub.add_parser('showcase'); show.add_argument('--output', type=Path, default=Path('run/showcase.json'))
    sub.add_parser('validate')
    rep = sub.add_parser('reproduce'); rep.add_argument('--dataset', choices=['lung', 'breast'], required=True);rep.add_argument('--output-dir',type=Path,required=True)
    args = ap.parse_args()
    if args.command == 'verify-manifest':manifest()
    elif args.command == 'showcase':showcase(args.output)
    elif args.command == 'validate':
        for name in ['test_chain.py', 'test_chain_nonflat.py', 'test_unit_gap.py', 'test_transport.py', 'test_atlas_binding.py', 'test_chain_contract.py']:
            run('tests/' + name)
    else:reproduce(args.dataset,args.output_dir)

if __name__ == '__main__':
    main()
