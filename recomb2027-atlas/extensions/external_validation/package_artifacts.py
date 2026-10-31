"""Package prepared evidence; keep native sampler outputs in the full evidence archive."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile


def package(analysis, destination):
    summary = json.loads((Path(__file__).with_name('results') / 'proseg_summary.json').read_text())
    assert len(summary) == 32 and sum(r['returncode'] == 0 for r in summary) == 31
    destination.mkdir(parents=True, exist_ok=True)
    manifest = {'scope': 'Prepared source crops, exact spectra and witnesses, saved Proseg core assignments, run records and native-output audit digests. Full native Proseg datasets and large context inputs are omitted.', 'archives': [], 'files': {}}
    for donor, prefix, count in [('LUNG', 'roi', 2), ('IDC', 'idc', 3), ('ILC', 'ilc', 3)]:
        for index in range(1, count + 1):
            region = analysis / donor / f'{prefix}_{index:02d}'
            paths = [p for sub in ['data', 'results'] for p in (region / sub).rglob('*') if p.is_file()]
            for name in ['run_status.json', 'context_identity_map.json']:
                if (region / name).exists():
                    paths.append(region / name)
            if index == 1 and (region.parent / 'selection.json').exists():
                paths.append(region.parent / 'selection.json')
            for setting in sorted((region / 'proseg').iterdir()):
                for name in ['run.json', 'run.log', 'audit.json', 'output_provenance.json',
                             'core_assignments.json.gz', 'counts.mtx.gz', 'cell-metadata.csv.gz', 'gene-metadata.csv.gz']:
                    path = setting / name
                    if path.exists():
                        paths.append(path)
            target = destination / f'{donor}_{region.name}.zip'
            with zipfile.ZipFile(target, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
                for path in sorted(paths):
                    name = str(path.relative_to(analysis))
                    data = path.read_bytes()
                    entry = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                    entry.compress_type = zipfile.ZIP_DEFLATED
                    entry.external_attr = 0o644 << 16
                    archive.writestr(entry, data, compresslevel=9)
                    manifest['files'][name] = hashlib.sha256(data).hexdigest()
            manifest['archives'].append({'name': target.name, 'bytes': target.stat().st_size,
                                        'sha256': hashlib.file_digest(target.open('rb'), 'sha256').hexdigest()})
            print(target.name, target.stat().st_size, flush=True)
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    package(args.analysis, args.output)
