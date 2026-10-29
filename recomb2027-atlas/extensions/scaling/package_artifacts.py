"""Package a completed, reconciled crop with deterministic ZIP metadata."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

from restore_artifacts import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared-root', type=Path, required=True)
    parser.add_argument('--size', choices=['400', '600'], required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--part-bytes', type=int, default=4 * 1024 * 1024)
    args = parser.parse_args()
    if args.part_bytes < 1:
        raise ValueError('Positive part size required')
    run = args.prepared_root
    measurement = json.loads((run / 'measurement.json').read_text())
    reconciliation = json.loads((run / 'results/output_reconciliation.json').read_text())
    if measurement['status'] != 'passed' or reconciliation['status'] != 'passed':
        raise ValueError('Package only completed measurements with passing output reconciliation')
    files = {}
    for directory in ['data', 'docs', 'results', 'previous_attempt']:
        for path in (run / directory).rglob('*'):
            if path.is_file():
                files[path.relative_to(run).as_posix()] = path
    for path in run.iterdir():
        if path.is_file() and path.suffix in {'.json', '.log'} and path.name != 'output_storage_note.json':
            files[path.name] = path
    extension = Path(__file__).parent
    for name in ['selection_protocol.json', 'source_provenance.json', 'environment.json']:
        files['provenance/' + name] = extension / name
    files_sha256 = {name: digest(path) for name, path in sorted(files.items())}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    filename = f'geocert-central{args.size}-checked.zip'
    target = args.output_dir / filename
    partial = target.with_suffix('.zip.tmp')
    with zipfile.ZipFile(partial, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
        for name, path in sorted(files.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            with path.open('rb') as source, bundle.open(info, 'w', force_zip64=True) as output:
                shutil.copyfileobj(source, output, length=8 * 1024 * 1024)
    partial.replace(target)
    with zipfile.ZipFile(target) as bundle:
        for name, expected in files_sha256.items():
            h = hashlib.sha256()
            with bundle.open(name) as stream:
                for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                    h.update(block)
            if h.hexdigest() != expected:
                raise ValueError('Packaged file mismatch: ' + name)
    manifest = json.loads(args.manifest.read_text()) if args.manifest.exists() else {'format_version': 1, 'artifacts': {}}
    parts = []
    with target.open('rb') as source:
        index = 0
        while block := source.read(args.part_bytes):
            index += 1
            name = filename + f'.part{index:02d}'
            part = args.output_dir / name
            temporary = part.with_suffix(part.suffix + '.tmp')
            temporary.write_bytes(block)
            temporary.replace(part)
            parts.append({'filename': name, 'bytes': len(block), 'sha256': hashlib.sha256(block).hexdigest()})
    manifest['artifacts'][args.size] = {'filename': filename, 'bytes': target.stat().st_size,
                                      'sha256': digest(target), 'files_sha256': files_sha256,
                                      'parts': parts,
                                      'status': 'completed_and_reconciled', 'count_method': measurement.get('count_method', 'chain'),
                                      'output_audit_script_sha256': digest(extension / 'audit_outputs.py')}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest['artifacts'][args.size].items() if k != 'files_sha256'}), flush=True)


if __name__ == '__main__':
    main()
