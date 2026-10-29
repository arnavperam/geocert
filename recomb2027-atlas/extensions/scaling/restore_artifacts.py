"""Restore hash-pinned scaling inputs and checked results without recomputation."""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import shutil
import tempfile
import zipfile


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-dir', type=Path, default=Path(__file__).with_name('artifacts'))
    parser.add_argument('--manifest', type=Path, default=Path(__file__).with_name('artifact_manifest.json'))
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--sizes', nargs='+', choices=['400', '600'], default=['400'])
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError('Restore only into an empty output directory')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with ExitStack() as resources:
        for size in args.sizes:
            record = manifest['artifacts'][size]
            archive = args.archive_dir / record['filename']
            if not archive.exists():
                staging_root = Path(resources.enter_context(tempfile.TemporaryDirectory(prefix='geocert-scaling-restore-')))
                archive = staging_root / record['filename']
                with archive.open('wb') as output:
                    for part in record['parts']:
                        source = args.archive_dir / part['filename']
                        if source.stat().st_size != part['bytes'] or digest(source) != part['sha256']:
                            raise ValueError('Archive part integrity mismatch: ' + source.name)
                        with source.open('rb') as stream:
                            shutil.copyfileobj(stream, output, length=8 * 1024 * 1024)
            if archive.stat().st_size != record['bytes'] or digest(archive) != record['sha256']:
                raise ValueError('Archive integrity mismatch: ' + archive.name)
            staging = args.output_dir / f'.central{size}.restore'
            staging.mkdir()
            with zipfile.ZipFile(archive) as bundle:
                names = set()
                for member in bundle.infolist():
                    path = PurePosixPath(member.filename)
                    mode = member.external_attr >> 16
                    if (path.is_absolute() or '..' in path.parts or not path.parts
                        or member.filename in names or stat.S_ISLNK(mode) or member.is_dir()):
                        raise ValueError('Unsafe or duplicate archive member: ' + member.filename)
                    if member.filename not in record['files_sha256']:
                        raise ValueError('Unlisted archive member: ' + member.filename)
                    names.add(member.filename)
                    bundle.extract(member, staging)
                if names != set(record['files_sha256']):
                    raise ValueError('Archive member list mismatch')
            for name, expected in record['files_sha256'].items():
                if digest(staging / name) != expected:
                    raise ValueError('Restored file integrity mismatch: ' + name)
            target = args.output_dir / f'central{size}'
            staging.rename(target)
            print(json.dumps({'size_um': int(size), 'directory': str(target), 'files': len(record['files_sha256']),
                              'archive_sha256': record['sha256'], 'status': 'restored_and_verified'}), flush=True)


if __name__ == '__main__':
    main()
