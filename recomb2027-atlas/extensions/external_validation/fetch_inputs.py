"""Retrieve the selected archive members and verify their recorded SHA-256 values."""
from pathlib import Path
import argparse
import hashlib
import json
import os

from remotezip import RemoteZip


def fetch(donor, output):
    catalog = json.loads(Path(__file__).with_name('sources.json').read_text())
    entry = catalog[donor]
    destination = output / donor
    destination.mkdir(parents=True, exist_ok=True)
    for member in entry['members']:
        path = destination / member['name']
        if not path.exists():
            temporary = path.with_name(path.name + '.partial')
            with RemoteZip(entry['source_archive_url'], timeout=180) as archive:
                with archive.open(member['name']) as source, temporary.open('wb') as target:
                    while block := source.read(8 * 1024 * 1024):
                        target.write(block)
            os.replace(temporary, path)
        with path.open('rb') as source:
            digest = hashlib.file_digest(source, 'sha256').hexdigest()
        if path.stat().st_size != member['size'] or digest != member['sha256']:
            raise ValueError(f'Source member differs from the published input: {path}')
        print(donor, member['name'], digest, flush=True)
    (destination / 'provenance.json').write_text(json.dumps(entry, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--donors', nargs='+', choices=['IDC', 'ILC', 'LUNG'], default=['IDC', 'ILC'])
    args = parser.parse_args()
    for donor in args.donors:
        fetch(donor, args.output)
