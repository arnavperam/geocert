"""Restore and checksum the published analysis inputs, spectra and core assignments."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile


def restore(output):
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / 'artifacts/manifest.json').read_text())
    output.mkdir(parents=True, exist_ok=True)
    for record in manifest['archives']:
        archive = root / 'artifacts' / record['name']
        with archive.open('rb') as source:
            assert hashlib.file_digest(source, 'sha256').hexdigest() == record['sha256']
        with zipfile.ZipFile(archive) as source:
            for entry in source.infolist():
                target = (output / entry.filename).resolve()
                assert target.is_relative_to(output.resolve())
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(source.read(entry))
        print(record['name'], 'verified and restored', flush=True)
    for name, digest in manifest['files'].items():
        with (output / name).open('rb') as source:
            assert hashlib.file_digest(source, 'sha256').hexdigest() == digest, name
    print(len(manifest['files']), 'restored file hashes verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    restore(parser.parse_args().output)
