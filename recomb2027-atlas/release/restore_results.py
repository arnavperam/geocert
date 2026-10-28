#!/usr/bin/env python3
"""Restore checksum-verified large records and shared source aliases.

No downloads are made. Existing files are verified rather than overwritten.
"""
from __future__ import annotations
from pathlib import Path
import gzip
import hashlib
import json
import os
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def safe(relative: str) -> Path:
    p = (ROOT / relative).resolve()
    if not p.is_relative_to(ROOT) or Path(relative).is_absolute():
        raise ValueError(f"Path leaves release: {relative}")
    return p


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    index = ROOT / 'release/compressed_results.json'
    if not index.exists():
        print(json.dumps({'status': 'working tree already expanded', 'restored': 0}))
        return
    records = json.loads(index.read_text())
    restored = 0
    for row in records['compressed']:
        packed, target = safe(row['compressed_path']), safe(row['path'])
        if sha(packed) != row['compressed_sha256']:
            raise ValueError(f"Compressed checksum mismatch: {row['compressed_path']}")
        if target.exists():
            if sha(target) != row['sha256']:
                raise ValueError(f"Existing expanded record differs: {row['path']}")
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(prefix='.restore-', dir=target.parent)
        try:
            with os.fdopen(fd, 'wb') as out, gzip.open(packed, 'rb') as source:
                shutil.copyfileobj(source, out, 1024 * 1024)
            temp = Path(name)
            if temp.stat().st_size != row['bytes'] or sha(temp) != row['sha256']:
                raise ValueError(f"Expanded checksum mismatch: {row['path']}")
            temp.replace(target)
        finally:
            Path(name).unlink(missing_ok=True)
        restored += 1
    aliases = 0
    for row in records['aliases']:
        source, target = safe(row['target']), safe(row['path'])
        # Sources are internal, already checksum-verified release members. A copy
        # avoids platform-specific symlink behavior after ZIP extraction.
        members = list(source.rglob('*')) if source.is_dir() else [source]
        for item in members:
            if not item.is_file():
                continue
            destination = target / item.relative_to(source) if source.is_dir() else target
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                if sha(destination) != sha(item):
                    raise ValueError(f"Existing source alias differs: {destination}")
            else:
                shutil.copyfile(item, destination)
                aliases += 1
    print(json.dumps({'status': 'passed', 'restored': restored, 'shared_input_files': aliases}))


if __name__ == '__main__':
    main()
