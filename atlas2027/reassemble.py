#!/usr/bin/env python3
"""Verify and reassemble the complete Atlas archive using only the standard library."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import tempfile
import zipfile


def leaf(value: object, label: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", value) or value in {".", ".."}:
        raise ValueError(f"invalid {label}")
    return value


def digest(value: object) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError("invalid SHA-256 digest")
    return value


def size(value: object) -> int:
    if type(value) is not int or value < 0:
        raise ValueError("invalid byte size")
    return value


def reassemble(manifest_path: Path, output: Path | None = None) -> tuple[Path, str]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    name = leaf(manifest["archive"], "archive filename")
    root = leaf(manifest["root"], "archive root")
    expected_hash, expected_size = digest(manifest["sha256"]), size(manifest["size"])
    parts = manifest["parts"]
    if not isinstance(parts, list) or not parts:
        raise ValueError("manifest must contain a nonempty part list")
    seen = set()
    for part in parts:
        part_name = leaf(part["name"], "part filename")
        if part_name in seen:
            raise ValueError("duplicate part filename")
        seen.add(part_name)
        digest(part["sha256"])
        size(part["size"])
    if sum(p["size"] for p in parts) != expected_size:
        raise ValueError("part sizes do not sum to archive size")
    output = output or manifest_path.parent.parent / name
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_symlink():
        raise ValueError("output must not be a symbolic link")
    handle, temporary_name = tempfile.mkstemp(prefix=output.name + ".", suffix=".partial", dir=output.parent)
    temporary = Path(temporary_name)
    total_hash, total_size = hashlib.sha256(), 0
    try:
        with os.fdopen(handle, "wb") as dest:
            for part in parts:
                source = manifest_path.parent / part["name"]
                if source.is_symlink() or not source.is_file():
                    raise ValueError(f"part is missing or is a link: {part['name']}")
                part_hash, part_size = hashlib.sha256(), 0
                with source.open("rb") as src:
                    for block in iter(lambda: src.read(1024 * 1024), b""):
                        dest.write(block)
                        part_hash.update(block)
                        total_hash.update(block)
                        part_size += len(block)
                        total_size += len(block)
                if part_size != part["size"] or part_hash.hexdigest() != part["sha256"]:
                    raise ValueError(f"part checksum or size mismatch: {part['name']}")
        if total_size != expected_size or total_hash.hexdigest() != expected_hash:
            raise ValueError("complete archive checksum or size mismatch")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return output, root


def extract(archive: Path, destination: Path, root: str) -> None:
    if destination.is_symlink() or (destination.exists() and (not destination.is_dir() or any(destination.iterdir()))):
        raise ValueError("extraction destination must be absent or an empty directory")
    with zipfile.ZipFile(archive) as zipped:
        members = zipped.infolist()
        names = set()
        for info in members:
            name = info.filename
            path = PurePosixPath(name)
            mode = info.external_attr >> 16
            if (not name or "\\" in name or ":" in name or name.startswith("/")
                    or ".." in name.split("/") or "." in name.split("/")
                    or not path.parts or path.parts[0] != root
                    or (stat.S_IFMT(mode) not in {0, stat.S_IFREG, stat.S_IFDIR})):
                raise ValueError(f"unsafe archive member: {name}")
            canonical = str(path)
            if canonical in names:
                raise ValueError(f"duplicate archive member: {name}")
            names.add(canonical)
        if zipped.testzip() is not None:
            raise ValueError("archive member failed CRC verification")
        destination.mkdir(parents=True, exist_ok=True)
        resolved_destination = destination.resolve()
        for info in members:
            target = destination / info.filename
            if not target.resolve().is_relative_to(resolved_destination):
                raise ValueError("archive member leaves extraction directory")
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with zipped.open(info) as src, target.open("xb") as dest:
                    shutil.copyfileobj(src, dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).resolve().parent / "artifact_parts/manifest.json")
    parser.add_argument("--output", type=Path, help="assembled ZIP path; defaults beside this script")
    parser.add_argument("--extract", type=Path, help="optionally extract into an absent or empty directory")
    args = parser.parse_args()
    archive, root = reassemble(args.manifest, args.output)
    if args.extract:
        extract(archive, args.extract, root)
    print(json.dumps({"status": "verified", "archive": str(archive), "extracted_to": str(args.extract) if args.extract else None}))


if __name__ == "__main__":
    main()
