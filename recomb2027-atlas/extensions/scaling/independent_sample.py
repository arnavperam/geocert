#!/usr/bin/env python3
"""Preselect multicopy queries, enumerate independently, and compare spectra.

Selection depends only on the input atlas. The reference solver enumerates target
assignments using scalar closure and tent integration; it imports no production
chain, cut, geometry, envelope or integration modules. A completed sample check
does not establish correctness of unselected queries.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib
import json
from pathlib import Path
import signal
import sys
import time


RULE = ("Enumerate every (patch,gene) with 2 through 6 copies inclusive; "
        "order by SHA256 of UTF8 roi:patch:gene and take the first 16. "
        "Selection uses only the input atlas.")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    with Path(path).open() as stream:
        return json.load(stream)


def write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


def candidates(atlas):
    rows = []
    for patch in atlas["patches"]:
        for gene, copies in Counter(p["gene"] for p in patch["points"]).items():
            if 2 <= copies <= 6:
                key = f'{atlas["roi"]}:{patch["id"]}:{gene}'
                rows.append({"key": key,
                             "selection_sha256": hashlib.sha256(key.encode()).hexdigest(),
                             "patch": patch["id"], "gene": gene, "copies": copies,
                             "total_patch_points": len(patch["points"])})
    return sorted(rows, key=lambda r: r["selection_sha256"])


def checked_selection(atlas_path, selection_path):
    atlas, selection = read(atlas_path), read(selection_path)
    if digest(atlas_path) != selection["atlas_file_sha256"]:
        raise ValueError("input atlas digest changed")
    if atlas["atlas_id"] != selection["atlas_id"] or atlas["roi"] != selection["roi"]:
        raise ValueError("atlas identity changed")
    eligible = candidates(atlas)
    if len(eligible) != selection["candidate_count"] or eligible[:16] != selection["selected"]:
        raise ValueError("selection is not the prescribed first 16 hash-ranked queries")
    return atlas, selection


def select(args):
    atlas = read(args.atlas)
    eligible = candidates(atlas)
    write(args.output, {"selection_rule": RULE,
                        "selected_at_utc": datetime.now(timezone.utc).isoformat(),
                        "atlas_file": Path(args.atlas).name,
                        "atlas_file_sha256": digest(args.atlas),
                        "atlas_id": atlas["atlas_id"], "roi": atlas["roi"],
                        "candidate_count": len(eligible), "selected": eligible[:16]})


def enumerate_reference(args):
    atlas, selection = checked_selection(args.atlas, args.selection)
    patches = {p["id"]: p for p in atlas["patches"]}
    sys.path.insert(0, str(Path(args.base_root).resolve()))
    scalar = importlib.import_module("atlas.scalar_reference")
    result = {k: v for k, v in selection.items() if k != "selected"}
    result.update(reference_solver="atlas.scalar_reference.spectrum; no production chain/cut imports",
                  reference_source_sha256=digest(scalar.__file__),
                  selection_sha256=digest(args.selection),
                  per_query_limit_seconds=args.timeout_seconds, results=[])

    class Deadline(Exception):
        pass

    def alarm(_signum, _frame):
        raise Deadline(f"{args.timeout_seconds}-second query deadline")

    previous = signal.signal(signal.SIGALRM, alarm)
    try:
        for query in selection["selected"]:
            row = dict(query)
            started = time.monotonic()
            signal.alarm(args.timeout_seconds)
            try:
                truth, checked = scalar.spectrum(patches[query["patch"]], query["gene"])
                row.update(status="computed", assignments_checked=checked,
                           minima={str(q): str(cost) for q, cost in sorted(truth.items())})
            except Deadline as exc:
                row.update(status="resource_limit", error=str(exc))
            except Exception as exc:
                row.update(status="error", error=repr(exc))
            finally:
                signal.alarm(0)
            row["seconds"] = time.monotonic() - started
            result["results"].append(row)
            write(args.output, result)
    finally:
        signal.signal(signal.SIGALRM, previous)
    result["counts"] = {status: sum(r["status"] == status for r in result["results"])
                        for status in ("computed", "resource_limit", "error")}
    write(args.output, result)
    print(json.dumps(result["counts"]))


def compare(args):
    atlas, selection = checked_selection(args.atlas, args.selection)
    reference = read(args.reference)
    if reference["atlas_file_sha256"] != selection["atlas_file_sha256"]:
        raise ValueError("reference atlas binding differs")
    records = reference["results"]
    if len(records) != len(selection["selected"]):
        raise ValueError("reference omitted a selected query")
    wanted = {(atlas["roi"], r["patch"], r["gene"]) for r in selection["selected"]}
    indexed = {}
    for row in read(args.spectra):
        key = row["roi"], row["patch"], row["gene"]
        if key in wanted:
            if key in indexed:
                raise ValueError(f"duplicate production query {key}")
            indexed[key] = row
    comparisons = []
    for query, reference_row in zip(selection["selected"], records):
        if any(reference_row[k] != v for k, v in query.items()):
            raise ValueError("reference query identity or selection changed")
        key = atlas["roi"], query["patch"], query["gene"]
        out = {"key": query["key"], "status": reference_row["status"]}
        if reference_row["status"] == "computed":
            actual = indexed.get(key)
            if actual is None:
                out.update(status="mismatch", reason="production query missing")
            elif actual["copies"] != query["copies"]:
                out.update(status="mismatch", reason="production target multiplicity differs")
            else:
                truth = {int(q): Fraction(v) for q, v in reference_row["minima"].items()}
                got = {int(q): Fraction(v) for q, v in actual["minima"].items()}
                out.update(status="matched" if truth == got else "mismatch",
                           count_states=len(truth),
                           assignments_checked=reference_row["assignments_checked"])
        comparisons.append(out)
    statuses = Counter(r["status"] for r in comparisons)
    result = {"status": "passed" if statuses["matched"] == len(comparisons) else "incomplete_or_failed",
              "scope": "Exact rational count spectra for 16 input-selected multicopy queries only",
              "atlas_id": atlas["atlas_id"], "atlas_file_sha256": digest(args.atlas),
              "selection_sha256": digest(args.selection), "reference_sha256": digest(args.reference),
              "production_spectra_sha256": digest(args.spectra),
              "counts": dict(statuses), "queries": comparisons}
    write(args.output, result)
    print(json.dumps({"status": result["status"], "counts": result["counts"]}))
    return 0 if result["status"] == "passed" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("select", "enumerate", "compare"):
        command = sub.add_parser(name)
        command.add_argument("--atlas", required=True)
        command.add_argument("--output", required=True)
        if name != "select":
            command.add_argument("--selection", required=True)
        if name == "enumerate":
            command.add_argument("--base-root", required=True)
            command.add_argument("--timeout-seconds", type=int, default=15)
        if name == "compare":
            command.add_argument("--reference", required=True)
            command.add_argument("--spectra", required=True)
    args = parser.parse_args()
    return {"select": select, "enumerate": enumerate_reference, "compare": compare}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
