# Immutable scientific source archive

The browsable files provide the exact implementation and protocols. The checked parts in `artifact_parts/` contain `GeoCert_RECOMB_1.1.0_Source.zip`, whose top-level directory is `recomb2027-atlas/`. This archive contains the scientific code, original inputs, compressed exact results, retained Proseg outputs, protocols and abstract. Its SHA-256 and every part's SHA-256 and byte count are recorded in `artifact_parts/manifest.json`.

The scientific archive intentionally excludes manuscript PDFs and LaTeX. The paper cites the immutable commit publishing this scientific source. The subsequent publication commit adds the final PDFs and LaTeX under `recomb2027-atlas/manuscript/`; its source code and scientific results remain identical. The separate complete evaluation ZIP includes those final documents as well. These are distinct archive identities and should retain their separate checksums.

From a clone or complete download of this publication branch, run:

```sh
python recomb2027-atlas/reassemble.py --output GeoCert_RECOMB_1.1.0_Source.zip --extract atlas-replay
cd atlas-replay/recomb2027-atlas
python atlas_cli.py verify-manifest
python release/restore_results.py
python release/check_interval_hashes.py
python atlas_cli.py validate
python analysis/reach/test_forced_reduction.py
python atlas_cli.py showcase --output run/showcase.json
python atlas_cli.py reproduce --dataset lung --output-dir run/lung
python atlas_cli.py reproduce --dataset breast --output-dir run/breast
python analysis/reach/reproduce.py --dataset lung --width 1 --output-dir run/lung-width1
python analysis/reach/check_outputs.py
python analysis/reach/audit_proseg.py
```

Use Python 3.12 or later with assertions enabled. All commands above use the Python standard library; optional source extraction and figure generation have separate analysis dependencies. Allow approximately 1.2 GB for restored records, plus fresh replay outputs.

Reassembly verifies each part and the complete archive before installing the ZIP. Optional extraction requires an absent or empty destination and accepts only regular files and directories inside the declared package root. Omitting `--extract` produces only the verified ZIP. Existing unrelated extraction contents are preserved by rejecting a nonempty destination.

The read-only CI workflow tests Python 3.12 and 3.13. It checks public/archive scientific-source equality, all archived file and interval identities, algorithm and source-contract tests, both primary datasets, the strict one-micrometre lung family, all recorded strict-width compositions and interval rows, and imported Proseg arithmetic. It does not execute Proseg, acquire new data, rerun all eight width extractions, or repeat the historical timing studies. Independent scalar checks cover a declared sample of local spectra; their coverage is distinguished from the complete global composition check.

The manifest structure is:

```json
{
  "archive": "GeoCert_RECOMB_1.1.0_Source.zip",
  "root": "recomb2027-atlas",
  "sha256": "<complete archive SHA-256>",
  "size": 123,
  "parts": [{"name": "part000", "sha256": "<part SHA-256>", "size": 123}]
}
```

The actual manifest contains computed values. The part list specifies assembly order. This transport format uses no external download service and makes no claim that an uncreated version tag or conference submission exists.
