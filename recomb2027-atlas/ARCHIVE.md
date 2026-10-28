# Scientific source archive

The browsable code and protocols accompany the scientific archive cited in the paper at [commit `6581f628`](https://github.com/arnavperam/geocert/tree/6581f628a5ac8fb83629f9e0dee0ffa4d685375d/recomb2027-atlas). The archive contains the implementation, prepared inputs, compressed exact results, Proseg execution records and protocols. Its top-level directory is `recomb2027-atlas/`.

| Archive | Value |
|---|---|
| Filename | `GeoCert_RECOMB_1.1.0_Source.zip` |
| Size | 74,715,174 bytes |
| Parts | 18 |
| SHA-256 | `1b08829e7e80a1e966a1bc39920dc4509af99f79ffd11dadfaa71b8338ca59d5` |

[artifact_parts/manifest.json](artifact_parts/manifest.json) records the assembly order, byte count and SHA-256 of each part. All parts are stored in this repository; reassembly requires no external service. [ARCHIVE_SHA256SUMS.txt](ARCHIVE_SHA256SUMS.txt) lists the files inside the reconstructed archive. It is the archive manifest, not a manifest of the browsable checkout; extraction restores it under the name `SHA256SUMS.txt` for the verification command below.

## Extract and verify

Run from the repository root using Python 3.12 or 3.13:

```sh
python3 recomb2027-atlas/reassemble.py --output GeoCert_RECOMB_1.1.0_Source.zip --extract atlas-replay
cd atlas-replay/recomb2027-atlas
python3 atlas_cli.py verify-manifest
python3 release/restore_results.py
python3 release/check_interval_hashes.py
python3 atlas_cli.py validate
python3 atlas_cli.py showcase --output run/showcase.json
```

Keep Python assertions enabled. These commands need only the standard library. Allow about 1.2 GB for the expanded records, plus space for fresh replay outputs. The [reproduction instructions](README.md) cover the primary datasets, width study, geometry and figures.

Reassembly checks every part and the complete ZIP. Extraction requires an absent or empty destination, checks ZIP integrity, and accepts only regular files and directories within the package root. Omit `--extract` to produce only the verified ZIP. Result restoration checks both compressed and expanded file hashes.

## Paper files and CI

The scientific archive is fixed and excludes the PDFs and LaTeX. Current [paper files](manuscript/) are stored separately in the repository; document corrections leave the archived code, inputs and numerical results unchanged.

The [CI workflow](../.github/workflows/atlas-final.yml) checks Python 3.12 and 3.13. It compares the browsable scientific code with the archive, checks file and interval identities, runs algorithm and input-contract tests, recomputes both primary datasets and the one-micrometre lung family, checks all recorded width-study compositions and interval rows, and audits Proseg transcript joins and count arithmetic. Independent scalar checks cover a deterministic sample of local spectra. Proseg execution, new data acquisition, all-width geometry extraction and historical timing measurements are outside this CI replay.
