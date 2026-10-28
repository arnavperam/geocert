# GeoCert Atlas

Exact cell-by-gene count intervals under a shared segmentation-deformation budget.

Companion to **Certified Expression Intervals under Shared Segmentation Uncertainty**, by Arnav Peram. Version **1.1.0** combines continuous assignment-cost reduction, an exact count-chain dynamic program, and compatible boundary patches. It includes the physical-width reach study and a reanalysis of supplied Proseg outputs.

A result is conditional on the recorded atlas: patch seams, uncovered boundaries, junction neighborhoods and nuclear owners stay fixed. The method measures sensitivity of retained observed counts; it does not estimate corrected anatomical assignments.

## Obtain the scientific inputs

From a public repository checkout, reconstruct the checked scientific archive before running the solver:

```sh
python recomb2027-atlas/reassemble.py --extract atlas-replay
cd atlas-replay/recomb2027-atlas
```

If this directory already contains `data/`, `external/` and `reach/`, the archive has been extracted and you can proceed directly. The public `ARCHIVE.md` explains the separate scientific archive and final-document package.

## Run the example

The prepared-input solver uses the Python standard library. Agent 1 recorded a CPython 3.13.5 replay. The final evaluation revision is checked again locally with Python 3.12 and in the dedicated Python 3.12/3.13 workflow. Historical measurements are retained separately from these verification runs.

```sh
python3 -m venv .venv
. .venv/bin/activate
python atlas_cli.py verify-manifest
python release/restore_results.py
python atlas_cli.py showcase --output run/showcase.json
python atlas_cli.py validate
```

On Windows, use `.venv\Scripts\python.exe`. Keep assertions enabled; do not use `-O`.

Large records are gzip-compressed to keep individual repository files manageable. Restoration checks the compressed and expanded SHA-256 identities and reconstructs shared input paths without platform-specific symlinks. It makes no network request. Allow approximately 1.2 GB for the expanded research records.

The example is cell 81, LDHB, in breast format-test crop 1. Source count 22 has exact interval `[19,23]` at 0.25 µm² and `[17,23]` at 0.5 µm². Count 17 is attained at approximately 0.484040265 µm² by two compatible patch changes. The JSON gives exact rational area, five transferred LDHB observations and one collateral KRT8. The format-test dataset does not support biological conclusions from this example.

## Reproduce the primary count tables

```sh
python atlas_cli.py reproduce --dataset lung --output-dir run/lung
python atlas_cli.py reproduce --dataset breast --output-dir run/breast
python release/check_interval_hashes.py
```

Each command recomputes the spectra and compares all exact minima and CSV rows with the supplied records. It uses the count chain for multi-copy questions and closure for single-copy questions.

| Collection | Complete cells | Patches | Cell–gene spectra | Original interval rows |
|---|---:|---:|---:|---:|
| Lung crops | 153 | 544 | 3,304 | 26,432 |
| Breast format-test crops | 467 | 1,708 | 22,976 | 183,808 |

The primary eight-budget marginal comparisons tighten **4/46 lung** and **203/1,026 breast** cell–gene questions whose marginal width can reach at least two. Repeated budget rows are not biological replicates. Sparse zero-entry semantics are recorded in the matrix indexes.

## Reproduce the reach study

The protocol in `docs/PREREGISTRATION.md` was fixed before computing the new width outcomes. It uses already available crops, not a new biological validation cohort. Required physical support halfwidths are 0.5, 1, 2 and 3 µm; accepted displacement is at most three-quarters of support. Each setting constructs a different atlas without fallback to smaller widths.

Recompute one complete prepared family:

```sh
python analysis/reach/reproduce.py --dataset lung --width 1 --output-dir run/lung-width1
```

Use either dataset and one of `0.5`, `1`, `2`, `3`. Breast width 3 has much larger exact records and took about 20 minutes in the original single-run environment. This is a resource observation, not a universal runtime guarantee or repeated benchmark.

Check every supplied global spectrum and all 1,807,967 new interval rows, including a deterministic independent scalar sample:

```sh
python analysis/reach/check_outputs.py
python analysis/reach/test_forced_reduction.py
python analysis/reach/audit_proseg.py
```

The new study allows displacement caps up to approximately 2.23 µm. Wider supports reduce accepted boundary coverage and aggregate alternative-reference count reach under this extractor. All widths and failures of model reach remain included.

## Tables and figures

Raster analysis and figures need the optional analysis environment:

```sh
python -m pip install -r requirements-analysis.txt
```

Run the following in a working copy; report-generation scripts replace derived summaries, not the underlying exact spectra.

| Output | Command |
|---|---|
| Primary input/coverage table | `python analysis/coverage.py` |
| Original exact count verification | `python analysis/check_global_counts.py` and `GEOCERT_DATA_ROOT=external python analysis/check_global_counts.py` |
| Additive-cost baseline | `python analysis/additive_baseline.py` and `GEOCERT_DATA_ROOT=external python analysis/additive_baseline.py` |
| All-width statistics and fixed-key reference reach | `python analysis/reach/report.py` |
| Complete supplemental budget tables | `python analysis/reach/make_tables.py` |
| Source-mask and attaining-boundary figures | `python analysis/make_figures.py` |
| Reach-versus-coverage and budget figures | `python analysis/reach/make_figures.py` |
| Geometry checks, original lung atlas | `python analysis/verify_atlas.py` |
| Geometry checks, original breast atlas | `GEOCERT_DATA_ROOT=external python analysis/verify_atlas.py` |
| Source-count re-extraction/check | `python analysis/check_source_inputs.py` |

For a new-width geometry check, set `GEOCERT_DATA_ROOT=reach/lung/r1`, or the corresponding recorded setting. For full width extraction, use `GEOCERT_DATA_ROOT=reach/lung/r1 GEOCERT_SUPPORT_UM=1 python analysis/reach/build_width_atlas.py`; the shared source paths must first be restored. These POSIX environment-variable examples can be set through the shell's environment facility on other platforms.

Build both PDFs with `pdflatex` twice from `manuscript/`. LaTeX and fonts are standard system dependencies, not distributed font files.

## Proseg and additional data

`results/proseg/` preserves the supplied Proseg 3.2.0 records. One adapted no-factorization execution per clipped lung crop completed after default factorization failed on the first crop. We independently checked transcript joins, observed coordinates, source-cell mapping and both sparse count matrices, then evaluated the returned counts under each Atlas family. These are reanalyzed existing runs, not five new runs or an accuracy benchmark. Missing records and crop-edge exclusions are retained in fixed comparison denominators.

`acquisition/` and `release/acquire.yml` provide a prospective real-breast/additional-lung acquisition and five-run Proseg workflow. Its ZIP-range helper was tested against controlled responses, but the external download and five-run experiment were **not executed in this release**. No results from that proposed study appear in the paper.

## Layout

- `atlas/`: scalar geometry, count-chain solver, gluing maps and independent scalar check.
- `geocert_expression/`: preserved continuum/cut primitives and query routines.
- `data/`, `external/data/`: reference crops and primary atlases.
- `reach/`: all eight width-specific atlases, spectra and witnesses.
- `results/`: original experiments, supplied Proseg records and complete reach summaries.
- `analysis/`, `tests/`: reproduction and correctness checks.
- `manuscript/`: paper, supplement, LaTeX and figures.
- `docs/`: assumptions, evaluation, author explanations and data attribution.
- `release/`: restoration, interval checks, proposed CI and publication instructions.

## Release and rights

This is the **Atlas 1.1.0 publication candidate**, separate from public v1.0.0 and the corridor-only draft PR #3. The canonical evaluation branch is `submission/atlas-1.1.0`, under `recomb2027-atlas/`. The manuscript cites an immutable scientific source commit; the later document commit supplies the rebuilt PDFs. This is an evaluation candidate, not an assertion that a stable `v1.1.0` tag or conference submission exists. See `ARCHIVE.md` in the public source folder for the checked scientific archive; the complete review package also includes both PDFs and their LaTeX sources.

Software is MIT licensed. Source-derived 10x data retain their CC BY 4.0 attribution; the trimmed breast data retain their format-test restriction. Proseg's upstream rights are separate; its binary is not included. See `NOTICE.md` and `docs/DATA_PROVENANCE.md`.

Generative AI assisted mathematics, implementation, experiments, verification and drafting. The author is responsible for the submission, claims and ability to defend them. Complete-family anatomical calibration, molecular ground truth, independent-donor biological validation and validated downstream inference are outside the completed evidence.

## Final evaluation checks

The point-count solver rejects unsupported free location intervals, multisite groups and overlapping owner thresholds. Both primary and strict-width certification bind the atlas to the prepared observations, source files and completed rational geometry audit. These are reproducibility checks; they do not replace running the geometry verifier for a newly constructed atlas. Existing point-data results remain unchanged. Regression records are in `results/final_guard_checks/`.
