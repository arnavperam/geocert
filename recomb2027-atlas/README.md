# GeoCert Atlas

Exact cell-by-gene count intervals under a shared segmentation-deformation budget.

Code and data for **Certified Expression Intervals under Shared Segmentation Uncertainty**, by Arnav Peram. Version 1.1.0 combines continuous assignment costs, an exact count-chain dynamic program and compatible boundary patches.

[Paper](manuscript/main.pdf) · [Supplement](manuscript/supplement.pdf) · [Fixed scientific source](https://github.com/arnavperam/geocert/tree/6581f628a5ac8fb83629f9e0dee0ffa4d685375d/recomb2027-atlas)

A result is conditional on the recorded atlas: patch seams, uncovered boundaries, junction neighborhoods and nuclear owners stay fixed. The method measures sensitivity of retained observed counts; it does not estimate corrected anatomical assignments.

## Obtain the scientific inputs

From the repository root, reconstruct the checked scientific archive before running the solver:

```sh
python recomb2027-atlas/reassemble.py --extract atlas-replay
cd atlas-replay/recomb2027-atlas
```

If this directory already contains `data/`, `external/` and `reach/`, the archive has been extracted and you can proceed directly. See [ARCHIVE.md](ARCHIVE.md) for the archive identity and extraction checks.

## Run the example

The prepared-input solver uses only the Python standard library and is checked on Python 3.12 and 3.13.

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

Each `reproduce` command recomputes the spectra and compares all exact minima and CSV rows with the archived records. It uses the count chain for multi-copy questions and closure for single-copy questions. The final command checks the recorded interval-table hashes.

| Collection | Complete cells | Patches | Cell–gene spectra | Original interval rows |
|---|---:|---:|---:|---:|
| Lung crops | 153 | 544 | 3,304 | 26,432 |
| Breast format-test crops | 467 | 1,708 | 22,976 | 183,808 |

The primary eight-budget marginal comparisons tighten **4/46 lung** and **203/1,026 breast** cell–gene questions whose marginal width can reach at least two. Repeated budget rows are not biological replicates. Sparse zero-entry semantics are recorded in the matrix indexes.

## Reproduce the reach study

The [width protocol](docs/PREREGISTRATION.md) was fixed before computing the outcomes and uses the same lung and breast format-test crops. Required physical support halfwidths are 0.5, 1, 2 and 3 µm; accepted displacement is at most three-quarters of support. Each setting constructs a different atlas without fallback to smaller widths.

Recompute one complete prepared family:

```sh
python analysis/reach/reproduce.py --dataset lung --width 1 --output-dir run/lung-width1
```

Use either dataset and one of `0.5`, `1`, `2`, `3`. Breast width 3 has larger exact records; its original single run took about 20 minutes. Runtime depends on the environment.

Check all 106,351 recorded global spectra and 1,807,967 width-study interval rows, with an independent scalar check of a deterministic local sample:

```sh
python analysis/reach/check_outputs.py
python analysis/reach/test_forced_reduction.py
python analysis/reach/audit_proseg.py
```

The study allows displacement caps up to approximately 2.23 µm. Wider supports reduce accepted boundary coverage and aggregate alternative-reference count reach under this extractor.

## Finite union of the five atlases

The current paper also evaluates these five complete atlas families as a finite union. Each count minimum selects one member, preserving its global partition and physical-area budget. The original primary results and the four width-specific results remain available separately.

Run the extension from the repository checkout, with the extracted scientific package as its input:

```sh
python recomb2027-atlas/extensions/atlas_union/union.py --data-root atlas-replay/recomb2027-atlas --output-dir run/atlas-union
python recomb2027-atlas/extensions/atlas_union/independent_verify.py --data-root atlas-replay/recomb2027-atlas --union-output run/atlas-union --output-dir run/atlas-union-check
```

The [extension README](extensions/atlas_union/README.md) defines the output census, exact count-membership checks and witness command. [Recorded spectra and interval tables](extensions/atlas_union/results/) contain 29,212 questions and 496,604 budget rows. The independent verification report binds the recorded output files by checksum.

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

The manuscript sources are in the repository checkout, outside the scientific archive. From `recomb2027-atlas/manuscript/`, run `pdflatex main.tex` twice and `pdflatex supplement.tex` twice. LaTeX is a separate system dependency.

## Proseg and additional data

`results/proseg/` preserves the Proseg 3.2.0 execution records. One adapted no-factorization execution per clipped lung crop completed after default factorization failed on the first crop. We independently checked transcript joins, observed coordinates, source-cell mapping and both sparse count matrices, then evaluated the returned counts under each Atlas family. This comparison measures model reach using the existing executions; it does not measure anatomical accuracy. Missing records and crop-edge exclusions are retained in fixed comparison denominators.

`acquisition/` and `release/acquire.yml` contain scripts for a separate data-acquisition and five-run Proseg study. That study was not executed and contributes no results to the paper.

## Layout

- `atlas/`: scalar geometry, count-chain solver, gluing maps and independent scalar check.
- `geocert_expression/`: preserved continuum/cut primitives and query routines.
- `data/`, `external/data/`: reference crops and primary atlases.
- `reach/`: all eight width-specific atlases, spectra and witnesses.
- `results/`: primary experiments, Proseg records and complete reach summaries.
- `analysis/`, `tests/`: reproduction and correctness checks.
- `manuscript/`: paper, supplement, LaTeX and figures.
- `docs/`: assumptions, evaluation, author explanations and data attribution.
- `release/`: restoration and interval checks.

## Input checks

The point-count solver rejects unsupported free location intervals, multisite groups and overlapping owner thresholds. Primary and width-study certification bind the atlas to its prepared observations, source files and completed rational geometry audit. Run the geometry verifier for any newly constructed atlas. Regression records are in `results/final_guard_checks/`.

## License and disclosure

Software in this package is [MIT licensed](LICENSE). Source-derived 10x data retain their CC BY 4.0 attribution, and the trimmed breast data retain their format-test restriction. Proseg's upstream rights are separate; its binary is not included. See [NOTICE.md](NOTICE.md) and [DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md).

Generative AI assisted mathematics, implementation, experiments, verification and drafting. The author is responsible for the work and its claims.

The evidence is limited to the declared boundary families and recorded crops. Anatomical calibration, molecular ground truth, independent-donor biological validation and downstream inference remain unvalidated.
