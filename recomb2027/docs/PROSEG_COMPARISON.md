# Executed Proseg development-crop audit

Proseg 3.2.0 was executed on the two existing Xenium lung development crops on
2026-10-28. This adds an actual external-method execution and a reproducible
comparison of assignments. It does **not** add a donor, known transcript ownership,
or held-out evaluation of the corridor uncertainty model.

## Executed protocol and adaptation

The official Proseg repository is <https://github.com/dcjones/proseg>. Its documented
Xenium input preset accepts transcript CSV or Parquet. The executable was obtained
from the Bioconda 3.2.0 Linux package; its package and executable SHA256 values are
in `results/proseg/provenance.json`.

The initial unmodified factorization run on roi_01 failed with `EmptyInput` in
Proseg feature selection. The failure reproduced for Parquet and an equivalent
CSV. The retained log is `results/proseg/roi_01_default_failure.log`. We then
disabled factorization using the documented `--no-factorization` option for BOTH
crops. This adaptation was selected because of the runtime failure, before
comparing assignment results. All other model/sampler settings remain upstream
defaults: 200 burn-in iterations, 200 sampling iterations, 50 hill-climb steps;
two CPU threads were requested. Neither diffusion nor source reassignment was
disabled. This is an adapted small-input execution, not a benchmark of Proseg's
default full-tissue workflow.

The inputs are prefiltered development-crop tables from a single lung donor,
already supplied with GeoCert. Crop edges truncate cells and remove neighboring
context. The runner converts their rows to CSV without adding, removing, or
changing transcript labels. Proseg itself excludes some input rows. Every excluded
ID is retained in the corresponding run record. Default initialization uses
nucleus-associated assignments and can omit source cells lacking retained nuclear
transcripts. A single stochastic run was retained per crop; Proseg exposes no seed
argument in this version, so a repeat can change numerical outcomes.

## Descriptive results

The source and Proseg owners are joined by exact transcript ID. Proseg cell indices
are mapped through `original_cell_id`; missing assignment values are background.
Genes and observed XY coordinates are checked before comparisons. Counts below
use only the identical retained molecules in each method. Source labels are a
comparison reference, **not truth**.

| Quantity | roi_01 | roi_02 |
|---|---:|---:|
| Supplied molecules | 5,847 | 1,523 |
| Proseg-retained molecules | 5,774 | 1,482 |
| Excluded molecules | 73 | 41 |
| Source cells represented among retained molecules | 106 | 85 |
| Proseg output cells | 93 | 72 |
| Same owner, including shared unassigned state | 2,532 | 717 |
| Different owner | 3,242 | 765 |
| Proseg unassigned/background | 2,484 | 1,335 |
| Active source-or-output cell/gene pairs | 3,373 | 616 |
| Changed cell/gene counts | 2,794 | 590 |
| Sum of absolute cell/gene count changes | 3,720 | 777 |
| Executable wall time (seconds) | 14.14 | 10.67 |

The large number of background assignments, particularly in roi_02, reinforces
that these clipped, sparse inputs are insufficient for a strong field benchmark.
Neither agreement nor disagreement measures correctness. We report no AUROC,
error reduction, posterior calibration, interval coverage, or biological effect
from these runs. The timing covers the Proseg invocation on these inputs, excludes
installation and adapter work, and provides no whole-slide scaling claim.

## Reproduction

Install Proseg 3.2.0 separately, for example with a compatible Rust toolchain:

```sh
cargo install --locked proseg --version 3.2.0
python -m pip install pandas pyarrow
python experiments/proseg_comparison.py --binary proseg \
  --crop data/source_crops/roi_01 --output run/proseg/roi_01
python experiments/proseg_comparison.py --binary proseg \
  --crop data/source_crops/roi_02 --output run/proseg/roi_02
python tests/test_proseg.py
```

The runner refuses an existing output directory, a missing executable or a
different version. It saves the exact executable command, binary/input/output
hashes, stdout/stderr, exit code, elapsed time, all exclusions, and comparison
results. A process or validation failure is marked `failed`, never completed.
Nine validator tests cover successful joins, exclusions, duplicate/unknown IDs,
gene/coordinate mismatches, unknown cells, and a missing executable.

All generated CSV/MatrixMarket and SpatialData outputs from the two completed runs
are retained under `results/proseg/`; each SpatialData directory is losslessly
archived as `proseg-output.zarr.zip` (extract before reading). The large third-party executable is excluded
from this project's deliverables. Installation uses Proseg's own distribution and
license. The 72-crop RNA2seg external study was not rerun.

## Remaining decisive study

Acquire complete independent specimens with aligned reference measurements,
freeze extraction and uncertainty settings without using outcomes, run segmenters
on adequate tissue context, and evaluate conditional count width together with
geometric containment and independently adjudicated ownership. These two runs
make the execution barrier concrete; they cannot supply that missing evidence.
