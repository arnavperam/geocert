# Reproduction

## Exact workflow

Create an isolated Python 3.13 environment (the new extension was executed with
3.13.5). The package is declared compatible with Python>=3.12; this candidate has
not run a new3.12 CI job. The inherited public core had earlier3.12 checks.

```sh
python3.13 -m venv .venv
. .venv/bin/activate
pip install --no-index --no-deps dist/geocert_expression-0.2.0a1-py3-none-any.whl
geocert-expression demo --output run/demo.json
python validation/replay.py --full --cpp --output run/replay.json
```

The full replay verifies the manifest, creates a disposable copy, regenerates all
84 native models, and compares every interval spectrum and summary to the frozen
records. It runs complete small-instance, calibration, probability, composition,
resource-budget and joint-count tests. With `--cpp`, g++ and Boost headers must be
available; the additional flow test verifies 100 exhaustive graphs and 252 native
objectives. No NumPy, pandas, SciPy or NetworkX is required for this exact replay.
Original recorded timings are not overwritten. Omit `--full` for the smaller
correctness suite. Assertions must remain enabled.

## Source ingestion and raster analysis

These use optional scientific dependencies and are recorded separately:

```sh
pip install -r requirements-analysis.txt
python experiments/prepare_counts.py
python tests/test_raster_adapter.py
python experiments/raster_relaxation.py
```

The first command uses the included restricted Parquet reader, an inherited shared
dependency. It is not a new standard-engine ingestion cross-check. The original
source crops and SHA-256 identities are retained. A separate source-matrix query
file carries fixed out-of-corridor counts explicitly.

## Engineering measurement

```sh
python experiments/benchmark_streaming.py
```

This compiles the streaming accumulator, checks it against Python on smaller
examples and measures three processes at each supplied synthetic size. Its time
excludes construction of candidate neighborhoods and actual segmentation. It is
not an end-to-end whole-slide performance test. The report preserves the original
measurements and describes any replay separately.

## Direct use and exports

The CLI accepts physical models as chart JSON, not arbitrary uninterpreted cell
polygons. The optional raster adapter accepts a labeled NumPy array with explicit
physical axis order, spacing and origin. A SpatialData transform is not inferred.
The returned conservative radius-based intervals use a different uncertainty
family from the exact corridor area intervals.

`results/expression_matrix_intervals.csv` and its index describe84 distinct
conditional matrices. Do not merge overlapping corridor rows into one global
certificate. The min-plus composer requires independently justified, disjoint
editable modules and compatible fixed cores.
