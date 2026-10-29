# Larger contiguous computation

Two central squares of the original breast format-test dataset were processed
from native masks through exact count intervals. Their locations were fixed
from the full raster extent before transcript or count outcomes were examined.

| Region (µm) | Exact backend | Complete cells | Input transcripts | Ribbons | Completed stages (s) | All attempts (s) | Peak stage RSS (MiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 399.925 × 399.925 | Chain | 1,363 | 226,983 | 5,102 | 1,353.03 | 1,358.53 | 2,488.77 |
| 599.888 × 599.888 | Enumeration | 3,171 | 470,293 | 11,748 | 1,783.07 | 2,213.11 | 2,251.76 |

| Region | Cell–gene spectra | Exact interval rows | Checked local witnesses | Checked cell-state traces | Incident-edge coverage |
| --- | ---: | ---: | ---: | ---: | ---: |
| 400 µm | 69,278 | 554,224 | 39,519 | 95,253 | 64.94% |
| 600 µm | 152,925 | 1,223,400 | 85,335 | 208,555 | 66.67% |

Each result covers the eight original budgets. Independent rational geometry
validation checked every retained transcript, all accepted ribbon pixels and
their nuclear constraints. A separate output audit reconciled the complete
local and cell–gene census, witness counts, scalar min-plus composition, cell
traces and every exported interval. All checks passed. The completed results
have no unresolved local queries. The 600 µm workload had one earlier global
timeout; its cost is included in the all-attempts column and recorded below.

This expands the contiguous computational footprint within the same 10x breast
format-test specimen. It supplies neither another biological donor nor a
whole-slide measurement. The two squares overlap each other and parts of the
original crops, so their cell counts are not independent sample totals.

## Timing and scope

The 400 µm run used the unchanged primary-atlas source. Extraction took
152.46 s, independent rational geometry validation 199.66 s, and exact count
spectra with witnesses and serialization 1,000.91 s. The reported pipeline
time sums these three fresh-process wall times, including startup and imports.
Source preparation (12.44 s), the coverage audit (1.02 s), final file readback
and independent output reconciliation (19.65 s) are separate from that sum.
The 600 µm completed stages took 423.93 s, 480.19 s and 878.95 s respectively;
preparation took 20.74 s, coverage 2.38 s and independent reconciliation 33.69 s.
The all-attempts column includes the measured wrapper wall time, including
coverage and final readback, and includes the interrupted attempt for 600 µm.
It excludes source preparation, independent reconciliation and time between
attempts. Peak RSS is the largest individual stage process, not aggregate
machine memory.

Measurements used Linux and Python 3.12.14 on a shared container with an
8-CPU quota. The stages ran sequentially, with numerical-library thread counts
set to one and no CPU-affinity pinning. Other work ran concurrently. Exact
package versions are in `environment.json`. The extractor retains its original
global candidate ordering and overlap checks; these measurements establish
completed workloads, not a linear or whole-slide scaling law.

The pre-count amendment for the 600 µm region permits the existing exact
enumeration backend only if its validated atlas has at most 16 copies of any
local gene, within the fixed 65,536-pattern limit. That decision is saved
before count outputs are computed. Different backends are identified
explicitly; their runtimes are not used as a scaling-speedup comparison.

The first 600 µm attempt reached its 1,321-second whole-run allowance during
count computation after extraction and rational validation had passed. One
additional count-only attempt was then recorded, preserving the same atlas,
backend and per-query limits, with a 1,200-second allowance. The failed attempt
is retained. The retry completed in 892.04 s including count computation,
coverage and final readback. Adding the first attempt (1,321.07 s) gives 2,213.11 s
(36.89 min) of cumulative measured attempt time. The sum of successful
scientific stages alone is 1,783.07 s (29.72 min); it excludes the interrupted
work.

## Restore checked results

`artifact_manifest.json` pins the archives, their 4 MiB parts and every
contained file. Restoration reassembles the parts before extraction. Each
archive includes the prepared masks and observations, compiled atlas, exact
results, geometry audit, resource records and output reconciliation.

```sh
python extensions/scaling/restore_artifacts.py \
  --archive-dir extensions/scaling/artifacts \
  --output-dir /tmp/geocert-checked --sizes 400 600

python extensions/scaling/audit_outputs.py \
  --prepared-root /tmp/geocert-checked/central400

python extensions/scaling/audit_outputs.py \
  --prepared-root /tmp/geocert-checked/central600
```

The destination must be empty. Restoration uses only the Python standard
library and checks archive hashes, member paths and every restored file hash.
The saved-output audit reconstructs the complete expected question census,
checks source counts and atlas identities, structural owner bits and witness
transcript counts, independently repeats scalar min-plus composition, and
checks every endpoint and trace cost. Local geometric optimality is supplied
by the exact local solver and its geometric verifier.

## Reproduce from the original source

The [public 10x breast two-FOV source](https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip)
is 379,018,013 bytes; SHA-256 is
`cc1e987b06aa748a6b24d3d6f51fc0d6765daa4836c483f14ec4f1bd18b1779b`.
The downloader verifies this identity and extracts all seven source members
needed by the adapter and independent source audit: cell masks, transcripts,
cells, experiment metadata, gene panel, cell boundaries and nuclear boundaries.

Install the restored scientific source's dependencies plus `pyarrow` and
`zarr>=2.18,<3`. The measured versions are recorded in `environment.json`.
`--frozen-code-root` refers to the original scientific source restored according
to the repository's main README.

```sh
python extensions/scaling/download_source.py --output-dir /path/to/source

python extensions/scaling/prepare_contiguous.py \
  --source-dir /path/to/source/breast_source \
  --source-provenance /path/to/source/breast_archive_provenance.json \
  --selection-protocol extensions/scaling/selection_protocol.json \
  --frozen-code-root /path/to/restored/recomb2027-atlas \
  --side-um 400 --output-dir /path/to/fresh/central400

python extensions/scaling/run_contiguous.py \
  --frozen-code-root /path/to/restored/recomb2027-atlas \
  --prepared-root /path/to/fresh/central400 \
  --count-method chain --total-timeout-seconds 2400

python extensions/scaling/audit_outputs.py \
  --prepared-root /path/to/fresh/central400
```

Use `--side-um 600` to prepare the second fixed region. Its actual side is
599.8875 µm (2,823 pixels); the 400 µm region has 1,882 pixels per side. Both use
the exact reciprocal of the original transform, rather than the rounded
nominal pixel scale. The adapter preserves transcript IDs and coordinates,
checks the full cell-ID mapping, filters Q20 measured genes, and retains the
vendor owner alongside the reference mask owner. Both prepared regions have
zero owner disagreements.

After preparing `central600`, the first measured execution used:

```sh
python extensions/scaling/run_contiguous.py \
  --frozen-code-root /path/to/restored/recomb2027-atlas \
  --prepared-root /path/to/fresh/central600 \
  --count-method enumeration --total-timeout-seconds 1321
```

If that run reaches its wall-time limit after geometric validation, this
command reproduces the additional attempt. It verifies the prior failure,
preserves its logs, reuses the same hash-bound atlas and records the amendment
before starting a fresh count computation:

```sh
python extensions/scaling/run_contiguous.py \
  --frozen-code-root /path/to/restored/recomb2027-atlas \
  --reuse-validated-input /path/to/fresh/central600 \
  --prepared-root /path/to/fresh/central600-retry \
  --count-method enumeration --total-timeout-seconds 1200

python extensions/scaling/audit_outputs.py \
  --prepared-root /path/to/fresh/central600-retry
```

The source audit imports no adapter code and independently compares the full
raw-source census with prepared masks, records, CSV rows and count dictionaries:

```sh
python extensions/scaling/independent_source_audit.py \
  --archive /path/to/source/Xenium_V1_human_Breast_2fov_outs.zip \
  --source-dir /path/to/source/breast_source \
  --selection-protocol extensions/scaling/selection_protocol.json \
  --prepared-root /path/to/central400 /path/to/central600 \
  --output /path/to/source-audit.json
```

The extent-only region figure can be reproduced with
`python extensions/scaling/plot_regions.py --source-dir /path/to/source/breast_source`.
It reads masks and frozen crop coordinates, not transcript or count outcomes.

The independent scalar checks cover 16 hash-selected multicopy questions per
region; all 32 match exactly. The 400 µm sample was frozen before production.
The 600 µm sample applies the same previously fixed input-only hash rule after
production, before inspecting its spectra; it is not described as a frozen
pre-production sample. The checker enumerates assignments using scalar closure
and tent integration without importing the
production chain, cut or geometric integration routines. The saved selection,
reference and comparison are in `results/central400_sample_*.json` and
`results/central600_sample_*.json`.

```sh
python extensions/scaling/independent_sample.py compare \
  --atlas /tmp/geocert-checked/central400/data/atlas_breast_central_400.json \
  --selection extensions/scaling/results/central400_sample_selection.json \
  --reference extensions/scaling/results/central400_sample_reference.json \
  --spectra /tmp/geocert-checked/central400/results/patch_spectra.json \
  --output /tmp/geocert-sample-check.json
```

## Earlier partial workload replay

`benchmark.py` requests three fresh sequential repetitions of each increasing
prefix of the original four 200 µm crops: lung 1, lung 2, breast 1, breast 2.
It compares every regenerated atlas and exact minimum against the frozen
results. This earlier replay adds no tissue area.

| Completed prefix | Repeats | Cells | Ribbons | Median pipeline time | Peak stage RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Lung 1 | 3 | 79 | 276 | 17.70 s | 103.89 MiB |
| Lung 1 + 2 | 3 | 153 | 544 | 38.31 s | 104.89 MiB |

A three-crop replay failed CSV-completeness validation and is excluded from
accepted timings. The four-crop replay was not run. This replay series is
partial and is separate from the completed contiguous experiment.

```sh
python extensions/scaling/benchmark.py \
  --data-root /path/to/restored/recomb2027-atlas \
  --output-dir /path/to/fresh/workload-replay --repeats 3
```

All scripts require Python assertions to remain enabled; do not use `-O`.
