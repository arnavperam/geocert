# Reproducing the exact workflow

## Environment

The reviewed release was executed locally with CPython 3.13.5. A separate supplied audit reports the same showcase, tests, and contracted workload under CPython 3.12.14. These versions are also the basis for the 3.12/3.13 CI matrix; a workflow file is not itself evidence of a successful run.

No external Python packages are required by `geocert.py`, `core/`, or `tests/`. Historical numerical comparisons used additional scientific libraries and a different environment. They are not rerun by this repository's commands.

## Commands and expected output

```sh
python geocert.py verify-manifest
python geocert.py showcase --output run/showcase.json
python tests/test_scc.py
python tests/test_prefix_bound.py
python geocert.py benchmark --repetitions 1 --scc --output-dir run/comparison
```

The showcase's exact deformation is `143227474810907125420487/40850916262281216000000`; its exact objective at penalty 1 is `-61027106500498954579513/40850916262281216000000`. The benchmark must report 84 models, 252 queries, and objective hash `66f6a364012ce85e587c9d71a3005c57d9b7a845d83b2513f9cdce5018e8017f` for each representation.

The recorded SCC suite covers 200 graphs, 1,523 feasible-cut checks, and 72 optional-objective comparisons. The prefix suite covers 168 arrangements and 1,260 assignment checks. These finite tests support implementation correctness; they are not independent proof review or biological validation.

Use a new output directory for a new benchmark. The benchmark launches each representation in a fresh process and verifies agreement of exact objectives. It includes compilation, optimization, reconstruction and checking, but excludes input parsing/imports. Graph reuse is compile-once, not warm-started flow. Runtime depends on the machine. The reviewed local workload took seconds to tens of seconds, with historical total process memory around 205 MiB including additional imports; those are observations, not resource guarantees.

## Recorded results

`results/summary.json` is the unchanged RC2 summary. `results/model_timings.csv` reports per-model medians across the three archived repetitions for each method/contraction setting. These are derived from the 12 recorded runs, not a fresh timing experiment. Per-model medians must not be summed to reproduce the median aggregate runtime; the latter is retained separately in the unchanged summary. `results/run_metadata.json` retains each original file's SHA-256, run total, objective hash, and peak process RSS. `results/environment.json` is the original measured environment.

The paper's main comparison uses equal hard-implication SCC contraction. The full 84-model denominator includes 48 models with optional ownership and 36 fully fixed models. The latter favor the height implementation because of its cached boundary; this subgroup is not omitted. There is no universal solver-superiority or anatomy-accuracy claim.

## Artifact boundaries

This public repository contains the current exact workflow and its derived native inputs. It is not the complete historical development archive or a fresh raw-Parquet extraction. The frozen RC2 archive is identified by hash in `docs/RELEASE_PROVENANCE.json`; numerical studies, source-mask crops and previous exploratory implementations remain in that separately distributed artifact. The exact algorithms, CLI, test code, and three data files included here are byte-identical to RC2.

The checksum manifest covers committed release files, excluding itself. Run outputs are intentionally ignored. Do not overwrite recorded scientific outputs when comparing new measurements; store new results under `run/`.
