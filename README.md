# GeoCert-ST

Exact joint deformation costs for continuous RNA-ownership sensitivity.

Companion code for **Exact Joint Deformation Costs for Continuous RNA-Ownership Sensitivity**, by Arnav Peram. The method represents the minimum integrated deformation for every feasible ownership vector in a local ordered-cell corridor. It includes a direct owner-layer construction and an exact critical-mesh height reference.

## Run the native example

Python 3.12 or 3.13 is sufficient. The current exact workflow uses only the Python standard library; there are no packages to install.

```sh
git clone https://github.com/arnavperam/geocert.git
cd geocert
python3 -m venv .venv
. .venv/bin/activate
python geocert.py verify-manifest
python geocert.py showcase --output run/showcase.json
```

On Windows, use `.venv\Scripts\python.exe` in place of `python`. Run without `-O`; scientific assertions must remain enabled.

The example is corridor `4037abe9fea38e` at deformation penalty 1. Both exact representations return five RNA reassignments, deformation approximately 3.506102 square micrometres, and joint-model objective approximately -1.493898. The JSON output contains exact rational values, a continuous boundary, and a flow/cut certificate.

## Reproduce the representation comparison

```sh
python geocert.py benchmark --repetitions 1 --scc --output-dir run/comparison
python tests/test_scc.py
python tests/test_prefix_bound.py
```

Each representation evaluates 84 corridor/setting models and 252 optional objectives. Expected objective hash:

```text
66f6a364012ce85e587c9d71a3005c57d9b7a845d83b2513f9cdce5018e8017f
```

Use `--repetitions 3` for repeated fresh-process measurements. Omit `--scc` to run the uncontracted ablation. Timings and tied choices may vary; exact objective values must agree. See [reproduction notes](docs/REPRODUCIBILITY.md) for scope and recorded environments.

## Repository contents

- `core/`: exact geometry, direct and height constructions, flow, contraction, and witness checks.
- `data/`: the 42 source-derived rational corridors, expected objectives, and separately enumerated showcase costs.
- `tests/`: bounded correctness tests.
- `results/`: reported summaries, per-model timing records, and validation outputs.
- `docs/`: mathematical notes, data attribution, and release provenance.

## Scope

This is a research implementation for model-conditioned assignment sensitivity, not an estimator of corrected anatomical ownership. The inputs are selected local corridors from one development donor. Boundaries are ordered graphs, and the geometric uncertainty settings are not anatomically calibrated. The theorem covers nonadditive geometric deformation and unary molecular costs; it does not automatically solve arbitrary higher-order molecular objectives or hard-budget problems.

Version 1.0.0 packages the reviewed RC2 exact implementation without changing its algorithms or scientific inputs. The [provenance record](docs/RELEASE_PROVENANCE.json) identifies the frozen predecessor and the documentation-only correction. Historical numerical workflows and the large development archive are not included in this compact source release.

## Citation and rights

Citation metadata is in `CITATION.cff`. The manuscript is not represented as accepted or published. Data attribution and modification details are in [DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md). No software license has been selected; see [NOTICE.md](NOTICE.md).
