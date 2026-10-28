# GeoCert: exact joint deformation in ordered corridors

An exact finite graph represents the minimum continuous boundary deformation for every feasible transcript-assignment vector in a specified ordered corridor. Hard-budget gene-count extrema are a secondary, explicitly corridor-conditioned query.

**Scope is part of the answer:** observations outside the selected corridor remain fixed. Overlapping corridors define different uncertainty families; their intervals must never be presented as one jointly varying whole-cell certificate. The native marginal-count comparison finds zero baseline improvements over the singleton enclosure and only five improvements in relaxed geometry. The package does not claim demonstrated biological utility from these marginal intervals.

This candidate extends the reviewed GeoCert continuous-boundary engine. It does not estimate a corrected segmentation. It conditions on the supplied corridor, retained transcripts, fixed nuclei, and geometry. Native budgets in this release are sensitivity settings, not anatomically calibrated error probabilities.

## Complete results and replay

This browsable source tree omits bulky saved arrays. Download and extract [the complete source/results archive](GeoCert_RECOMB_Revision_0.2.0a2.zip) for the full 318-file artifact, hash manifest, all native frontiers, scoped matrices and raw pair certificates. Run the commands below from the extracted `GeoCert_RECOMB_Revision` directory. The ZIP SHA256 is `75ff637722c52df066d014a33262b89bb16f94d9b806f1b6d1fa195b8608f122`.

## Run an example

Python 3.12 or newer is required. The exact query workflow uses the standard library.

```sh
python -m geocert_expression.cli demo --output run/demo.json
python -m geocert_expression.cli certify geocert_expression/examples/native_count.json --output run/count.json
```

A wheel is included in `dist/`:

```sh
python -m pip install --no-deps dist/geocert_expression-0.2.0a2-py3-none-any.whl
geocert-expression demo --output demo.json
```

The example queries PECAM1 in the left cell of one source-derived corridor. Its interval is `[0,0]` at deformation 1 µm² and `[0,1]` at 2 and 4 µm². These are bounds on measured counts, with observations outside the corridor held fixed.

## Guarantees

`auto` computes the complete least-deformation spectrum for the scored transcript groups when the source is admissible and the pattern limit permits it. It keeps every other RNA free to move as geometric collateral. `bounds` uses exact rational min-cut objectives to produce Lagrangian outer bounds and feasible witnesses. The bounds are sound for every nonnegative multiplier but can be loose. A resource limit is not a preservation result.

The interval spans its exact minimum and maximum when `status` is `exact`; intermediate integers and endpoint combinations need not be attainable. A gene absent from the observation set is not proved biologically absent. See `docs/THEORY.md` for the scope of each guarantee.

```sh
python -m geocert_expression.cli certify input.json --method bounds --output bounds.json
python tests/test_budget.py
python tests/test_joint.py
python tests/test_binding.py
python validation/check_native_intervals.py
```

## Native count table

`results/expression_intervals.csv` contains 14,384 intervals from the unchanged 42 development corridors, two geometry settings, and eight physical budgets. The cell, gene, corridor, fixed exterior count and model are explicit. `results/native/` contains the exact score spectra and witnesses. The same donor and some transcripts recur across queries; these are not biological replicates.

```sh
python experiments/native_expression.py
python experiments/raster_relaxation.py
```

The second command requires NumPy and uses a **different** uncertainty model: maximum source-label transport displacement, followed by an independent-candidate relaxation. Its intervals apply to the measured RNAs in the supplied crops. They are not exact connected-boundary optima.

## Calibration and downstream use

`calibrate` implements the split-conformal order statistic for held-out independent region/slide scores, including infinite model-miss scores. Calibration/test labels must not be used to fit the family. If the result is unbounded, no finite geometric certificate is obtained. This release tests calibration only in simulation.

Fixed linear contrasts and joint expression-presence predicates are supported by the Python modules. Marginal positive upper bounds do not establish co-expression, ligand–receptor signaling, or statistical differential expression. No false-discovery-rate guarantee is inferred from interval arithmetic.

## Optional C++ components

Two ordinary reference backends are included, not new flow algorithms:

```sh
g++ -std=c++17 -O2 cpp/count_bounds.cpp -o run/count_bounds
g++ -std=c++17 -O2 cpp/exact_flow.cpp -o run/exact_flow
python tests/test_cpp_flow.py
```

The flow backend uses the system Boost multiprecision headers and returns a certificate checked by Python integers. It uses fresh Dinic flows; dynamic residual reuse and a whole-slide flow benchmark are not implemented. The count accumulator processes already constructed candidate lists; its million-record throughput excludes geometric candidate extraction.

## Reproduction

`python experiments/prepare_counts.py` reconstructs the included observed counts using NumPy, pandas and the inherited restricted Parquet reader. That shared ingestion dependency is disclosed. `python validation/replay.py --full --cpp --output run/replay.json` re-executes the scientific checks and full native census in a disposable copy, then compares the results. The C++ stage needs a C++17 compiler and Boost headers; omit `--cpp` for Python-only replay. Raw-data and raster checks are separate.

The self-contained mathematical proof is in `manuscript/main.pdf` and `manuscript/supplement.pdf`. Historical development notes are retained for provenance. The scientific inputs and engine provenance are in `docs/PREDECESSOR.json` and `docs/DATA_PROVENANCE.md`. This revision is published on the `recomb-theory-revision` review branch in the `recomb2027/` directory. Earlier releases remain preserved. Software is MIT licensed; data attribution remains separate in `NOTICE.md`.

The complete sparse count-matrix export is
`results/expression_matrix_intervals.csv`; read its companion
`expression_matrix_index.json` before treating absent entries as zeros. Fixed
positive entries are included without counting them as new optimization tasks.
Different corridor keys identify different conditional matrices.

Additional APIs are documented in the scientific supplement:

- `resources.resource_spectrum`: exact count spectra with an additional all-RNA change cap.
- `joint.outcome_spectrum`: jointly attainable counts, rather than Cartesian combinations of marginal endpoints.
- `evaluation.evaluate_aligned`: common evaluation for supplied aligned reference owners and method outputs; it does not run other segmenters.

The new native table is development-only. The reported calibration is a generated
holdout experiment, not calibrated confidence for those native intervals. See
`docs/VALIDATION_STATUS.md` for the completed and missing evidence.

## Revision checks

The exhaustive pair-assignment audit is `experiments/joint_assignment_audit.py`, independently checked by `experiments/joint_assignment_audit_verify.py`. It tests joint assignment costs, a different estimand from marginal per-gene count intervals. Read `docs/joint_assignment_audit.md` for all denominators. Scope safeguards and the overlap audit document the different conditional families; a union allows one corridor to vary at a time, never all corridors simultaneously.

Simulation is software validation only. Proseg execution on the two development crops is documented separately; agreement with source assignments is not ownership truth. The project still has one native development donor and lacks held-out biological calibration.

The original exact-deformation experiment records are retained under `results/original_joint_deformation/`. Earlier 72-crop raster results belong to a separate uncertainty model and are not reused as validation of this continuum theorem.
