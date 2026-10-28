# Minimal correctness guards carried forward

Agent 1's 1.1.0 input omitted the prior point-input contract fix and prepared-source binding. These are now retained in the evaluation candidate, with no new algorithm or experiment.

## Changes

- `atlas/chain_spectrum.py`: reject free horizontal or longitudinal location boxes, free multisite owner groups, and overlapping free-owner thresholds. True point transcripts remain supported in all three existing location modes. Fixed structural regions remain allowed.
- `atlas/binding.py`: exact atlas digest, completed saved geometry-audit identity, source-file hashes, coordinate transform, source count reconstruction, unique patch and transcript ownership, source labels/genes/coordinates and physical Jacobian checks.
- `analysis/certify_cells.py` and `analysis/reach/certify_width.py`: require those checks before producing count output; reject Python optimized mode because downstream geometric witness checks use assertions. The strict-width driver binds the original collection source files through the existing restored shared-input paths.
- `atlas_cli.py`: adds the three regression scripts to `validate`.
- `tests/test_atlas_binding.py`, `tests/test_chain_contract.py`, `tests/test_width_binding.py`: regression coverage.

## Executed checks

All 4 primary and 16 strict-width prepared atlases pass source/audit binding. The tests reject 17 invalid primary binding cases, 3 invalid strict-width entry/source cases, and 7 unsupported free location or multisite cases. Three supported point modes match separate fixed-owner minima.

The original chain checks pass 201 spectrum comparisons; 80 nonflat models pass 1,644 exhaustive owner comparisons. The strict-width forced-bit preprocessing passes 160 spectrum comparisons on 80 models, including 773 independently enumerated target assignments and 287 checked pruned witnesses. No full eight-width optimization sweep was rerun by this guard task.

All original Python files outside the listed owned changes remain byte-identical to the supplied archive. Geometry extraction, forced-bit preprocessing, spectrum recurrence, budget grid and scientific result tables are unchanged. The original validation records were preserved, with new execution records in `results/final_guard_checks/` and the three new guard JSON reports in `results/`.

A saved geometry audit is a reproducibility attestation, not a cryptographic defense against someone who forges both inputs and audit files. New geometry still requires the independent rational geometry validator. These checks establish the declared input contract; they do not establish anatomical truth or peer-review acceptance.

The first copy-based guard attempt retained old modification times and was reverted by workspace synchronization; failing regressions exposed that before any release. The final fresh-timestamp writes pass and have hashes recorded in `results/final_guard_checks/summary.json`.
