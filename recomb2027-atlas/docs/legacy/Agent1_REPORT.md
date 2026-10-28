# Atlas consolidation and physical-reach study

## Decision

Retain **Certified Expression Intervals under Shared Segmentation Uncertainty** and the Atlas model. The uploaded 0.2.0a2 revision has stronger release assets and actual Proseg records, but its corridor-local intervals and exponential target enumeration omit the Atlas composition and count-chain results. The final manuscript combines the stronger research structure with the useful audited release work. Simulated reviews are assessments, not acceptance evidence.

## Adopted and preserved

The MIT license is copied byte-for-byte from the supplied licensed source. Proseg inputs, exclusions, outputs, commands and the failed default run are preserved. The main theorem explicitly separates logical hard arcs from numerical hard capacities chosen after molecular costs. The general coupled proof is retained in the supplement; the main paper proves the scalar case used by the Atlas and count-chain experiments.

All 31 inherited `atlas/` and `geocert_expression/` Python files remain byte-identical to the Atlas predecessor. The four original atlas files, both primary interval tables and both primary global spectra are unchanged. New study drivers live in `analysis/reach/`. The native pair census reported by the other draft largely repeats earlier work; it is not promoted as a fresh contribution.

## Corrections to the critique

The original eight-budget denominators reproduce: exact intervals tighten 4 of 46 eligible lung cell–gene questions and 203 of 1,026 breast questions, affecting 4/153 and 139/467 cells. Eligibility means marginal width at least two for at least one tested budget. These are conditional computational questions, not biological replicates.

Non-singleton totals are constant from area 0.25 onward, but endpoints are not all constant: 61 breast intervals widen from 0.25 to 0.5 and three from 0.5 to 1. Last endpoint-saturation areas are approximately 0.240525 and 0.573939 µm². The previous inference that the budget never binds above 0.25 is too strong.

The primary ribbon permits displacement at most `3*p*B/(4*sqrt(2))`, approximately 0.4508 µm for `B=4`. This limitation is now explicit.

## Completed frozen width study

The protocol was saved before the new extraction/count outcomes. It fixes required physical support halfwidths 0.5, 1, 2 and 3 µm, with no narrower fallback, and a 17-budget grid. All original arc priorities, seam rules, slope limits, nuclear owners and observations are unchanged. Each width yields a separately extracted atlas, not a nested family. This is an outcome-frozen analysis on already seen crops, not held-out biological validation.

| Collection / required support | Patches | Editable edge coverage | Displacement cap (µm) | Last endpoint threshold (µm²) |
|---|---:|---:|---:|---:|
| Lung, adaptive primary | 544 | 68.3% | 0.451 | 0.241 |
| Lung, 0.5 | 314 | 45.1% | 0.366 | 0.154 |
| Lung, 1 | 249 | 33.8% | 0.733 | 1.498 |
| Lung, 2 | 185 | 23.6% | 1.493 | 3.090 |
| Lung, 3 | 135 | 16.9% | 2.226 | 4.878 |
| Breast format test, adaptive primary | 1,708 | 65.8% | 0.451 | 0.574 |
| Breast format test, 0.5 | 994 | 44.6% | 0.366 | 0.376 |
| Breast format test, 1 | 803 | 34.6% | 0.733 | 1.281 |
| Breast format test, 2 | 609 | 24.7% | 1.493 | 5.096 |
| Breast format test, 3 | 481 | 17.8% | 2.226 | 13.783 |

All eight new computations completed: 106,351 global count spectra and 1,807,967 interval rows. The largest patch contains 136 molecules; the maximum local gene multiplicity is 17. Larger reach is real within selected patches, but rejected wider ribbons freeze more boundary portions.

The original spectra were separately projected to the same fine budget grid for a fair comparison. Every original row remains identical. Finer sampling finds 12/46 tightened lung and 361/1,026 breast questions; these secondary numbers do not replace the original 4/46 and 203/1,026. Full per-width/per-budget denominators are supplied in `results/reach/` and the supplement.

## Proseg records and reference reach

The other agent ran Proseg 3.2.0 once per clipped lung crop after the default-factorization run failed on the first crop. We rechecked those executions, not reran or replicated Proseg. Exact transcript IDs, genes, original owners, observed XY values and the source-cell map agree. Independent reconstruction verifies all 1,443 and 67 nonzero output count entries. The outputs contain 5,774 and 1,482 transcripts, excluding 73 and 41; 2,484 and 1,335 output transcripts are background. These cropped, adapted runs are poor inputs for an accuracy claim and are treated only as model-reach diagnostics.

All widths use identical eligible comparison keys. The nuclear reference has a five-micrometre crop-edge guard. Missing-reference RNA induces the same source and possible-destination key exclusions in every family.

| Reference, fixed changed-count denominator | Adaptive primary: attainable | Required support 3: attainable |
|---|---:|---:|
| Lung Proseg, 2,798 | 185 (6.6%) | 107 (3.8%) |
| Lung nuclear reference, 1,114 | 184 (16.5%) | 111 (10.0%) |
| Breast nuclear reference, 10,965 | 1,773 (16.2%) | 1,333 (12.2%) |

These are unbounded-budget count attainability rates. Every finite-budget reach curve is retained as well. Complete reference assignment vectors are ruled out by fixed RNA or unavailable destinations in every family. Entrywise inclusion would not prove joint matrix realizability. Wider support therefore does not fix the family-adequacy issue; coverage loss dominates aggregate reach here. The result is retained rather than selecting a width that gives a more favorable story.

## Validation

Exact rational geometry verification completed for all eight new atlases. Source side labels, all point memberships, nuclear constraints and disjoint interiors were checked. Computation directly checks every local attaining boundary.

A separate composition routine rechecked all 106,351 global spectra, 131,376 witness traces and 1,807,967 CSV rows. A deterministic independent scalar sample checked 246 local multi-copy spectra and 1,628 target assignments with zero failures or resource limits. It does not independently enumerate every new local optimum. The new forced-bit preprocessing passed 160 spectrum comparisons on 80 models against unpruned DP and separate scalar enumeration, with 773 assignments and 287 witnesses.

The prospective ZIP range helper passed a controlled two-member/five-request test. That is a helper test, not a successful live download. The final frozen-package replay is recorded separately and must not be confused with recomputation of the full eight-width study or inherited timing censuses.

## Paper and release

The paper retains the Atlas title and self-contained scalar proof, gluing and count-chain arguments. The full coupled extension is in the typeset supplement. The title page is unnumbered. Figure 1 now includes observed transcripts, a physical scale and a legend; a reach/coverage plot shows every width. All budget tables, the original null results, explicit AI contribution and format-test qualification remain.

Version 1.1.0 is a local MIT-licensed publication candidate. The supplied PR #3 and its GitHub checks are public, but they do not publish or verify the new Atlas. The available connector exposes no write action. Publication instructions and CI are included; no public Atlas tag, merge, conference registration or submission is claimed.

## Remaining scientific limits

No real Janesick breast dataset, additional lung region, new donor, five-run Proseg comparison, anatomical calibration or independent molecular truth was acquired. Prospective acquisition code is included without experimental claims. The new study quantifies, rather than resolves, fixed-seam family inadequacy. A specialist must still judge the exact reduction's novelty and significance. This is a bounded theory-and-methods submission candidate, not a demonstrated biological correction tool or guaranteed acceptance.
