# Atlas experiments

This report summarizes the prepared-data experiments for **Certified Expression Intervals under Shared Segmentation Uncertainty**. Counts and deformation costs are exact within each recorded atlas. The uncertainty families fix patch seams, uncovered boundaries, junction neighborhoods and nuclear owners.

## Primary count intervals

The primary eight-budget analysis tightens the marginal interval for 4 of 46 eligible lung cell–gene questions and 203 of 1,026 breast questions, affecting 4/153 and 139/467 cells. Eligibility requires a marginal width of at least two for at least one tested budget. These are computational questions; repeated budgets are not biological replicates.

Non-singleton totals are constant from 0.25 µm² onward, while some endpoints still change: 61 breast intervals widen from 0.25 to 0.5 and three from 0.5 to 1. The last endpoint-saturation areas are approximately 0.240525 and 0.573939 µm² for lung and breast. The primary ribbon permits displacement at most `3*p*B/(4*sqrt(2))`, approximately 0.4508 µm for `B=4`.

## Physical-width study

The [protocol](PREREGISTRATION.md) fixes required support halfwidths of 0.5, 1, 2 and 3 µm, with no fallback to narrower ribbons, and a 17-budget grid. It was recorded before these extraction and count outcomes. Arc priorities, seam rules, slope limits, nuclear owners and observations match the primary analysis. Each width produces a separate atlas; the families are not nested. The crops had already been examined, so this study is not held-out biological validation.

| Collection / required support | Patches | Editable edge coverage | Displacement cap (µm) | Last endpoint threshold (µm²) |
|---|---:|---:|---:|---:|
| Lung, adaptive primary | 544 | 68.3% | 0.451 | 0.241 |
| Lung, 0.5 | 314 | 45.1% | 0.366 | 0.154 |
| Lung, 1 | 249 | 33.8% | 0.733 | 1.498 |
| Lung, 2 | 185 | 23.5% | 1.493 | 3.090 |
| Lung, 3 | 135 | 16.9% | 2.226 | 4.878 |
| Breast format test, adaptive primary | 1,708 | 65.8% | 0.451 | 0.574 |
| Breast format test, 0.5 | 994 | 44.6% | 0.366 | 0.376 |
| Breast format test, 1 | 803 | 34.6% | 0.733 | 1.281 |
| Breast format test, 2 | 609 | 24.7% | 1.493 | 5.096 |
| Breast format test, 3 | 481 | 17.8% | 2.226 | 13.783 |

All eight computations completed, producing 106,351 global count spectra and 1,807,967 interval rows. The largest patch contains 136 molecules, and the maximum local gene multiplicity is 17. Wider supports allow larger displacements within accepted patches but leave more boundary portions fixed.

The original spectra were also projected to the fine budget grid. Every original row remains identical. Finer sampling finds 12/46 tightened lung questions and 361/1,026 breast questions; these secondary counts are separate from the original eight-budget results. Complete denominators by width and budget are in `results/reach/` within the reconstructed scientific archive and in the [supplement](../manuscript/supplement.pdf).

## Alternative-reference reach

Proseg 3.2.0 completed one adapted execution per clipped lung crop after its default factorization failed on the first crop. The audit rechecked the stored executions: transcript IDs, genes, original owners, observed coordinates and the source-cell map agree. Independent reconstruction verifies all 1,443 and 67 nonzero output count entries. The outputs contain 5,774 and 1,482 transcripts, excluding 73 and 41; 2,484 and 1,335 output transcripts are assigned to background. These runs provide model-reach diagnostics rather than an ownership-accuracy benchmark.

Every width uses the same eligible comparison keys. The nuclear reference has a five-micrometre crop-edge guard. Molecules missing from a reference induce the same source and possible-destination key exclusions across all families.

| Reference, fixed changed-count denominator | Adaptive primary: attainable | Required support 3: attainable |
|---|---:|---:|
| Lung Proseg, 2,798 | 185 (6.6%) | 107 (3.8%) |
| Lung nuclear reference, 1,114 | 184 (16.5%) | 111 (10.0%) |
| Breast nuclear reference, 10,965 | 1,773 (16.2%) | 1,333 (12.2%) |

These rates allow an unbounded deformation budget. Finite-budget curves are also recorded. Fixed molecules or unavailable destinations rule out the complete reference assignment vector in every family. Entrywise count inclusion alone would not establish simultaneous attainment of the reference matrix. Under this extractor, the loss of coverage at larger widths outweighs their greater local displacement range.

## Computational checks

Rational geometry checks cover the eight width-specific atlases: source-side labels, all point memberships, nuclear constraints and disjoint interiors. Every local attaining boundary is checked during computation.

A separate composition routine verifies all 106,351 global spectra, 131,376 witness traces and 1,807,967 interval rows. A deterministic independent scalar sample checks 246 local multi-copy spectra and 1,628 target assignments, with zero failures or resource limits. This sample does not enumerate every local optimum. Forced-bit preprocessing passes 160 spectrum comparisons on 80 models against the unpruned dynamic program and separate scalar enumeration, including 773 assignments and 287 witnesses.

Input guards reject unsupported free location intervals, multisite groups and overlapping owner thresholds. Source binding checks the atlas, observations and completed geometry audit before both primary and strict-width certification. The four primary atlas files, both primary interval tables and both primary global spectra remain unchanged. See [validation status](VALIDATION_STATUS.md) for the completed replay scope and CI records.

## Scope

The breast inputs are the vendor's trimmed format-test data. The experiments provide no additional donor, published breast cohort, anatomical calibration, independent molecular truth or downstream biological validation. Fixed seams and junctions limit the represented changes. The Proseg diagnostics use two cropped, adapted executions and do not establish replication or accuracy.

Prospective acquisition and five-run Proseg scripts are included separately. Their ZIP range helper passed a controlled two-member/five-request test; the proposed live acquisition and five-run experiment were not executed and contribute no paper results. Original full-study records and later verification runs are distinguished in the stored provenance.
