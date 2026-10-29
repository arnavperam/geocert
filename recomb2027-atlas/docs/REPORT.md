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

These rates allow an unbounded deformation budget. Finite-budget curves are also recorded. Fixed molecules or unavailable destinations rule out the complete reference assignment vector in every family. Entrywise count inclusion alone would not establish simultaneous attainment of the reference matrix. When a single fixed-width atlas replaces the adaptive primary atlas, its reduced coverage outweighs its greater local displacement range on these comparison keys. The finite union below retains the alternatives from all five atlases.

## Finite union of complete atlases

The [finite-family extension](../extensions/atlas_union/README.md) combines the primary atlas and four fixed-width families. All members share the reference partition, observations and physical-area units. For each count, it takes the least exact cost among members and retains that member's complete attaining trace. Each witness uses one entire atlas; ribbons from different members are never mixed. Primary feasible states are preserved.

The [recorded spectra and intervals](../extensions/atlas_union/results/) contain 29,212 questions, 42,123 finite count states and 496,604 rows on the common 17-budget grid. The following fractions distinguish the original question census from the expanded census, which also includes genes editable only in another member.

| Collection and question census | Questions | Primary invariant at 4 µm² | Union invariant at 4 µm² |
|---|---:|---:|---:|
| Lung, original primary questions | 3,304 | 2,780 (84.1%) | 2,524 (76.4%) |
| Breast, original primary questions | 22,976 | 17,717 (77.1%) | 15,143 (65.9%) |
| Lung, complete union questions | 3,766 | — | 2,905 (77.1%) |
| Breast, complete union questions | 25,446 | — | 17,160 (67.4%) |

On the original questions, 298 lung and 3,595 breast intervals widen at 4 µm². Across the 17-budget grid, exact geometry tightens the member-wise independent-molecule envelope for 29/102 eligible lung and 791/2,076 eligible breast questions in the complete union census. On the original census alone, these values are 28/101 and 783/2,058. Eligibility requires outer width at least two at some tested budget. Complete per-budget counts are in the [budget summary](../extensions/atlas_union/results/budget_summary.csv).

| Unbounded-budget reference reach | Common changed-count keys | Primary attainable | Union attainable |
|---|---:|---:|---:|
| Lung Proseg | 2,798 | 185 (6.6%) | 315 (11.3%) |
| Lung nuclear reference | 1,114 | 184 (16.5%) | 299 (26.8%) |
| Breast nuclear reference | 10,965 | 1,773 (16.2%) | 2,942 (26.8%) |

These [reference comparisons](../extensions/atlas_union/results/reference_summary.csv) use exactly the earlier comparison keys and exclusions. They test membership in the discrete count spectrum, including any holes, rather than containment in an endpoint interval. Gains reflect alternative patch selection as well as displacement range. Coverage remains member-specific: the largest covered fraction and largest displacement cap can belong to different atlases. Fixed junctions remain part of every member, and countwise reach does not imply simultaneous realization of a complete reference matrix.

## Larger contiguous computational regions

The [contiguous-region experiment](../extensions/scaling/README.md) uses two central squares selected from the full breast format-test raster before extracting transcripts or computing counts. Their nominal sides are 400 and 600 µm; pixel-aligned sides are 399.925 and 599.888 µm. The same primary adaptive atlas rule produces the following completed results.

| Quantity | Central 400 µm region | Central 600 µm region |
|---|---:|---:|
| Crop-contained connected reference cells | 1,363 | 3,171 |
| Retained Q20 measured-gene transcripts | 226,983 | 470,293 |
| Disjoint ribbons | 5,102 | 11,748 |
| Transcripts within ribbons | 29,011 | 62,108 |
| Editable incident edges | 137,516 / 211,770 (64.94%) | 318,946 / 478,423 (66.67%) |
| Cell–gene spectra | 69,278 | 152,925 |
| Exact interval rows on eight budgets | 554,224 | 1,223,400 |
| Exact count backend | Chain | Assignment enumeration |

For the 400 µm region, atlas extraction from masks, rational geometry validation and chain count/witness computation take 152.46, 199.66 and 1,000.91 seconds, respectively: 22.55 minutes combined. Peak individual-stage RSS is 2,488.77 MiB (2.43 GiB).

For the 600 µm region, the existing exact enumeration backend was selected from input multiplicities: every local gene has at most seven copies. An initial attempt reached the global runtime deadline during count computation. The completed count retry reuses the validated atlas and its geometry audit. Completed extraction, geometry and final count stages total 1,783.07 seconds (29.72 minutes); cumulative measured wall time across both attempts is **2,213.11 seconds (36.89 minutes)**, including the interrupted attempt. Peak completed-stage RSS is 2,251.76 MiB (2.20 GiB). The preserved deadline record is distinct from the final result, which contains every required query and has no incomplete local spectrum.

Source preparation, coverage accounting and independent output reconciliation are measured separately. The [measurement record](../extensions/scaling/results/contiguous_summary.json) identifies stage times, both 600 µm attempts and the backend choice. These are completed workloads on a shared Linux host, with serialization included; the two backends and interrupted attempt preclude treating them as a controlled chain-scaling comparison.

The regions overlap each other and earlier crops within the same format-test specimen; their cell and RNA totals are reported separately. They measure larger contiguous workloads under the primary family, independently of the five-family union comparison. The [source audit](../extensions/scaling/results/source_audit.json) verifies both prepared populations against the original archive.

## Computational checks

Rational geometry checks cover the eight width-specific atlases: source-side labels, all point memberships, nuclear constraints and disjoint interiors. Every local attaining boundary is checked during computation.

A separate composition routine verifies all 106,351 global spectra, 131,376 witness traces and 1,807,967 interval rows. A deterministic independent scalar sample checks 246 local multi-copy spectra and 1,628 target assignments, with zero failures or resource limits. This sample does not enumerate every local optimum. Forced-bit preprocessing passes 160 spectrum comparisons on 80 models against the unpruned dynamic program and separate scalar enumeration, including 773 assignments and 287 witnesses.

Input guards reject unsupported free location intervals, multisite groups and overlapping owner thresholds. Source binding checks the atlas, observations and completed geometry audit before both primary and strict-width certification. The four primary atlas files, both primary interval tables and both primary global spectra remain unchanged. See [validation status](VALIDATION_STATUS.md) for the completed replay scope and CI records.

The union checks 76 checksum-pinned scientific inputs before merging. Its [independent verification](../extensions/atlas_union/results/independent_verification.json) reconstructs all 29,212 spectra, 42,123 count states and 496,604 interval rows; checks all 2,254,727 component rows on the common grid; and validates 39,191 explicit traces, 2,932 fixed-count completions and 18,762 distinct continuous local witnesses. All 68 population/budget summaries and three unbounded reference comparisons agree. Five invalid boundary mutations are rejected.

For both contiguous regions, the independent source audit verifies the original archive and seven extracted members, then reconstructs exact rational raster coordinates and the complete Q20 measured-gene census. Prepared masks, transcript identities, genes, coordinates, quality values, owners and source counts agree with the original data. The saved-output audit reconstructs every expected local and cell–gene question and independently repeats min-plus composition. For the 400 µm region it reconciles 39,519 recorded local witness counts/costs, 95,253 cell-state traces and all 554,224 interval rows; for 600 µm, the corresponding counts are 85,335, 208,555 and 1,223,400. Continuous witness checks and unsampled local optimality use the unchanged exact solver. Both completed datasets and the preserved retry records are checksum-bound in the [artifact manifest](../extensions/scaling/artifact_manifest.json).

For the 400 µm region, a separate [scalar check](../extensions/scaling/results/central400_sample_verification.json) matches all 16 input-selected multicopy spectra: 92 target assignments and 38 finite count states, with no errors or resource limits. Candidates have two through six target copies; the first 16 SHA-256-ranked query keys contain two through four copies. The reference enumerator uses scalar closure and tent integration without importing the production chain, cut or geometric integration routines. This sample complements the complete output audit.

The [600 µm scalar check](../extensions/scaling/results/central600_sample_verification.json) also matches all 16 spectra: 152 target assignments and 33 finite count states, with no errors or resource limits. Its selected queries contain two through six target copies. The same previously established input-only hash rule was applied after the 600 µm count run completed, before inspecting its production spectra; the 400 µm selection was recorded before production spectra. Both samples are limited to their selected local queries.

An [earlier workload replay](../extensions/scaling/results/partial_replay_summary.json) completes three repeats each for one and two original lung crops. The first three-crop attempt fails output reconciliation and is excluded from accepted timings; the four-crop workload has no completed measurement. These partial replay records are separate from the completed contiguous-region experiment.

## Scope

The breast inputs are the vendor's trimmed format-test data. The experiments provide no additional donor, published breast cohort, anatomical calibration, independent molecular truth or downstream biological validation. Fixed seams and junctions limit the represented changes. The Proseg diagnostics use two cropped, adapted executions and do not establish replication or accuracy.

The complete breast format-test archive was acquired and checked for the contiguous-region experiment. The separate proposed five-run Proseg experiment remains unexecuted and contributes no paper results. Its prospective acquisition helper has only the previously recorded controlled two-member/five-request test. Original study records and later verification runs remain distinguished in the stored provenance.
