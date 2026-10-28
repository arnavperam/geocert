# Comparison audit: Atlas and corridor revision

This audit reads the supplied tables without changing their values. Reproduce with `python analysis/audit_revision_csv.py` from the package root. `results/comparison_csv_audit.json` records input SHA256s, every budget, endpoint changes, denominators and the LDHB example.

## Overall assessment

Atlas has the stronger methodological/application structure: a physically disjoint compatible atlas yields one declared count family per reference cell, and its scalar count-chain is a polynomial hard-budget solver for explicitly represented point transcripts. The corridor revision supplies useful scope safeguards, a cleaner public release, actual Proseg execution records and additional independent replay. Subject to source/theory verification, the best combined paper uses Atlas as its scientific base and adopts the release and diagnostic work from the corridor revision.

These are fixed-seam, topology-preserving sensitivities. They are not unconstrained whole-cell uncertainty sets, and neither version demonstrates anatomical calibration, independent transcript truth or downstream statistical error control.

## Independently reproduced supplied CSV census

| Quantity | Lung | Breast format-test |
|---|---:|---:|
| Interval rows | 26,432 | 183,808 |
| Distinct cell/gene questions | 3,304 | 22,976 |
| Status `exact` | 26,432 | 183,808 |
| Rows strictly narrower than marginal outer interval | 4 | 228 |
| Questions narrower at any tested budget | 4 | 203 |
| Questions with maximal tested marginal width at least 2 | 46 | 1,026 |
| Tightened fraction within that targeted denominator | 8.70% | 19.79% |
| Tightened fraction of all questions | 0.12% | 0.88% |
| Cells appearing in table | 148 | 467 |
| Cells with any tightened question | 4 | 139 |

The report states 153 complete lung cells, with five represented as implicit zero rows; this explains the 148 cells actually appearing in the CSV. Breast 139/467 is 29.76%. The targeted denominator is scientifically relevant, but is not the uniquely correct denominator. Retain both it and the all-question denominator. It is defined here by the maximum marginal interval width over the declared eight-budget grid, not inferred biological susceptibility.

## Important correction to the pasted evaluation: budget saturation

A constant number of non-singleton intervals does not imply constant intervals. The breast CSV contains 61 questions whose endpoints change from 0.25 to 0.5 square micrometres and three additional changes from 0.5 to 1. The cited LDHB example itself changes from [19,23] to [17,23] between 0.25 and 0.5. Thus the statement that "above 0.25 the budget never binds" is false for breast.

| Budget (square micrometres) | Lung nonsingletons | Lung sum of widths | Lung tightened rows | Breast nonsingletons | Breast sum of widths | Breast tightened rows |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0.05 | 366 | 380 | 0 | 3,684 | 4,279 | 84 |
| 0.1 | 429 | 451 | 4 | 4,414 | 5,339 | 80 |
| 0.25 | 524 | 572 | 0 | 5,259 | 6,754 | 61 |
| 0.5 | 524 | 572 | 0 | 5,259 | 6,820 | 3 |
| 1 | 524 | 572 | 0 | 5,259 | 6,823 | 0 |
| 2 | 524 | 572 | 0 | 5,259 | 6,823 | 0 |
| 4 | 524 | 572 | 0 | 5,259 | 6,823 | 0 |

Correct language: on the evaluated budget grid, lung interval endpoints plateau from 0.25 onward; breast endpoints plateau from 1 onward. Both non-singleton counts plateau from 0.25. This is an endpoint-census statement, not proof that every full assignment or every untested budget is unconstrained. All observed tightening occurs at budgets at most 0.5.

## Physical reach

The supplied theory specifies

`Psi(t,h) = (p/2) (sx*(t+r(t)+h), sy*(t-r(t)-h))`.

The Euclidean displacement of a reference interface point under its fixed-t map is therefore `p*abs(h)/sqrt(2)`. The recorded band `abs(h) <= 3B/4` gives a maximum displacement `3*p*B/(4*sqrt(2))`. For the largest attempted halfwidth B=4 and the supplied lung pixel size 0.21250001266598778 micrometres, the cap is about 0.4507806 micrometres. Smaller accepted B produces smaller caps. This is a displacement bound under the chart map; the transverse direction need not be the local Euclidean normal along a raster corner.

Narrow reach is a genuine limitation. However, it has not been isolated as the sole reason nuclear-reference counts lie outside the family. Fixed seams, junctions, immutable labels, topology and exterior constraints can also exclude alternatives. A width/coverage/reach study can measure this dependence; causal attribution to width alone is currently unsupported.

## Correction to the corridor revision's novelty framing

The existing original supplement at the earlier repository commit3d44852a58582f350231bb838abaa7943192cd66, `manuscript/supplement.tex`, section "Joint feasibility and the native minimum-thickness obstruction" (around lines204–216), already reports:

- 80 baseline and652 largest-family pair questions;
- 22 and85 jointly infeasible;
- 58 and567 jointly feasible;
- 23 and144 with joint cost below the sum of singleton thresholds.

It also reports the finite two-target scalar bound matches56/58 baseline and567/567 largest-family costs, and own-boundary slope explains22/22 and84/85 infeasible pairs. Earlier lines state that all732 pair queries match the source-clipped shortcut.

Therefore the recent pair census is a reproduction and expanded diagnostic comparison. Its independently checked integer-flow certificates and fixed-budget admission analysis are useful verification/analysis, but the80/652,22/85 and23/144 results cannot be called newly discovered evidence. The main combined paper should not claim them as a new scientific contribution.

## Other cautions for merging

- The breast data remain vendor format-test material. They provide external computational transfer, not an independent biological cohort. Avoid hiding that limitation through an abstract edit.
- The CSVs contain no additive-inner columns, so the additive LDHB result and68 additive-miss rows require the source result files; they are not independently verified by these attached tables alone.
- Actual Proseg runs are an asset, but are adapted no-factorization runs on clipped same-donor inputs, with many background assignments. They can support an explicitly qualified reach diagnostic after matching source/mask identities and molecular populations, not an accuracy claim.
- A frozen width sensitivity protocol should predeclare values before running the new sweep, but must be called a post-review sensitivity study, not retrospective preregistration of the earlier results.
- Updating a semantic version to1.1.0 does not establish stability or validation. Versioning should reflect the software release status rather than presentation optics.
