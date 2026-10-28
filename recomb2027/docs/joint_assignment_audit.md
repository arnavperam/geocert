# Exhaustive pairwise assignment audit

The count-endpoint null result remains unchanged. This separate, post-review experiment asks the question addressed by the original theorem: can two explicitly specified transcript reassignments occur together, and how much shared interface deformation do they require?

## Fixed universe and denominators

Every retained corridor is evaluated at both inherited geometry settings. For each model, include every individually feasible move of an observed transcript to a non-source cell. Pair all moves involving distinct transcript groups. Other transcripts are free collateral variables. There is no filtering by gene, desired effect, pair distance or cost outcome. The protocol is written before pair-area outcomes; the work is a post-review development analysis, not a preregistered biological study.

The eight inherited area budgets are 0, 0.05, 0.1, 0.25, 0.5, 1, 2 and 4 square micrometres. A model with fewer than two eligible distinct groups contributes no pair query but remains in the denominator of 42 models per setting.

| Setting | Models with pairs / all | Individual moves | Pair queries | Jointly infeasible | Jointly feasible | Joint cost above max singleton | Joint cost below sum singleton |
|---|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 8/42 | 48 | 80 | 22 | 58 | 51/58 | 23/58 |
| Maximal | 26/42 | 171 | 652 | 85 | 567 | 519/567 | 144/567 |

Seven of 42 baseline models contain an incompatible pair; seven contain a feasible pair with shared deformation savings. In maximal geometry these counts are 17 and 21, respectively. Models and queries share cells and are not independent biological replicates.

## Budget decisions

At an area budget of 1 square micrometre:

| Setting | Both singleton costs within budget | Joint pair within budget | Marginal false admissions | Sum-singletons plus common feasibility admits | Attainable pairs missed by additive test |
|---|---:|---:|---:|---:|---:|
| Baseline | 58 | 40 | 18 | 39 | 1 |
| Maximal | 435 | 315 | 120 | 304 | 11 |

The complete eight-budget table, all 732 pair records, all 219 move costs, and every attaining boundary or infeasibility chain are retained in `results/joint_assignment_audit/`. Reporting the budget-1 row does not define a new favorable evaluation grid; every declared grid point remains available. At budgets 1 and 4 the baseline additive test each misses one pair, while it misses none at the six other grid points.

## What the comparisons mean

Let c(A) be the minimum swept-interface area for pins A, allowing collateral reassignments. For individually feasible pins A and B, max(c(A),c(B)) is a sound lower bound on c(A union B). It is an outer feasibility test, not evidence that a common boundary exists. Incompatibility gives c(A union B)=infinity.

When A union B is feasible in this source-admissible implication family, c(A union B) is at most c(A)+c(B). Forward and backward closure of the union of pin sets equals the respective unions of their closures. Hence every positive source-relative lower or upper excursion forced by the combined query appears in at least one single-query envelope. The combined excursion in either direction is the pointwise maximum of the two nonnegative excursions. Integrating max(u,v) <= u+v proves the area upper bound. Feasibility is essential: incompatible upward and downward requirements cannot be combined into a boundary. Thus the additive comparison here explicitly checks common feasibility; a raw sum without that check has no such guarantee.

The inequalities max(c(A),c(B)) <= c(A union B) <= c(A)+c(B) are checked with rational arithmetic on every feasible pair. A strict upper inequality measures shared deformation; a strict lower inequality measures the failure of independent singleton budgets to control total deformation. These elementary inequalities are not claimed as a new theorem.

The primary algorithm uses source-clipped closure plus exact envelope realization. A separate replay compiles the event graph and solves all 732 pinned area objectives by integer-capacity minimum cut. Its full flow certificates are supplied. The script checks max-flow/min-cut equality and agrees with every saved feasibility decision and minimum area.

## What this does and does not establish

This audit exhibits a nonseparable joint-assignment problem on the original development geometry. It supports centering the paper on the ownerwise exact deformation representation. It does not establish useful marginal gene-count improvements, calibrated geometry, error prediction, a coherent whole-cell uncertainty family or a biological conclusion. Only 1 baseline pair and 22 maximal pairs involve the same gene; the overwhelming majority of tested interactions cross genes. The existing count-endpoint null result is therefore compatible with these interactions and must remain in the paper.

The source CSV was independently recounted: 14,384 exact interval rows, 1,798 corridor/setting/cell/gene questions, of which 1,584 involve one local transcript, 166 involve two and 48 involve three. Width frequencies are 12,597 zero, 1,679 one, 94 two and 14 three. This confirms the critique's principal count statistics.

## Complexity claim correction

The inherited weak-hardness reduction uses binary-encoded bundle multiplicities. It supplies no hardness theorem for explicitly listed unit molecules or for arbitrary planar cell tessellations. The finite penalty-gap example with weights 2,3,3 can, however, be expanded into eight explicitly listed co-located molecules. It remains a valid small counterexample to exact budget recovery by penalization; the encoding limitation applies to the asymptotic hardness inference.

## Replay

Run `python experiments/joint_assignment_audit.py` followed by `python experiments/joint_assignment_audit_verify.py` from the project root. Neither command accesses external data or runs the 72-crop study. Output directories are confined to this audit. The protocol and pair table have recorded SHA256 digests in the results.
