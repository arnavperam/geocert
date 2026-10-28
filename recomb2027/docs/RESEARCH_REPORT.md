# Executed extension: certified expression intervals

This release changes the output from a deformation optimizer to a cell-by-gene
count certificate under an explicit hard physical budget. It is an author-review
research candidate. The broad biological program in the request has not been
completed, and no acceptance, perfection or historical-priority claim is made.

## Implemented outputs

`geocert_expression` provides exact rational Lagrangian count bounds, attained
boundary witnesses, explicit gaps, exact target-specific count spectra, joint
count predicates, two-resource physical/RNA budgets, and all-penalized-optimum
residual DAGs. Target projection keeps collateral genes free and enforces one
shared geometry. It does not simply add independent transcript vulnerabilities.

The output schema distinguishes `exact`, `bounded`, `empty_uncertainty_set`,
`model_infeasible` and `feasibility_unresolved`. An incomplete target enumeration
cannot be called exact. Persisted cut queries can be bound to the original
physical model, score, multiplier, owner bits and rational boundary.

## Native expression results

The new count extraction checks all 7,370 retained transcript records against the
source masks and all 464 repeated corridor records against their source identities.
The unchanged 42 corridors are evaluated under baseline and the inherited enlarged
geometry. There are 1,798 cell/gene questions at eight budgets, producing **14,384
exact intervals**. The solver checks 25,172 penalized cuts. The maximum local gene
multiplicity is three; the exact projection examines 2,358 target assignments with
1,076 unique minimum geometric completions.

The sampled multiplier bounds and their attained witnesses close **14,371**
intervals on their own. Exact target projection closes the remaining **13**.
The outer endpoints alone match truth in 14,379 cases; equality to an independently
computed answer is not the same as a self-contained primal/dual closure. Five
sampled outer bounds are loose, while eight more need an attaining endpoint.
A larger multiplier set might close some of these; they are not called intrinsic
duality gaps.

A separately implemented checker uses the old scanning compiler, enumerates
5,182 full target pin patterns, solves 1,346 distinct pinned-cut queries and
reconstructs 263 unique minimum boundaries. It verifies all 1,798 spectra,
14,384 interval rows and 1,176 saved root flow certificates without discrepancies.
Geometry primitives and ingestion are still shared dependencies.

The sparse full-matrix export contains 32,480 rows. Of these, 18,096 are unchanged
positive counts outside the edited gene set; they are algebraic copies, not new
optimizations. Absent entries are exact zero only within the recorded 335 observed
genes, retained molecule set and listed corridor/cell/setting/budget. Different
corridors remain different conditional matrices.

At deformation budget 1 square micrometre, 80 of 899 baseline intervals and 265 of 899
enlarged-family intervals are non-singleton. At budget 4 the counts are 94 and 307.
This is conditional sensitivity in one development donor, not a prevalence estimate
of actual expression error. Counts outside one corridor remain fixed; crop-border
cells are not asserted to have complete whole-cell profiles.

### Two concrete outputs

For PECAM1 in cell 153124 of corridor 4037abe9fea38e, the observed count is 0 and the
exact interval is [0,0] at budget 1, then [0,1] at budgets 2 and 4. The delivered witness
attains count 1 and is independently checked against the physical geometry.

For the previously identified corridor c8f8188430d092 under the enlarged family
at budget 4, MDM2 in cell 153965 and TFPI in cell 153959 each have interval [0,1].
They cannot both become positive together under the original minimum thickness.
Relaxing only that additional minimum changes the joint predicate from `never` to
`variable`. This is a joint count exclusion, not a cell-state or signaling finding.

### Strong simple-baseline result

The post-census attribution retains every original question. A sound marginal
owner-domain envelope is already tight except for **five enlarged-family**
questions; none tightens at baseline. A feasibility-matched sum of singleton
costs misses no count endpoint on this census. The latter is an INNER
approximation, not generally a safe exclusion certificate. These results show
that most native single-gene queries are easy. We do not claim a broad empirical
advantage over simple marginal analysis.

## Hard budgets and theoretical scope

A constructed continuous example has three RNA bundles with count gains 2,3,3 and
separated minimum areas 4,9,9. At budget 10, the exact maximum is 3, while the best
nonnegative penalty bound is 4. An explicit mixture proves the persistent gap;
exact target projection returns 3. This prevents the central overclaim that one
penalized cut automatically solves a hard-budget count query.

The source-feasible projection is exact in a number of assignments exponential
only in the scored groups. A source-clipped completion also minimizes all-RNA
source-owner-change weight, giving an exact extension for two hard budgets.
Min-plus composition is exact only for independently gluable patches with disjoint
observations, fixed compatible traces and connected cores. No native overlapping
corridors were silently assembled into a whole tissue.

A succinct weighted-count knapsack embedding establishes weak hardness in a
restricted disjoint-patch model. It is not a new planar three-label theorem or
hardness for explicitly listed unit RNA. Parametric cut nesting is established
only under monotone-terminal conditions; a middle-cell counterexample is supplied.
All-optimal residual DAGs refer to penalized optima, not all budget-feasible labels.
The proofs and limitations are in `THEORY.md`.

## Calibration and probabilistic bridge

A training-only synthetic fit uses 30 generated paired-boundary units; 99 calibration
units determine a 90% split-conformal threshold. On 200 independently generated
test regions,184 true owner vectors and 185 complete count matrices are covered.
The simultaneous matrix coverage is 92.5%, with a Wilson 95% interval
[0.8799648,0.9540251], and mean marginal width 0.71. The source point matrix is
correct in 97 regions. A misspecification guard retaining 20 infinite calibration
scores correctly refuses a finite geometric threshold.

This is simulation-only calibration with known generated owners, not a membrane
image/segmenter benchmark or biological confidence validation. The transcript-level
AUROC is 0.8655 for exact fragility and 0.8612 for distance to the flat source;
within-region dependence and the small difference preclude a superiority claim.

A separate enumerated four-state profile-Gibbs experiment compares full-state
Gumbel sampling with independent unary perturb-and-MAP. Across three 2,500-draw
replicates, full-state sampling has total-variation discrepancy 0.0044--0.0144;
unary perturbations have 0.0732--0.0784. This is a deliberately adverse example,
not a posterior benchmark against ResolVI. It demonstrates why randomized cuts
are not automatically exact or calibrated posterior samples. Profiling over
boundaries is also different from marginalizing their volume.

## Whole-raster relaxation and engineering

A distinct, conservative model permits an arbitrary2D/3D label topology. Candidate
owners are all reference labels within a maximum label-transport radius, enlarged
for RNA localization error. Dropping coherent-map constraints gives sparse outer
count intervals, with an optional global RNA-change cap. The radius and the
corridor's summed-area budget are not interchangeable.

All 7,370 crop transcripts were processed at radii 0,0.25,0.5,1 micrometre. Outside-crop
windows include all known labels plus an unknown-exterior label; label0 denotes
reference background/unassigned space, not an additional biological cell. The
method models retained crop molecules only. The 3D adapter was checked on generated
label arrays, not new 3D tissue.

The optional C++ accumulator matches 12,668 small-test interval entries. Three fresh
processes at one million PRECOMPUTED candidate records have median 0.341 seconds
for reading, accumulating and writing counts and maximum recorded RSS 6.5MiB.
This omits neighborhood construction, images and segmentation; it is not
whole-slide-in-minutes evidence.

The optional Boost-multiprecision C++ Dinic backend passes 100 graph optima against
exhaustive cuts, including 10 tests beyond 64-bit capacities, and all 252 archived
native optional objectives with independent flow and physical checks. No dynamic
residual reuse, BK/pseudoflow comparison or general speed advantage is claimed.
The exact Python package needs only the standard library; the raster/source
analysis adapters additionally use NumPy/pandas. A wheel and CLI are supplied.

## Validation inventory

-60 small physical models:4,350 full owner vectors,183 score minima,398 budget
 enclosures and 600 exact penalized cuts; no discrepancies.
-20 joint models:185 joint-outcome minima and 100 all-optimal marginal-domain
 comparisons; no discrepancies.
-20 two-resource models:80 complete spectra and 480 rational bound checks; none fail.
-120 candidate-label cases:2,880 budget/key checks and 49,450 independent feasible
 assignments; none fail.
-120 module compositions:12,924 Cartesian vectors and 3,040 budget checks; none fail.
-120 extreme-cell nesting comparisons plus an explicit nonnested middle-cell case.
-908 generated 2D/3D membership checks; eight malformed query/certificate mutations
 rejected; source-inconsistent universal-region cases distinguish empty from bounded.

These are internal checks, not an external proof review. The frozen extraction
record states exactly what was rerun and which scientific outputs matched.

## What was not completed

No new donor was acquired; no reporter, sex-mismatch or xenograft ownership truth
was measured; no deployed segmentation tool was run against common held-out truth;
no real parameter calibration, DE/FDR study or biological claim reanalysis occurred.
The library includes an aligned-output evaluation harness and a conditional
familywise-error gate, but does not produce valid DE p-values by assertion.
There is no implemented exact general-tessellation or arbitrary 3D continuum solver,
validated SpatialData/napari integration, or 10^8-transcript end-to-end benchmark.
`VALIDATION_STATUS.md` maps every requested pillar to the evidence actually present.

The previous ownerwise result and event integrator are retained. This extension
uses classical duality, convolution, conformal calibration and residual-cut
principles. Its useful output is a correctly qualified expression certificate;
these facts do not establish historical novelty, flawlessness or RECOMB acceptance.
