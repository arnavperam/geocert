# Certified expression intervals under shared boundary uncertainty

## 1. Observations, family and the object being certified

A measured RNA group i has an integer gene multiplicity vector w_ig and one cell
owner c_i. Single observed molecules have multiplicity one. Counts are

    N_kg(c) = n_fixed,kg + sum_i w_ig 1[c_i = k].

The fixed contribution is explicitly supplied; for the native experiment it
contains retained crop observations outside one editable corridor. All RNA inside
the corridor, including other genes, remains part of the shared geometric model.
We do not certify molecules that were not detected or genes outside the assay.

Let F be an explicitly declared geometric family and D*(c) its least deformation
realizing c, with infinity for infeasible owners. In the implemented continuous
family, cell boundaries are ordered Lipschitz graphs with bands, fixed endpoints,
minimum/maximum thickness and immutable nuclear regions. D is the sum of the
integrated absolute movements of its interfaces, in square micrometres. In a
multicell corridor this counts a location more than once if it crosses multiple
interfaces; it is not necessarily tissue area changing labels counted once.

For delta >= 0 define C_delta = {c: D*(c) <= delta}. The exact endpoints are

    L_kg(delta) = min_(c in C_delta) N_kg(c),
    U_kg(delta) = max_(c in C_delta) N_kg(c).

An empty family is an error/model result, not a zero-width biological certificate.
An interval denotes extrema; intermediate integer values need not all be
attainable. Endpoints for different entries need not be attainable together.

## 2. The inherited representation and its proof obligation

The previous GeoCert engine gives a finite rational-capacity graph that preserves
D*(c) for EVERY feasible owner vector before selecting molecular costs. The
continuous proof propagates all lower and upper ownership requirements through
Lipschitz and thickness inequalities. Its envelope families L,U and the source
f0 satisfy the same geometric constraints, which are closed under componentwise
minimum and maximum. Thus clamp(f0,L,U) is jointly feasible and pointwise nearest
to f0. Exact integration of its upper/lower excursion envelopes gives nonnegative
coverage coefficients, encoded by standard OR helpers. This is why arbitrary
unary molecular costs can be added without changing the feasible family.

The integration implementation uses the preserved event-driven routine. Its
ordering events and lazy exact integration change construction work, not the
ownerwise function. The current extension does not claim to invent maximum flow,
Lagrange duality, conformal prediction, coverage graphs, or min-plus convolution.
The complete earlier theory is retained in PREDECESSOR_THEORY.md. Graph proof
validity, physical compilation and physical witness validity are separate checks.

## 3. Sound hard-budget bounds

Let S(c)=a+sum_i s_i(c_i) be any integer-valued additive score. The fixed offset a
can include unchanged observations. For each rational mu >= 0 an exact cut gives

    m_plus(mu)  = min_c [mu D*(c) - sum_i s_i(c_i)],
    m_minus(mu) = min_c [mu D*(c) + sum_i s_i(c_i)].

THEOREM 1. Every c in C_delta satisfies

    a + ceil(m_minus(mu) - mu delta) <= S(c)
       <= a + floor(mu delta - m_plus(mu)).

PROOF. The definition of a minimum gives mu D*(c)-sum s_i(c_i) >= m_plus.
Rearrange, use D*(c)<=delta and mu>=0, and round downward because S is integer.
Apply the same argument to -S for the lower endpoint. Each multiplier is valid
independently of how it was selected. The maximum lower bound and minimum upper
bound over any finite multiplier set remain sound. QED.

An attaining boundary within budget gives a primal score. The software reports

    L_bound <= min S <= L_witness,
    U_witness <= max S <= U_bound.

Only equality at both endpoints earns status `exact`. Otherwise it reports
`bounded`, not preservation. A source inconsistent with universal RNA regions is
handled by an initial exact minimum-deformation solve; the source-clipping
shortcut below is unavailable. If even the minimum cost exceeds delta, the
uncertainty set is empty. All rounding and endpoint comparisons use rationals.

### A gap that no penalty can close

Three independent bundled groups can add counts 2,3,3 to one cell. Their separated
slope-one tents have exact areas 4,9,9. At delta=10 the maximum is 3: any pair costs
at least 13. At mu=1/3, the Lagrangian upper bound is 4. It cannot be improved by
another multiplier: a mixture with weight 2/3 on (D,S)=(4,2) and 1/3 on (22,8)
has expected (10,4); every affine penalty upper bound is at least that mixture's
score. Thus a complete penalty path does not recover every budgeted optimum.
This is a standard unsupported-outcome phenomenon, instantiated and executed in
an exact continuous ownership family. The exact refinement returns 3.

## 4. Exact projection onto the transcripts that affect the question

Assume the source realizes all source-owner and structural requirements at zero
deformation. Geometric source feasibility alone does not suffice. Pin some owner
labels. In the implication graph, propagate required true bits forward and
required false bits backward. A collision proves infeasibility. Otherwise retain
the source on every unforced bit. This is a feasible completion: mixed
forced/unforced implication violations are excluded by propagation, and an
entirely unforced implication holds at the source.

Every newly activated positive displacement requirement in this completion is
forced by the pins in every feasible completion. All other source-active
requirements have zero source-relative excursion. The canonical envelope of the
source-clipped completion therefore has pointwise no larger movement than every
other completion. This proves minimum deformation for those pins, including
necessary collateral RNA changes.

THEOREM 2. For an integer score with r owner-dependent groups, enumerate their
at most (m+1)^r assignments and apply the minimum completion above. Recording

    d(v) = min { D*(c) : S(c)=v }

for each attainable score preserves every hard-budget extreme. The cost is
O((m+1)^r T_completion), with cheaper pruning by individual feasible domains.
Only target groups are enumerated; the remaining molecular owners are eliminated
by an exact geometric completion, not fixed arbitrarily.

PROOF. Every complete feasible owner vector restricts to one enumerated target
assignment. Its cheapest completion has the same score and no larger cost.
Conversely every accepted completion is a full feasible vector attaining its
stored score and cost. Minimizing over completions gives d(v) in both directions.
Consequently min/max {v:d(v)<=delta} are the exact endpoints. QED.

This is a query-local exact algorithm, not an FPT theorem in physical budget
alone. Binary-encoded multiplicities are retained in the score, not expanded by
assumption. Resource caps yield an incomplete spectrum; its feasible samples can
improve witnesses but cannot establish missing outcomes. The code refuses to
present an incomplete spectrum as exact.

The same projection computes a joint vector of counts. A predicate can then be
certified `always`, `never`, or `variable` over all budget-feasible joint outcomes.
Marginal upper endpoints alone do not establish co-expression or a ligand--receptor
expression conjunction. The native MDM2/TFPI example illustrates this distinction.
Neither detection nor its geometric robustness proves communication or cell type.

### A second all-RNA reassignment budget

Give every group its full nonnegative integer reassignment weight w_i (zero for
structural regions). Source clipping also minimizes C(c)=sum_i w_i[c_i!=c_i^0]
for each pinned target pattern: every changed group has a forced bit different
from its source, so that group changes in every completion; every other group
retains its source. Thus the same completion minimizes D and C simultaneously.
For a fixed RNA cap K, discard a target pattern only when its least C exceeds K,
then minimize D among remaining patterns with each score. This yields exact
two-budget spectra under the same source-admissibility premise.

Without enumeration, every pair mu,nu>=0 gives the valid upper bound

    max_(D<=delta,C<=K) S <= a+floor(mu delta+nu K-min_c[mu D+nu C-(S-a)]),

and the symmetric lower bound. C is unary, so the existing exact cut computes
the penalized minimum. Tests compare both this bound and exact two-budget spectra
against complete owner enumeration. The two-budget API is in resources.py;
native interval tables use the physical budget alone and are not relabeled as
experiments with an additional RNA cap.

## 5. Composition over nonoverlapping editable regions

Suppose a reference partition has connected fixed cell cores. Editable patches
are pairwise interior-disjoint, have fixed traces on their exterior boundaries,
contain disjoint molecular groups, and every local cell piece attaches to its
appropriate fixed core. Local changes must leave each such piece connected.
Under these conditions, any collection of feasible local alternatives glues into
a global partition: the patches do not overlap, their labels agree with the fixed
exterior, and every piece attaches to a connected core. Counts and deformation
add. These conditions need independent verification; a list of overlapping charts
is not a gluing certificate.

THEOREM 3. Let d_j(v) be a complete local score spectrum on patch j. The global
least deformation for total score v is the min-plus convolution

    d_global(v) = min_(v_1+...+v_p=v) sum_j d_j(v_j).

PROOF. Every global alternative restricts to local alternatives whose costs bound
it below by the right side. Conversely local attaining alternatives for a
minimizing tuple glue by the stated conditions and attain its total cost. QED.

For bounded explicit count range W the ordinary dynamic program costs O(p W^2)
arithmetic operations (often less with sparse spectra). This is a known
convolution algorithm applied to a geometrically justified product. It is
pseudo-polynomial for binary-encoded W. The 42 native overlapping corridors are
NOT combined this way in the reported application.

### Restricted hardness with explicit encoding

Even disjoint two-cell patches can encode 0--1 knapsack if transcript multiplicity
is binary-encoded. For item (a_i,v_i), use a flat source on [-a_i,a_i], slope1/a_i,
band[-2,2], a point at transverse0.9/longitudinal0, and clearance0.1. Moving its
bundle into a common receiving cell forces a height-one tent with area a_i,
attained exactly. Let the bundle contribute integer count v_i. Disjoint patches
can attach to the connected receiving core and separate leaf-cell cores. A budget
B and target count V then give exactly the knapsack decision. All rational input
coordinates and multiplicities have polynomial binary length. Verification is
polynomial, yielding weak NP-hardness for this succinct restricted model.

This does not establish hardness for explicitly listed unit RNA, general planar
continuous segmentation, or three fixed labels on a planar graph. In particular,
three-terminal hardness on arbitrary graphs must not be transferred to planar
multiway cut: planar multiway cut is polynomial for each fixed terminal number.
The exact frontier for general tessellations, star-shaped cells and general 3D
layered surfaces remains open in this project; no unsupported dichotomy is stated.

## 6. A different broad geometry: bounded label transport

Let c0(x) be any 2D/3D reference label field. Suppose a common homeomorphism phi
satisfies ||phi-id||_infinity<=r, and actual ownership is c0(phi^{-1}(x)). Observed
molecule x may also have localization error at most e in that norm. Every possible
owner lies among reference labels intersecting the cube centered at x of radius
r+e. The raster adapter computes a conservative set from every intersecting voxel,
using physical spacing and array-axis coordinates. It does not require ordered
corridors, fixed label count, or a supplied tissue adjacency graph.

Discarding the shared transport yields an OUTER family of independent candidate
assignments. For one cell/gene let n be the source count, a the number of source
RNAs that can leave, and b the number of other RNAs that can enter. At a global
unit-RNA reassignment cap K, exact marginal extrema in this relaxation are

    [n-min(K,a), n+min(K,b)].

PROOF. At most a source molecules and K total changes can lower the count;
selecting any min(K,a) allowed departures attains the bound independently.
The analogous construction selects incoming molecules. QED.

These are sound outer intervals for any coherent family contained in the
candidate model; their endpoint assignments need not define a coherent boundary.
Streaming accumulation costs O(N+E) where E is total candidate-list length.
Candidate extraction, coordinate transformation, image loading and actual
segmentation inference are separate costs. Radius r is NOT the corridor summed-
area budget delta. Neither metric may be substituted for the other.

When a candidate cube leaves a crop, all known labels and an unknown-exterior
label are retained. Only the measured crop molecules are modeled. RNA outside the
crop is not invented, and whole-cell completeness of crop-border cells is not
asserted. The NumPy adapter supports both 2D and 3D, but is not a tested SpatialData
transform adapter or whole-slide benchmark.

## 7. Finite-sample calibration and simultaneous count coverage

Fit the uncertainty family using training units only. Units are independent
regions/donors/slides as warranted by the design, not correlated transcripts.
For a paired calibration unit (X,Y), let Y be an aligned complete reference-owner
vector and define A(X,Y)=D*_X(Y), infinity if Y is not realizable. Disagreement
with fixed exterior owners or structural labels must also yield infinity; it
cannot be omitted merely because the implementation edits only a local patch. This score can
be smaller than the actual boundary's movement; it calibrates OWNERSHIP coverage.

For n exchangeable calibration units and one exchangeable test unit, take rank
ceil((n+1)(1-alpha)) among the calibration scores, padded by infinity. Let q be
that value. If q is finite, use the budgeted ownership family at q; if it is
infinite, return an unrestricted assignment set or abstain from finite geometric
claims. Incompatible calibration units must never be dropped.

THEOREM 4. Conditional on the independently fitted family, the resulting set
contains the test owner vector with probability at least1-alpha. All sound
count intervals propagated from that set simultaneously contain every test
cell/gene count with at least the same probability.

PROOF. Exchangeability makes the rank of the test nonconformity score uniform
under random tie-breaking and conservative without tie-breaking. The padded
order statistic therefore accepts it with probability at least1-alpha. On that
single event, every count is a deterministic function of the same included owner
vector and belongs to its sound extrema. No union bound over cell/gene entries
is needed. QED.

This is a standard split-conformal rank argument plus deterministic propagation,
not a new general calibration theorem. It is marginal over calibration and test
units, not conditional coverage for each donor, each score or every frozen
calibration draw. Paired segmentations calibrate agreement with their reference;
only a justified ownership standard can support claims about biological truth.
The current native budget settings are uncalibrated. The executed holdout uses
known simulated assignments only.

## 8. Statistical and probabilistic boundaries

Fixed-weight count contrasts remain additive and can be bounded directly. Their
sign stability is not a DE significance result. Suppose, separately, true-data
null p-values are valid, hypotheses are fixed in advance, and a calibrated family
covers the entire true count matrix with probability>=1-alpha. If each computed
p_upper bounds its test's p-value over the WHOLE family, Bonferroni cutoff
(beta-alpha)/M controls unconditional FWER<=beta for beta>alpha: outside the
coverage event costs at mostalpha, and within it any false rejection implies a
true null p-value at most that cutoff. The library implements this gate, NOT a
p-value generator, and reports no native DE/FDR experiment.

A proper finite profile model can define P(c) proportional to exp(-lambda D*(c)
plus molecular utility). This is not automatically the owner marginal of a prior
on continuous boundaries: integration over all realizing boundaries generally
contributes owner-dependent volume factors that minimization discards. A continuum
Laplace-type prior additionally needs a reference measure and normalization.
Configuration-level independent Gumbel perturbations yield the usual exact
categorical identity in ideal arithmetic. Independent unary perturbations
usually define a different law. The delivered counterexample measures that
bias rather than presenting randomized mincuts as calibrated posterior samples.

## 9. Parametric flow and all penalized optima

For a leftmost-cell nonnegative count reward, dividing a positive area penalty
shows the form D-t N_left. Count reward enters only source-terminal capacities,
monotonically in t; the remaining graph is fixed. The classical monotone-terminal
conditions therefore permit nested canonical mincuts. The analogous orientation
works for the rightmost cell. Middle-cell membership is a difference of two
threshold bits; opposite-sign terms change both terminal directions, so the same
condition does not follow. A two-RNA example has nonnested owner bits. The code
uses repeated cuts, not an implemented single parametric-flow sweep.

For any solved penalized objective, the residual SCC condensation represents ALL
minimum cuts as forward-closed component sets containing source and excluding
sink. Decoding owner thresholds recovers the corresponding optimal assignments.
Marginal labels computed from reachability are not a Cartesian product. This is
the classical Picard--Queyranne result, used explicitly; it characterizes optimal
penalized assignments, not all hard-budget-feasible labels.

## 10. Verification and contribution boundary

The public API distinguishes exact, bounded, empty and unresolved results.
Persisted cut proofs must be checked against the requested model, score, multiplier
and graph before reuse. The binding checker verifies capacities, integer flow,
cut equality, owner bits and rational physical witnesses. Binding shares the
production compiler; the native upper extrema are separately checked by the
older scanning compiler and full target-owner pinning. These are different code
paths within one project, not independent theorem review.

The implemented contribution is the transfer of the exact shared geometry into
query-specific hard-budget expression certificates, together with conservative
broader-domain relaxations and an auditable calibration interface. Lagrangian
bounds, target enumeration, min-plus composition, conformal ranks and residual
cuts have prior art. Their correct use does not establish historical priority,
whole-tissue anatomical calibration or RECOMB acceptance.
