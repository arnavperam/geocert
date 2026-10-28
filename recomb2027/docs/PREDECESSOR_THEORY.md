# Exact coupled continuous boundaries, uncertain RNA locations, and optimal ownership domains

Derived and implemented 26 October 2026. These are internal mathematical arguments, not external peer review, proof-assistant verification, or established bibliographic priority. Interacting-surface cuts, Lipschitz envelopes, difference constraints, coverage energies, and residual minimum-cut structure have substantial prior art. The new-in-this-project result is the particular continuum-to-ownership representation and its precisely defined uncertainty semantics.

## 1. Model and what is being optimized

There are m ordered graph boundaries f_a(t), a=0,...,m-1, on a common compact interval I. They partition a rectangular corridor into m+1 labeled cells. Fixed exterior rails lie beyond all boundary bands. Every boundary is s_a-Lipschitz with its own positive rational s_a, lies in a specified constant band [l_a,u_a], and has the source's fixed endpoint values. Adjacent differences obey

    0 < d_min[a] <= f_(a+1)(t)-f_a(t) <= d_max[a],

where the upper bound may be absent. The rational piecewise-linear source f0 obeys all these constraints. This ensures nonempty connected strips. It does NOT describe arbitrary tissue tessellations, holes, triple junctions, nuclei in the moving corridor, 3D surfaces, or changes in cell topology. Such generality cannot be obtained merely by relabeling this coordinate chart.

Site i belongs to a group g(i); a group has one shared cell label c_g. A group may contain multiple sites. A location region is an axis-aligned rectangle X_i x I_i with rational endpoints. Positive horizontal clearance rho is a deterministic constraint, not an estimated probability or a Euclidean error radius. The source label must agree with the observed rectangle center; it need not be valid for the whole rectangle. In particular the robust feasible family can be empty even though some possible true placements are valid.

The primary exact mode requires ONE boundary family and ONE owner per group that work for EVERY point in every site rectangle. Define threshold bits

    b_(g,a) = 1[c_g <= a],    b_(g,a) <= b_(g,a+1).

A lower (left-owner) constraint has value ell_i=xhi_i+rho throughout I_i. An upper (right-owner) constraint has value u_i=xlo_i-rho throughout I_i. All boundaries to the right of an assigned cell get lower constraints, all boundaries to its left get upper constraints.

For arbitrary rational unary cell costs c_g(c), and lambda>=0, the primary problem is

    min { lambda D(f) + sum_g c_g(c_g) },
    D(f) = sum_a integral_I |f_a(t)-f0_a(t)| dt.

D is SUMMED INTERFACE SWEPT AREA, also an ordinal cell-label transport cost. It is not in general the tissue area whose cell owner changes counted once. Section 8 establishes that distinction. The result is a penalized optimization theorem, not automatic exactness for hard RNA/area budgets, top-k loss, or sums of failed comparisons.

In a multi-cell region, a pair's total RNA need not remain fixed if RNA moves to a third cell. Therefore the predecessor's two-cell determinant linearization cannot silently be reused. Valid unary objectives include fixed-coefficient regional count contrasts or fixed-denominator quantities; arbitrary re-normalized pairwise ratios require a different solver.

## 2. Difference closure and heterogeneous slope propagation

Put a directed edge a->a+1 of weight d_min[a], and a+1->a of weight -d_max[a] when present. A weight d means f_b(t)>=f_a(t)+d. Let D_ab be the maximum weight of any directed a-to-b path, D_aa=0, and D_ab=-infinity when no path exists. The valid source excludes positive cycles. Floyd-Warshall in the max-plus convention computes D in O(m^3) arithmetic operations.

Define the propagated kernel

    K_ab(r) = max_c { D_ac + D_cb - s_c r },    r>=0,                 (1)

omitting terms with missing paths. This is NOT generally D_ab-s_b r or D_ab-s_a r. A smoother neighboring boundary can constrain an otherwise rougher boundary through thickness constraints.

### Lemma 1: exact path closure

Every admissible family satisfies

    f_b(v) >= f_a(u) + K_ab(|u-v|).                                 (2)

Moreover, (1) is the strongest bound obtained by any finite composition of the declared Lipschitz and same-coordinate difference inequalities between these endpoint variables.

**Proof.** For a particular c, first transfer from a to c at u using D_ac, move from u to v on c at rate s_c, and transfer from c to b at v using D_cb. This proves each term in (1). Conversely take any composite path. Select a surface c of minimum slope among surfaces on which that path moves longitudinally. The vertical steps before and after a visit to c have total weights at most D_ac and D_cb. The summed absolute longitudinal travel is at least |u-v|, and all its rates are at least s_c. Therefore this path's implied bound is at most D_ac+D_cb-s_c|u-v|. A path with no longitudinal travel is covered at r=0. This proves exact closure. QED.

The kernel obeys the composition inequality

    K_ab(r1+r2) >= K_ac(r1)+K_cb(r2).

Indeed concatenate paths attaining the right side with same-direction travel; Lemma 1 bounds the concatenation by the left. It is nonincreasing, and K_ab(0)=D_ab. For an output boundary b, any term with slope s_c>s_b is dominated by the direct b term, since D_ac+D_cb<=D_ab. Similarly rates larger than s_a are dominated. Thus each relevant output cone has slope at most the destination boundary's allowed slope. Kernels with equal rates or everywhere dominated offset/rate pairs can be safely removed; this is an implementation simplification, not a distinct theorem of originality.

### Concrete failure of independent slopes

Let f1=f0+1 exactly on [0,4], both with constant source, s0=2 and s1=1/4. A left-owner obligation f0(1/2)>=3/5 passes the left surface's own endpoint bound of 1. But f1's fixed endpoint and smaller slope imply f0(1/2)<=1/8. The joint assignment is infeasible. Equation (1) gives K_00(r)=-(1/4)r, correctly transmitting the stricter bound. The construction is executed in tests/test_kernels.py.

## 3. Exact ownership feasibility

A selected lower obligation (surface a, interval I_i, level ell_i) and selected upper obligation (surface b, interval I_j, level u_j) conflict exactly when

    ell_i + K_ab(dist(I_i,I_j)) > u_j.                             (3)

If no path exists between a and b there is no such derived inequality. Since each K is the maximum of at most m affine terms, checking (3) is finite rational arithmetic.

Unary rules compare each owner-side obligation to all propagated fixed bands and endpoints. For example a lower interval on a is impossible if, for any b and kernel term (d,s), ell_i+d exceeds b's upper band, or exceeds f0_b(e)+s dist(I_i,{e}) at either fixed endpoint e. Upper obligations are checked symmetrically. A conflicting owner fixes its bit to the opposite value. Conflicting fixes mean an empty family, not a biological robustness certificate.

For every conflict (3), add the implication

    b_(g(i),a) => b_(g(j),b).

Add threshold nesting and any user pins. Same-bit conflicts are vacuous because one bit cannot select both polarities simultaneously.

### Theorem 1: closure-system characterization

The threshold bits correspond to a realizable universal-region ownership pattern if and only if these unary rules, pairwise implications, and nesting constraints hold.

**Necessity.** Use (2) on every pair of interval positions. The minimum separation of two closed intervals maximizes the nonincreasing propagation bound, yielding (3). Bands, endpoints and nesting are necessary by definition.

**Sufficiency.** Construct, for each boundary b, lower and upper envelopes from every selected lower/upper obligation, every fixed endpoint, and the fixed bands:

    L_b(t) = max{ ell_i + K_ab(dist(t,I_i)), propagated fixed data },
    U_b(t) = min{ u_j - K_ba(dist(t,I_j)), propagated fixed data }.

For any lower origin a and upper origin c, the sum of their path kernels through b is at most K_ac(dist(I_i,I_j)), by composition, the triangle inequality and monotonicity. Every lower component is therefore below every upper component by the pairwise/unary conditions. Hence L_b<=U_b. The propagated kernels ensure L and U individually satisfy all difference constraints. Their output slopes are at most s_b, and endpoint intervals collapse to the prescribed values. Their band bounds follow by propagation. Thus L and U are feasible boundary families. This proves existence. QED.

Absolute owner threshold vectors form an implication closure system. This is not a downward-closed assertion about source-relative flips.

## 4. Exact minimum deformation and the finite cut

### Theorem 2: componentwise clamp is globally feasible and pointwise nearest

For any feasible ownership pattern,

    f*_b(t)=min(U_b(t),max(L_b(t),f0_b(t)))                           (4)

is feasible and minimizes D among all its continuous realizations.

**Proof.** Families satisfying all lower difference inequalities are closed under componentwise minimum and maximum. For example, if f_b>=f_a+d and g_b>=g_a+d, then min(f_b,g_b)>=min(f_a,g_a)+d, and likewise for max. All other constraints, including each slope bound, are also preserved by real-valued min/max. L,U and f0 obey the geometric constraints, so (4) does too. It lies between L and U, giving ownership feasibility. Any realizing boundary lies between these envelopes, and (4) is the nearest point to f0 at each t and b. Summing/integrating the pointwise inequality proves minimum D. QED.

**Source-region caveat.** f0 is geometrically feasible but may not satisfy its whole uncertain boxes. Area costs must therefore include all selected ownership constraints, not only those called flips by the nominal source. The implementation does this; omitting non-flip cones would be incorrect.

For a lower obligation and kernel term, the source-relative excursion is

    [ell_i+d-s dist(t,I_i)-f0_b(t)]_+,

activated by its threshold bit. Upper obligations give the analogous nonnegative downward excursions activated by the bit's complement. Fixed geometric data are already met by f0 and cannot push it outward. For feasible ownership patterns, upward and downward excursions on the same boundary do not overlap positively; their maxima integrate to the exact clamped deformation. Sum across boundaries.

On each rational interval where positive excursion functions have a fixed order h1>=...>=hk>0, with h_(k+1)=0, their active maximum is

    sum_j (h_j-h_(j+1)) OR{active variables among the j highest functions}.

Integrate each affine difference exactly. All resulting rational weights are nonnegative. A positive OR layer is represented by a helper with bit=>helper hard arcs and a helper-to-sink cost. A complement OR uses the reversed helper implications and source-to-helper cost. These are classical cut representations. Duplicate cones carrying one variable remain one activation; equal variable sets can be merged. Nested predecessor helpers reduce the arcs per distinct layer to a constant.

A cell cost is encoded exactly by

    c_g(m) + sum_(a=0)^(m-1) (c_g(a)-c_g(a+1)) b_(g,a).

Negative unary coefficients are shifted by an explicit constant; hard feasibility arcs are added. Multiplying all rational capacities by their denominator LCM produces an integer network. A hard capacity greater than the sum of every finite capacity prevents invalid cuts whenever any feasible owner pattern exists. A minimum cut using hard capacity identifies model infeasibility.

### Theorem 3: exact finite optimization

For rational data, arbitrary unary cell costs, and lambda>=0, one finite integer minimum cut solves the joint universal-region ownership/continuous-boundary problem exactly and reconstructs (4). No raster sampling or molecular-owner enumeration is required.

The lower/upper construction, analytic layer weights and unary telescoping prove equality for every feasible owner pattern, hence also after optimization. A maximum-flow/minimum-cut equality certifies the resulting finite-network optimum. It does not independently certify that arbitrary code compiled the physical problem correctly: continuum equivalence remains a theorem/compiler obligation, tested separately below.

### Polynomial size, not a practical scaling promise

Each surface-pair kernel has at most m surviving affine terms. For N sites and m boundaries there are at most O(N m^3) propagated ownership cones across all destinations. On destination b let M_b be its cone count and J_b its source knots. Flat-top cone pair intersections have constant complexity; cone/source zeros have O(J_b) complexity. A full arrangement has O(M_b^2+M_b J_b+J_b+1) slabs, with at most M_b layers per slab. The readable implementation explicitly constructs this polynomial arrangement. All arithmetic is rational; denominator/numerator bit lengths and integer-capacity scaling have polynomial bit costs. Maximum flow has a separate polynomial complexity.

This reference compiler is more expensive than the predecessor's single-point sparse routine. A box's flat-top source-relative support can be disconnected when the source oscillates, so that earlier single-support proof cannot be reused silently. No sparse bound depending only on nearby point count is claimed here.

## 5. Ownership across every optimum

Given an optimal residual network, a source-side minimum cut is exactly a residual-forward-closed set containing source and excluding sink. This is the classical all-minimum-cut characterization (Picard–Queyranne), not a new principle.

For a proposed group owner c, require its threshold-one vertices and the source to be inside, and its threshold-zero vertices and the sink outside. The label occurs in some optimum iff residual reachability from all required-inside vertices avoids every required-outside vertex. Otherwise no such minimum cut exists. Repeat for each group/cell after one maximum flow.

A singleton marginal domain establishes the same owner in every optimum of that penalized model. Multiple allowed marginal labels do NOT imply every combination of those labels is feasible or optimal. Nor does an optimality domain describe every boundary under a hard uncertainty budget. The tests enumerate all minimizers independently on small instances.

## 6. Uncertainty quantifiers and safe bounds

Let z be any fixed placement of the site points inside the stated rectangles. Define F_z as its feasible geometric/owner set. Define F_forall using the simultaneous all-position constraints above. Always F_forall is a subset of F_z.

There is also an implemented OUTER relaxation for possible placements. At rectangle midpoint t0, longitudinal half-width eta, any point assigned to the left of boundary a must satisfy

    f_a(t0) >= xlo+rho-s_a eta,

and any point assigned to the right must satisfy

    f_a(t0) <= xhi-rho+s_a eta.

This follows by applying s_a-Lipschitz continuity from an actual position to the midpoint. These point constraints are necessary for every possible placement, so F_z is contained in F_outer. They are not sufficient for one shared location realizing all boundary constraints. An outer solution must NOT be labeled an existential-location witness.

For the same objective and any fixed z,

    min(F_outer) <= min(F_z) <= min(F_forall),                       (5)

with infinity when appropriate. The nominal-center optimum is also bracketed. For a maximization the directions reverse. This is a deterministic inclusion guarantee, not probabilistic calibration. An empty universal family can coexist with a feasible true placement and does not establish that an RNA cannot belong to a cell.

## 7. Independence failures and effective regularity examples

Two source boundaries 0 and 1 over [0,10] have slope1 and minimum gap3/5. At t=5, RNA sites x=2/5 and3/5 are both originally in the middle cell, with rho=1/20. Their costs prefer the left and right exterior cells respectively. Independent binary optimization can choose both transfers, leaving gap1/10, which violates the joint constraint. The exact joint optimum transfers one site and maintains the gap. Its objective is -159/200 and D=41/200; independently optimizing gives the infeasible objective -319/200.

There are two equivalent joint optima. Their marginal owner domains are {0,1} and {1,2}; the joint combination (0,2) is infeasible. This supplies an explicit test of why individual alternative owners do not form an admissible Cartesian product.

Widening the locations by horizontal radius3/100 and longitudinal radius1/10 gives objectives -561/625 (possible-location outer relaxation), -159/200 (nominal), and -407/625 (universal). These are not three estimates of one posterior probability; they answer different quantified problems.

## 8. Deformation cost semantics

Let c_f(x,t) be the ordered cell index. Since both old and new boundaries are ordered, at a fixed (x,t) their threshold vectors are monotone and disagree in exactly |c_f-c_f0| positions. Integrating those binary differences gives

    D(f) = integral |c_f(x,t)-c_f0(x,t)| dx dt.

If A_change counts a point once when its owner changes, then

    A_change <= D <= m A_change.

This distinction is measurable, not cosmetic. The supplied three-boundary example has D=63/8 but A_change=97/16. A point crossing several interfaces contributes several units to D. The finite-cut theorem optimizes D, not A_change; the latter is only recomputed by a separate geometric checker. If swept regions do not overlap (for example appropriately disjoint bands), equality holds.

## 9. Novelty and claim boundary

Multiple interacting surface cuts were developed long before this project (Li et al.,2006), and irregular/subvoxel multiple-surface methods also exist (Shah et al.,2019). Heterogeneous smoothness, minimum/maximum distances, graph cuts, Lipschitz extensions, and residual minimum-cut domains cannot be advertised as inventions. The candidate technical distinction is an observation-derived, exact finite representation of the stated continuous multi-cell region-ownership objective, including the heterogeneous propagation closure and explicit universal/possible-location bounds. Search has not established historical priority. This result could still be judged a specialization of established mathematics.

The executable evidence is synthetic/designed and predecessor regression, not end-to-end native tissue validation. No new donor, calibrated location confidence region, true anatomy, clinical decision, public code release, or conference submission is claimed. The restricted corridor geometry must remain explicit in any title or abstract.
