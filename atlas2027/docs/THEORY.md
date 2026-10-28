# Exact expression bounds on a glued boundary atlas

## 1. The computational object

Fix a reference labelled partition, measured point transcripts with unit gene contributions, immutable structural regions, and a geometrically specified family. For an assignment vector c, D*(c) is the minimum physical deformation of any allowed partition realizing c, or infinity if none exists. A count interval is an extremum over D*(c) <= delta. It describes the observed count under that one family; it does not estimate the true membrane, missing transcripts, or posterior probability.

The local continuous-to-finite ownerwise reduction remains the engine. For ordered boundary graphs with heterogeneous constant slope bounds, constant bands and thickness inequalities, fixed endpoints, and interval ownership constraints, propagate lower and upper requirements using

K_ab(r) = max_c {Delta_ac + Delta_cb - s_c r}.

Here Delta is maximum-path closure of same-coordinate difference constraints. Pairwise lower/upper conflicts and unary endpoint/band conflicts exactly characterize feasible owners. Their propagated envelopes L,U satisfy the original difference inequalities and individual slope bounds. The reference f0 does too. Componentwise min/max closure makes min(U,max(L,f0)) jointly feasible and pointwise nearest. Its integrated displacement is D*(c), for every feasible c. Integrating ordered positive excursion layers gives a weighted OR expansion and a rational cut representation before molecular unary costs are selected.

These ingredients have established foundations: Lipschitz extensions, difference-potential lattices, graph-representable higher-order energies and surface cuts. The specialized identity preserves the complete minimum realization cost, not merely continuous feasibility or one optimized query.

## 2. A globally compatible family on a fixed tessellation

Let P_i be finitely many closed ribbon patches with disjoint interiors. Each has an injective chart Psi_i:[a_i,b_i] x [-B_i,B_i] -> P_i. The reference interface is v=0, its two halves carry exactly the stated two labels, and both rails and endpoint traces are fixed. Every proposed interface v=h_i(t) is continuous, h_i(a_i)=h_i(b_i)=0 and |h_i|<=q B_i for a common q<1. All observations belong to at most one patch; unassigned observations and unselected interfaces retain their source labels.

**Gluing theorem.** Every tuple of individually admissible patch interfaces is induced by one global homeomorphism of the reference partition which is identity outside the patch interiors. Consequently reference connectedness, holes, component counts, and fixed junction incidences are preserved; no independently solved interfaces compete for the same physical region.

**Constructive proof.** At fixed t define T_i(t,v)=(t,H_i(t,v)) by

H_i(t,v) = -B_i + (v+B_i)(h_i(t)+B_i)/B_i, for v<=0;
H_i(t,v) = h_i(t) + v(B_i-h_i(t))/B_i, for v>=0.

Both slopes in v are strictly positive and bounded away from zero. The branches meet at h_i(t), fix +/-B_i, and have continuous explicit inverses. At t=a_i,b_i the map is identity. Therefore Psi_i T_i Psi_i^{-1} is a homeomorphism of P_i fixing its entire boundary. Extend by identity outside. Patch intersections lie only on fixed boundaries, so the finite maps and their inverses glue continuously. Their global composition is a homeomorphism. Ownership constraints additionally preserve each nuclear region's label; material nuclear coordinates are not inferred.

The evaluated atlas uses the exact source raster interface, without replacing it by row-center interpolation. A monotone edge path has t=sx*x+sy*y and u=sx*x-sy*y in pixel coordinates, sx,sy in {-1,1}. Write its exact graph u=r_i(t). With pixel size p, the chart is

Psi_i(t,h) = p/2 * (sx*(t+r_i(t)+h), sy*(t-r_i(t)-h)),

plus the recorded coordinate offset. Its absolute Jacobian is p^2/2 on each affine segment, and it remains an injective piecewise-affine chart across source corners. Thus physical changed-label area is exactly

D_i(h_i) = (p^2/2) integral |h_i(t)| dt.

Since patch interiors are disjoint, physical area counted once is D=sum_i D_i. This differs from summing overlapping corridor areas. Curved/sheared source arcs have exact physical conversion here; a general chart needs its own Jacobian, and ordinary area integration cannot silently ignore it.

The atlas is restrictive: unselected reference boundaries, all patch seams, junction neighborhoods, crop-border cells and prescribed nuclear labels are fixed. The guarantee concerns complete retained counts for crop-contained reference cells under this family, not freely varying whole perimeters. Selection and coverage are part of every exported atlas identity.

## 3. Count-specific exact critical-height chain

Assume one boundary on a compact rational interval, a rational piecewise-affine reference, positive constant slope bound, a constant band, fixed endpoints, any fixed structural interval constraints, and one point per free RNA group. Point locations may coincide, but each RNA is an explicitly represented unit observation. Free interval-location groups require additional memory and are outside this algorithm.

Collect all affine pieces of all lower and upper propagated observation cones, fixed geometric envelopes, and the source, before choosing owners. Include piece endpoints, source and observation knots, and every pairwise crossing in one common mesh T. Candidate heights V_j are function values at t_j, within forced envelopes. Remove candidates violating fixed constraints or any point's clearance gap. Feasible canonical realizations only use the remaining candidates. Crossings outside the forced envelopes cannot be a switch in a canonical realization; omitting them is safe.

**Count-spectrum theorem.** Let q be a nonnegative integer score obtained by summing per-point owner scores, and W its variable score range. For every q, the exact least continuous area at score q is obtained by a height/count shortest-path dynamic program on this common mesh. If H_j=|V_j|, its transition work is O((W+1) sum_j H_j), beyond polynomial critical-mesh construction. With unit molecule counts W is at most the number of relevant input molecules; the algorithm is polynomial in the explicit input size and rational encoding. With binary-compressed multiplicities it is only pseudo-polynomial in W.

**Proof of geometric exactness.** For any feasible complete assignment, the canonical minimizing boundary is made of candidate pieces. Between consecutive events its selected piece is affine, and displacement from the reference has constant sign. Its trapezoidal nodal cost equals its continuous integral, and its owners/count are retained. Conversely any candidate-height path satisfying adjacent slope inequalities and nodal fixed-region/ownership constraints interpolates to a continuous feasible boundary: the mesh contains every constraint endpoint, and every intervening constraint is linear. The trapezoidal absolute-displacement cost is at least that interpolant's exact area, which is at least the minimum for its owners. The retained canonical minimizer attains equality, including after restricting to a count. Therefore minimum nodal cost at each count equals the continuous minimum.

**DP recurrence.** Assign every RNA exactly once, at its point's mesh node. A permitted height determines its unique owner because clearance excludes the intervening gap. Let z_j(h) be the total local score, w_j half the sum of adjacent mesh widths, and a_j(h)=w_j|h-f0(t_j)|. With F_j(h,q) the least prefix cost ending at h with score q,

F_j(h,q)=a_j(h)+min_{u in V_(j-1): |u-h|<=s(t_j-t_(j-1))} F_(j-1)(u,q-z_j(h)).

Initialization includes fixed groups and the first node. The usual shortest-path induction proves exactness. For each predecessor count, the valid sorted predecessor heights form a sliding interval as h increases. A monotone queue gives all their minima in O(H_(j-1)+H_j) operations; no pairwise height-arc expansion is necessary. Sum over counts and nodes. Backpointers give attaining boundaries, independently checked in continuous coordinates.

If P is the total number of explicit affine candidate pieces, the naive mesh has O(P^2+P) events and at most O(P) heights per event. Preparation is polynomial, and the DP is conservatively O((W+1)P^3) operations. Memory for all backpointers is O((W+1)sum H_j); value-only rolling memory can be smaller. Rational operand bit lengths and witness serialization must be accounted for separately. Resource limits return an incomplete outcome, never an exact spectrum.

The resource dimension and shortest-path recurrence are standard. The contribution here is their exact continuous realization via a complete critical mesh and physical atlas, rather than treating an arbitrary sampled mesh as exact. The algorithm does not claim to solve all multi-boundary hard-budget problems.

## 4. One cell interval under one budget

For cell k and gene g, let n0 be its retained source count. For each patch i incident to that cell, compute S_i(d), the least physical area attaining count change d in that patch, while allowing all other genes to move consistently. Patches with no relevant transcript have zero-cost source states. With min-plus convolution,

S_kg = S_1 star ... star S_r,    (A star B)(q)=min_u (A(u)+B(q-u)).

The exact interval is min/max {n0+d : S_kg(d)<=delta}. Local source counts are subtracted only to form changes, so the fixed contribution is counted once. Every endpoint trace chooses one local attaining state per patch; the gluing theorem yields a valid global partition of that total area. Conversely every allowed global partition decomposes into those local states with at least the sum of their minima. This proves both inequalities and exactness.

Convolution is classical. The substantive conditions are unique transcript ownership by patches, compatible fixed traces, disjoint physical support, a cost-preserving chart, and one shared global budget. The input atlas verifier checks them. Shared cell identity alone is insufficient. Marginal extrema for different cells/genes need not occur together, although every feasible full count matrix is contained componentwise in the interval table.

For bounded unit counts the total state range is at most the number of relevant observed molecules. Finite min-plus convolution is polynomial; a simple implementation takes O(r W^2). Together with the scalar count-chain theorem this yields exact per-cell/gene hard-budget spectra for arbitrary fixed reference adjacency equipped with such an atlas. This does not solve arbitrary topology-changing deformations or independently movable junctions.

## 5. Penalty bounds, explicit unit-count gap, and target projection

For any multiplier mu>=0, m(mu)=min_c {mu D*(c)-N(c)} gives U(delta)<=floor(mu delta-m(mu)). A feasible boundary gives the opposite bound. Applying the same argument to -N gives a lower endpoint. This is weak duality, not a new duality theorem, and equality need not hold.

A gap already occurs with eight distinct unit-RNA points. Use source zero on [0,30], slope1, band[-4,4], clearance rho=1/10. At centers5,15,25, use respectively2,3,3 points. Point j in a cluster has t=center+j/100 and x=height-j/100-rho, where heights are2,3,3. Within a cluster, clearance and slope force common owners; the leading cone dominates the other cones. Clusters' positive supports are disjoint. The exact count/area spectrum is {0:0,2:4,3:9,5:13,6:18,8:22}. All coordinates are distinct; no binary-encoded multiplicities are used. At delta10 the maximum count is3, while the best penalty bound is4. For mu<=1/3, the count8 state gives B(mu)>=8-12mu>=4; for mu>=1/3 the count2 state gives B(mu)>=2+6mu>=4. At1/3 equality holds. The count chain recovers the unsupported count3 state.

For sparse targets, a separate exact fast path enumerates only target owners and uses source-clipped implication closure to optimally complete collateral owners. This requires the source to satisfy every source-owner and structural obligation at zero cost. The count-chain theorem does not require that stronger premise. In the measured hybrid workflow, single-copy targets use the closure fast path and all multi-copy targets use the count chain. Complete enumerated spectra and an independent scalar implementation are retained as verification baselines, not as a hidden dependency of the polynomial algorithm.

## 6. Comparator semantics

The marginal-owner outer envelope lets every individually affordable transcript move at once. It is sound but can overspend the shared budget or combine incompatible moves. The feasibility-matched additive-cost approximation retains all allowed complete assignments but charges the sum of individual minimum transfer areas, including collateral transfers. Under source admissibility this exceeds or equals the true shared area. Its feasible-budget set is a subset, so its count interval is an inner approximation and can omit valid endpoints. These two baselines answer different approximation questions and must not be described interchangeably.

A deterministic five-micrometre nearest-nuclear-pixel segmentation is an additional image-only reference. Its disagreement with a count interval tests compatibility with the atlas, not interval coverage of biological truth. Countwise inclusion is weaker than feasibility of its entire assignment vector. No Proseg, Baysor or ResolVI run is represented by this baseline.

## 7. Evidence and limits

The accompanying report gives all fixed-rule development and external-format-test denominators, exact witness and source checks, count-chain checks, and the independent integration/Cartesian product validation. These test correctness under declared geometry. No theorem implies that actual membrane errors stay inside that geometry. The external vendor test data are explicitly unsuitable for drawing biological conclusions. The restricted exact reduction, budget algorithm and topological composition should be assessed as methodological results, not as a clinical or corrected-cell claim.
