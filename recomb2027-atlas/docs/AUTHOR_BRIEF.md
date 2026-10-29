# Author briefing

**Paper:** Certified Expression Intervals under Shared Segmentation Uncertainty.

## What the paper claims

For a declared geometric family, the method preserves the least continuous area needed for every feasible transcript-assignment vector. A finite critical-height/count dynamic program solves all hard-budget count levels for scalar point-RNA patches. Disjoint patches with fixed traces glue into a single topology-preserving partition, and their spectra compose under one area budget. The output is exact for that family and retained observations; it is not a posterior, a segmentation correction, or evidence that actual errors obey the family.

The main text proves the scalar case used in the atlas. The supplement retains the full multiple-boundary extension, including heterogeneous slopes and thickness propagation. The single-interface application does not demonstrate a native benefit from thickness coupling.

## The gluing theorem — worked example

Take one chart with longitudinal coordinate t in [0,2] and transverse rails at -2 and 2. Let h(t)=min(t,2-t)/2. It is zero at both ends and reaches 1/2 at t=1.

At t=1, the lower-half map has slope (2+1/2)/2=5/4 and the upper-half map has slope (2-1/2)/2=3/4. Both are positive. The source interface v=0 maps to v=1/2; the rails remain -2 and 2. For example, v=-1 maps to -3/4, and v=1 maps to 5/4. The two branches agree at the interface and have continuous inverses. At both longitudinal endpoints the entire transverse map is identity.

Thus a chart-conjugated map fixes every point of the ribbon boundary. Two such maps on disjoint interiors agree wherever their supports touch; extending by identity outside gives one homeomorphism. Disjointness and fixed traces are essential. Shared cell identifiers alone do not prove gluing.

The chart-area change is the triangular integral 1/2. With pixel scale p=1 and the raster chart's Jacobian p²/2, physical changed-label area is 1/4. The exact Jacobian is necessary even when the reference staircase turns. Reference connectivity, holes and junction incidences are preserved, but fixed seams limit which errors can be represented.

## Lemma 5 — why the critical mesh is exact

Use reference f0=0 on [0,2], slope bound 1, and a transcript at (t,x)=(1,1/5), with clearance 1/20. Reassigning it requires f(1)>=1/4. The nearest feasible boundary is max(0,1/4-|t-1|), with knots 3/4, 1, 5/4 and area 1/16.

A mesh containing only 0,1,2 would charge trapezoidal area 1/4 to that changed assignment. The critical mesh includes all cone/source intersections, so it includes 3/4 and 5/4 and retains the canonical minimizer exactly.

The proof needs two directions. Every feasible nodal interpolation is a feasible continuous boundary, and its trapezoidal absolute-displacement cost upper-bounds its continuous area. Conversely, for every assignment the canonical minimum is available on the critical mesh and has constant-sign affine displacement on each interval, so its trapezoidal cost is exact. Therefore the nodal minimum, including after a count restriction, equals the continuous minimum. The integral of every possible interpolant is not claimed exact.

## Theorem 6 — count state and a worked hard budget

Put two point transcripts at (t,x)=(1,1/5) and (5/4,1/5), on a zero source over [0,3], slope 1 and clearance 1/20. Both begin on the right. Each can move alone at area 1/16. Both can move together at area 7/64, because their tents overlap.

The exact left-count spectrum is {0:0, 1:1/16, 2:7/64}. At budget 1/10, count 2 is impossible although both singletons are individually affordable. At budget 11/100, count 2 is possible although the additive singleton sum 1/8 exceeds that budget. These are the outer-bound and additive-inner failure modes.

In the dynamic program, F_j(h,q) stores the least prefix area ending at candidate height h with accumulated count q. A point is counted once at its longitudinal mesh node. A predecessor height u is allowed when |u-h| is within the slope times the interval width. For each count, sorted feasible predecessors form a sliding interval; a monotone deque gives all transition minima in linear work in the two candidate-list sizes. Backpointers recover an attaining boundary.

The transition bound O((W+1) sum H_j) excludes critical-mesh construction and rational operand costs. Unit transcripts bound W by explicit input size. Compressed multiplicities instead give pseudo-polynomial dependence. The algorithm is scalar and point-based; arbitrary coupled multi-boundary counts need additional state.

## Ownerwise identity and source admissibility

The source and propagated envelopes satisfy a min/max-closed geometric family. Their clamp is therefore jointly feasible and nearest to the source at every coordinate. Integration preserves the entire minimum-deformation function before molecular unary terms are chosen. Standard weighted-coverage helpers encode that function; those gadgets and lattice closure are prior work.

Logical hard arcs are part of the representation. Their numerical capacity must be set above the sum of finite capacities after the particular unary objective and multiplier are known.

The source-clipping shortcut needs the source to realize every source-owner and structural requirement at zero cost. Merely satisfying slopes, bands and endpoints is insufficient. The general representation and count-chain argument have weaker source assumptions.

## What the experiments mean

- Original eight-budget primary analysis: exact geometry tightens 4/46 eligible lung and 203/1,026 eligible breast cell–gene questions. Eligibility means marginal width at least two at some tested budget. The complete query and row denominators remain available.
- The breast inputs are 10x format-test data, not a biological validation set.
- Constant numbers of non-singleton intervals do not imply fixed endpoints. The original breast table has 61 further endpoint changes after 0.25 and three after 0.5 square micrometres.
- Width sensitivity uses separately extracted atlases. They need not be nested: wider supports reject more arcs. The analysis measures the reach/coverage trade-off of this extraction rule, rather than estimating plausible membrane error.
- The Proseg comparison uses one recorded run per clipped lung crop, with factorization disabled after the initial failure. Transcript joins and count matrices are checked independently. These runs provide a reference reassignment, rather than an ownership-accuracy test.
- Entrywise count reach does not imply the complete alternative assignment vector or its count matrix is simultaneously realizable.

## Four questions to answer aloud

1. Why is the clamp feasible? Because the propagated envelopes and source obey the same min/max-closed geometry; selected ownership bounds are also preserved.
2. Why do singleton costs lose information? Several changes can share area, or be mutually incompatible, so neither a sum nor a Cartesian product describes the joint family.
3. Why can source detail change weights without enlarging the prefix bound? A common source changes positivity but cancels from relative cone order.
4. Which parts are new claims versus known machinery? The specialized cost-preserving continuum reduction, exact critical-mesh count realization, and compatible physical composition are the results to assess. Cuts, coverage gadgets, shortest-path resource states and min-plus convolution are established foundations.

AI contributed to mathematical development, algorithms, implementation, checking and writing. The named author must review, understand and take responsibility for the proofs, data interpretation and submission. A successful replay is not an independent human proof review or an acceptance guarantee.

## Finite-family extension

The five atlases share one reference partition, observed molecules and physical-area cost. A permitted partition belongs to one of those complete families. For each count q, its minimum union cost is therefore the minimum of the five member costs at q. An endpoint inherits the selected member's complete partition witness and atlas identity.

The primary family is a member of the union, so every primary alternative remains available. Wider intervals are expected: at area 4 µm² on the original question census, invariant counts change from 84.1% to 76.4% in lung and from 77.1% to 65.9% in breast. Additional source-zero questions expand the full union census to 29,212; report that denominator separately.

The comparison checks whether the exact reference count occurs in the spectrum, not merely between its endpoints. A spectrum can have missing interior integers. On fixed comparison keys, union reach is 315/2,798 for lung Proseg, 299/1,114 for lung nuclear expansion and 2,942/10,965 for breast nuclear expansion.

Be able to answer these questions:

1. Which atlas supplies this endpoint? The witness records one complete member and the compatible local traces within it.
2. Can two endpoints from different members be combined? Each is a separate alternative; a joint question requires a joint predicate or vector spectrum within each member before taking the finite union.
3. Does the union allow the widest movement over all primary edges? Each member retains its own coverage and displacement range. There is no single atlas with both maximum coverage and maximum displacement.
4. What certifies the component optima? They are the original exact, checksum-bound spectra. The new verifier independently reconstructs the union and checks its selected continuous witnesses.

## Larger contiguous computation

The central 400 µm square was fixed from the original raster extent before transcript or count outcomes. It has 1,363 complete cells and 226,983 retained molecules; the unchanged primary-atlas pipeline returns 69,278 spectra and 554,224 exact rows. Atlas extraction, rational geometry validation and counts/witness output took 22.55 minutes in one measured run, with a largest-stage RSS of 2.43 GiB.

This result establishes a completed larger contiguous workload. The region overlaps earlier crops and remains the same vendor format-test specimen. Its cell count is not an independent cohort total. Per-ribbon independence alone does not establish linear whole-slide runtime: extraction ordering, overlap checks, geometric validation, local count work and output size also matter.

The saved-output audit independently reconstructs the query census, min-plus compositions, witness counts/costs and every interval. Continuous witness checks and local optimality come from the exact solver; sixteen hash-selected multicopy queries additionally match independent scalar enumeration. The source audit separately checks every retained record against the original vendor files.

The concentric 600 µm square has 3,171 complete cells and 470,293 retained molecules; exact enumeration returns 152,925 spectra and 1,223,400 interval rows. It uses the existing enumeration backend, selected after geometry validation established at most seven copies per local gene. The first attempt reached its global deadline during counts. A separately recorded count-only retry completed; cumulative measured attempt time is 36.89 minutes. Completed stages alone total 29.72 minutes, which excludes the interrupted work. These are separate measurements with different backends, rather than a runtime-scaling ratio.

Sixteen additional 600 µm queries match independent scalar enumeration. Their identities were generated using the pre-existing input-only hash rule after production completed but before inspecting its spectra; the 400 µm sample was generated before production. Both complete output censuses and all interval rows pass the separate audits.
