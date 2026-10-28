# Event-driven exact integration of ownership costs

This note extends the existing cost-preserving reduction. It changes neither the
continuous family nor the molecular objective. Intersection arrangements, adjacent
sorting, and lazy accumulation are standard techniques. The contribution here is
an exact integration algorithm and work bound for the derived activation costs.
It is not a literature-priority finding or a new maximum-flow algorithm.

## 1. The integration problem

After fixed-owner substitution, one boundary and one displacement direction have
M unforced functions

    g_i(t) = alpha_i - s_i dist(t, [a_i,b_i]),

on a common interval, activation labels v_i (which may repeat), and a continuous
piecewise-affine reference r with J knots. The reference includes any mandatory
excursion and therefore can have more knots than the original source. Mandatory
area is an already extracted constant. The remaining cost for active labels X is

    F(X) = integral max(0, max_{i: v_i in X} (g_i(t)-r(t))) dt.

Each g_i has at most three affine pieces. The objective is to construct exact
rational weights W_S >= 0 with

    F(X) = sum_S W_S 1[S intersects X].

No owner vector is enumerated. Both polarities are treated separately by the
existing exact feasibility and cost-preservation theorem.

## 2. Include the reference in the order

Sort the M cones and the reference by decreasing value on an open interval with
no piece breakpoint or intersection. Let the reference occupy position q. For an
edge k between adjacent positions k and k+1, k<q, let S_k be the activation labels
in the first k+1 positions. Then, pointwise,

    max(r, max active g_i) - r
      = sum_{k<q} (h_k-h_{k+1}) 1[S_k intersects X].

This is the usual telescoping identity. Every gap is nonnegative. Every S_k is a
prefix of the cone order; the reference has no activation label. Identical labels
are included once, even if several cones use the same molecular bit.

Between events, the adjacent gap has an affine expression u*t+v. Store its prefix,
coefficients, and start time, but do not integrate it at every intervening event.
When that record changes at time t, add

    (u*(a+t)+2*v)*(t-a)/2

to the prefix weight, where a was its start time. The expression integrates the
whole unchanged gap exactly, including any unrelated events inside [a,t].

## 3. Which events require work?

The event list contains domain endpoints, cone plateau endpoints, reference knots,
cone-cone intersections, and cone-reference intersections. Coincidence intervals
require their piece endpoints, not infinitely many equality events. An endpoint's
right-hand order is determined by value, right slope, and a fixed index tie rule.

Three local changes suffice:

* A cone slope change modifies only its two adjacent gaps.
* A reference slope change modifies only its adjacent gaps; only the gap above
  the reference can contribute a nonzero layer.
* An adjacent swap modifies at most three adjacent gaps and one full prefix.

At a multiple event, update the affected affine pieces, then process adjacent
inversions until the exact right-hand order is restored. Each swap decreases its
inversion count. Enqueue neighboring comparisons after each swap; a single pass
through affected functions is insufficient in general. Each dirty gap is flushed
using its saved old record and then installed with the final new order. Zero-time
intermediate swaps contribute no area.

For a swap at position k, all prefixes except S_k are unchanged. Recompute S_k as
the prefix ending at k-1 union the new activation label at k. A repeated label
may leave S_k unchanged. Prefixes below the reference are tracked for future
updates but receive no weight until they lie above it.

## 4. Exactness theorem

**Theorem.** The event-driven integrator returns the same nonzero rational
coverage coefficients as exhaustive analytic-slab integration, and hence the same
cost for every activation vector.

**Proof.** Between consecutive events all affine pieces and their total order are
fixed. The telescoping identity therefore specifies exactly the active prefix and
affine coefficient on every adjacent gap. Initialization installs these records.
At an event, a gap can change only through a piece change of one endpoint, an
adjacent ordering change, a prefix change, or movement across the reference.
The local update rules include all these positions. Unaffected records retain the
same prefix and affine expression across the event. Flushing a record integrates
that same expression from its creation to replacement; splitting this integral
at omitted unrelated events gives exactly the slab-by-slab sum. The final flush
covers the last interval. Every time belongs to some such record interval, apart
from finitely many zero-measure events. Aggregating by prefix therefore gives
the identical coefficient for every nonzero OR term. QED.

Helper predecessors are recorded whenever a new prefix arises. If its last new
label is v, its predecessor is S without v, itself a prefix at that occurrence.
The construction retains necessary zero-weight predecessors, then removes unused
zero-weight nodes. Predecessor cardinality decreases strictly, preventing cycles.
Each helper represents the intended OR irrespective of which valid predecessor
was recorded first. This preserves exact ownerwise costs with repeated labels.

The old geometric compiler, owner feasibility, mandatory-envelope construction,
flow solver, and continuous witness checker are unchanged. Substituting this
integrator in each direction consequently preserves D*(c) for every feasible c,
not only the optimized values tested experimentally.

## 5. Work bound and its limits

Each cone pair has O(1) affine pieces and order changes. Cones versus a J-knot
reference create O(MJ) crossings and segment endpoints. Thus the event list and
its affected-item records have size

    E = O(M^2 + MJ + J + 1).

Events can be generated by pairwise intersection calculations and indexed source
lookup, then sorted in O(E log(M+J+2)) rational-arithmetic operations. The initial
order and prefixes take O(M log(M+1)) work. Across all events there are O(E)
changes of affine coefficients and adjacent swaps: a swap is an inversion between
the orders immediately before and after an event, hence corresponds to an actual
pair-order change (including the fixed tie convention). Queue updates and gap
flushes are constant in number per piece change or swap. Therefore total work is

    O((M^2 + MJ + J + 1) log(M+J+2))

rational-arithmetic and prefix-bitset operations. Bitsets have width equal to the
number of owner bits; copying, hashing, and union are not unit-cost operations
when that width grows. Rational coefficient sizes and event coordinates have
the same polynomial encoding obligations as in the original reduction.

This replaces scanning and sorting up to M positive layers on every event
interval. It is a bound for the **integration stage**, given the fixed reference.
It is not a bound for the complete pipeline independent of feasibility closure,
mandatory-envelope construction, structural rectangles, flow, or reconstruction.
The implemented event generator still considers every cone pair; a sparse regime
can therefore have substantial pair-generation overhead. All limits remain
explicit rather than inferred from favorable timing ratios.

The prefix-helper bound remains O(M^2+M), independent of J: the reference changes
where the positive prefix ends, not the relative cone order. Reference detail
still enters integration work and coefficient encoding.

## 6. The quadratic positive-coverage expansion is worst-case sharp

For n point sites on one boundary, take zero source, slope one, clearance 1/10,
centers t_i=i, required flipped height a=n+1, domain [-a,n+a+1], and a band
containing [-a-1,a+1]. Every site starts on the right at x=a-1/10. Any binary owner
pattern is feasible: sites are one unit apart, exceeding the conflict clearance
2/10, and the endpoints leave ample space for the tents. The n positive cones
are equal-height tents a-|t-i|.

Every nonempty contiguous interval of indices is a strict top-k set on a
nonzero-length interval near its midpoint: the farthest included center is
strictly nearer than the nearest excluded center. Its corresponding adjacent
height gap is positive there. Thus every one of the n(n+1)/2 contiguous-index
sets has a strictly positive integrated coverage coefficient. A top-k set of
nearest points on a line is contiguous, so these are exactly the nonzero sets.

The OR functions 1[S intersects X] for nonempty S are linearly independent over
all binary X. To see this, if their weighted sum vanishes for all X, evaluation
at the full set gives the sum of all weights zero. Evaluating at the complement
of T then shows the sum of weights indexed by subsets of T is zero for every T.
Induction on |T| (equivalently, subset inversion) makes every weight zero.
Therefore this family needs n(n+1)/2 nonzero terms in this positive-OR basis.

This proves sharpness of the canonical **coverage expansion**. It does not prove
a quadratic lower bound on every possible auxiliary cut graph, on mixed
representations, or on the work needed to optimize a single molecular objective.
The independent tests recover all these sets and their areas for n=1,...,9.

## 7. Evidence and prior-work boundary

The implementation first passed randomized coefficient equality and independent
selected-envelope integration, then physical-model checks and the complete
unchanged native optional workload. Benchmark inputs and source indices are the
same for the scanning and event methods. Native and designed results are reported
separately. The development data provide no external biological validation.

This is a new implementation and derived specialization within the project.
Ordering sweeps, analytic integration and standard graph gadgets are established.
Neither correctness tests nor a measured improvement proves bibliographic priority.
