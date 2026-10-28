# Aggregated prefix families and hard-implication contraction

27 October 2026 (Phoenix). This note develops the refinement suggested in the two supplied reviews. The companion edit-pack ZIP mentioned by the second review was not attached; the proof and tests below were written from its stated argument. Neither ordinary SCC contraction nor adjacent-swap counting is claimed as a new general technique.

## Prefix-family graph-size proposition
For one boundary and one polarity let g_i(t)=alpha_i-s_i dist(t,I_i), i=1,...,M, and let r be the common piecewise-affine reference. The layer functions are [g_i-r]_+ (the other polarity changes the common reference's sign). The reference can already include mandatory excursions after pin substitution.

**Proposition.** O(M^2+M) distinct activation-prefix helper vertices and helper arcs suffice after aggregation, independently of the number of reference knots.

**Proof.** Each g_i has at most three affine pieces. A pair therefore has only a constant number of relative-order changes, including the start/end of coincidence intervals. Fix a global index tie rule on open intervals. The entire sequence of orderings has O(M^2) inversions across event transitions. Replace any simultaneous permutation by adjacent swaps for counting. One swap changes only the prefix between the swapped positions; there are at most M initial prefixes. Hence O(M^2+M) cone-index prefixes occur. A source crossing changes the positive cutoff but not the relative order. Every positive layer is one of these full prefixes. Replacing cone indices by repeated activation labels can only identify sets, not increase their number.

For any newly appearing label set S, take its first appearance in one ordering, ending at the last new label v. The previous distinct-label prefix P=S\{v} is in the same family. Include P even when its integrated weight is zero. A helper for S depends on P's helper and v; its cardinality decreases along predecessor links, so no cycle arises. Two hard arcs and one weighted terminal arc suffice; complement polarity reverses the implications. Sum exactly the weights of every repeated S. This yields the claimed helper bound. QED.

Counting RNA and structural obligation rectangles after exact union merging as N, molecular/structural groups as G and boundaries as m, with M_b=O(Nm^2), the total graph can have

    V = O(Gm + sum_b(M_b^2+M_b))
    E = O(Gm + N^2 m^2 + sum_b(M_b^2+M_b)).

Thus the auxiliary-size contribution is O(N^2 m^5), rather than charging O(M_b) raw layers on every analytic slab. This is an aggregated graph-size bound. The reference compiler still scans slabs/layers; its runtime does NOT automatically improve to this bound. Source knots affect integration, capacity bit lengths, input size and reconstruction. The argument makes no lower-bound claim against every compressed height representation.

**Executed support.** 168 polarity-arrangements and 1,260 exact assignment-integral comparisons passed. Tests include simultaneous events, identical cones, repeated activation labels, mandatory-style common references and up to 41 reference knots. The truth routine integrates selected maxima separately; production integration supplies only the compared layers. Finite tests support implementation, not the general proof or literature priority.

## Hard-SCC cost preservation
In a directed hard-implication graph every hard-feasible binary assignment is constant on each strongly connected component. Contract those components, discard internal arcs, and sum parallel capacities. Quotient cuts and original hard-feasible cuts correspond by expansion, with equal finite cost. Fixing any feasible owner vector preserves this correspondence, so minimizing auxiliaries preserves its ownerwise minimum too. The partition depends on hard constraints, not the unary objective; it can be reused while those constraints remain unchanged. Late pins require revalidation/recompilation. A hard component containing both source and sink denotes infeasibility, not a valid merged-terminal instance.

The delivered implementation preserves a conventional graph certificate in the terminal-conflict case by retaining its original infeasible network. Every optimized quotient is lifted, its original cut value recomputed, and its continuous boundary checked in original coordinates. Parallel-arc aggregation is applied identically to both constructions.

This standard simplification substantially shrinks the height reference; its effect is measured rather than assumed absent. It preserves the comparison's objective but does not remove cost of constructing the original critical mesh.
