# Author briefing: exact expression intervals

## The paper to defend

**Certified Expression Intervals under Shared Segmentation Uncertainty.**
The main evaluated model has one interface per disjoint patch. Its count-chain dynamic program computes exact hard-budget spectra, and the gluing theorem makes patch combinations one valid partition. The supplement gives the more general coupled-interface cost theorem. Do not claim that thickness coupling is responsible for the atlas results.

## Gluing: why separate changes coexist

For chart halfwidth B and displacement h(t), map the transverse coordinate v by

    T(t,v) = -B + (v+B)(h(t)+B)/B       for -B <= v <= 0,
             h(t) + v(B-h(t))/B        for 0 <= v <= B.

Both pieces meet at h(t). Their slopes are positive when |h|<B; the two rails remain fixed. When h=0 at the longitudinal endpoints, the entire outer trace is fixed. Disjoint supports therefore extend by the identity to a common homeomorphism of the reference partition.

Worked example: B=2 and h=1/2 give transverse slopes 1.25 and 0.75. The middle interface moves but no fold or overlap is introduced. For h(t)=min(t,2-t)/2 on [0,2], the chart area is 1/2. With pixel size p=0.2 micrometres and Jacobian p^2/2, its physical changed area is 0.01 square micrometres. Two disjoint copies cost 0.02, not two independent budgets of 0.02.

## Critical heights (Lemma 5)

Candidate events contain source knots, constraint endpoints and cone intersections before labels are selected. A canonical minimum is affine with displacement of one sign between events. Its trapezoidal integral is exact. Other interpolants may have larger trapezoidal costs; the retained canonical optimum gives the equality of minima.

Worked example: source zero on [0,2], slope one, and a required height 1/4 at t=1. Its attaining tent has knots 3/4,1,5/4 and area 1/16. A mesh containing only 0,1,2 charges 1/4 instead. The critical mesh includes the missing corners. We preserve minimum cost for each assignment, not the area of every interpolant.

## Hard budgets (Theorem 6)

At each critical node, a state records candidate height and accumulated gene count. Slope-compatible predecessors form a contiguous candidate range. A range-minimum sweep performs the transition for every count; local transcripts are counted once. Backtracking gives an attaining boundary. Bounds are polynomial in explicitly listed molecules, not in an arbitrarily large binary-compressed multiplicity.

Worked example: eight distinct unit transcripts form compatible groups of sizes 2,3,3 with costs 4,9,9. Attainable count/minimum-area pairs are (0,0),(2,4),(3,9),(5,13),(6,18),(8,22). Budget 10 admits maximum count 3. Every penalty bound is at least 4: mix (4,2) with weight 2/3 and (22,8) with weight 1/3. This is why a complete penalty path is not a complete hard-budget solution.

## Explain the biological output

Cell 81 in the breast format-test crop has 22 retained LDHB molecules. Its interval is [19,23] at area 0.25 and [17,23] at 0.5 square micrometres. Count 17 is attained at approximately 0.484040265 by moving five LDHB molecules through two disjoint patches, plus a collateral KRT8. This is a conditional measurement statement, not corrected ownership or a biological LDHB finding.

## Interpret the width study

The physical support halfwidth is pB/sqrt(2); the movement ceiling is 3pB/(4sqrt(2)). B=4 and p approximately 0.2125 permit at most 0.451 micrometres, while B=20 permits at most 2.254. Wider proposed ribbons overlap more often, so their accepted families are not necessarily nested. Report every frozen width, coverage and changed-count denominator. A fixed non-singleton indicator is weaker than fixed endpoints: 61 breast intervals still change after budget 0.25.

## Two comparisons, two conclusions

Joint versus singleton-additive geometry tests information lost by summing individual costs. Direct versus critical-height cut implementations optimize the same energy and test representation cost. Proseg/nuclear-reference reach instead asks whether alternative *counts* can be represented in the atlas. It is not an accuracy score, and separate reachable entries need not jointly realize a whole alternative matrix.

## Boundaries of the evidence

The old pair census is not new. The two Proseg runs were supplied by the other agent, adapted after default factorization failed, and use clipped context with high background assignments. They were re-ingested and checked, not rerun as five replicates. The breast dataset is vendor format-test data. No new true ownership, donor calibration or statistical differential-expression result is established. Know the assumptions and every denominator before registration or presenting the work.
