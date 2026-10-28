# Small-crop Proseg amendment

The locked Proseg 3.2.0 run and repeated diagnostic fail on roi_01 in featureselection.rs with EmptyInput before segmentation. Both initial roi_02 runs complete. Original commands, statuses, binary hashes, and returned data remain preserved.

Before evaluating count compatibility, repeat the same fixed-XY and default-coordinate settings on BOTH crops with the upstream documented --no-factorization option. This bypasses expression-feature factorization without modifying Proseg's source. Preserve all four initial outcomes and all amendment outcomes; do not replace a failure with an unlabelled success. Sampling defaults, two threads, observed inputs, and geometric-certificate settings are unchanged.

Compare two distinct readouts: (1) reported assignments, including explicit transcript filtering/background exclusion, and (2) consensus polygons evaluated at every original retained coordinate. The second readout isolates a geometric point-assignment comparison but is not identical to Proseg's three-dimensional or diffusion-aware output. Preserve unmatchable cells and out-of-family changes. No output is anatomical truth, no parameter is fitted to improve coverage, and one stochastic run does not estimate performance variation.

The disjoint-patch protocol is also applied unchanged to the previously locked public breast crops on evaluation/breast-transfer. This is an external-tissue transportability check, not an independent blinded biological ground-truth study: those crops were already processed by a separate branch of this project. All extraction failures and the frozen boundary fraction remain reported.
