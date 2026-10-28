# Cell-count readout audit (frozen before outcomes)

Question: under the unchanged local geometric families, which cell/gene readouts can change, and how much of each source-cell count lies outside the modeled corridor?

Inputs: all 42 published corridors, baseline and largest pre-existing settings. No new corridor, gene panel, geometry setting, or penalty is selected for a favorable result. Only source cells whose mask does not touch the supplied crop edge are eligible for complete count accounting. A source cell crossing the crop boundary is excluded and recorded.

For each eligible corridor/cell and every gene represented by at least one RNA anywhere in that corridor, count source RNA in the complete retained Q20 measured-gene crop; subtract that cell's in-corridor RNA to obtain an immutable out-of-corridor contribution. Independently minimize and maximize the resulting cell count over the joint feasible family. Use zero deformation penalty for exact feasible count extrema. Reconstruct one minimum-deformation boundary for each optimal count. RNA outside the corridor remains fixed. Report all outcomes including invariant readouts.

Additionally evaluate the original three deformation penalties (1/4, 1, 4) for count decrease and increase on every eligible baseline readout whose feasible count interval is nonconstant. These are penalized sensitivity optima, not hard-area-budget extrema. Selection on nonconstant intervals is an explicit efficiency filter; it does not select genes by an observed positive comparison result.

Verification: independently recompute full counts from selected owners plus fixed contributions; verify every boundary using the original-coordinate checker; check source inclusion and count conservation; compare an outcome-blind subset with the exact critical-height representation. Retain failures. For endpoint comparisons, arbitrary owner ties are allowed, but gene counts must agree. No classification accuracy, true assignment, new donor, statistical confidence, or global stitched segmentation is inferred.

This is a new analysis of existing development data, not a held-out validation. Generic ability to optimize linear readouts follows from the existing theorem and is not a new algorithmic claim.

Implementation clarification before solving: choose each count extremum first and minimum deformation second using the exact scalarization `D +/- M*count`, with integer `M > |I| sum(boundary band widths)`. This is equivalent to a zero-penalty count extremum with an exact deformation tie-break and avoids falsely interpreting an arbitrary count-optimal witness as the least-deforming one. Hash-select 1/25 of readouts, independently of outcomes, for two-sided critical-height verification in both settings.
