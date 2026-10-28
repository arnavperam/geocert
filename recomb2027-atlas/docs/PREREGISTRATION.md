# Frozen reach-analysis plan

Recorded before the new ribbon-width extraction and count computations. This is an outcome-frozen sensitivity analysis on already inspected inputs, not a registered external validation study. The inherited atlas and all of its outcomes remain the primary reference.

## Questions

Does increasing the physical support width of a ribbon increase the set of representable count changes, and how much boundary coverage is lost to geometric exclusion and disjointness? Do the supplied Proseg assignments lie in this model?

## Fixed design

Use all four included crops, all complete connected reference cells, and all retained measured RNA. Keep the original geometry-only candidate arcs, ordering, endpoint trimming, slope 1, clearance, nuclear restrictions and fraction 3/4 of support width for the displacement band. For each physical support halfwidth 0.5, 1, 2 and 3 micrometres, extract a separate atlas. Convert to chart units B by flooring 4*sqrt(2)*r/p and dividing by 4; report the actual resulting support radius pB/sqrt(2) and displacement cap 3pB/(4sqrt(2)). A required-width run has no narrower fallback. A candidate that cannot fit is rejected. Report all four widths, including empty or unsuccessful settings. These separately selected atlases need not be nested. The original adaptive-width atlas is reported separately.

Use area budgets 0, 0.005, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.375, 0.5, 0.75, 1, 2, 4, 8 and 16 square micrometres. Compute complete spectra with the existing exact count chain; one-copy targets may use the proved closure shortcut. Keep original 60-second local limits. Resource-limited tasks are incomplete, not exact or preserved.

Report accepted patches, reference-edge coverage, RNA inside patches, non-singleton counts, summed interval widths, endpoint changes between budgets, and tightening among questions whose marginal width is at least two. Keep all eligible cells and fixed entries.

## Alternative references

Reuse the original deterministic 5-micrometre nuclear-expansion labels with the same crop-edge guard. In addition ingest both supplied Proseg 3.2.0 runs, mapping output cells through original cell IDs and joining transcript IDs. Verify observed XY/gene consistency and keep exclusions. These are one adapted stochastic run per clipped lung crop, with factorization disabled after the recorded failure. They are an execution/reach diagnostic, not accuracy, calibration, a five-run experiment, or a fresh Proseg execution.

For Proseg, compare a fixed set of complete-cell/gene entries on retained output transcripts. Exclude every key whose full count or editable contribution can involve a Proseg-excluded transcript in any tested atlas. Report this exclusion set. No excluded molecule is assigned an invented Proseg owner. For both references, countwise inclusion is weaker than simultaneous feasibility of the complete matrix. Report unbounded-budget representability and budget-dependent reach; do not interpolate an unsupported new geometry.

## No outcome-driven changes

No cell, gene, width or unfavorable case is omitted to improve success rates. Any necessary implementation repair is recorded separately. The analysis does not enlarge the reference data, calibrate physical parameters, establish true membranes, or promote format-test breast input to biological validation.
