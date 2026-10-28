# Frozen reach-study protocol

This continuation evaluates the supplied Atlas model and incorporates the other
agent's existing Proseg output. It is a post-review sensitivity study, not a new
anatomically calibrated family or an independent-donor study.

## Locked before any new width-specific outcomes

Use all four supplied crops: roi_01, roi_02, breast_01, breast_02. Keep the original
connected/non-border cell selection, 128-edge run limit, fixed end edges, candidate
ordering, orientation rules, zero displacement endpoints, slope one, clearance,
structural nuclear constraints, exact pixel checks, and disjoint-ribbon checks.

Width schedules in chart units (each proposal tries them in order):
- cap4: 4,2,1,1/2 (original reference)
- cap8: 8,4,2,1,1/2
- cap16: 16,8,4,2,1,1/2
- cap20: 20,16,8,4,2,1,1/2

At pixel size p, B gives physical ribbon halfwidth p B / sqrt(2); the allowed
transverse displacement is at most 3 p B / (4 sqrt(2)). A 0.2125-micrometre pixel
and B=20 therefore allow a halfwidth about 3 micrometres and displacement about
2.25 micrometres. Report realized widths too; fallback and endpoints can restrict
reach. Each schedule regenerates a geometry-only atlas. Families need not be nested.

Budgets in square micrometres:
0, .01, .025, .05, .075, .1, .15, .2, .25, .35, .5, .75, 1, 2, 4, 8, 16.
All budgets use the same complete count spectra; never use a timeout as exactness.

Report all proposed/accepted patches, editable unique interface length, editable
RNA, allowed displacement, counts and widths at every budget, questions with
marginal width at least two at any listed budget, all tightened questions, and all
cells with a tightened gene. Retain the all-question and all-row denominators.
Non-singleton indicator saturation is distinct from complete endpoint saturation.

The 5-micrometre nearest-nuclear reference is identical to the prior experiment.
For Proseg, import the single recorded no-factorization run per lung crop. Verify
raw input/output checksums, transcript identity, original observed XY, gene, cell
mapping and excluded IDs. Use common output RNA for count comparison. Keep all RNA
as geometric constraints; subtract fixed excluded source counts. Exclude a cell/gene
comparison if an excluded transcript is editable at that key. Report exclusions
and changed-count reach, not accuracy. Retain background/absent destinations and
fixed-exterior failures. No crop-level run is relabeled as a haloed/full-slide run.

Compute exact spectra by target enumeration for small target multiplicity; use
the independently tested count chain if necessary. Cross-check against the other
algorithm for the full new lung spectra and a predetermined ordered breast sample.
No successful claim will depend on choosing a width or example after seeing it.
The original breast data remain explicitly labeled vendor format-test data.

The manuscript will use the validated one-boundary application, keeping the
coupled multi-boundary theorem with its separate evidence in the supplement;
no unvalidated three-label atlas will be described as completed.
