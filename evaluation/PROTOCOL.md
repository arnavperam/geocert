# External corridor transfer check

Protocol fixed before extraction or optimization of this dataset.

## Purpose and scope

Test whether the existing corridor extraction and exact sensitivity workflow transfer to a previously unused public breast-tissue example. This is a technical transfer check, not segmentation accuracy, anatomical calibration, or an independent-patient efficacy study. The vendor supplies this two-FOV subset to demonstrate file formats; its aggregate metadata must not be interpreted as a population sample.

Source: 10x Genomics, Xenium Human Breast Gene Expression, multimodal segmentation, XOA 2.0 two-FOV example. Documentation: https://www.10xgenomics.com/support/software/xenium-onboard-analysis/2.0/resources/xenium-example-data . Archive: https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip . Published MD5: 4632914eca973a1d532231ea646e10cc.

## Selection fixed in advance

Use two 200-by-200 micrometre crops. From the supplied cell-centroid table, choose the cell nearest the coordinatewise 25th-percentile point and the cell nearest the coordinatewise 75th-percentile point. Break distance ties by cell identifier. Center one crop at each selected centroid. If the crops overlap, retain both and report overlap; do not search for replacements. Selection uses geometry only, before RNA-level outcomes.

Retain measured-gene RNA with Q score at least 20; identify gene features from the vendor feature table or panel. Use native cell and nuclear raster masks, their recorded coordinate transform, and the cell-mask owner of each nuclear pixel. Record disagreements between transcript-assignment fields and raster-derived ownership; do not silently alter inputs to improve agreement.

Apply the existing extraction rules at both cardinal orientations, window heights 12 and 24 native rows, and two or three ordered cells. Keep the original 3--24 RNA eligibility rule, 0.001 micrometre clearance, 0.5 micrometre band padding, source-admissible slopes max(1, maximum source slope), and minimum/maximum gaps at one-half/twice the source gap extrema. Preserve all intersecting nuclear regions. Select nonoverlapping eligible charts within each crop/setting in lexicographic SHA-256 order, at most 12 per crop/setting. Record every rejection category, all eligible charts, and all selected charts. Zero eligible charts is a valid result, not grounds to change the screen.

## Questions

For every selected chart, retain the baseline and the previously defined enlarged family (double slopes and padding, endpoint extension of half the original length). Evaluate penalties 1/4, 1 and 4 for the existing RNA-reassignment objective. Compare exact joint cost with the feasibility-matched singleton-cost surrogate, with favorable tie handling for the surrogate. Record positive and zero regret, failures, and exact witness checks separately. This remains an energy ablation, not a comparison to a segmentation method.

Evaluate each measured gene represented in a selected chart as a cell-count readout where source-cell count accounting is complete inside the crop. Outside-chart transcripts keep their source assignments. Report count extrema and all exclusions. No marker panel, positive example, geometry, or penalty is chosen from the resulting answers.

## Verification and failure policy

Use the frozen v1.0.0 exact core. Check all returned boundaries at the original coordinates and verify integer flow/cut certificates. Record dependency versions and hashes. No result is called exact after a timeout. Allow 20 minutes for extraction and 15 seconds per exact optional question; retain censored cases. Do not replace cells or charts because their results are uninteresting. Implementation repairs may correct parsing only, with the change and its effect recorded. Changing scientific selection or parameters requires a separately labeled amendment.

This protocol is internal prospective experimental bookkeeping, not a registered clinical or statistical study. No new biological-generalization claim follows from this check alone.
