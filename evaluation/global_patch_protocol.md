# Global patch evaluation protocol

Recorded before the new patch construction and Proseg outputs are inspected.

The scientific target is one observed cell-by-gene interval under one globally compatible, explicitly restricted patch-deformation family. Existing overlapping corridor outputs remain historical local certificates; they will not be averaged or relabelled as whole-cell bounds.

## Data and construction

Use both existing public Xenium lung crops without selecting additional successful cases. Retain every supplied Q20 measured-gene transcript and all cell/nuclear mask records. Exclude crop-border or disconnected source cells only from claims about complete observed cell profiles, retaining their records as geometric neighbours. Scan all two-cell shared interfaces. Protect junctions and require disjoint patch interiors, fixed traces on patch boundaries, unique transcript membership, and preserved nuclear regions. Construct patches from geometry without consulting gene identity or certificate outcomes. Report the frozen complement and extraction failures. Gene-count budgets are 0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4 square micrometres. Model parameters are sensitivity settings, not biological calibration.

For each cell and gene, convolve all incident patch count/cost spectra with a single budget, while counting every outside-patch transcript exactly once. No min/max aggregation across inconsistent corridor models is permitted. Joint-cell results must also preserve transcript conservation. Incomplete target enumeration must yield an explicitly bounded/unresolved output, never an exact label.

## Published-method diagnostic

Run unmodified Proseg 3.2.0 on each complete supplied crop. The primary comparable run uses its Xenium preset, two-dimensional coordinates, one voxel layer, and disables transcript diffusion, so the geometric readout uses measured positions. A secondary run uses the Xenium preset with default dimensionality and diffusion. Use the published sampling defaults and two threads. Retain process status, version, binary hash, commands, input hashes, logs and output tables. Proseg is stochastic; one run per crop/configuration is a feasibility/comparison diagnostic, not a replicated accuracy ranking. Do not substitute a surrogate implementation.

The primary readout is the minimum certificate budget needed to contain Proseg's retained observed counts, with unmatchable cell identities, unsupported transcript changes, and unreachable counts reported explicitly. An independent segmenter is a comparator, not anatomical truth. Do not use its output to retune the patch family before reporting this comparison. Standard Proseg output may remove background or move RNA, which are distinct from the certificate's fixed-observation geometry.

## Claims and release

This branch does not alter main or v1.0.0. No conference submission or software licence grant is implied. New results are development evidence unless a separate external dataset and a locked evaluation are actually completed. The paper will contain the ownerwise reduction and its proof as well as the count application. AI contributions to mathematical development, implementation, experiments and writing must be disclosed accurately.
