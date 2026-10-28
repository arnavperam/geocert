# GeoCert: compatible cell-count certificates and an exact hard-budget algorithm

## Executed advance

This continuation addresses the current critique's two mathematical gaps: separate overlapping corridor models did not provide one cell-level uncertainty family, and the expression version omitted the original cost-representation proof. The merged manuscript restores that proof and adds a disjoint, topology-preserving patch atlas. It also implements a new exact critical-height/count dynamic program for scalar patches with explicitly represented point RNAs. The dynamic program returns all count minima, including unsupported hard-budget outcomes, without enumerating every RNA assignment.

The result is a conditional sensitivity system, not the completed biological-validation program in the earlier request. No acceptance probability or historical-firstness claim is made.

## One geometric family and one count per cell/gene/budget

A reference raster partition is preserved exactly. Geometry-only extraction follows every eligible shared pixel interface, splits it into monotone arcs, reserves fixed seams, and accepts interior-disjoint physical ribbons. Source counts include all retained molecules in crop-contained connected reference cells; molecules outside ribbons remain fixed. Background is label zero, not a biological cell.

The chart is a shear and rotation with constant physical Jacobian p²/2. A continuous transverse map sends the reference interface to an admissible graph while fixing every ribbon boundary. These maps glue into a homeomorphism of the complete reference partition. The proof preserves topology rather than relying on an untested assumption about connected residual cell cores. Changed-label area is counted once across patches.

The fixed-seam qualification matters. Uncovered boundary portions, junction neighborhoods, crop-border cells and structural labels stay fixed. The atlas does not represent all realistic membrane errors. All results retain an atlas identity, coverage record and coordinate/cost definitions.

## Complete data and count results

| Quantity | Lung development crops | Breast format-test crops |
|---|---:|---:|
| Retained observed RNA | 7,370 | 81,182 |
| Crop-contained connected reference cells | 153 | 467 |
| Accepted disjoint patches | 544 | 1,708 |
| Editable distinct RNA | 745 | 8,262 |
| Editable/all incident reference edges | 15,821/23,164 | 52,912/80,463 |
| Boundary-length fraction | 68.3% | 65.8% |
| Exact cell/gene spectra | 3,304 | 22,976 |
| Eight-budget interval rows | 26,432 | 183,808 |
| Largest interval width | 3 | 10 |
| Rows tightening a strong marginal outer bound | 4 | 228 |
| Rows whose endpoint the additive-cost inner model misses | 2 | 68 |

The 620 cells are reference objects, not independent biological replicates. Five lung cells have no positive or editable retained count and are represented by implicit zero entries under the matrix index. The sparse export covers every positive-source or editable cell/gene entry; it does not infer biological absence or unmeasured expression.

Most endpoints equal the simple baselines. The outer marginal model can spend the budget separately on each RNA and combine incompatible changes. The additive inner model retains feasibility but overcharges shared deformation. They are different approximation contracts, not two state-of-the-art segmenters.

### Two-patch expression example

In breast format-test crop 1, reference cell 81 has 22 retained LDHB transcripts. At 0.25 µm² its exact count interval is [19,23], whereas the marginal outer bound is [17,23]. At 0.5 µm² the exact interval is [17,23], whereas the additive inner model gives [18,23].

The minimum deformation attaining 17 is exactly

1977348833190267596699567429 / 4085091626228121600000000000 µm²

or approximately 0.484040265 µm². The witness transfers three LDHB molecules through one patch and two through a second, together with one collateral KRT8 molecule. It respects all other patches at their source states and spends one common budget. This is a post-hoc illustration from the full census, not a biological finding about LDHB or corrected segmentation.

## New exact hard-budget class

For one continuous boundary with a rational piecewise-affine reference, point-valued free RNA groups and fixed structural intervals, a common critical mesh contains a canonical minimum for every assignment. Nodal area is exact at those minima. A height/count shortest path therefore yields the exact minimum area at each count.

If W is the variable count range and H_j the number of candidate heights at node j, a monotone queue computes the dynamic-programming transitions in O((W+1) sum_j H_j) arithmetic operations after the mesh and local node scores are formed. Preparation, score formation, witness output/checking and rational bit lengths are separate. For explicitly listed unit molecules, W is bounded by the relevant input count. Binary-compressed multiplicities would instead give a pseudo-polynomial dependence; free RNA-location intervals and multiple moving boundaries are outside this chain algorithm's contract.

Min-plus composition of patch spectra is exact under the verified gluing conditions and produces one interval per complete reference cell under one atlas. The original min-cut representation remains the engine for arbitrary unary penalized queries and sound Lagrangian bounds; it is not falsely claimed to solve every hard-budget problem in one cut.

### A duality gap with explicit unit molecules

An eight-point example, with no succinct bundles or binary multiplicities, has exact count/area spectrum {0:0, 2:4, 3:9, 5:13, 6:18, 8:22}. All eight point locations are distinct. At area 10 the exact maximum is 3, while the best nonnegative-penalty upper bound is 4. The count chain recovers the unsupported count-three state. This replaces the weaker bundle-only illustration; it is a gap example, not a planar NP-hardness theorem.

### Complete validation of the new algorithm

All 944 multicopy local gene spectra from the two data collections match separately enumerated exact minima. A second full end-to-end calculation uses the chain for every multicopy query and the source-feasible closure fast path for single-copy queries. All 26,280 global spectra and all 210,240 CSV interval rows match the independently checked enumerated workflow exactly. Every produced witness is directly checked; tied witness selection is not required to be identical.

Additional tests cover 201 scalar spectra, 80 nonflat models with 1,644 owner vectors, positive mandatory deformation, malformed/free-region rejection, and the explicit unit-count gap. Standard resource-state dynamic programming and critical-height constructions are credited; the paper's claim is the exact continuous realization and compatible budget computation, not invention of dynamic programming.

## External data and comparison limits

A completed GitHub workflow retained a downloadable breast acquisition artifact. It was retrieved, hash-verified and re-ingested. The same atlas rules were applied without tuning to its results. The vendor explicitly describes this two-FOV dataset as trimmed data for testing output formats, unsuitable for biological conclusions. It is therefore an external-data transfer test, not ground truth, a clinical cohort or validated biological replication.

An actual image-only alternative was computed: nearest labelled nuclear pixel within 5 µm, background otherwise. Original RNA positions are fixed, and a crop-border guard prevents using absent outside nuclei. Of 1,128 changed lung cell/gene counts, 185 are attainable within the atlas; for the breast test the corresponding numbers are 1,818/11,151. All four full retained alternative assignment vectors are infeasible. This is an adverse model-reach result. It cannot be converted into a confidence calibration or a claim that the alternative is wrong.

Proseg was not executed. Its repository/release and installation routes were inspected, but executable acquisition/build was not available through the current environment. The GitHub connector exposed read/download actions but no create, dispatch or push operation. No implemented nuclear-distance rule is labelled Proseg, Baysor or ResolVI. No competing method's accuracy, AUROC, FDR or runtime was invented.

## Broader biology requested earlier

Reporter/mosaic ownership truth, xenograft species truth, held-out biological donors, membrane-calibrated parameter estimation, statistical DE/FDR validation, a published biological claim reanalysis, and whole-slide 10^8-transcript deployment are not completed. The previous synthetic conformal and probabilistic demonstrations remain development checks, not replacements for those experiments. They are not the merged paper's headline.

The current marker-pair census supplies exact co-detection predicates, but its results do not establish cell type, leakage, doublets or signaling. Simple marginal intervals already determine all tested predicate classifications.

## Independent checking and reproducibility

Source ingestion checked 88,552 retained records against prepared identities/counts and original masks. Lung ingestion uses the inherited restricted Parquet decoder; breast ingestion uses ordinary CSV/gzip. This shared-dependency limitation is explicit.

A separate rational scalar implementation imports no production implication, envelope, count-chain, integrated-layer or flow routine. It reconstructs 7,774 local spectra and checks global spectra against 36,635 Cartesian combinations. It checks 33,675 global attaining count witnesses and all 210,240 unique interval rows. Geometry checks cover 2,252 disjoint patches, 137,466 half-ribbon/pixel intersections, 2,093 relevant patch-pair intersections and 49,341 nuclear fragments. These are internal independent-code-path checks, not independent human proof review.

Twelve complete fresh-process constructed runs measure all count levels for 16,32,64,128 point RNAs. Median total times are 0.465,1.574,7.138,42.206 seconds. The 128-RNA run visits 605,632 states and returns all 129 levels with checked witnesses. The family has a known analytical solution, so this is a larger operating-range demonstration rather than a general solver-superiority claim.

The new core primitives from the preceding expression release remain byte-identical (25 Python files). New atlas and count-chain programs have their own hashes. The final frozen-package replay and document preflight are recorded separately; stored timing censuses are not silently counted as rerun.

## Submission scope

This is a stronger, self-contained exact-methods candidate with a sound global model and a new explicit-count hard-budget algorithm. It is not the full ground-truth-validated “perfection” program. Its novelty and significance require specialist assessment; internal checking cannot guarantee acceptance.

No conference submission, abstract registration, public push of this extension or license selection occurred. The current public base repository does not yet contain this local atlas release. Public documentation, an actual source revision, author-approved licensing and accurate assistance disclosure remain release tasks.
