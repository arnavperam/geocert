# GeoCert: biological replication and segmentation context

31 October 2026

The shared-boundary effect replicates in two biological breast cancer samples.
It remains a modest part of the overall count uncertainty, concentrated at small
area budgets. A separate context experiment shows that the very high background
fractions in the original small lung Proseg runs depend on the surrounding input.
Neither experiment changes the atlas geometry or supplies anatomical ground truth.

## Biological donor replication

We used the IDC and ILC samples from the public 10x
[FFPE Human Breast with Pre-designed Panel dataset](https://www.10xgenomics.com/datasets/ffpe-human-breast-with-pre-designed-panel-1-standard).
These are biological tumor samples from two donors, processed with XOA 1.0.2.
The native reference segmentation expands DAPI nuclei by 15 µm.

For each donor, three 200 µm regions were selected using cell-centroid coordinates
before querying expression. All six regions were retained. They are disjoint and
contain 1,518 distinct complete cells and 189,552 retained Q20 transcripts.
The donor behind the older breast format-test dataset is unspecified; we do not
count it as another independent breast donor.

| Result across three regions | Ductal donor | Lobular donor |
|---|---:|---:|
| Complete cells | 659 | 859 |
| Cell–gene questions | 32,360 | 36,409 |
| Eligible questions: outer width reaches at least two | 1,271 | 829 |
| Questions tightened by shared geometry | 176 (13.8%) | 97 (11.7%) |
| Tightened as a fraction of all questions | 0.54% | 0.27% |
| Questions with an endpoint missed by additive costs | 58 | 37 |
| Invariant at 4 µm² within the primary atlas | 79.6% | 83.5% |

The invariance rates are conditional on the unchanged primary atlas: displacement
is capped at about 0.45 µm and about 76.5% of eligible shared edges are editable.
The results do not certify invariance under unrestricted segmentation errors.
Questions are positive or editable cell–gene pairs, and repeated budgets are not
independent samples.

The exact intervals tighten the matched independent-molecule bounds at budgets
0.05, 0.1 and 0.25 µm². On the tested grid, those bounds equal the exact intervals
from 0.5 µm² onward. Additive costs also miss attainable endpoints at smaller
budgets. The effect therefore replicates without becoming the dominant source
of count uncertainty.

All 550,152 interval rows have complete exact spectra. Independent scalar
calculations reproduce the local minima; Cartesian composition and explicit
boundary checks reproduce the global counts and areas. Native masks and original
Parquet transcript records were reconciled independently, with zero source-owner
discrepancies in the six retained cores. A conservative distance-only comparison
also contains every exact interval, but it relaxes several geometric constraints
and is not used as the sole evidence for a joint effect.

## Proseg with surrounding tissue

The six breast cores and two original lung cores were each run with 200 µm and
800 µm inputs. Default factorization and no factorization were kept separate,
giving 32 planned settings. Thirty-one completed successfully; default factorization
on the first small lung core failed again and is retained as a failed setting. All use Proseg 3.2.0, two threads and the upstream
200 burn-in, 200 sampling and 50 hill-climbing iterations.

In the original lung cores, matched no-factorization runs give:

| Core | Background, 200 µm input | Background, 800 µm input |
|---|---:|---:|
| Lung 1 | 42.5% | 19.4% |
| Lung 2 | 89.8% | 62.1% |

Each comparison uses the same core transcript identifiers across completed
settings. Missing records, source-cell identity mapping and background assignment
semantics are accounted for explicitly. This avoids comparing the whole large
window with the small core or silently changing denominators.

The full setting inventory, retained/input transcript totals, background fractions,
changed-count denominators, countwise reach and transfers outside editable ribbons
are in `results/proseg_summary.csv` and the supplement. Each setting is a single
stochastic run. The differences are evidence of sensitivity to context, not a
precise causal estimate or a segmentation-accuracy ranking.

Every completed full assignment vector lies outside the primary atlas. This is
different from an individual cell–gene count falling inside its interval: separate
counts can be attainable using different boundaries. The new comparisons do not
erase the fixed-seam or coverage-versus-reach limitations.

## Manual whole-cell reference

The public SPATCH files were downloaded and inspected. Its
[paper](https://www.nature.com/articles/s41467-025-64292-3) describes the manual
annotations as nuclear boundaries. Nuclear labels cannot be substituted for
whole-cell boundaries when assigning cytoplasmic transcripts.

This part is **not completed as a whole-cell accuracy or calibration experiment**.
The specific limitation is annotation semantics, not access or a failed download.
The inspected files and decision are recorded in `manual_reference_eligibility.json`.
Obtaining suitable expert whole-cell annotations remains necessary for that claim.

## Manuscript changes

- Added the two-donor replication to the main Results and a full protocol, source
  figure and tables to the supplement.
- Added the matched lung context finding and all Proseg settings to the supplement.
- Removed the repeated LDHB numbers from the Discussion; retained the shared-cost lesson.
- Kept the registered title and abstract, theorem statements, proofs and numbered
  equations unchanged.
- Retained the ten-page body, separate title page, public source link and disclosure.

The paper now has evidence from biological breast samples in addition to the
lung and format-test inputs. The geometry remains narrow, and no result here
establishes biological calibration or whole-cell assignment accuracy.

## Reproduction

`README.md` gives the commands for source acquisition, fresh analysis and replay
of the saved evidence. `sources.json` records retrieved archive-member hashes;
`artifacts/manifest.json` binds the prepared inputs, spectra, witnesses and Proseg
core assignments. The original count and geometry implementations are unchanged.

The public compact artifacts support exact replay of the reported tables and
checks on prepared data. Complete native Proseg outputs are preserved separately
in the evidence archive; their reconciliation with output count matrices was
performed before extracting core assignments. Repeating the Proseg sampler is
stochastic and need not give the same assignments.
