# Replication on biological breast tissue

This analysis applies the existing primary atlas to six predefined regions from
two public breast cancer donors, and checks how the Proseg comparison changes
when the input includes surrounding tissue. The geometry and count algorithms
are unchanged. `protocol.json` records the selection rule before count queries.

## Results

| | Ductal donor | Lobular donor |
|---|---:|---:|
| Complete cells in three regions | 659 | 859 |
| Cell–gene questions | 32,360 | 36,409 |
| Questions with outer width at least two | 1,271 | 829 |
| Questions tightened by shared geometry | 176 (13.8%) | 97 (11.7%) |
| Questions with an additive missed endpoint | 58 | 37 |

Questions are counted once across the original eight budgets. Tightening affects
0.54% and 0.27% of all questions, respectively, and occurs at the smaller budgets.
The exact and independent-molecule intervals agree from 0.5 µm² onward on this
grid. The new samples replicate the effect; they do not make it a large fraction
of all expression counts.

The primary atlas covers about 76.5% of eligible shared raster edges here and
allows at most about 0.45 µm of displacement. Counts outside that family remain
fixed. All 550,152 interval rows have complete spectra and verified attaining
boundaries. Independent source reads reconcile the native masks, transcript
identities and counts. The six cores are disjoint, with 1,518 distinct complete
cells and 189,552 retained transcripts.

`results/proseg_summary.csv` reports every core/context and factorization setting,
including the failed default run on the first small lung crop. In the original
lung cores, matched no-factorization runs assign 42.5% versus 19.4%, and 89.8%
versus 62.1%, of the same core transcripts to background when input width grows
from 200 to 800 µm. Each setting is one stochastic run. These are sensitivity
comparisons, not reference truth. Countwise reach and complete-assignment
representability are reported separately.

## Replay the published evidence

From the repository root, using Python 3.12 or later:

```bash
python3 -m pip install numpy==2.3.5
python3 recomb2027-atlas/extensions/external_validation/restore_artifacts.py --output work/external-validation
python3 recomb2027-atlas/extensions/external_validation/replay.py --analysis work/external-validation
```

The archives contain prepared source crops, native masks, exact spectra,
boundary witnesses, comparison tables, Proseg core assignments and run/audit
records. Restoration verifies every saved file hash. Replay independently checks
the geometry and exact spectra, reconstructs witnesses, and recomputes all core
comparisons and donor summaries.

This replay uses saved Proseg assignments. It does not rerun the stochastic
sampler or repeat the original audit against complete native Proseg outputs.
That audit was performed before extracting the core assignments. Its results
and native-output hashes are preserved. An optional Proseg CSV export was
truncated in one completed run; the complete native Parquet and count matrix
passed reconciliation and were used throughout the analysis.

## Repeat from the vendor data

The input catalog in `sources.json` gives the official dataset pages, archive
URLs and hashes of the retrieved members. Whole archives were not downloaded
or assigned invented hashes. The breast sources are the IDC and ILC samples
from the 10x **FFPE Human Breast with Pre-designed Panel** release (XOA 1.0.2).
Their cell masks use DAPI nuclei and 15 µm expansion. The two donors are distinct
from the lung source; independence from the earlier breast format-test donor
is unknown.

Install the versions recorded in `requirements.txt`, then:

```bash
python3 recomb2027-atlas/extensions/external_validation/fetch_inputs.py --output work/raw --donors IDC ILC
python3 recomb2027-atlas/extensions/external_validation/prepare.py --raw work/raw/IDC --output work/new-analysis
python3 recomb2027-atlas/extensions/external_validation/prepare.py --raw work/raw/ILC --output work/new-analysis
python3 recomb2027-atlas/extensions/external_validation/run.py work/new-analysis/IDC/idc_01 work/new-analysis/IDC/idc_02 work/new-analysis/IDC/idc_03 work/new-analysis/ILC/ilc_01 work/new-analysis/ILC/ilc_02 work/new-analysis/ILC/ilc_03
python3 recomb2027-atlas/extensions/external_validation/audit_inputs.py --raw work/raw/IDC --analysis work/new-analysis
python3 recomb2027-atlas/extensions/external_validation/audit_inputs.py --raw work/raw/ILC --analysis work/new-analysis
python3 recomb2027-atlas/extensions/external_validation/distance_baseline.py work/new-analysis/IDC/idc_01 work/new-analysis/IDC/idc_02 work/new-analysis/IDC/idc_03 work/new-analysis/ILC/ilc_01 work/new-analysis/ILC/ilc_02 work/new-analysis/ILC/ilc_03
```

Selection uses source cell centroids, with no expression-based filtering. All
three selected cores per donor are retained. The adapter preserves original
physical coordinates, IDs, raster masks and quality filtering.

Use Proseg 3.2.0 for the context study:

```bash
python3 recomb2027-atlas/extensions/external_validation/proseg_runs.py --binary /path/to/proseg work/new-analysis/IDC/idc_01 work/new-analysis/IDC/idc_02 work/new-analysis/IDC/idc_03 work/new-analysis/ILC/ilc_01 work/new-analysis/ILC/ilc_02 work/new-analysis/ILC/ilc_03
python3 recomb2027-atlas/extensions/external_validation/evaluate_proseg.py work/new-analysis/IDC/idc_01 work/new-analysis/IDC/idc_02 work/new-analysis/IDC/idc_03 work/new-analysis/ILC/ilc_01 work/new-analysis/ILC/ilc_02 work/new-analysis/ILC/ilc_03
```

The runner uses two threads per job and the upstream 200 burn-in, 200 sampling
and 50 hill-climbing iterations. Default and no-factorization settings are kept
separate. Fresh stochastic runs need not reproduce the recorded assignments.
Use fresh output directories; the runner deliberately does not treat an
unverified existing directory as a successful run.

For the original lung cores, restore the original scientific package with
`recomb2027-atlas/reassemble.py`, fetch `LUNG` with `fetch_inputs.py`, and run
`prepare_lung_context.py --raw work/raw/LUNG --old /path/to/restored/recomb2027-atlas --output work/new-analysis/LUNG`.
Then use the same Proseg runner and evaluator on `roi_01` and `roi_02`.

## Manual-reference assessment

The SPATCH release was accessible and its Xenium annotation files were inspected.
Its publication identifies the annotations as nuclear boundaries. They therefore
cannot establish whole-cell transcript-assignment truth. No empirical calibration
or whole-cell accuracy result is claimed from those labels.
`manual_reference_eligibility.json` records the source and decision.

## Files and licensing

- `protocol.json`, `lung_context_protocol.json`: fixed selection and comparison rules.
- `results/`: donor, region, budget and Proseg summaries.
- `artifacts/`: compact inputs, spectra, witnesses and core assignments, with hashes.
- `sources.json`, `environment.json`: acquisition and software provenance.
- `REPORT.md`: interpretation and remaining limits.

Code follows the repository's MIT license. Source data and derived crops retain
10x Genomics attribution under CC BY 4.0. This analysis adds biological donor
replication and diagnoses crop-context sensitivity. It does not change the
fixed-seam geometry or provide expert whole-cell ground truth.
