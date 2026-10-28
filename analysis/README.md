# Assignment-sensitivity analyses

These scripts add cell-count readouts and a fixed-rule transfer check to the unchanged `v1.0.0` exact core. They do not estimate corrected segmentation or evaluate annotation accuracy.

## Development cell-count audit

Run from the repository root with Python 3.12 or 3.13, without `-O`:

```sh
python analysis/readout_inputs.py
python analysis/audit_cell_counts.py --output run/readout_results.json
python analysis/check_readouts.py --input run/readout_results.json
python analysis/audit_shared_cost.py --output run/shared_cost.json
```

These commands use only the standard library. `source_counts.json` records the retained source-mask counts and crop-edge exclusions. The query universe includes every gene represented in a corridor and each source cell that passes the crop-edge rule. RNA outside the corridor keeps its source assignment. There are 750 readouts per geometric setting: 72 have nonconstant count extrema under the baseline and 246 under the enlarged family. This is a new analysis of the existing development sample, not independent biological validation.

The shared-cost script repeats all 252 optional questions with identical feasibility and favorable tie handling for the additive model. Expected strict losses: 1/126 baseline and 16/126 enlarged-family. It is an energy ablation, not a comparison with a segmentation method.

## Fixed-rule breast transfer check

The selection and analysis rules are recorded in [`evaluation/PROTOCOL.md`](../evaluation/PROTOCOL.md), committed before the dataset was extracted. The 10x Genomics dataset is a trimmed two-field format-test example; it is explicitly unsuitable for drawing biological conclusions. This evaluation tests technical transfer to previously unused inputs.

To reconstruct the input crops:

```sh
python -m pip install numpy==2.2.6 pandas==2.2.3 pyarrow==20.0.0 zarr==2.18.7 numcodecs==0.15.1
curl --fail --location --output /tmp/geocert-breast-2fov.zip https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip
python evaluation/inspect_external.py
python evaluation/run_external.py
mkdir -p analysis/external_inputs results/external
cp -r run/external/breast_01 run/external/breast_02 analysis/external_inputs/
```

The download is 379 MB. Its published MD5 and computed SHA-256 are checked before extraction. Two crop centers are selected using only cell-centroid quantiles; no failed or uninteresting crop is replaced.

```sh
python analysis/screen_external.py --height 24 --cells 3
python analysis/screen_external.py --height 12 --cells 3
python analysis/screen_external.py --height 24 --cells 2
python analysis/screen_external.py --height 12 --cells 2
python analysis/analyze_external.py
python analysis/check_external.py
```

The screen retains 29 corridors. Expected results are 174 completed optional questions, with 2/87 baseline and 6/87 enlarged-family strict additive losses. All 166 null comparisons remain. There are 716 cell–gene readouts per setting; 108 baseline and 197 enlarged-family readouts have differing extrema. The saved analysis has no resource-limit failures.

`results/external/` retains all eligibility counts, all eligible and selected charts, optional answers, count extrema, and witness checks. Repeated records across settings are not independent samples. Parameters are sensitivity settings, not calibrated anatomical uncertainty.

## Verification scope

The witness check recomputes original-coordinate constraints and cell counts. A fixed hash subset is solved by the exact height reference. These paths share model definitions; they are not independent theorem verification. The GitHub workflow reruns extraction and every selected query in a separate environment. The original core release and its checksum manifest are unchanged.

The vendor source is CC BY 4.0. See `docs/DATA_PROVENANCE.md` and the dataset link in the protocol. Software licensing remains an author decision; public visibility is not an open-source license.
