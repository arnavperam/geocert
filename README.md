# GeoCert

Exact expression intervals under a shared segmentation-deformation budget.

Companion software for **Certified Expression Intervals under Shared Segmentation Uncertainty**, by Arnav Peram. A geometry-only atlas selects disjoint editable pieces of reference cell boundaries. Exact count spectra and min-plus composition allocate one common area budget across those pieces, producing one interval per retained cell–gene count. Uncovered interfaces, junctions and seams remain fixed.

## Quick start

The prepared-input exact workflow uses the Python standard library. Python 3.12 and 3.13 are supported; the release record names the environments actually tested.

```sh
git clone --branch atlas-finalization https://github.com/arnavperam/geocert.git
cd geocert
python3 -m venv .venv
. .venv/bin/activate
python fetch_inputs.py
python atlas_cli.py verify-manifest
python atlas_cli.py showcase --output run/showcase.json
python atlas_cli.py validate
```

On Windows use `.venv\Scripts\python.exe`. Assertions must remain enabled: do not run with `-O`.

The showcase uses the vendor breast format-test crop. A retained LDHB count of 22 has interval `[19,23]` at 0.25 square micrometres and `[17,23]` at 0.5. Its count-17 boundary moves five LDHB molecules through two compatible patches, plus one collateral KRT8 molecule, at minimum area approximately 0.484040265 square micrometres. The output contains rational costs and the attaining boundaries.

## Reproduce the paper

Run from the repository root. These commands preserve the stored baseline results by writing to fresh directories.

| Result | Command |
|---|---|
| Baseline lung intervals and exact spectra | `python atlas_cli.py reproduce --dataset lung --output-dir run/lung` |
| Baseline breast-test intervals and exact spectra | `python atlas_cli.py reproduce --dataset breast --output-dir run/breast` |
| Scalar and topology checks | `python atlas_cli.py validate` |
| Independent local-spectrum/global-count checks | `python analysis/check_global_counts.py` |
| Same check on breast test inputs | `GEOCERT_DATA_ROOT=external python analysis/check_global_counts.py` |
| Four-width, two-dataset study | `python analysis/run_reach_study.py --output-dir run/reach` |
| Width replay versus frozen scientific hashes | `python analysis/verify_reach_replay.py --input-dir run/reach` |
| Imported Proseg ID, coordinate and count checks | `python analysis/proseg_reach.py` |
| Tables from all recorded width outcomes | `python analysis/summarize_reach.py` |
| Native witness and scaling figures | `python analysis/make_figures.py` |
| Source-mask, reach and budget figures | `python analysis/make_final_figures.py` |

The last four commands read the supplied `reach_study/` records. In a source-only checkout, run the width study first with `--output-dir reach_study`. Source extraction, raster checking and figures require the analysis environment:

```sh
python -m pip install -r requirements-analysis.txt
python analysis/check_source_inputs.py
python analysis/verify_atlas.py
GEOCERT_DATA_ROOT=external python analysis/verify_atlas.py
```

The original eight-budget tables contain 26,432 lung and 183,808 breast-test rows. The wider-family study uses a separate frozen 17-budget grid. `PREREGISTRATION.md` specifies its width schedules, sample retention, budget grid, reference comparisons and failure handling. A wider schedule may select different patches, so its uncertainty family need not contain the narrower one.

## Algorithm and files

- `atlas/`: rational geometry, topology-preserving maps, critical-height/count dynamic programming, and an independent scalar reference.
- `geocert_expression/`: exact assignment-cost representation and source-clipped completion.
- `data/`, `external/data/`: masks, observed molecules, fixed nuclei, prepared atlases and attribution.
- `data/proseg/`: imported Proseg 3.2.0 execution records and raw outputs.
- `results/`, `external/results/`, `reach_study/`: complete spectra, witnesses, intervals, reported summaries and replay identities.
- `analysis/`, `tests/`: reproduction programs and correctness tests.
- `manuscript/`: LaTeX sources and figures.

The hard-budget count-chain algorithm is exact for a single boundary and explicitly represented point transcripts. The more general coupled-boundary theorem supports unary penalized objectives. The implementation does not extend the scalar count-chain complexity claim to arbitrary multi-boundary or higher-order queries.

## Interpretation

An interval is an extremum over a declared geometric family, not a confidence interval or corrected segmentation. Complete-cell retained counts include fixed contributions outside the editable atlas. Different interval endpoints need not occur together. The supplied source masks do not establish anatomical truth.

The breast inputs are the vendor's trimmed format-test dataset and support computational transfer tests only. The Proseg outputs are two previously executed, adapted runs on clipped lung crops, not new executions or a held-out accuracy benchmark. Excluded transcripts, cell mapping and background assignments remain explicit. Wider-ribbon reach is measured against common changed-count questions, with separate coverage denominators.

## License and citation

Software is under the [MIT license](LICENSE). Data retain their separate original terms in [DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md). Citation metadata is in `CITATION.cff`. Version 1.1.0 identifies this source release; no conference acceptance is implied.

The manuscript discloses AI assistance in mathematical development, algorithms, code, experiments, checking and writing. The author is responsible for the final scientific claims.
