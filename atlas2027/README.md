# GeoCert: expression intervals on a compatible boundary atlas

Research implementation accompanying **Certified Expression Intervals under Shared Segmentation Uncertainty**, by Arnav Peram.

The input is a reference cell segmentation, measured RNA coordinates, fixed nuclear regions, and a declared geometric uncertainty family. The output is one count interval per reference cell, gene and physical-area budget. All editable patches belong to one disjoint atlas, share one budget and glue into a topology-preserving partition. Unselected boundaries and patch seams stay fixed.

This release adds a critical-height/count dynamic program for exact hard-budget spectra with explicitly listed point RNAs. It retains the earlier exact continuous assignment-cost reduction. It does not infer which segmentation is anatomically correct.

## Run

Python 3.12 or later is required by the included code. The inherited clean replay used Python 3.13.5; this merged candidate is checked again with Python 3.12 and CI. The prepared-input exact workflow uses only the Python standard library.

```sh
python3 -m venv .venv
. .venv/bin/activate
python atlas_cli.py verify-manifest
python atlas_cli.py showcase --output run/showcase.json
python atlas_cli.py validate
```

On Windows, substitute `.venv\Scripts\python.exe` for `python`. Keep assertions enabled; do not pass `-O`.

The showcase checks a recorded exact spectrum and its continuous witness. In breast format-test crop 1, cell 81 has an LDHB source count of 22. Its exact interval is `[19,23]` at area 0.25 µm² and `[17,23]` at area 0.5 µm². Count 17 has minimum area approximately 0.484040265 µm². The witness transfers five LDHB RNAs through two compatible patches and one collateral KRT8 RNA. The test dataset does not support a biological interpretation of that example.

## Recompute the full count tables

```sh
python atlas_cli.py reproduce --dataset lung --output-dir run/lung
python atlas_cli.py reproduce --dataset breast --output-dir run/breast
```

The commands compute all multicopy spectra with the exact count chain and use implication closure for single-copy targets. They compare every rational minimum area and every interval row with the stored independent-reference-checked results. A fresh output directory preserves the stored records.

Expected outputs:

| Input collection | Complete reference cells | Local patches | Cell–gene spectra | Eight-budget interval rows |
|---|---:|---:|---:|---:|
| Lung crops | 153 | 544 | 3,304 | 26,432 |
| Breast format-test crops | 467 | 1,708 | 22,976 | 183,808 |

Five eligible lung cells have no retained positive or editable gene entry. Their entries are implicit zero within the gene universe recorded in `matrix_index.json`. All counts refer to detected retained molecules, not biological abundance or unmeasured genes.

For the independent scalar/Cartesian check:

```sh
python analysis/check_global_counts.py
GEOCERT_DATA_ROOT=external python analysis/check_global_counts.py
```

These scripts regenerate their small validation summaries; they do not alter spectra or inputs. A copied working directory is appropriate for replay. The scalar reference does not import the production compiler, geometry envelopes, count chain or optimizer. It shares prepared input extraction.

## Geometry and source-data checks

The exact count solver does not require scientific libraries. Source extraction, image geometry indexing, nuclear-distance analysis and plotting do:

```sh
python -m pip install -r requirements-analysis.txt
python analysis/check_source_inputs.py
python analysis/verify_atlas.py
GEOCERT_DATA_ROOT=external python analysis/verify_atlas.py
```

`analysis/build_atlas.py` regenerates the lung atlas from the supplied mask and transcript inputs. For the external collection, set `GEOCERT_DATA_ROOT=external`. The extraction rule is geometry-driven, with the protocol and coverage denominators in `docs/` and `results/`.

## Files

- `atlas/`: physical charts and gluing maps, exact count-chain algorithm, independent scalar reference.
- `geocert_expression/`: preserved exact continuum primitives and earlier count APIs.
- `data/`, `external/data/`: original crops, prepared observations and rational atlas definitions.
- `results/`, `external/results/`: complete spectra, witnesses, interval tables and baseline comparisons.
- `tests/`: exhaustive, nonflat, unit-molecule gap and topology-map checks.
- `analysis/`: full evaluation and figure commands.
- `manuscript/`: LaTeX, figures, paper and supplement.

The repeated 12-process scaling experiment is in `results/chain_benchmarks/`. Run `python analysis/run_chain_benchmarks.py` in a working copy to repeat it. The special equal-tent family has an analytical spectrum; these measurements are not a comparison against every specialized method.

## Interpretation

The complete-cell count table uses a restricted *fixed-seam* uncertainty family, not unrestricted deformation of every cell perimeter. Roughly two-thirds of incident boundary length is editable in these data. The remaining portions, junctions, crop-border cells and nuclear owners are fixed. Exported atlas identifiers specify this conditioning. Marginal endpoints of different matrix entries need not be jointly attainable.

The breast data are the vendor's trimmed format-test data and are used only for computational transfer testing. Neither collection supplies true transcript ownership. A deterministic five-micrometre nuclear-expansion reference is implemented and often lies outside the atlas; it is not Proseg or a biological accuracy standard. Proseg 3.2.0 was executed on the two lung crops using the documented no-factorization adaptation. Raw outputs and exclusions are retained in results/proseg; reference reach uses shared retained molecules and original observed coordinates, not Proseg-diffused positions. Baysor and ResolVI were not executed.

The merged candidate is `0.3.0a2`, on the public `atlas-merged-revision` review branch under `atlas2027/`; it is not a stable release tag. Public v1.0.0 is an earlier exact-workflow release, not this extension. See `docs/PROVENANCE.json` for input identities, `docs/REPORT.md` for results and `docs/RESPONSE.md` for the critique response.

## Rights and disclosure

Source-data attribution and permitted use are documented in `docs/DATA_PROVENANCE.md`. Project software is MIT licensed; `NOTICE.md` preserves separate dataset and dependency rights.

AI assistance contributed to mathematics, algorithms, code, experiments, verification and writing. The listed author is responsible for the submission and for reviewing and defending its claims. This package is not a conference acceptance or biological-validation claim.

## Comparison and merge

The Atlas is the scientific base because its disjoint supports give one compatible fixed-seam count family and its scalar count-chain algorithm handles explicit unit molecules in polynomial time. The preceding corridor revision contributes actual Proseg outputs, MIT licensing and release checks. Its old pair-census totals are not a new discovery; the original supplement already contained them. See docs/COMPARISON_AUDIT.md.

Both denominators matter: joint geometry changes 203/1,026 breast questions with marginal width at least two, but 203/22,976 questions overall. Constant numbers of non-singleton intervals do not imply fixed endpoints: on the supplied grid, lung endpoints plateau from 0.25 square micrometres, breast endpoints from 1. The breast data remain explicitly vendor format-test data.

`python analysis/audit_revision_csv.py` recomputes these denominator and saturation checks. New width/reach results are in results/revision_reach; their protocol was fixed before execution. The study uses lung crops only, and external-method reach is not accuracy or anatomical calibration.

## Complete public artifact

The public branch stores the complete archive in checked parts; see `ARCHIVE.md` beside the repository source and run `python atlas2027/reassemble.py --extract atlas-replay` from the repository root. This extracted package already includes every prepared input and result. The scientific source revision is `1cdbb9e0b4ea386c210fd231372900b2844dddba`. Main and supplement sources compile from `manuscript/`.
