# Complete reproducibility archive

The browsable files contain the merged source. The complete package also includes source crops, prepared atlas inputs, exact spectra, all interval rows, saved boundary witnesses, Proseg outputs, the manuscript and supplement.

Download or clone this branch with every file in `atlas2027/artifact_parts/`, then run from the repository root:

```sh
python atlas2027/reassemble.py --output GeoCert_Atlas_Merged_0.3.0a2.zip --extract atlas-replay
cd atlas-replay/GeoCert_Atlas_Merged_0.3.0a2
python atlas_cli.py verify-manifest
python atlas_cli.py validate
python atlas_cli.py reproduce --dataset lung --output-dir run/lung
python atlas_cli.py reproduce --dataset breast --output-dir run/breast
```

Python 3.12 or later is required. Reassembly uses the Python standard library and verifies each part and the assembled archive against `artifact_parts/manifest.json`. Extraction accepts only regular files and directories within the declared package root and requires an absent or empty destination. Omitting `--extract` produces just the verified ZIP.

The exact prepared-input workflow uses the standard library. For source-geometry analysis and the saved width-study checks, install `requirements-analysis.txt`. CI repeats the package checks, algorithm tests, every baseline spectrum and interval, and the saved width-study mapping and witness checks under Python 3.12 and 3.13. CI does not rerun stochastic Proseg inference or the historical timing experiment. See the package README and validation status for the scope of each experiment.

Scientific source revision: `1cdbb9e0b4ea386c210fd231372900b2844dddba`. The publication commit adds the complete archive, manuscript and automation to that source revision.
