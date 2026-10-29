# GeoCert Atlas

Exact cell–gene count intervals under a shared segmentation-deformation budget.

Code and data for **Certified Expression Intervals under Shared Segmentation Uncertainty**, by Arnav Peram.

[Paper](recomb2027-atlas/manuscript/main.pdf) · [Supplement](recomb2027-atlas/manuscript/supplement.pdf) · [Reproduction instructions](recomb2027-atlas/README.md) · [Scientific source cited in the paper](https://github.com/arnavperam/geocert/tree/6581f628a5ac8fb83629f9e0dee0ffa4d685375d/recomb2027-atlas)

GeoCert computes the smallest and largest observed transcript counts attainable within a declared boundary family and total deformation budget. Exact rational costs and attaining boundaries accompany the endpoints. The family holds patch seams, junction neighborhoods, uncovered boundaries and nuclear owners fixed; the intervals describe sensitivity within that family rather than anatomical ground truth.

## Run an example

Use Python 3.12 or 3.13. The prepared-input solver needs only the standard library. Allow about 1.2 GB for the restored records.

```sh
git clone https://github.com/arnavperam/geocert.git
cd geocert
python3 recomb2027-atlas/reassemble.py --extract atlas-replay
cd atlas-replay/recomb2027-atlas
python3 atlas_cli.py verify-manifest
python3 release/restore_results.py
python3 atlas_cli.py showcase --output run/showcase.json
```

The example returns `[19, 23]` at 0.25 µm² and `[17, 23]` at 0.5 µm² for LDHB in cell 81 of breast format-test crop 1, together with a boundary attaining count 17. Keep Python assertions enabled. The [archive notes](recomb2027-atlas/ARCHIVE.md) describe checksum verification and extraction.

## Compare complete atlas families

The finite-family extension evaluates the primary atlas and four wider-support alternatives. Each endpoint comes from one complete member atlas and retains its boundary witness. The [extension instructions](recomb2027-atlas/extensions/atlas_union/README.md) reproduce all 496,604 interval rows and run a separate verifier.

On the same comparison keys, the union attains 315 of 2,798 Proseg count changes, compared with 185 for the primary atlas. At 4 µm², 76.4% of the original lung count questions and 65.9% of the breast questions remain invariant under the expanded family. Both the primary and union results are reported in the paper.

## Reproduce the paper

[Atlas 1.1.0](recomb2027-atlas/) contains the implementation, tests, protocols and checked input/results archive. Its [README](recomb2027-atlas/README.md) gives commands for the primary count tables, the four-width study, geometry checks and figures. [CI](https://github.com/arnavperam/geocert/actions/workflows/atlas-final.yml) checks the exact replay on Python 3.12 and 3.13.

The study uses lung crops and breast format-test crops. The reported comparison to Proseg measures count changes attainable within the model; it is not an ownership-accuracy benchmark. Dataset attribution and exclusions are recorded in [DATA_PROVENANCE.md](recomb2027-atlas/docs/DATA_PROVENANCE.md).

The earlier corridor implementation remains in the repository root and is preserved at [v1.0.0](https://github.com/arnavperam/geocert/tree/v1.0.0).

## Citation and license

Citation metadata is in [CITATION.cff](CITATION.cff). The software in `recomb2027-atlas/` is [MIT licensed](recomb2027-atlas/LICENSE). Data have separate attribution and reuse terms; see the [Atlas notice](recomb2027-atlas/NOTICE.md). The earlier root implementation retains its [original notice](NOTICE.md).
