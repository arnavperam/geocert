# Data provenance and use

## Lung development crops

The supplied mask and transcript crops derive from 10x Genomics, **Preview Data: FFPE Human Lung Cancer with Xenium Multimodal Cell Segmentation**, released 5 February 2024, one donor, CC BY 4.0:

https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard

The two original cropped mask/Parquet/metadata sets and the previously prepared count inputs are retained. This continuation rebuilds boundary patches from the original pixel-edge partition. It does not claim a new donor, a corrected segmentation or true ownership.

## External breast format test

Source: 10x Genomics **Example Xenium Datasets, Xenium Onboard Analysis 2.0**, two-FOV human breast, multimodal segmentation, CC BY 4.0.

https://www.10xgenomics.com/support/software/xenium-onboard-analysis/2.0/resources/xenium-example-data

Original archive:

https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip

Original archive SHA-256: `cc1e987b06aa748a6b24d3d6f51fc0d6765daa4836c483f14ec4f1bd18b1779b`.

The vendor says these trimmed data demonstrate output formats and are not intended for biological conclusions. That restriction on interpretation is retained. The two cropped sources were acquired in an existing GitHub workflow; this continuation recovered its output artifact, independently re-ingested CSV/gzip records and applied the current fixed atlas rules. Prior workflow analysis exists, so this is not described as pristine, prospectively preregistered biological validation.

Artifact workflow: `arnavperam/geocert`, run `40630851940`, commit `d677bd96768d665983b095aaa63e0117fb943c89`, artifact `12791482673`.

Artifact SHA-256: `f53bee7d82a385a320c14a180049e66bf6f60963dda8d7a75574fda90c28febf`.

Modifications consist of previously defined 200-µm crops, Q20 measured-gene filtering, count/coordinate extraction, geometry-only patch construction, exact spectra and analysis outputs. Original coordinates and source mask labels are retained. These data rights are separate from the as-yet-unselected software license.
