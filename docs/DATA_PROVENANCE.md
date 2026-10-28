# Native data and conditioning

The inputs derive from two development crops of the public 10x Genomics dataset **Preview Data: FFPE Human Lung Cancer with Xenium Multimodal Cell Segmentation**, 5 February 2024, one donor, published under CC BY 4.0:
https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard

The selected 42 corridors contain 464 repeated RNA records / 387 unique transcripts and 51 source cells; the retained source crops have 7,370 RNA records. There are eight three-cell and 34 two-cell corridors. These repeated queries are not independent biological samples. The archived full artifact contains input provenance, extraction, exclusion and replay records. The current clean workflow uses unchanged source-derived rational JSON, without downloading new tissue or reingesting raw Parquet. Fixed nuclear regions retain universal ownership semantics. Corridor selection, candidate cells, endpoints, slopes, bands, gaps and clearance are conditioning assumptions. The largest geometric family is a declared stress setting, not calibrated anatomy.
