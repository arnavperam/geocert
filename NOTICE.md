# Rights and attribution

The Atlas software in `recomb2027-atlas/` is distributed under the MIT License in [recomb2027-atlas/LICENSE](recomb2027-atlas/LICENSE). Its source-derived data and Proseg execution records have separate attribution described in the [Atlas notice](recomb2027-atlas/NOTICE.md) and [data provenance](recomb2027-atlas/docs/DATA_PROVENANCE.md).

## Earlier corridor implementation

The root-level `core/`, `geocert.py` and accompanying corridor artifacts are the earlier v1.0.0 implementation. No software license was selected for that release; this update leaves those terms unchanged. Contact the repository owner concerning reuse of that earlier implementation.

The source-derived native data originate from the 10x Genomics dataset identified in `docs/DATA_PROVENANCE.md`, released under CC BY 4.0. That attribution and license apply to the source data; no ownership of third-party data is claimed here.

The implementation uses Python's standard library and does not vendor NumPy, SciPy, NetworkX, font files, or a third-party surface-segmentation package. Algorithmic foundations and their citations are documented in the manuscript and mathematical notes.
