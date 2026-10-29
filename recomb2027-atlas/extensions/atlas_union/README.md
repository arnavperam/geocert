# Finite union of certified atlases

This extension evaluates the primary atlas and the four fixed-width atlases as five alternative complete deformation families. All five use the same reference partition, retained transcripts and physical-area units. Each attaining boundary selects one whole member; patches from different members are never combined.

For count `q`, the union spectrum is the smallest member cost at `q`. An unattainable count remains absent from the spectrum. Endpoint intervals summarize these feasible counts and may contain unattainable interior integers.

## Run

First extract the pinned scientific package and restore its compressed records using the repository's reproduction instructions. Keep the output directory separate from that immutable source directory.

```sh
python recomb2027-atlas/extensions/atlas_union/test_union.py
python recomb2027-atlas/extensions/atlas_union/union.py \
  --data-root /path/to/extracted/recomb2027-atlas \
  --output-dir /path/to/union-results
python recomb2027-atlas/extensions/atlas_union/independent_verify.py \
  --data-root /path/to/extracted/recomb2027-atlas \
  --union-output /path/to/union-results \
  --output-dir /path/to/independent-results
```

The extension uses Python's standard library and the companion's existing exact source-binding routines. `frozen_inputs.json` pins the 76 consumed scientific files to the published 1.1.0 source. It rejects altered or missing records before calculating any interval. This CLI evaluates those fixed component families; newly generated atlases require their own certification and an explicitly revised input contract.

Before union composition, it checks source-file hashes, physical calibration, atlas identities, complete observation records, source counts, complete component question censuses, every global component composition and attaining trace, and every member interval on the common 17-budget grid. The local exact spectra and geometry audits are the previously certified, checksum-bound inputs. The separate verifier reconstructs union results without importing the extension implementation.

## Outputs

| File | Contents |
|---|---|
| `union_spectra.json.gz` | Exact rational area for every attainable count, with member identity and a complete component trace |
| `intervals.csv.gz` | All 17 budgets, endpoint costs and member identities, primary-atlas endpoints and union marginal bounds |
| `budget_summary.csv` | Invariance, interval width, tightening and changes from the primary atlas |
| `reference_counts.csv.gz` | Exact membership of each reference count on the unchanged comparison keys |
| `reference_summary.csv` | Finite-budget and unbounded-budget reference-count reach |
| `input_manifest.json` | Input checksums, budget grid, member settings and union model identities |
| `summary.json` | Complete census, input-binding checks and output checksums |

`baseline_questions` retains the original positive-or-editable primary census. `union_questions` includes every source-positive or editable-in-any-member question. An absent member row becomes a fixed source count only after the complete component census proves that no relevant editable transcript exists there. Missing expected spectra are errors.

The marginal comparison takes the smallest lower bound and largest upper bound among the five member-wise independent-molecule relaxations. It therefore refers to the same five-family union and the same budget. No budget is allocated across different member atlases.

Reference reach uses the existing common comparison keys and missing-transcript exclusions. It tests membership in the full count spectrum, rather than containment in the endpoint interval. Additional countwise reach does not establish simultaneous realization of an entire reference assignment matrix.

## Inspect an attaining boundary

```sh
python recomb2027-atlas/extensions/atlas_union/witness.py \
  --data-root /path/to/extracted/recomb2027-atlas \
  --results-dir /path/to/union-results \
  --dataset lung --roi roi_01 --cell 51003 --gene CAVIN2 --count 0 \
  --output /path/to/witness.json
```

This example removes the source CAVIN2 count of one in lung cell 51003 using the 0.5 µm support member; the primary atlas cannot attain count zero. Other recorded keys and attainable counts are listed in `union_spectra.json.gz`. The command verifies the minimum across all five members, loads the selected component trace, checks its continuous patch boundaries and area, and writes the attaining boundary with its member atlas. All other patches in that member remain at their source boundary. Different query endpoints may select different members.

## Recorded results

The union contains 3,766 lung and 25,446 breast questions and produces 496,604 exact interval rows on the common 17-budget grid. On the original primary question census at area 4 µm², invariant counts are 2,524/3,304 (76.4%) in lung and 15,143/22,976 (65.9%) in breast. The complete union census gives 2,905/3,766 (77.1%) and 17,160/25,446 (67.4%); the denominators differ because additional source-zero questions can become editable.

| Unbounded-budget comparison | Primary attainable | Union attainable | Common changed-count denominator |
|---|---:|---:|---:|
| Lung Proseg | 185 | 315 | 2,798 |
| Lung nuclear expansion | 184 | 299 | 1,114 |
| Breast nuclear expansion | 1,773 | 2,942 | 10,965 |

These comparisons retain all original exclusions and use exact count membership. The union includes alternative patch selections as well as differing displacement ranges. It preserves fixed tri-junctions and the other constraints of each component family.
