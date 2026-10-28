# Lung-only ribbon-width and reference-count reach study

The Atlas count-chain and disjoint-gluing construction are the basis of this
study. Existing Proseg 3.2.0 outputs from both lung crops are integrated as an
additional assignment reference. No new Proseg fit or donor was added here.
Neither Proseg nor nuclear expansion is treated as truth.

## Frozen scope and populations

Before the new width runs, `REVISION_REACH_PROTOCOL.json` fixed caps B=4,12,28
chart units, fallback widths 4,2,1,1/2, the original geometry-only candidate
ordering and a 16-budget grid between 0 and 4 square micrometres. The two supplied
lung crops are the entire width-study population. Breast results are not
recomputed at alternative widths. Every rebuilt atlas repeats the rational
pixel-containment, nuclear-obligation and patch-disjointness checks.

A physical chart moves a point by p|h|/sqrt(2). The displacement band is
|h|<=3B/4, with p=1048576/4934475 micrometres. The caps therefore correspond to
maximum point-displacement bounds 0.4507806, 1.3523418 and 3.1554642 micrometres.
Each band cap occurs in at least one accepted patch; most patches use narrower
fallback widths. At B=12, 132 of 544 patches use the top width; at B=28, only 41 of 543 do. The other accepted patches use B<=4. These are band-implied displacement bounds, not measured anatomical errors or solved attainable maximum movements after nuclear and slope constraints.

The atlases are rebuilt greedily at each cap. They are not nested uncertainty
sets. Wider proposals may change orientation/patch selection and reduce reachable
count changes; monotonicity across these three models is not guaranteed.

Proseg comparisons use exactly its 7,256 retained transcript IDs: 5,774 and 1,482
in the two crops. Its 114 exclusions are retained explicitly. Gene, source-mask
identity and observed coordinate joins are checked. Proseg's diffused coordinates
are ignored; the original observed positions define the geometric model. Proseg
cell indices are mapped through original_cell_id and the original mask label.

Nuclear comparisons use the existing five-micrometre nearest-nucleus rule and
five-micrometre crop-edge guard. Each reference therefore has its own fixed
shared molecule population. In both cases excluded molecules receive zero count
weight but remain as free collateral RNA in the original geometric constraints.
The compared source and alternative counts always use the same molecule IDs.

## Complete results

All counts below aggregate the same two lung crops. Boundary coverage uses all
23,164 distinct reference edges incident to complete reference cells, including
frozen interfaces; the narrower eligible-interface denominator is also retained
in per-crop records.

| Largest accepted band-implied displacement cap (µm) | Accepted patches / 782 candidates | Editable RNA | Editable boundary | Nuclear changed counts attainable | Proseg changed counts attainable |
|---:|---:|---:|---:|---:|---:|
| 0.4508 | 544 | 745 | 15,821/23,164 (68.300%) | 186/1,134 (16.402%) | 185/2,798 (6.612%) |
| 1.3523 | 544 | 1,288 | 15,730/23,164 (67.907%) | 223/1,134 (19.665%) | 241/2,798 (8.613%) |
| 3.1555 | 543 | 1,220 | 15,765/23,164 (68.058%) | 190/1,134 (16.755%) | 200/2,798 (7.148%) |

“Attainable” means the alternative count belongs to the exact cell/gene score
spectrum with no additional finite area cap. It is stricter than merely lying
between the two marginal endpoints. It does not mean that all alternative counts
can coexist in one partition. The stored 16-budget tables further constrain the
same spectra by one physical-area budget. Reach, non-singleton frequencies and
marginal-width-at-least-two denominators are retained for every budget and crop.

The 1.35-µm setting reaches more reference count changes, but the 3.16-µm setting
falls back. Boundary coverage stays near 68%. No setting was chosen as a new
primary model on the basis of these outcomes. Increasing a width cap alone does
not establish a realistic uncertainty family; most reference count changes remain
outside every tested atlas, even with unlimited area.

## Why the nuclear denominator is now 1,134, not 1,128

The earlier nuclear analysis dropped a cell/gene query whenever any associated
editable RNA failed its crop-edge guard. The revised calculation scores only
shared eligible IDs and retains excluded RNA as unscored geometric collateral.
It can therefore include those queries exactly. Six changed-count questions are
added; five are unattainable and one (roi_01, cell 153965, C7) is attainable.
All shared old/new query source counts, alternative counts and attainability
classifications agree. The old 185/1,128 becomes 186/1,134 for this explicitly
expanded query population, not by changing a count target or a geometric setting.
The six records are listed in `nuclear_denominator_audit.json`.

## Verification and reproducibility

Run from the Atlas package root:

```sh
python analysis/revision_reach.py --cap 4
python analysis/revision_reach.py --cap 12
python analysis/revision_reach.py --cap 28
python analysis/revision_reach.py --summarize
python analysis/revision_reach_verify.py
```

Each run retains generated atlas definitions, independent geometry-validation
records, full reference-count spectra, every finite-budget row, and local and
composed witnesses. The verification replay directly checks 8,270 local witness
counts/geometry and 25,610 global count/area traces; all pass. Computation uses
exact rational count-chain optimization for multicopy queries and source-feasible
closure for singleton queries. All RNA, including zero-scored collateral, remains
in the geometric problem. A timeout would terminate the run rather than produce
an exact-looking unresolved result.

The reused Proseg execution itself is documented in `PROSEG_COMPARISON.md`.
`execution_record.json` discloses the package-relative path amendment and the
all-incident-edge plotting denominator. Initial execution summaries and hashes
are preserved under `results/revision_reach/initial_execution/`. The intermediate
baseline protocol bytes were not retained; its historical hash is preserved and
that run is not described as a byte-identical source-frozen execution. A final
replay binds all three cap runs to the final point-input contract and one current
protocol, comparing every saved numerical spectrum, CSV and local witness with
the initial runs. Scientific width/budget/population
settings were unchanged after the new outcomes were observed. This is a
post-review sensitivity analysis, not an independent preregistered experiment.

The figures `reach_width.pdf` and `reach_budget.pdf` show the complete tested
width settings and budget grid. They are descriptive model-reach figures,
not accuracy, confidence, calibration, or biological validation.
