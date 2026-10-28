# Conditional families and safe aggregation

Every numeric endpoint is conditional on a declared uncertainty family. A cell ID,
gene, and budget do not identify that family. The same cell can have different
exact endpoints when a different boundary segment is permitted to move.

## Corrected public output

The `certify` and `demo` commands use schema `geocert-expression-0.3`. Each interval
includes its full scope, a stable family ID derived from the supplied chart and
geometry settings, a query identifying its region/cell/gene/budget/metric, and a
model/score/offset binding. Every corridor output explicitly sets:

- `kind: single_corridor_conditional`
- `exterior_policy: all_observations_outside_corridor_fixed`
- `whole_cell_boundary_uncertainty: false`
- `simultaneous_corridor_deformation: false`

A family ID changes if the chart geometry, setting, ownership mode, or localization
radius changes. Direct Python spectra and enclosures carry the physical model key,
score key, and fixed offset. A direct score spectrum is scoped to the supplied
geometric model; it does not gain whole-cell meaning by adding a fixed exterior
count. References to exactness refer to optimization within that family.

## Reproduced overlap audit

`python validation/audit_scope.py` reindexes the existing results without solving
the native models again. It independently reproduces the reviewer's counts:

| Quantity | Result |
|---|---:|
| Original interval rows | 14,384 |
| Cell/gene/setting/budget keys | 11,392 |
| Keys with different endpoints across corridors | 178 |
| Baseline / maximal differences | 74 / 104 |
| Region-qualified cells | 51 |
| Cells appearing in several corridors | 19 |
| Maximum corridors per cell | 7 |

For cell 153131, CD8A, budget 1 square micrometre, corridor 02d9a53ef82968
returns [2,2] and corridor 4c4f466415dd22 returns [2,3]. Both are conditional
answers to different questions. They must not be presented as a single invariant
whole-cell certificate.

The original numeric results are preserved for provenance. Use
`expression_intervals_scoped.csv` and `expression_matrix_intervals_scoped.csv` as
the explicit replacements. `corridor_scope_disagreements.csv` lists every one of
the 178 keys and its contributing corridors. `corridor_scope_audit.json` binds the
source files by SHA256.

## A safe additional question: one corridor may vary

The new `union_intervals` function and `union` CLI accept fully scoped intervals
with a shared explicit `reference_id`. They require equal source counts and equal
region, cell, gene, budget, and budget metric, and reject duplicate family IDs,
unresolved inputs, reversed bounds, or mixed references.

For a declared union of families, the exact lower endpoint is the minimum of the
member lower endpoints and the exact upper endpoint is the maximum of the member
upper endpoints, when every member result is exact. Each endpoint keeps the
identity of the family attaining it. This is a set-union identity, not a new
optimization theorem. For the CD8A example, the union interval is [2,3]. An
intersection would wrongly exclude an admitted assignment.

`corridor_union_intervals.csv` provides one answer for each cell/gene/setting/budget
under the union of all listed corridor families for that cell. The 16,953 stored
rows are sparse over an explicitly indexed 273,360-row domain. The index lists
member corridors and the measured gene universe. Its exact-zero rule applies
only within that declared domain; unknown genes and unmeasured transcripts remain
unresolved. Each endpoint names its attaining corridor for lookup in the original
local frontier. Source geometry witnesses fixed genes when no local frontier is
needed.

This union permits a choice of one corridor family at a time. It does not allow
several corridors to move simultaneously and does not cover all possible
whole-cell perturbations. Endpoint witnesses for different genes remain separate
optimizations and need not coexist.

## Geometric composition is blocked until gluing can be checked

The previous `compose_spectra` function accepted arbitrary local spectra and put
its disjointness assumption only in a docstring. It now raises an error. Local
cost spectra do not contain enough information to verify disjoint edit support,
common boundary traces, additive area, and disjoint transcript accounting.

`convolve_score_tables` preserves the algebraic min-plus utility and returns an
unbound score table. Its output is not a geometric certificate. The composition
unit test uses analytically independent toy modules; this does not establish
independence of the native corridors.

A future geometric composition interface must verify the edit supports and
boundary compatibility in a common coordinate frame before returning a global
certificate. A caller-supplied Boolean flag would not establish those conditions.

## Verification

`python tests/test_scope.py` checks exact union endpoints against 30 explicit
finite sets, rejects 22 unsafe operations, exercises the native CLI, verifies
family identities change with geometry, and checks equivalent rational radii
share an identity. Existing budget, resource, and composition tests exercise the
unchanged numerical algorithms. No native census or 72-crop external run is
needed for this correction.
