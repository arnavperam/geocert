# Validation status

The scientific implementation, prepared inputs and results are frozen at [6581f628](https://github.com/arnavperam/geocert/tree/6581f628a5ac8fb83629f9e0dee0ffa4d685375d/recomb2027-atlas). Subsequent paper and documentation edits preserve those files.

The completed CI runs on document commit `d804f70e0296d13ed694457ab92d067b16ad8b24` validate that scientific source:

- [Atlas exact replay, Python 3.12 and 3.13](https://github.com/arnavperam/geocert/actions/runs/40779065551): passed.
- [Legacy corridor checks, Python 3.12 and 3.13](https://github.com/arnavperam/geocert/actions/runs/40779065656): passed.

These links identify completed runs, rather than the status of a later documentation commit. The Atlas workflow verifies archive parts and file identities, matches 71 exposed scientific source/environment files to the archive, restores records, runs input and algorithm checks, recomputes both primary collections and the lung support-1 family, checks all saved width compositions and interval rows, and audits Proseg transcript and count-matrix identities.

| Check or evidence | Scope |
|---|---|
| Continuous assignment-cost function | Scalar proof in the main paper; coupled proof in the supplement |
| Compatible partition | Explicit homeomorphism, gluing argument and rational chart-area conversion |
| Hard-budget count spectra | Scalar point-transcript count chain and composition |
| Primary count replay | 26,280 spectra and 210,240 interval rows recomputed |
| Four required physical widths per dataset | Eight recorded families and 1,807,967 interval rows |
| Width-study geometry | Accepted ribbons, side labels, nuclei, memberships and disjoint interiors checked |
| Width-study exact output | Every local witness checked during computation; 106,351 global spectra, 131,376 traces and every interval row rechecked |
| Independent scalar check | Deterministic sample of 246 local spectra and 1,628 assignments; zero failures or resource limits |
| Input binding | All 20 prepared atlases accepted; 20 inconsistent bindings and seven unsupported free-input cases rejected |
| Proseg outputs | One adapted execution per lung crop; transcript joins, coordinates and both sparse count matrices checked |
| Archive availability | Public archive parts downloaded without authentication and matched to the recorded SHA-256 |

## Evidence limits

Widths and budgets are sensitivity settings without anatomical calibration. The study has no independent transcript-ownership truth, additional biological donor or validated downstream inference. Breast material is a vendor format-test collection. Proseg replication and ownership accuracy were not assessed. Three-label atlas ribbons, moving junctions, topology changes and whole-slide deployment lie outside the completed experiments.

The full eight-width optimization study was run once in the recorded environment. CI replays both primary collections and one width family, then checks every saved width composition and interval row. It does not rerun the entire optimization study or independently enumerate every local optimum. These computational checks are distinct from specialist review of the mathematical proofs.

Historical development and review notes are retained in [history/](history/) and [legacy/](legacy/). Their release instructions and status statements describe the versions at which they were written.
