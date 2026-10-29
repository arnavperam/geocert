# Validation status

The primary and fixed-width scientific implementation, prepared inputs and results remain frozen at [6581f628](https://github.com/arnavperam/geocert/tree/6581f628a5ac8fb83629f9e0dee0ffa4d685375d/recomb2027-atlas). The finite-union and contiguous-region extensions retain that implementation and its original results, with their additional code and records under [extensions/](../extensions/).

The completed CI runs on document commit `d804f70e0296d13ed694457ab92d067b16ad8b24` validate that scientific source:

- [Atlas exact replay, Python 3.12 and 3.13](https://github.com/arnavperam/geocert/actions/runs/40779065551): passed.
- [Legacy corridor checks, Python 3.12 and 3.13](https://github.com/arnavperam/geocert/actions/runs/40779065656): passed.

The finite-union source at [9fb63a2](https://github.com/arnavperam/geocert/tree/9fb63a22f42a191ffc54c415d5f53d526ee91182/recomb2027-atlas) also has a completed [Atlas exact replay on Python 3.12 and 3.13](https://github.com/arnavperam/geocert/actions/runs/40796267306): passed. That run additionally checks the union input contracts, computes the union and independently verifies its spectra, intervals and witnesses.

These links identify completed runs, rather than the status of a later commit. The Atlas workflow verifies archive parts and file identities, matches 71 exposed scientific source/environment files to the archive, restores records, runs input and algorithm checks, recomputes both primary collections and the lung support-1 family, checks all saved width compositions and interval rows, and audits Proseg transcript and count-matrix identities. The larger-region measurements and additional checks below are separately recorded; the linked union CI run does not claim to rerun them.

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
| Finite-union input contracts | 76 pinned scientific files; common source partition, observations, units and complete component question censuses checked |
| Finite-union exact output | 29,212 spectra, 42,123 finite states and 496,604 rows independently reconstructed; all 2,254,727 component rows checked |
| Finite-union witnesses | 39,191 explicit traces, 2,932 justified fixed-count completions and 18,762 distinct continuous local witnesses checked; five invalid mutations rejected |
| Larger-region source integrity | Full original breast archive and seven members checked; complete prepared 400 and 600 µm populations agree with exact source pixels, masks, transcript records and counts |
| Completed contiguous 400 µm run | Primary adaptive atlas, 1,363 cells, 5,102 ribbons, 69,278 spectra and 554,224 exact rows; no geometry failures or incomplete local queries |
| Completed contiguous 600 µm run | Primary adaptive atlas, 3,171 cells, 11,748 ribbons, 152,925 spectra and 1,223,400 exact rows; initial global deadline and completed retry preserved |
| Contiguous-region output audit | Complete expected local/cell census, zero-cost source states and independent min-plus composition; 95,253 and 208,555 state traces and every interval row checked for 400 and 600 µm |
| Contiguous-region scalar samples | 16 spectra per region match separate scalar enumeration: 92 assignments / 38 states at 400 µm and 152 assignments / 33 states at 600 µm; no errors or resource limits |

## Extension records and denominators

The [union verification report](../extensions/atlas_union/results/independent_verification.json) binds the recorded output files by checksum. Its [budget table](../extensions/atlas_union/results/budget_summary.csv) separates the original primary questions from the complete union census. At 4 µm², union invariance on the original census is 2,524/3,304 in lung and 15,143/22,976 in breast; the unchanged primary fractions are 2,780/3,304 and 17,717/22,976. On the expanded union census, invariance is 2,905/3,766 and 17,160/25,446. Each endpoint selects one complete member atlas.

The [reference-reach record](../extensions/atlas_union/results/reference_summary.csv) uses the same comparison keys and exclusions as the width study. Unbounded-budget reach increases from 185 to 315 of 2,798 lung Proseg changes, from 184 to 299 of 1,114 lung nuclear-reference changes, and from 1,773 to 2,942 of 10,965 breast nuclear-reference changes. These are exact discrete count-membership checks, rather than interval-containment tests or claims about a jointly attainable count matrix.

The [contiguous-region measurements](../extensions/scaling/results/contiguous_summary.json) and [independent source audit](../extensions/scaling/results/source_audit.json) identify both completed regions and their checks. Separate scalar comparisons verify local spectra for [400 µm](../extensions/scaling/results/central400_sample_verification.json) and [600 µm](../extensions/scaling/results/central600_sample_verification.json). The [reproduction instructions](../extensions/scaling/README.md) restore both source/output artifacts against the [manifest](../extensions/scaling/artifact_manifest.json). The 600 µm archive contains all 38 recorded files, including the initial interrupted attempt and final completed measurement.

The contiguous-region output audit independently checks composition and exported witness counts/costs. Continuous geometry and unsampled local optimality use the unchanged exact solver. Each scalar sample supplies a separate local-optimality check for its 16 selected queries only. The 400 µm selection was recorded before production spectra. The 600 µm selection applies the same pre-existing input-only hash rule after the count run completed and before its production spectra were inspected; it is not a preproduction-frozen selection.

## Evidence limits

Widths and budgets are sensitivity settings without anatomical calibration. The study has no independent transcript-ownership truth, additional biological donor or validated downstream inference. Breast material is a vendor format-test collection. Proseg replication and ownership accuracy were not assessed. Three-label atlas ribbons, moving junctions, topology changes and whole-slide deployment lie outside the completed experiments.

The full eight-width optimization study was run once in the recorded environment. CI replays both primary collections and one width family, then checks every saved width composition and interval row. It does not rerun the entire optimization study or independently enumerate every local optimum. These computational checks are distinct from specialist review of the mathematical proofs.

Both contiguous regions use the same breast format-test specimen, overlap one another and earlier crops, and are not pooled as independent samples. The 400 µm region uses the chain backend; its 22.55-minute pipeline measurement includes extraction, geometry validation, count computation and serialization, with peak stage RSS 2.43 GiB.

The 600 µm region uses the existing exact assignment-enumeration backend, selected from input multiplicities before count outcomes. An initial global deadline interrupted count computation; a retry reused the validated atlas and completed all required questions. The accepted extraction, geometry and final count stages total 29.72 minutes. Cumulative measured time across both attempts is **36.89 minutes**, with peak completed-stage RSS 2.20 GiB. The final run has no incomplete local spectra; the earlier deadline interruption remains explicit. Different backends and the retry mean these two measurements do not define a controlled chain-scaling curve. Source preparation, coverage accounting and output reconciliation are separate, and the host is shared.

The [partial original-crop replay](../extensions/scaling/results/partial_replay_summary.json) preserves a separate failed three-crop output check and reports no accepted timing for that failed run.

Historical development and review notes are retained in [history/](history/) and [legacy/](legacy/). Their release instructions and status statements describe the versions at which they were written.
