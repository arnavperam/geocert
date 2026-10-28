# Reproduce a claim

The exact core is unchanged from release v1.0.0. This revision adds a reader-facing
manuscript, gene-count analysis, source-count verification and missing artifact
access. Each command has a distinct scope.

| Claim | Command | Evidence |
|---|---|---|
| Exact five-transfer witness | `python geocert.py showcase --output run/showcase.json` | Exact rational costs, continuous boundaries, flow certificate |
| Direct versus height objectives | `python geocert.py benchmark --repetitions 1 --scc --output-dir run/comparison` | 252 objectives per representation; expected hash in README |
| All 26 feasible showcase costs | `python experiments/independent_showcase.py . run/independent_showcase.json` | Separate rational scalar integration and complete enumeration; shared primitive ingestion |
| Complete additive ablation and gene-count consequences | `python experiments/readout_audit.py --output-dir run/readout` | 756 fresh optimizations, 1,008 archived objective/regret checks, all 252 readout rows |
| Derived counts bound to masks and transcripts | `python experiments/build_readout_inputs.py --source-root . --output run/readout_inputs.json` | 7,370 crop records and 464 repeated corridor records; requires NumPy and PyArrow |
| SCC and prefix correctness tests | `python tests/test_scc.py` and `python tests/test_prefix_bound.py` | Bounded exact implementation tests |

The independent-showcase command requires the output parent directory to exist;
the preceding showcase command creates `run/`.

## Full historical companion

The complete, unchanged publication companion is stored in two binary parts under
`artifacts/`. Download and
join them with the standard-library helper:

```sh
python tools/fetch_companion.py run/publication_companion.zip
```

The reconstructed archive must have SHA-256
`c786d48ce8925c5e5fe45dd4a96a7782edc5ba61377d2f833e9ddd0c9c53ba60`.
It contains `GeoCert_ST_Publication_Release/archive/reviewed_research_archive.zip`,
which preserves historical numerical, extraction and exploratory studies. The
helper never extracts or executes the archive. The revision's scientific
supplement focuses on the current theory and experiments; historical release
records remain available here rather than occupying the scientific narrative.

## Interpret the new census correctly

47 baseline and 88 relaxed-family joint answers change a gene-detection indicator
relative to the source. These are **not** counts of failures of the additive
ablation. The strict-regret counts remain 1 and 16; detection calls differ between
the returned joint/additive answers in 1 and 13 of those cases. Detection means a
retained Q20 transcript count greater than zero, not biological expression truth.
Other crop assignments stay fixed, each corridor is evaluated separately, and
alternative tied optima outside the unique showcase may give different readouts.
