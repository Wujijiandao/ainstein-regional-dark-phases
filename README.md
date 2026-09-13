# AInstein Regional Dark Phases

Reproducibility code for the manuscript **“Regional gravitational phases and a geometric active-charge law for dark-sector phenomenology.”**

## Current status

The default branch tracks the **PRD science-closure candidate v2.2.0**. The earlier Communications Physics snapshot remains permanently preserved by the GitHub tag/release **v1.2.0** and is not rewritten.

The repository is intentionally split into two reproducibility layers:

- the top-level `src/`, `tests/`, `results/`, `figures/`, and `docs/` directories preserve the public v1.2 baseline calculations;
- `current/` contains the additional nonlinear particle-mesh, descendant-history, control-state, causal-response, nonlinear-1PI, stochastic-conversion, and phase–geometry feedback calculations used by the PRD v2.2.0 manuscript.

This layout keeps the historical release reproducible while making the current manuscript additions explicit instead of silently overwriting the original submission snapshot.

## Scientific scope

AInstein tests whether part of the phenomenology usually assigned separately to dark matter and dark energy can be organized as different **regional constitutive states of one coarse-grained matter–geometry sector**. It is an effective-theory programme, not a completed microscopic unification.

The current code covers:

- tri-stable phase structure and interface/nucleation diagnostics;
- regional expansion/geometry control variables and exact spherical benchmarks;
- non-recursive basin antichain bookkeeping;
- Gaussian, Zel'dovich, and nonlinear particle-mesh obstruction tests;
- descendant-anchored basin histories and selector-ceiling audits;
- leaf-constrained active-charge geometry and non-spherical field benchmarks;
- branch-parity, causal phase kinetics, and passive-memory tests;
- retarded/Kubo response identification;
- nonlinear 1PI reconstruction and stochastic-conversion benchmarks;
- reciprocal and nonreciprocal phase–geometry feedback audits;
- shared-control/global-competition and global-budget identifiability audits;
- the SPARC radial-acceleration consistency diagnostic.

The numerical structure-formation calculations are deliberately **stress/obstruction tests**, not precision late-time cosmological simulations. The PM-derived regional variables are operational proxies and are not claimed to be direct relativistic Buchert measurements.

## Reproduce

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
pytest -q
```

For the v2.2.0 clean GitHub candidate this source-available suite reports **67 passed**. The larger reviewer archive reports **47 passed** with frozen history tables; the two counts overlap and are not additive.

The root pytest configuration runs the historical v1.2 suite plus all current source-only deterministic suites. The two PM-history-dependent branch-parity/causal-kinetics suites require the larger frozen descendant-history tables supplied in the submission-matched reviewer archive and are documented separately rather than silently vendored as multi-megabyte CSV files. For analysis-by-analysis commands and scope notes, see:

- `docs/REPRODUCIBILITY.md`
- `current/README_REPRODUCIBILITY.txt`

The public repository does not redistribute the third-party SPARC `RAR.mrt` table. Fetch it from its public source with:

```bash
python scripts/fetch_sparc.py
python src/analysis_emg.py
```

## Repository layout

- `src/` — baseline analytical/numerical modules preserved from v1.2
- `tests/` — baseline deterministic regression tests
- `results/` — baseline frozen text/CSV outputs
- `figures/` — baseline reproducible manuscript figures
- `current/` — PRD v2.2.0 manuscript-facing additions
- `scripts/` — data-fetch helpers
- `docs/` — reproducibility and scope notes
- `metadata/` — release/data/code metadata

## Versioning policy

- **v1.2.0**: immutable historical Communications Physics submission snapshot.
- **v2.1.1**: PRD referee-hardened checkpoint.
- **v2.2.0**: current bounded PRD science-closure candidate on `main`.

A tag or release must never be force-moved to a different commit. New scientific revisions receive a new version.

## Citation and archive

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22042564.svg)](https://doi.org/10.5281/zenodo.22042564)

- Repository: https://github.com/Wujijiandao/ainstein-regional-dark-phases
- Stable Zenodo project DOI: `10.5281/zenodo.22042564`
- ORCID: https://orcid.org/0009-0000-3121-7972

The stable project DOI identifies the software archive across versions. `CITATION.cff`, `.zenodo.json`, and `codemeta.json` describe the current default-branch snapshot; the exact historical v1.2 metadata remain recoverable from tag `v1.2.0`.

## License

MIT License for original source code. Third-party data remain subject to their original terms and are not relicensed by this repository.
