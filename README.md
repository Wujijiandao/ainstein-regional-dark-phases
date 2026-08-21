# AInstein Regional Dark Phases

Reproducibility code for the manuscript **“Regional gravitational phases and a geometric active-charge law for dark-sector phenomenology.”**

## Scope

This repository reproduces the analytical and numerical diagnostics used in the manuscript. The framework tests whether dark-matter-like and dark-energy-like phenomenology can be organized as regional phases of one coarse-grained gravitational sector. It is an effective theory, not a microscopic completion.

The repository contains deterministic modules for:

- minimal tri-stable phase potential and hysteresis;
- spherical regional state coordinate `Upsilon`;
- two-barrier first-passage bookkeeping;
- Planck-normalized power-spectrum timing audits;
- matter–geometry equality selector;
- hierarchical basin antichain and history-dependent occupation;
- conversion-sign and interface-coarsening tests;
- diffuse phase-wall and thin-interface constants;
- Hamiltonian-bias / finite-size matching;
- geometric active-charge and SPARC diagnostics;
- 3D Gaussian-field watershed stress test;
- history-dependent 3D Zel'dovich Lagrangian-basin stress test.

## Reproduce

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
pytest -q
python run_all.py
```

The GitHub-safe repository does not redistribute the third-party SPARC `RAR.mrt` table. To run the SPARC fit, fetch the public file using:

```bash
python scripts/fetch_sparc.py
python src/analysis_emg.py
```

## Frozen v1.2 checks

- regression tests: **30 passed**;
- 3D synthetic-basin stress test: finite-size-active volume fraction **0.67%–2.94%** across the frozen seed/persistence grid;
- Zel'dovich history benchmark at `D=1`: strict instantaneous **0.26%–3.89%**, absorbing-history **0.73%–4.85%**, kinematic upper envelope **1.50%–8.72%**;
- SPARC diagnostic (when the public table is present): `N=2693`, `a_E=1.097458984856508e-10 m s^-2`, `chi2/dof=1.610215393925522`, raw RMS `0.1329097576083664 dex`.

The 3D synthetic test is intentionally an obstruction test, not an N-body prediction. Its low active volume fraction shows that an instantaneous single-scale Gaussian watershed is insufficient to generate a volume-dominant V phase.

## Repository layout

- `src/` analysis modules
- `tests/` deterministic regression tests
- `scripts/` third-party data fetch helper
- `results/` frozen text/CSV outputs
- `figures/` reproducible manuscript figures
- `docs/` derivation and reproducibility notes
- `metadata/` GitHub and Zenodo release metadata

## Citation and archival release

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22042564.svg)](https://doi.org/10.5281/zenodo.22042564)

- GitHub repository: https://github.com/Wujijiandao/ainstein-regional-dark-phases
- Submission release: https://github.com/Wujijiandao/ainstein-regional-dark-phases/releases/tag/v1.2.0
- Stable software archive DOI: `10.5281/zenodo.22042564`
- Author ORCID: https://orcid.org/0009-0000-3121-7972

`CITATION.cff`, `.zenodo.json` and `codemeta.json` contain machine-readable metadata. The manuscript cites the stable project DOI together with the exact GitHub release tag used for the submission. Zenodo may additionally expose version-specific DOIs; the project-level DOI above is retained as the persistent software citation used in the manuscript.

## License

MIT License for original source code. Third-party data remain subject to their source terms and are not relicensed by this repository.
