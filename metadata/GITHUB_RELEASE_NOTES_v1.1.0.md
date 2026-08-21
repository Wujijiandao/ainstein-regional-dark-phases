# AInstein Regional Dark Phases v1.1.0

This release accompanies the revised Communications Physics submission manuscript.

## New in v1.1.0

- adds `analysis_zeldovich_basin_history.py`, a history-dependent 3D Lagrangian-basin Zel'dovich stress test;
- follows fixed watershed basins through quasi-linear deformation and computes regional volume, density, expansion and `Gamma` directly from the deformation Jacobian;
- freezes strict pre-shell-crossing, absorbing-history and kinematic-upper-envelope active-volume diagnostics;
- adds `fig_zeldovich_history.pdf` and corresponding regression tests;
- updates manuscript-facing terminology, ORCID and repository metadata.

## Key new obstruction result

At `D=1` over four fixed seeds and three persistence settings:

- strict instantaneous finite-size-active volume: **0.261%--3.891%**;
- history-absorbing active volume: **0.730%--4.849%**;
- kinematic upper envelope without the basin-wide single-stream requirement: **1.499%--8.722%**.

Thus minimal Zel'dovich history does not generate a volume-dominant expanding phase in the frozen benchmark. The test is not an N-body abundance prediction.

## Existing archive

Stable project archive DOI: 10.5281/zenodo.22042564. Zenodo may additionally expose release-specific version DOIs. This historical release note is retained for provenance.
