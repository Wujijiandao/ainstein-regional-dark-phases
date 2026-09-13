# AInstein Regional Dark Phases — PRD submission-matched public release v2.1.1

This archive is the release-ready reproducibility snapshot corresponding to the PRD referee-hardened manuscript candidate v2.1.1.

## Status

- `main` is organized to match the PRD v2.1.1 reproducibility state.
- Historical release `v1.2.0` remains immutable and recoverable by tag.
- The stable Zenodo project DOI remains 10.5281/zenodo.22042564; a new Zenodo version should be minted only when the online release is actually deposited.

## Scientific scope

The archive preserves the deterministic calculations used in the PRD manuscript: nonlinear PM/descendant histories, control-state identities, leaf-global field benchmarks, branch parity and causal kinetics, retarded/Kubo response tests, nonlinear 1PI/nucleation audits, and reciprocal/nonreciprocal feedback benchmarks.

The v2.1.1 manuscript hardening does not introduce new numerical results. It changes claim scope and terminology only, including the notation fix that reserves `r_0` for the 1PI quadratic baseline and uses `zeta_Sigma = E_0/epsilon_V` for the legacy interface normalization ratio.

## Reproducibility

Run the test suites from the contained project directories as documented in `README_REPRODUCIBILITY.txt`. The submission package also records a combined regression run.
