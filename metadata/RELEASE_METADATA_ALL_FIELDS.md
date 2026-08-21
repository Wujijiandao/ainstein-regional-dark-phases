# Release metadata — AInstein Regional Dark Phases v1.2.0

## GitHub release
- Repository: `Wujijiandao/ainstein-regional-dark-phases`
- Tag: `v1.2.0`
- Release title: `AInstein Regional Dark Phases v1.2.0`
- Release URL: https://github.com/Wujijiandao/ainstein-regional-dark-phases/releases/tag/v1.2.0
- Release date: 2026-08-21
- License: MIT

## Stable archival identifier
- Zenodo project DOI: `10.5281/zenodo.22042564`
- DOI URL: https://doi.org/10.5281/zenodo.22042564
- ORCID: `0009-0000-3121-7972`

The manuscript cites the stable Zenodo project DOI and the exact GitHub release tag. Zenodo may additionally assign a release-specific DOI; the project DOI above is retained as the long-lived software citation entry used in the manuscript.

## Release summary
Submission-hardening release following adversarial pre-submission review. No frozen scientific regression target has been changed. The release contains 30 deterministic tests, the static Gaussian-watershed obstruction test, the history-dependent Zel'dovich basin obstruction test, the conditional leaf-constrained active-charge calculation, phase/interface diagnostics, and the SPARC consistency check.

## GitHub release body
Use `metadata/GITHUB_RELEASE_NOTES_v1.2.0.md`.

## Before release
1. Commit the v1.2.0 metadata and code to `main`.
2. Run `pytest -q` and confirm `30 passed`.
3. Confirm the public SPARC helper works or run the self-contained archive privately.
4. Create tag/release `v1.2.0`.
5. Confirm the release URL resolves.
6. Confirm Zenodo GitHub integration archives the release under the project record.
7. Submit the manuscript only after the public GitHub release and DOI both resolve.
