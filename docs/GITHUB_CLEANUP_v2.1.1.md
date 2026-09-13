# GitHub cleanup audit — v2.1.1

The default branch had remained byte-for-byte aligned with the historical v1.2.0 submission snapshot while the manuscript advanced substantially. That state was scientifically confusing rather than corrupt: the old release was valid, but `main` no longer represented the current PRD manuscript.

The cleanup therefore follows a non-destructive policy:

1. preserve tag/release `v1.2.0` as an immutable historical snapshot;
2. retain the v1.2 baseline modules at repository root so old reproduction commands remain intelligible;
3. place PRD v2.1.1 follow-up source under `current/` instead of silently overwriting historical analyses;
4. update README, reproducibility documentation, package metadata, citation metadata, and changelog to identify v2.1.1 as the current default-branch state;
5. do not vendor multi-megabyte descendant-history intermediate CSV files into GitHub merely to duplicate the reviewer archive;
6. keep the full submission-matched reviewer archive as the exact heavy-data regression snapshot.

Before this cleanup, `docs/REPRODUCIBILITY.md` also contained stale bookkeeping (duplicate item numbering and an obsolete 22-test statement) even though the historical v1.2 suite had advanced to 30 tests. Those documentation residues are corrected on `main` without modifying the frozen `v1.2.0` tag.

Local validation of the cleaned source tree on 2026-09-12:

- `pytest -q`: 62 source-available tests passed;
- standalone PRD reviewer archive: 42/42 tests passed at packaging;
- historical v1.2 release: 30 deterministic tests frozen in that snapshot.

The test counts overlap conceptually; they are reported for reproducibility layers, not summed as independent scientific results.
