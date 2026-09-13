# AInstein Regional Dark Phases — PRD science-closure public release candidate v2.2.0

This repository snapshot corresponds to the bounded science-closure revision prepared before PRD submission.

## New relative to v2.1.1

- Adds the shared-control/global-competition audit under `current/research_v2.2/`.
- Generalizes the lowest-order passive reciprocal feedback sign restriction to arbitrary stable shared Gaussian controls with real reciprocal linear coupling.
- The exact Schur-complement elimination produces `-B^T K^{-1} B`, a negative-semidefinite induced kernel in the coupled regional observables; passive shared Gaussian mediation is therefore cooperative rather than a source of positive long-range competition.
- Adds a global-budget identifiability result: an odd conserved phase moment fixes only C/V imbalance, while an even phase moment can fix total converted fraction; the present one-sector scaffold conserves the material current, not `chi^2`.
- Keeps abundance predictive status explicitly open rather than fitting a target fraction.

## Repository policy

The root-level source/result tree preserves the immutable v1.2 public baseline. `current/` contains the post-v1.2 calculations used by the PRD manuscript. Large frozen descendant-history tables and reviewer-only intermediates remain in the submission-matched reproducibility archive and are intentionally not duplicated here.

This is a release candidate. It does not claim that GitHub tag/release v2.2.0 or a Zenodo version has already been published online.
