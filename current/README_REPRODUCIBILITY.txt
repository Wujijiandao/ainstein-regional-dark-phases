PRD v2.2.0 current-source reproducibility layer.

The GitHub default branch contains source code for all submission-matched follow-up modules. To keep the repository readable, multi-megabyte frozen descendant-history intermediate tables are not duplicated here; they remain in the submission-matched reviewer archive.

From repository root:
  pytest -q

runs the baseline v1.2 suite plus current source-only deterministic suites.

The GitHub source-available regression suite for this v2.2.0 candidate produced **67 passed** at packaging.

The standalone submission-matched reviewer archive additionally includes the large frozen descendant-history tables and runs the full follow-up suite, including `research_v2.2/tests`; that archive produced **47 passed** at packaging. These two counts overlap substantially and must not be added.

Public project repository: https://github.com/Wujijiandao/ainstein-regional-dark-phases
Stable Zenodo project record: https://doi.org/10.5281/zenodo.22042564
