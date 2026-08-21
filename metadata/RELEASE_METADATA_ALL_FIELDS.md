# GitHub + Zenodo 发布元数据总表 v1.0

本文件汇总公开代码仓库和 Zenodo 软件归档所需字段。机器可读文件已经放在 GitHub-ready 仓库根目录：`CITATION.cff`、`.zenodo.json`、`codemeta.json`。

## A. GitHub repository

- Repository name: `ainstein-regional-dark-phases`
- Visibility: Public（准备公开后）
- Description: `Reproducibility code for regional gravitational phase transitions, cosmic-web basin selectors, effective dark-sector emergence, and galaxy active-charge diagnostics.`
- Homepage / Website: 暂留空；Zenodo DOI 或论文 DOI 存在后再填稳定 landing page。
- License: MIT
- Default branch: `main`
- Release tag: `v1.0.0`
- Release title: `AInstein Regional Dark Phases v1.0.0`
- Pre-release: No
- Latest release: Yes

### GitHub topics

- `dark-matter`
- `dark-energy`
- `emergent-gravity`
- `cosmology`
- `cosmic-web`
- `phase-transition`
- `galaxy-dynamics`
- `sparc`
- `reproducible-research`
- `python`

### About text

> Reproducibility code for a regional gravitational phase framework linking cosmic-web organization, effective dark-sector phases, and capacity-like galactic active charge.

### Release body

使用仓库中的 `metadata/GITHUB_RELEASE_NOTES_v1.0.0.md`。

### Recommended repository settings

- Issues: On
- Discussions: Optional
- Wiki: Off initially
- Releases: On
- Branch protection for `main`: Optional before submission; recommended once v1.0.0 is archived.
- Zenodo integration: 在创建最终 GitHub Release **之前**启用仓库归档，确保 release 自动被 Zenodo 捕获。

## B. Zenodo software deposition

### Core fields

- Resource type: Software
- Title: `AInstein Regional Dark Phases: Reproducibility Code`
- Version: `1.0.0`
- Publication date: `2026-08-21` **仅在这确实是公开 release 日期时使用；否则改成真实日期**
- Access right: Open
- License: MIT
- Language: English (`eng`)

### Creator

- Name: `Zhang, Yuzhan`
- Given name: `Yuzhan`
- Family name: `Zhang`
- Affiliation: `Independent Researcher`
- ORCID: 留空，除非已有并核验真实 ORCID

### Description

> Deterministic reproducibility code for the manuscript “Regional gravitational phases as an effective framework for dark-sector unification”. The release reproduces regional phase-transition diagnostics, matter–geometry selector calculations, basin hierarchy and interface tests, geometric active-charge results, the SPARC radial-acceleration diagnostic, and a three-dimensional synthetic watershed stress test. The framework is an effective theory and the 3D Gaussian test is an obstruction test rather than an N-body prediction.

### Keywords

- dark matter
- dark energy
- emergent gravity
- phase transitions
- cosmic web
- galaxy dynamics
- SPARC
- reproducible research

### Contributors

当前没有必须添加的其他 contributor。不要把 AI 系统列为作者或 contributor。

### Related identifiers

等真实标识符存在后再填，不要提前杜撰：

- Journal article DOI → 软件记录与论文建立 `isSupplementTo` / Zenodo 当前 UI 中语义最接近的关系。
- Preprint DOI / arXiv identifier → 如以后确实公开，可作为 related identifier。
- GitHub repository URL → GitHub integration 通常自动记录来源；如手工填写，避免把源码 URL 错当 DOI。
- Previous / next software versions → 通过 Zenodo “New version” 管理，Zenodo 会自动建立版本关系。

### Funding / grants

- None declared for this work.
- 不填虚构 grant。

### Communities

- Optional。
- 只在真实符合某个 community 的 curator scope 时添加。

### Subjects / domain

如果 Zenodo UI 要求或允许，可使用：

- Physics
- Astrophysics and cosmology
- Gravitation
- Complex systems
- Computational physics

### Files

最优方式：由 GitHub release `v1.0.0` 经 Zenodo GitHub integration 自动归档。

不要把第三方 SPARC 原始数据打进公开 Zenodo 软件 release，除非明确确认其再分发许可。GitHub-safe 仓库通过 `scripts/fetch_sparc.py` 取得公开数据并记录预期哈希。

## C. Machine-readable metadata hierarchy

仓库同时提供：

1. `.zenodo.json` — Zenodo-specific metadata；若它与 `CITATION.cff` 同时存在，Zenodo 当前说明会优先采用 `.zenodo.json`。
2. `CITATION.cff` — GitHub 用于显示 citation suggestion；即使 Zenodo 有 `.zenodo.json` 也建议保留。
3. `codemeta.json` — 通用 software metadata / interoperability。
4. `LICENSE` — MIT license full text。

## D. Zenodo DOI strategy

推荐流程：

1. 在 Zenodo 连接 GitHub repository。
2. 检查 `.zenodo.json` 和 release metadata。
3. 在 GitHub 创建 `v1.0.0` release。
4. Zenodo 自动归档并 mint version DOI。
5. 记录 concept DOI + version DOI。
6. 稿件 Code availability 优先引用具体 `v1.0.0` version DOI；README 可同时显示 concept DOI 便于追踪后续版本。
7. 将真实 DOI 回填到稿件、README、CITATION/related identifier（如需要）。

## E. Citation sentence after DOI exists

软件 DOI 尚未生成前不要写假 DOI。生成后可使用：

> Zhang, Y. *AInstein Regional Dark Phases: Reproducibility Code*, version 1.0.0, Zenodo, DOI: [REAL DOI].

## F. GitHub release integrity

发布前执行：

```bash
python -m pytest -q
python run_all.py
sha256sum ...
```

当前冻结回归状态：**28 passed**。3D synthetic stress test 的 finite-size-active volume fraction 范围为 **0.006686–0.029367**。
