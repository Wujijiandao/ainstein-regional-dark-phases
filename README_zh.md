# AInstein 区域暗相框架

本仓库提供论文 **“Regional gravitational phases and a geometric active-charge law for dark-sector phenomenology”** 的可复现代码。

## 当前状态

默认分支 `main` 对应 **PRD 投稿候选版 v2.2.0**。此前 Communications Physics 投稿时使用的 **v1.2.0** 已由 GitHub tag/release 永久冻结，不改写历史快照。

仓库现在分两层：

- 根目录的 `src/`、`tests/`、`results/`、`figures/`、`docs/` 保留 v1.2 的公开基线计算；
- `current/` 保存 PRD v2.2.0 新增的 nonlinear PM、descendant history、control-state、causal/Kubo response、nonlinear 1PI、stochastic conversion 与 phase–geometry feedback 等计算。

这样既不会破坏旧版可复现性，也不会让默认分支继续假装停留在 v1.2。

## 研究范围

AInstein 检验的是：通常分别归因于暗物质与暗能量的部分引力现象，是否可以被组织为同一个粗粒化 matter–geometry sector 的不同**区域宏观构成态**。当前工作是可证伪的有效理论，不声称已经完成微观统一。

当前代码覆盖：三稳态/界面与成核、区域几何控制变量、非递归 basin 选择、Gaussian/Zel'dovich/nonlinear PM obstruction tests、descendant-anchored histories、leaf-constrained active charge、causal kinetics、retarded/Kubo response、nonlinear 1PI、stochastic conversion、phase–geometry feedback，以及 SPARC 一致性诊断。

PM 与 cosmic-web 数值部分是**压力测试/阻碍检验**，不是 precision late-time cosmology；PM 中的区域状态量是 operational proxy，不冒充直接的 relativistic Buchert measurement。

## 复现

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
pytest -q
```

根目录默认 `pytest -q` 运行历史 v1.2 测试和所有无需大型冻结中间表的当前确定性测试。依赖 multi-megabyte descendant-history CSV 的 branch-parity/causal-kinetics 完整回归仍保存在投稿匹配 reviewer archive 中，不在 GitHub 主仓库重复塞入大型中间表。更详细的分析入口和适用边界见 `docs/REPRODUCIBILITY.md` 与 `current/README_REPRODUCIBILITY.txt`。

SPARC 第三方 `RAR.mrt` 数据不随仓库再分发；使用 `python scripts/fetch_sparc.py` 从公开来源获取。

## 版本原则

- `v1.2.0`：历史 Communications Physics 投稿快照，永久冻结；
- `v2.2.0`：当前 PRD referee-hardened 投稿候选版。

已经发布的 tag/release 不强制移动；后续科学修改一律增加新版本号。

## 归档

稳定 Zenodo 项目 DOI：`10.5281/zenodo.22042564`  
ORCID：`0009-0000-3121-7972`

## 许可证

原创代码采用 MIT License。第三方数据继续遵守其原始许可条件。


### v2.2.0 新增

加入共享控制/global competition 与 global-budget identifiability 审计；不把目标 abundance 作为拟合输入。


## v2.2.0 GitHub 回归

默认仓库的 source-available 测试为 **67 passed**。submission-matched reviewer archive 的完整 follow-up 测试为 **47 passed**；两者有重叠，不能相加。
