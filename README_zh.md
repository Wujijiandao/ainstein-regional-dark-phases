# AInstein 区域暗相理论：可复现代码

本仓库对应稿件 **《Regional gravitational phases and a geometric active-charge law for dark-sector phenomenology》**，用于复现稿件中的解析常数、数值审计、星系诊断与三维合成 basin 压力测试。

核心研究问题是：暗物质式额外引力与暗能量式负压，是否可能是同一个粗粒化“物质—几何系统”在不同区域组织态下的两种宏观响应，而不必先验地把它们视为两种基本物质。

## 一键复现

```bash
python -m venv .venv
pip install -r requirements.txt
pytest -q
python run_all.py
```

GitHub-safe 版本不直接再分发第三方 SPARC `RAR.mrt`，可运行：

```bash
python scripts/fetch_sparc.py
python src/analysis_emg.py
```

## v1.2 冻结结果

- 回归测试：**30 passed**；
- 三维 Gaussian/watershed 压力测试：冻结的 seed/persistence 网格中，有限尺寸激活体积分数仅 **0.67%–2.94%**；
- SPARC 诊断：`N=2693`，`a_E=1.097458984856508e-10 m s^-2`，`chi2/dof=1.610215393925522`，raw RMS `0.1329097576083664 dex`。

三维结果是**否定性/压力测试**，不是 N-body 预测。它说明“瞬时、单尺度 Gaussian basin”不足以自动产生体积占优的 V 相，真实理论仍需要非线性 basin 层级、历史或不同的粗粒化规则。

## 许可与归档

原创代码采用 MIT License。`CITATION.cff` 用于 GitHub 引用展示，`.zenodo.json` 用于 Zenodo GitHub integration。Zenodo 发布后应把软件 DOI 回填到论文 Code Availability，并在论文 DOI 出现后建立软件—论文的 related identifier。


## v1.2 发布说明

新增 `analysis_zeldovich_basin_history.py`：固定 Lagrangian watershed basin 经 Zel'dovich 形变演化，直接由 Jacobian 计算区域体积、密度、膨胀率和 `Gamma`。在 `D=1` 的冻结基准中，严格 active 体积分数为 0.26%--3.89%，吸收历史为 0.73%--4.85%，去掉 basin-wide single-stream 条件后的运动学上限仍仅 1.50%--8.72%。这是一项 obstruction test，不是 N-body 暗能量丰度预测。

公开仓库：https://github.com/Wujijiandao/ainstein-regional-dark-phases

本轮投稿对应 GitHub Release：`v1.2.0`  
预期固定链接：https://github.com/Wujijiandao/ainstein-regional-dark-phases/releases/tag/v1.2.0

稳定 Zenodo 软件归档 DOI：`10.5281/zenodo.22042564`  
DOI 链接：https://doi.org/10.5281/zenodo.22042564

ORCID：https://orcid.org/0009-0000-3121-7972

论文与代码都使用上述稳定 DOI，并同时记录精确 GitHub release tag 以锁定本次投稿所对应的代码版本。Zenodo 仍可能显示各版本的单独 DOI；本项目在稿件中使用上述项目级稳定 DOI 作为长期引用入口。
