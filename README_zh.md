# AInstein 区域暗相理论：可复现代码

本仓库对应稿件 **《Regional gravitational phases as an effective framework for dark-sector unification》**，用于复现稿件中的解析常数、数值审计、星系诊断与三维合成 basin 压力测试。

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

## v1.0 冻结结果

- 回归测试：**28 passed**；
- 三维 Gaussian/watershed 压力测试：冻结的 seed/persistence 网格中，有限尺寸激活体积分数仅 **0.67%–2.94%**；
- SPARC 诊断：`N=2693`，`a_E=1.097458984856508e-10 m s^-2`，`chi2/dof=1.610215393925522`，raw RMS `0.1329097576083664 dex`。

三维结果是**否定性/压力测试**，不是 N-body 预测。它说明“瞬时、单尺度 Gaussian basin”不足以自动产生体积占优的 V 相，真实理论仍需要非线性 basin 层级、历史或不同的粗粒化规则。

## 许可与归档

原创代码采用 MIT License。`CITATION.cff` 用于 GitHub 引用展示，`.zenodo.json` 用于 Zenodo GitHub integration。Zenodo 发布后应把软件 DOI 回填到论文 Code Availability，并在论文 DOI 出现后建立软件—论文的 related identifier。
