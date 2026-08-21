# 下一刀：3D 非球形 synthetic basin 压力测试 v1.0

## 目标

测试当前 regional selector 在一个故意保守、可完全复现的三维非球形环境中会不会“过度激活” V 相。该测试不是 N-body 宇宙学预测，只是 selector stress test。

## 构造

- 周期盒：\(256\,h^{-1}{\rm Mpc}\)。
- 网格：\(80^3\)。
- 线性 Gaussian random field。
- transfer shape：Eisenstein–Hu-like。
- 精确归一化到 \(\sigma_8=0.811\)。
- Gaussian smoothing：\(R_{\rm smooth}=4\,h^{-1}{\rm Mpc}\)。
- watershed basin segmentation。
- persistence thresholds：\(0,0.05,0.10,0.20\) 倍 smoothed-field rms。
- 固定随机种子：11, 23, 47, 91。

每个 basin 计算：

1. basin volume；
2. effective radius \(R_{\rm eff}\)；
3. regional mean \(\delta_L\)；
4. spherical underdensity map 给出的 conditional \(\Gamma_\Omega\)；
5. finite-size activation condition

\[
\Gamma_\Omega\ge
1+\frac{2C_\sigma\xi_0}{\beta_\Gamma R_{\rm eff}},
\]

其中此处只为压力测试固定

\[
C_\sigma=0.396830116706,\quad
\xi_0=3.12938\,h^{-1}{\rm Mpc},\quad
\beta_\Gamma=1.
\]

## 数值核验

\[
\max|\sigma_8^{\rm measured}-0.811|=1.11\times10^{-16}.
\]

无 persistence pruning 时：

\[
\phi(\Gamma\ge1)=0.031426\text{--}0.053537.
\]

同时满足 finite-size boundary 的体积分数：

\[
\phi_{\rm active}=0.016977\text{--}0.029367.
\]

合并全部 persistence 设置：

\[
\boxed{
\phi_{\rm active}=0.006686\text{--}0.029367
}
\]

即约

\[
\boxed{0.67\%\text{--}2.94\%}.
\]

## 判决

这是一个**否定性结果**。

瞬时、单尺度、Gaussian watershed basin catalog 并不会自然产生 volume-dominant V phase。因此不能把此前 spherical first-passage 的数量级匹配直接提升为宇宙学验证。

当前最小结论变为：

\[
\boxed{
\text{若 V phase 真实存在，则 nonlinear basin hierarchy、history、merger/coarsening 或不同 coarse graining 至少有一项不可省略。}
}
\]

这一结果被放入 Supplementary Note 8，而不是主稿的成功证据。

## 对理论创新性的意义

这个结果不增加新的拟合能力，反而删除一个过于简单的实现。它使 regional phase hypothesis 更可证伪：未来真实 N-body / relativistic simulation 必须展示 basin history 如何把瞬时小 filling fraction 转化为合适的 late-time phase occupation；否则 V-branch cosmology 应被否定。

## 可复现入口

- `src/analysis_synthetic_3d_basins.py`
- `tests/test_repro.py`
- `results/synthetic_3d_basin_summary.txt`
- `results/synthetic_3d_basin_table.csv`

当前完整回归状态：**28 passed**。
