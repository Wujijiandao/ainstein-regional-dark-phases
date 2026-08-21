> **Historical note:** This document records the v0.12 milestone. The v1.0.0 release has 28 passing regression tests; any 27/27 count below is intentionally retained as historical provenance.

# AInstein / Regional Dark Phases — v0.12 Hamiltonian–Landau 驱动力匹配

## 1. 研究问题

v0.11 已把区域相场写成 Allen–Cahn 型动力学，并得到球形 daughter basin 的 sharp-interface 关系

\[
R_c=\frac{2\sigma_\Sigma}{\Delta f}.
\]

真正缺失的是区域 Hamiltonian 状态 `Gamma_Omega` 如何驱动 Landau 相势。直接写任意函数 `H(Gamma)` 会重新引入一个自由函数，因此 v0.12 改为匹配**bulk free-energy difference**，因为 sharp-interface dynamics 真正依赖的是 `Delta f`，而不是 `H` 本身。

## 2. 最小匹配闭合

在 V 分支定义

\[
\Delta w_V(H)=W_D(H)-W_V(H)>0,
\]

并测试

\[
\boxed{\Delta w_V(H_\Omega)=\beta_\Gamma(\Gamma_\Omega-1),\qquad \Gamma_\Omega\ge1.}
\]

这里 `beta_Gamma>0` 是唯一新增的无量纲归一化/耦合。该式不是 GR 定理，而是一个比任意 `H(Gamma)` 更严格的有效闭合。

立即得到：

- `Gamma=1` 自动对应 D–V coexistence；
- 近共存点
  \[
  -H-H_{\rm coex}=\frac{\beta_\Gamma}{\Delta\chi}(\Gamma-1)+O[(\Gamma-1)^2],
  \]
  其中 `Delta_chi=0.948007800871`；
- 任意函数自由度被压缩成 `beta_Gamma` 一个归一化。

## 3. 与相墙曲率合并：有限区域位移

sextic 相墙常数

\[
C_\sigma=0.396830116706,
\]

且

\[
\frac{R_c}{\xi_0}=\frac{2C_\sigma}{\Delta w}.
\]

代入匹配式后得到精确有限尺寸关系

\[
\boxed{
\frac{R_c}{\xi_0}=\frac{2C_\sigma}{\beta_\Gamma(\Gamma-1)}.
}
\]

等价地

\[
\boxed{
\Gamma_{\rm grow}(R)
=1+\frac{2C_\sigma}{\beta_\Gamma}\frac{\xi_0}{R}.
}
\]

物理意义：`Gamma=1` 只是无限大区域的共存条件；有限 basin 有界面曲率代价，因此必须继续进入 geometry/correlation dominated 区域才会自行增长。这是区域相场版的 Gibbs–Thomson / capillary finite-size shift。

## 4. 条件性尺度核验

沿用 v0.11 的独立条件性 benchmark：

- `d_10-90 = 5.298450046 h^-1 Mpc`；
- `xi0 = 3.12938181925 h^-1 Mpc`；
- hydraulic basin scale `R = 15.6 h^-1 Mpc`；
- `R/xi0 = 4.98501010776`。

仅取 `beta_Gamma=1` 作为 matched-normalization 示例，则

\[
\boxed{\Gamma_{\rm grow}=1.15920935289},
\]

并由完整 sextic 反解

\[
\boxed{H_V=-0.512592527206}.
\]

对应

\[
|H|-H_{\rm coex}=0.170904032188,
\qquad
\Delta w=0.159209352891.
\]

与 v0.11 从临界半径独立得到的数值完全一致。这只是 algebraic consistency check；`beta_Gamma=1` 不是微观推导。

## 5. 这一刀真正消掉了什么？

v0.11 的缺口是：

\[
H_\Omega(\Gamma_\Omega)=\text{任意未知函数}.
\]

v0.12 后变成：

\[
\boxed{
\{\beta_\Gamma,\xi_0\}
}
\]

两个物理量控制有限 basin 的驱动和曲率代价。也就是说，**函数自由度被消掉了，但微观尺度没有被消掉。**

下一步若要继续提高第一性原理程度，需要从粗粒化 matter–geometry dynamics 推出：

1. `beta_Gamma` 的归一化；
2. `xi0=sqrt(kappa/E0)` 的 Mpc 级来源；
3. `Gamma_Omega` 的协变或可操作替代；
4. basin lineage 上的噪声/成核/退相机制。

## 6. 可证伪预测雏形

若 `beta_Gamma xi0` 是统一的宏观尺度，则不同大小的 expanding basins 的真正增长阈值不相同：

\[
\Gamma_{\rm grow}-1\propto R^{-1}.
\]

因此小 basin 需要更强的 geometry/correlation dominance 才能进入稳定 V 相，大 basin 更接近 `Gamma=1` 即可转换。这是一个可在 basin-tree simulation 中直接检查的 size–state correlation，而不是一个自由 `w(z)`。

## 7. 数值与代码

新增：`src/analysis_bias_matching.py`

回归测试：27/27 passed。

主要 benchmark：

```text
Gamma_growth_threshold = 1.15920935289
matched_H_V = -0.512592527206
matched_Delta_w = 0.159209352891
R_c/xi backcheck = 4.98501010776
```
