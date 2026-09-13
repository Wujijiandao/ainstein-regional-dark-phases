# AInstein / Regional Dark Phases — v2.2 nonlocal/global competition 收束研究

日期：2026-09-13

## 1. 本轮问题与停止条件

v2.0 已证明：在最低阶局域、被动、互易 Gaussian phase--geometry feedback 中，消去稳定控制变量后只会得到 parent softening，不能自然把转换停在有限 abundance。剩余路线被压缩为：

1. higher-order saturating geometry response；
2. nonlocal/global competition among regions；
3. genuinely nonequilibrium nonreciprocal feedback。

本轮只攻击第 2 条，而且预先规定停止条件：如果最小 nonlocal/global completion 仍不能在不输入目标 abundance 的情况下产生自限，就把这一点作为新的 obstruction/identifiability result，并停止继续堆机制。这样 v2.2 可以作为 PRD 前最后一次有边界的科学扩展。

## 2. 文献边界：global constraint 可以限相分数，但“限”不等于“预测”

相场理论里存在两种必须区分的 nonlocality。

第一类是 **硬守恒约束**。Rubinstein--Sternberg 型 nonlocal Allen--Cahn 通过全局 Lagrange multiplier 强制质量/相体积守恒；现代工作也证明其 sharp-interface limit 对应 volume-preserving mean-curvature flow。这里总相量之所以固定，是因为它被预先指定为 conserved quantity。

第二类是 **显式长程竞争项**。Ohta--Kawasaki 模型含正的 nonlocal repulsive contribution，可阻止普通宏观相分离并形成有限尺度结构。这个正 nonlocal term 是新的物理成分，不是把一个普通稳定 Gaussian auxiliary control 消去后自动得到的结果。

宇宙学侧，Buchert/Wiegand 的 volume partition 也强调 regional domains 与 global averages 之间必须满足一致性关系；regional expansion 不能被当成彼此独立的小宇宙。但这些 partition identities 本身没有给出一个 conserved phase-order budget。

因此本轮真正的问题不是“global coupling 能不能写”，而是：

\[
\boxed{\text{现有 AInstein 守恒/几何结构是否自动产生具有正确符号的 global competition？}}
\]

## 3. Nonlocal reciprocal Gaussian sign theorem

设所有区域的相观测量组成向量

\[
q=q(\chi_1,\ldots,\chi_N),
\]

并引入任意有限维共享控制变量 \(y\)。考虑最一般的稳定实 Gaussian reciprocal benchmark：

\[
\mathcal F(q,y)=\mathcal F_0(q)
+\frac12(y-y_0)^T K (y-y_0)-y^TBq,
\qquad K=K^T\succ0.
\]

对 \(y\) 精确消元：

\[
y_*=y_0+K^{-1}Bq,
\]

得到

\[
\boxed{
\mathcal F_{\rm eff}(q)
=\mathcal F_0(q)-y_0^TBq-\frac12 q^T G q,
\qquad
G\equiv B^TK^{-1}B\succeq0.}
\]

所以 induced nonlocal kernel 是

\[
\boxed{-G\preceq0.}
\]

这给出一个比 v2.0 局域标量例子更一般的结论：

> **任何由稳定实 Gaussian 共享控制变量、通过 reciprocal 线性耦合产生的静态 nonlocal feedback，都不能自动生成一个正半定的 long-range competition penalty。**

它的效果是 cooperative，而不是 self-limiting。

- 若 \(q_i\propto\chi_i\)，负核偏好 daughter signs 的 collective alignment；
- 若 \(q_i\propto\chi_i^2\)，负核偏好更大的 collective conversion amplitude。

这不是所有 nonlocal theory 的 no-go。它只针对：real、stable、Gaussian、reciprocal、linear-source elimination 这一类。显式 repulsive kernels、hard constraints、non-Gaussian sectors、gauge/constraint structures、nonreciprocal dynamics 均不在定理范围内。

### 数值回归

100 组随机正定 \(K\) 与随机 \(B\) 上：

- 精确消元 identity 最大误差：\(2.13\times10^{-14}\)；
- \(G\) 最小特征值：\(-2.61\times10^{-15}\)，即数值精度内 \(G\succeq0\)；
- induced \(-G\) 最大特征值：\(2.61\times10^{-15}\)，即数值精度内 \(-G\preceq0\)。

这不是拟合结果，而是 Schur-complement 解析恒等式的 regression check。

## 4. 离散 converted-fraction 压力测试

为了避免把连续 sextic 的细节误认为结论来源，再做一个 branch-reduced audit。令每个区域只有 unconverted/converted 两种状态 \(s_i\in\{0,1\}\)，局域转换代价为 \(\Delta_i\)，等权体积分数为

\[
f=\frac1N\sum_i s_i.
\]

被动 Gaussian 共享 mediator 消元后的最小 rank-one 结构是

\[
E=\frac1N\sum_i\Delta_i s_i-\frac{\gamma}{2}f^2,
\qquad \gamma\ge0.
\]

固定转换数 \(k\) 时，最优集合是选择最小的 \(k\) 个 \(\Delta_i\)。不同 \(k\) 的能量关于 \(\gamma\) 是直线，而且 \(k\) 越大斜率越负。因此 global minimizer 的 converted fraction 随 \(\gamma\) **只能不减**。这是一种 avalanche/cooperation，而不是 conversion-capacity depletion。

四个冻结异质样本（仅作结构回归，不是 cosmological calibration）：

| seed | \(\gamma=0\) baseline | scan 末端 passive Gaussian | 同强度显式正 competition |
|---:|---:|---:|---:|
| 11 | 0.3150 | 1.0000 | 0.1275 |
| 23 | 0.3525 | 1.0000 | 0.1425 |
| 47 | 0.3350 | 1.0000 | 0.1350 |
| 91 | 0.3000 | 1.0000 | 0.1125 |

四个 passive scans 全部严格满足 converted fraction 非递减。

对照项

\[
+\frac{\kappa}{2}f^2
\]

确实可以压低转换分数，但这正是一个**额外的正 long-range competition ingredient**。它不能由上述稳定 reciprocal Gaussian elimination 自动得到，所以不能把 \(\kappa\) 事后调到目标 dark-energy abundance 后再称为预测。

## 5. Hard global budget：可以自限，但 abundance 可能只是换了名字

硬约束可以通过 Lagrange multiplier 完全阻止 runaway。关键问题是约束的 conserved quantity 从哪里来。

对 AInstein 的对称三分支 benchmark：

\[
\chi_D=0,\qquad
\chi_C=+\chi_*,\qquad
\chi_V=-\chi_*.
\]

设体积分数为 \(f_D,f_C,f_V\)，和为 1。

branch-odd global moment：

\[
m_1\equiv\langle\chi\rangle
=\chi_*(f_C-f_V).
\]

所以固定 \(m_1\) **只能固定 C/V imbalance**。总转换分数

\[
f_{\rm conv}=f_C+f_V
\]

仍可在

\[
\boxed{|m_1|/\chi_*\le f_{\rm conv}\le1}
\]

之间变化。

branch-even moment：

\[
m_2\equiv\langle\chi^2\rangle
=\chi_*^2(f_C+f_V)
\]

则会精确固定

\[
\boxed{f_{\rm conv}=m_2/\chi_*^2.}
\]

这说明 global constraint 的逻辑非常清楚：**要固定总 abundance，需要一个 branch-even conserved budget。**

但 v1.6 的 one-sector EFT 真正恒等保存的是 material current

\[
\nabla_\mu J^\mu\equiv0,
\qquad b=\sqrt{-J^2},
\]

而不是 \(\chi\) 或 \(\chi^2\)。目前没有推导

\[
\int\chi^2\,dV=\mathrm{const}.
\]

因此直接把 \(\langle\chi^2\rangle\) 固定成一个“conversion capacity”会是新假设。

一个未来 microscopic completion 仍可能通过 phase-specific carrier densities、global charge neutrality、或真正的 relativistic constraint 把 conserved carrier \(b\) 映射成一个有效 even budget。若能做到，lever-rule abundance 才会成为预测；当前版本尚未做到。

## 6. 与 Buchert volume partition 的关系

已有 multiscale backreaction work 明确要求 regional domains 组合后满足 global partition consistency，volume fractions 也参与 global backreaction。这支持 AInstein 不应把每个 basin 的 environment 当作独立 external field。

但要严格区分：

\[
\boxed{
\text{volume-weighted consistency identity}
\neq
\text{conserved phase-order budget}.}
\]

Buchert partition 可以约束 regional expansion/backreaction/curvature 如何组合成 global averages，却不会自动给出 \(\int\chi^2\) 或等价 converted-volume conservation。因此它提供了 global-coupling kinematics，但还没有提供 abundance-fixing thermodynamics。

## 7. 本轮得到的新 obstruction 与剩余路线

v2.0 的结论是：

\[
\text{local passive reciprocal Gaussian feedback}
\Rightarrow
\text{softening}.
\]

v2.2 把它推广为：

\[
\boxed{
\text{stable reciprocal Gaussian mediator, even if shared/nonlocal}
\Rightarrow
\text{negative-semidefinite induced kernel in control-observable space}.}
\]

所以“把 local control 改成 global background”本身不是 rescue。

而 hard global constraint 的结论是：

\[
\boxed{
\text{global budget can self-limit abundance, but does not predict it unless the budget is independently derived}.}
\]

因此剩余真正可物理化的路线被进一步缩成：

1. **derived non-Gaussian / explicitly repulsive global response**，其正 competition kernel 必须来自 microscopic/covariant physics，而不是 abundance fit；
2. **derived conserved even budget**，必须从 one-sector carrier 或 relativistic constraint 推出；
3. **genuinely nonequilibrium nonreciprocal feedback**，进入 SK / large-deviation sector；
4. higher-order local saturation 仍可存在，但同样必须由 response data 定义。

其中 1--3 已经明显超出“投稿前再修一小轮”的尺度。

## 8. 是否应该继续研究，还是收束投 PRD？

本轮之后，我建议**收束**。

理由不是理论已经完成，而是当前论文已经形成一条完整的、可审稿的科学逻辑：

- exact/conditional positive results：regional state identities、leaf active-charge theorem、one-sector existence scaffold、response/1PI identities；
- successive obstructions：static scaling、Gaussian/ZA/PM occupancy、selector ceiling、odd-only over-conversion、local reciprocal sign restriction；
- v2.2 新增：**shared/nonlocal reciprocal Gaussian sign obstruction + global-budget identifiability result**；
- open completion 被压缩到几类明确、不可再通过低成本调参绕过的物理机制。

继续做 nonreciprocal SK quasipotential、完整 relativistic budget 或 covariant phase-to-geometry backreaction，会变成真正的新研究项目，不再是当前 PRD 稿的合理 pre-submission extension。

因此 v2.2 的最佳定位是：

\[
\boxed{\textbf{最后一次科学收束版，而不是新一轮无限扩展。}}
\]

## 9. 本轮代码与 QA

新增：

- `research_v2.2/src/analysis_global_competition.py`
- `research_v2.2/tests/test_global_competition.py`
- `research_v2.2/results/global_competition/global_competition_summary.json`
- `research_v2.2/results/global_competition/global_competition_scan.csv`
- `research_v2.2/results/global_competition/global_competition_seed_summary.csv`
- `research_v2.2/results/global_competition/fig_global_competition.pdf/png`

专项测试：**5 passed**。

本轮所有数值均用于验证解析 sign/identifiability 结构，不被解释为 cosmological abundance fit。

## References added for this audit

- J. Rubinstein and P. Sternberg, “Nonlocal reaction--diffusion equations and nucleation,” *IMA Journal of Applied Mathematics* **48**, 249--264 (1992), doi:10.1093/imamat/48.3.249.
- T. Ohta and K. Kawasaki, “Equilibrium morphology of block copolymer melts,” *Macromolecules* **19**, 2621--2632 (1986), doi:10.1021/ma00164a028.
- A. Wiegand and T. Buchert, “Multiscale cosmology and structure-emerging dark energy: a plausibility analysis,” *Physical Review D* **82**, 023523 (2010), doi:10.1103/PhysRevD.82.023523.
