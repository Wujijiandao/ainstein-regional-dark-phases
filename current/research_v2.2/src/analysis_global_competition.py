from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research_v2.2' / 'results' / 'global_competition'
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1) Exact elimination theorem for a shared/nonlocal passive Gaussian control.
#
# F(q,y) = F0(q) + 1/2 (y-y0)^T K (y-y0) - y^T B q,
# K = K^T > 0.  Eliminating y gives
# F_eff(q) = F0(q) - y0^T B q - 1/2 q^T B^T K^{-1} B q.
# Hence the induced q-space kernel is -G, G=B^T K^{-1} B >= 0.
# A stable real reciprocal Gaussian mediator therefore supplies cooperative
# softening/alignment, not a positive-semidefinite competition penalty.
# ---------------------------------------------------------------------------

rng = np.random.default_rng(2200)
identity_errors = []
max_induced_eigs = []
min_G_eigs = []
for _ in range(100):
    ny = int(rng.integers(2, 7))
    nq = int(rng.integers(2, 9))
    A = rng.normal(size=(ny, ny))
    K = A.T @ A + 0.5*np.eye(ny)
    B = rng.normal(size=(ny, nq))
    q = rng.normal(size=nq)
    y0 = rng.normal(size=ny)
    ystar = y0 + np.linalg.solve(K, B @ q)
    F_controls = 0.5*(ystar-y0) @ K @ (ystar-y0) - ystar @ B @ q
    G = B.T @ np.linalg.solve(K, B)
    F_reduced = -y0 @ B @ q - 0.5*q @ G @ q
    identity_errors.append(abs(F_controls-F_reduced))
    min_G_eigs.append(float(np.min(np.linalg.eigvalsh(G))))
    max_induced_eigs.append(float(np.max(np.linalg.eigvalsh(-G))))

identity_max_error = float(max(identity_errors))
min_G_eigenvalue = float(min(min_G_eigs))
max_induced_kernel_eigenvalue = float(max(max_induced_eigs))

# Rank-one examples directly relevant to the AInstein branch variables.
# q_i = chi_i: induced term favours coherent daughter alignment.
# q_i = chi_i^2: induced term favours larger collective conversion amplitude.
n_demo = 8
g = np.linspace(0.5, 1.2, n_demo)
G_rank1 = np.outer(g, g)
rank1_nonzero_eigenvalue = float(np.max(np.linalg.eigvalsh(G_rank1)))
rank1_induced_max_eigenvalue = float(np.max(np.linalg.eigvalsh(-G_rank1)))

# ---------------------------------------------------------------------------
# 2) Exact global-budget identifiability for D/C/V = 0,+chi*,-chi*.
#
# m1=<chi>=chi*(f_C-f_V): fixes daughter imbalance only.
# m2=<chi^2>=chi^2*(f_C+f_V): fixes total converted fraction.
# The current one-sector scaffold conserves the material current b, not chi or
# chi^2, so an m2 budget is an additional constitutive/microscopic input unless
# a completion derives a map from the conserved carrier to this phase moment.
# ---------------------------------------------------------------------------

chi_star = 1.0
m1 = 0.20
# For fixed m1, f_tot can range from |m1|/chi* to 1 subject to nonnegative fractions.
f_total_min_from_m1 = abs(m1)/chi_star
f_total_max_from_m1 = 1.0
m2 = 0.40
f_total_from_m2 = m2/(chi_star**2)
fC_from_m1m2 = 0.5*(f_total_from_m2 + m1/chi_star)
fV_from_m1m2 = 0.5*(f_total_from_m2 - m1/chi_star)
fD_from_m1m2 = 1.0 - f_total_from_m2

# ---------------------------------------------------------------------------
# 3) Discrete branch audit: passive Gaussian global coupling is cooperative.
#
# s_i in {0,1} labels unconverted/converted for a branch-reduced audit.
# With equal weights and local conversion costs Delta_i,
# E(k;gamma)=sum_{j<=k} Delta_(j)/N - gamma/2*(k/N)^2.
# Because slopes in gamma become more negative as k increases, the globally
# minimizing converted fraction is nondecreasing with gamma.  The audit below
# checks this on frozen heterogeneous ensembles and contrasts a positive
# competition penalty +kappa f^2/2, which can self-limit but is a distinct
# physical ingredient rather than the Gaussian-elimination result.
# ---------------------------------------------------------------------------

def optimal_fraction(costs: np.ndarray, strength: float, sign: int) -> tuple[float, float]:
    """sign=-1 passive Gaussian attraction; sign=+1 explicit repulsive competition."""
    costs = np.sort(np.asarray(costs, float))
    N = len(costs)
    k = np.arange(N+1)
    f = k/N
    cumulative = np.r_[0.0, np.cumsum(costs)]/N
    E = cumulative + sign*0.5*strength*f*f
    j = int(np.argmin(E))
    return float(f[j]), float(E[j])

seeds = [11, 23, 47, 91]
gammas = np.linspace(0.0, 1.8, 91)
kappas = np.linspace(0.0, 1.8, 91)
rows = []
monotone_flags = []
for seed in seeds:
    rr = np.random.default_rng(seed)
    costs = rr.normal(loc=0.15, scale=0.35, size=400)
    f_passive = []
    f_repulsive = []
    for gam, kap in zip(gammas, kappas):
        fp, Ep = optimal_fraction(costs, gam, -1)
        fr, Er = optimal_fraction(costs, kap, +1)
        f_passive.append(fp)
        f_repulsive.append(fr)
        rows.append({'seed':seed, 'strength':float(gam), 'passive_gaussian_fraction':fp,
                     'explicit_repulsive_fraction':fr, 'passive_energy':Ep, 'repulsive_energy':Er})
    monotone_flags.append(bool(np.all(np.diff(f_passive) >= -1e-15)))

scan = pd.DataFrame(rows)
scan.to_csv(OUT/'global_competition_scan.csv', index=False)
passive_monotone_all = bool(all(monotone_flags))

summary_by_seed = scan.groupby('seed').agg(
    baseline_fraction=('passive_gaussian_fraction','first'),
    passive_final_fraction=('passive_gaussian_fraction','last'),
    repulsive_final_fraction=('explicit_repulsive_fraction','last')
).reset_index()
summary_by_seed.to_csv(OUT/'global_competition_seed_summary.csv', index=False)

# ---------------------------------------------------------------------------
# 4) A simple volume-partition consistency audit.
# A volume-weighted global mean of an additive regional observable is fixed by
# the regional fractions, but this is a bookkeeping identity rather than a new
# conserved phase charge.  For three states, normalization plus an odd moment
# still leaves total conversion undetermined.
# ---------------------------------------------------------------------------

fgrid = np.linspace(f_total_min_from_m1, 1.0, 161)
fC_grid = 0.5*(fgrid + m1/chi_star)
fV_grid = 0.5*(fgrid - m1/chi_star)
fD_grid = 1.0-fgrid
simplex_ok = bool(np.all((fC_grid >= -1e-15) & (fV_grid >= -1e-15) & (fD_grid >= -1e-15)))
odd_constraint_error = float(np.max(np.abs(chi_star*(fC_grid-fV_grid)-m1)))

# Figures
fig, ax = plt.subplots(figsize=(7.2, 4.6))
for seed in seeds:
    ss = scan[scan.seed == seed]
    ax.plot(ss.strength, ss.passive_gaussian_fraction, label=f'seed {seed}')
ax.set_xlabel(r'passive shared-control strength $\gamma$')
ax.set_ylabel('globally minimizing converted fraction')
ax.set_title('Passive Gaussian global feedback is cooperative, not self-limiting')
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_passive_global_cooperation.pdf')
fig.savefig(OUT/'fig_passive_global_cooperation.png', dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 4.6))
for seed in seeds:
    ss = scan[scan.seed == seed]
    ax.plot(ss.strength, ss.explicit_repulsive_fraction, label=f'seed {seed}')
ax.set_xlabel(r'explicit competition strength $\kappa$')
ax.set_ylabel('globally minimizing converted fraction')
ax.set_title('A positive nonlocal competition term can self-limit, but is extra physics')
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_explicit_global_competition.pdf')
fig.savefig(OUT/'fig_explicit_global_competition.png', dpi=180)
plt.close(fig)

summary = {
    'gaussian_elimination_identity_max_error': identity_max_error,
    'minimum_eigenvalue_of_G_should_be_nonnegative': min_G_eigenvalue,
    'maximum_eigenvalue_of_induced_minus_G_should_be_nonpositive': max_induced_kernel_eigenvalue,
    'rank1_G_nonzero_eigenvalue': rank1_nonzero_eigenvalue,
    'rank1_induced_max_eigenvalue': rank1_induced_max_eigenvalue,
    'odd_budget_m1': m1,
    'odd_budget_allowed_total_fraction_min': f_total_min_from_m1,
    'odd_budget_allowed_total_fraction_max': f_total_max_from_m1,
    'even_budget_m2': m2,
    'even_budget_fixed_total_fraction': f_total_from_m2,
    'm1_m2_example_fD': fD_from_m1m2,
    'm1_m2_example_fC': fC_from_m1m2,
    'm1_m2_example_fV': fV_from_m1m2,
    'odd_constraint_simplex_ok': simplex_ok,
    'odd_constraint_max_error': odd_constraint_error,
    'passive_fraction_nondecreasing_all_frozen_seeds': passive_monotone_all,
    'interpretation': (
        'Stable reciprocal Gaussian shared controls induce a negative-semidefinite q-space kernel and therefore '
        'do not generate positive long-range competition. Exact global constraints can self-limit phase volume, '
        'but the constrained budget is an input unless independently derived. In the symmetric D/C/V benchmark, '
        'a conserved odd moment fixes C-V imbalance but not total conversion; an even moment would fix total '
        'conversion, but no existing AInstein conservation law fixes chi^2.'
    )
}
(OUT/'global_competition_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')

print(json.dumps(summary, indent=2))
print('\nseed summary:')
print(summary_by_seed.to_string(index=False))

# Combined manuscript figure.
fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.2))
for seed in seeds:
    ss = scan[scan.seed == seed]
    axs[0].plot(ss.strength, ss.passive_gaussian_fraction, label=f'seed {seed}')
    axs[1].plot(ss.strength, ss.explicit_repulsive_fraction, label=f'seed {seed}')
axs[0].set_xlabel(r'passive shared-control strength $\gamma$')
axs[0].set_ylabel('globally minimizing converted fraction')
axs[0].set_title('(a) Gaussian elimination: cooperative sign')
axs[1].set_xlabel(r'explicit competition strength $\kappa$')
axs[1].set_ylabel('globally minimizing converted fraction')
axs[1].set_title('(b) Positive competition: distinct ingredient')
axs[1].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_global_competition.pdf')
fig.savefig(OUT/'fig_global_competition.png', dpi=180)
plt.close(fig)
