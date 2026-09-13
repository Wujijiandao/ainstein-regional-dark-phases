from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research_v1.7' / 'results' / 'causal_kinetics'
OUT.mkdir(parents=True, exist_ok=True)

H_COEX = 0.341688495018
H_SPIN = 0.928501997833
TAUS = np.array([0.0, 0.1, 0.3, 1.0, 3.0])  # e-fold memory times for audit only
GGRID = np.geomspace(0.3, 300.0, 241)

DATASETS = [
    (ROOT / 'research_v1.4' / 'results' / 'descendant_transport_n64_4seed_combined.csv', '64^3'),
    (ROOT / 'research_v1.4' / 'results' / 'descendant_transport_n80_s11.csv', '80^3'),
]


def pivot_history(path: Path):
    df = pd.read_csv(path)
    for seed, g in df.groupby('seed'):
        aa = np.array(sorted(g['a'].unique()), float)
        if not np.allclose(aa, [0.25, 0.50, 0.75, 1.00]):
            raise RuntimeError(f'unexpected snapshots: {aa}')
        p = g.pivot(index='final_basin_id', columns='a', values=['volume_fraction', 'Hloc_over_H'])
        w = p['volume_fraction'][1.0].to_numpy(float)
        w = w / np.sum(w)
        X = np.column_stack([p['Hloc_over_H'][a].to_numpy(float) - 1.0 for a in aa])
        yield int(seed), aa, w, X


def exp_memory_final(X: np.ndarray, aa: np.ndarray, tauN: float) -> np.ndarray:
    """Normalized passive causal exponential memory of X at a=1.

    This is a sensitivity kernel, not the claimed microscopic response law.
    tauN is measured in e-folds of N=ln a.
    """
    if tauN <= 0:
        return X[:, -1].copy()
    N = np.log(aa)
    grid = np.linspace(N[0], N[-1], 1201)
    Xi = np.vstack([np.interp(grid, N, row) for row in X])
    kw = np.exp(-(N[-1] - grid) / tauN)
    return np.trapezoid(Xi * kw, grid, axis=1) / np.trapezoid(kw, grid)


def weighted_req_coupling(M: np.ndarray, w: np.ndarray, target: float, threshold: float, sign: int) -> float:
    y = sign * M
    keep = y > 0
    if np.sum(w[keep]) + 1e-14 < target:
        return float('inf')
    req = threshold / y[keep]
    ww = w[keep]
    order = np.argsort(req)
    c = np.cumsum(ww[order])
    j = int(np.searchsorted(c, target, side='left'))
    return float(req[order[min(j, len(order)-1)]])


# --- PM forcing audit -------------------------------------------------------
trichotomy_rows = []
memory_rows = []
scan_rows = []
invert_rows = []
for path, label in DATASETS:
    for seed, aa, w, X in pivot_history(path):
        all_pos = np.all(X > 0, axis=1)
        all_neg = np.all(X < 0, axis=1)
        mixed = ~(all_pos | all_neg)
        trichotomy_rows.append({
            'dataset': label,
            'seed': seed,
            'persistent_positive_fraction': float(np.sum(w[all_pos])),
            'persistent_negative_fraction': float(np.sum(w[all_neg])),
            'sign_mixed_fraction': float(np.sum(w[mixed])),
            'persistent_one_sign_fraction': float(np.sum(w[all_pos | all_neg])),
        })
        for tau in TAUS:
            M = exp_memory_final(X, aa, float(tau))
            memory_rows.append({
                'dataset': label, 'seed': seed, 'tauN': float(tau),
                'memory_positive_fraction': float(np.sum(w[M > 0])),
                'memory_negative_fraction': float(np.sum(w[M < 0])),
                'memory_zero_fraction': float(np.sum(w[M == 0])),
                'memory_weighted_mean': float(np.sum(w * M)),
                'memory_weighted_abs_mean': float(np.sum(w * np.abs(M))),
            })
            for target in [0.10, 0.30, 0.50]:
                for threshold_name, threshold in [('coexistence', H_COEX), ('spinodal', H_SPIN)]:
                    invert_rows.append({
                        'dataset': label, 'seed': seed, 'tauN': float(tau),
                        'target_fraction': target, 'threshold': threshold_name,
                        'V_required_gX': weighted_req_coupling(M, w, target, threshold, +1),
                        'C_required_gX': weighted_req_coupling(M, w, target, threshold, -1),
                    })
            for gx in GGRID:
                scan_rows.append({
                    'dataset': label, 'seed': seed, 'tauN': float(tau), 'gX': float(gx),
                    'V_coexistence_fraction': float(np.sum(w[gx * M >= H_COEX])),
                    'C_coexistence_fraction': float(np.sum(w[-gx * M >= H_COEX])),
                    'V_spinodal_fraction': float(np.sum(w[gx * M >= H_SPIN])),
                    'C_spinodal_fraction': float(np.sum(w[-gx * M >= H_SPIN])),
                })

tri = pd.DataFrame(trichotomy_rows)
mem = pd.DataFrame(memory_rows)
scan = pd.DataFrame(scan_rows)
inv = pd.DataFrame(invert_rows)
tri.to_csv(OUT / 'persistent_sign_trichotomy.csv', index=False)
mem.to_csv(OUT / 'passive_memory_sign_audit.csv', index=False)
scan.to_csv(OUT / 'driver_coupling_scan.csv', index=False)
inv.to_csv(OUT / 'driver_coupling_inversion.csv', index=False)

# --- Causal relaxational scaffold checks -----------------------------------
# Equations: D_t chi = pi; tau D_t pi + pi = -Lambda * mu,
# mu = dW/dchi for a homogeneous static test.  The extended Lyapunov function
# Lext = W + tau*pi^2/(2 Lambda) obeys dLext/dt = -pi^2/Lambda.

def W0(x):
    return (4.0/3.0) * x**6 - 3.0*x**4 + 2.0*x**2


def dW0(x):
    return 8.0*x**5 - 12.0*x**3 + 4.0*x


def causal_rhs(state, tau=0.03, lam=1.0, H=0.12):
    chi, pi = state
    mu = dW0(chi) - H
    return np.array([pi, (-pi - lam*mu)/tau])


def rk4(y, dt, n, **kwargs):
    ys = np.empty((n+1, 2), float); ys[0] = y
    for i in range(n):
        z = ys[i]
        k1 = causal_rhs(z, **kwargs)
        k2 = causal_rhs(z + 0.5*dt*k1, **kwargs)
        k3 = causal_rhs(z + 0.5*dt*k2, **kwargs)
        k4 = causal_rhs(z + dt*k3, **kwargs)
        ys[i+1] = z + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0
    return ys

tau_test = 0.03
lam_test = 1.0
H_test = 0.12
dt = 2.0e-4
nstep = int(4.0/dt)
ys = rk4(np.array([0.0, 0.0]), dt, nstep, tau=tau_test, lam=lam_test, H=H_test)
t = np.arange(nstep+1)*dt
chi = ys[:,0]; pi = ys[:,1]
F = W0(chi) - H_test*chi
Lext = F + tau_test*pi*pi/(2.0*lam_test)
max_up = float(np.max(np.diff(Lext)))
# Integral balance check.
diss_int = np.concatenate([[0.0], np.cumsum(0.5*(pi[1:]**2 + pi[:-1]**2)/lam_test*dt)])
balance = Lext + diss_int
balance_error = float(np.max(np.abs(balance - balance[0])))

# Linearized dispersion and sign-preserving Green function in monotone causal sector.
m2 = 4.0  # W0''(0)
kappa = 0.02
lam = 0.5
tau = 0.02
c2 = lam*kappa/tau
disc0 = 1.0 - 4.0*tau*lam*m2
if disc0 < 0:
    raise RuntimeError('chosen audit point is not overdamped')
rp = (-1.0 + np.sqrt(disc0))/(2.0*tau)
rm = (-1.0 - np.sqrt(disc0))/(2.0*tau)
tg = np.linspace(1e-6, 5.0, 2000)
Green = lam/(tau*(rp-rm))*(np.exp(rp*tg)-np.exp(rm*tg))
green_min = float(np.min(Green))

# verify all Fourier modes in a representative k-range are stable; the asymptotic speed is sqrt(lam*kappa/tau)
max_real_root = -np.inf
for k in np.linspace(0, 100, 501):
    coeff = [tau, 1.0, lam*(m2+kappa*k*k)]
    roots = np.roots(coeff)
    max_real_root = max(max_real_root, float(np.max(np.real(roots))))

checks = {
    'H_coexistence': H_COEX,
    'H_spinodal': H_SPIN,
    'causal_relaxation_test': {
        'tau': tau_test, 'Lambda': lam_test, 'H': H_test,
        'max_single_step_extended_energy_increase': max_up,
        'integrated_energy_balance_max_error': balance_error,
    },
    'linearized_causal_test': {
        'm2': m2, 'kappa': kappa, 'Lambda': lam, 'tau': tau,
        'c_phase_squared': c2,
        'overdamped_discriminant_k0': disc0,
        'green_function_minimum': green_min,
        'max_real_part_of_roots': max_real_root,
    },
    'trichotomy_records': tri.to_dict(orient='records'),
}
(OUT / 'causal_kinetics_summary.json').write_text(json.dumps(checks, indent=2))

# --- Figures ---------------------------------------------------------------
# 1) robust trichotomy
plot_tri = tri.copy()
plot_tri['run'] = plot_tri.apply(lambda r: f"{r['dataset']} s{int(r['seed'])}", axis=1)
x = np.arange(len(plot_tri))
fig, ax = plt.subplots(figsize=(7.4,4.2))
bottom = np.zeros(len(plot_tri))
for col, lab in [('persistent_positive_fraction','persistent X>0'),('persistent_negative_fraction','persistent X<0'),('sign_mixed_fraction','sign-mixed')]:
    vals = plot_tri[col].to_numpy(float)
    ax.bar(x, vals, bottom=bottom, label=lab)
    bottom += vals
ax.set_xticks(x, plot_tri['run'], rotation=30, ha='right')
ax.set_ylabel('final-volume fraction')
ax.set_ylim(0,1.02)
ax.set_title('Persistent branch-odd sign trichotomy')
ax.legend(frameon=False, ncol=3, fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_causal_sign_trichotomy.pdf')
fig.savefig(OUT/'fig_causal_sign_trichotomy.png', dpi=180)
plt.close(fig)

# 2) passive memory sign robustness for 64^3 seeds
m64=mem[mem.dataset=='64^3']
fig,ax=plt.subplots(figsize=(6.8,4.2))
for seed,g in m64.groupby('seed'):
    ax.plot(g.tauN, g.memory_positive_fraction, marker='o', label=f'seed {int(seed)}')
ax.axhspan(tri[tri.dataset=='64^3'].persistent_positive_fraction.min(), tri[tri.dataset=='64^3'].persistent_positive_fraction.max(), alpha=0.12, label='kernel-independent persistent + support')
ax.set_xscale('symlog', linthresh=0.08)
ax.set_xlabel(r'exponential memory time $\tau_N$ [e-folds]')
ax.set_ylabel('positive final memory volume fraction')
ax.set_ylim(0.40,0.62)
ax.set_title('Branch sign is weakly sensitive to passive memory time')
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_passive_memory_sign.pdf')
fig.savefig(OUT/'fig_passive_memory_sign.png', dpi=180)
plt.close(fig)

# 3) coexistence/spinodal coupling scan, 64^3 min-max band, tau=0.3
s=scan[(scan.dataset=='64^3') & np.isclose(scan.tauN,0.3)]
agg=s.groupby('gX').agg(
    Vco_min=('V_coexistence_fraction','min'), Vco_max=('V_coexistence_fraction','max'), Vco_mean=('V_coexistence_fraction','mean'),
    Cco_min=('C_coexistence_fraction','min'), Cco_max=('C_coexistence_fraction','max'), Cco_mean=('C_coexistence_fraction','mean'),
    Vsp_mean=('V_spinodal_fraction','mean'), Csp_mean=('C_spinodal_fraction','mean')
).reset_index()
fig,ax=plt.subplots(figsize=(7.0,4.4))
ax.fill_between(agg.gX,agg.Vco_min,agg.Vco_max,alpha=0.15)
ax.plot(agg.gX,agg.Vco_mean,label='V: coexistence')
ax.fill_between(agg.gX,agg.Cco_min,agg.Cco_max,alpha=0.15)
ax.plot(agg.gX,agg.Cco_mean,label='C: coexistence')
ax.plot(agg.gX,agg.Vsp_mean,'--',label='V: spinodal')
ax.plot(agg.gX,agg.Csp_mean,'--',label='C: spinodal')
ax.set_xscale('log')
ax.set_xlabel(r'branch-odd coupling $g_X$ (fixed sextic normalization)')
ax.set_ylabel('final-volume fraction above threshold')
ax.set_ylim(0,0.62)
ax.set_title(r'Coupling identifiability audit ($\tau_N=0.3$)')
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_driver_coupling_identifiability.pdf')
fig.savefig(OUT/'fig_driver_coupling_identifiability.png', dpi=180)
plt.close(fig)

# 4) Lyapunov decay check
fig,ax=plt.subplots(figsize=(6.6,4.0))
ax.plot(t,Lext-Lext[-1])
ax.set_xlabel('dimensionless local time')
ax.set_ylabel(r'$\mathcal{L}(t)-\mathcal{L}(\infty)$')
ax.set_yscale('log')
ax.set_title('Causal relaxational scaffold: extended Lyapunov decay')
fig.tight_layout()
fig.savefig(OUT/'fig_causal_lyapunov_decay.pdf')
fig.savefig(OUT/'fig_causal_lyapunov_decay.png', dpi=180)
plt.close(fig)

print(json.dumps(checks['causal_relaxation_test'], indent=2))
print(json.dumps(checks['linearized_causal_test'], indent=2))
print(tri.to_string(index=False))
