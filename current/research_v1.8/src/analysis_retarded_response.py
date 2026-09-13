from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research_v1.8' / 'results' / 'retarded_response'
OUT.mkdir(parents=True, exist_ok=True)

# Deterministic internal benchmark. These are NOT cosmological fits.
M2 = 4.0
KAPPA = 0.02
LAMBDA = 0.5
TAU = 0.02
GAMMA = 1.0 / LAMBDA           # coefficient of -i omega in inverse response
INERTIA = TAU / LAMBDA         # coefficient of -omega^2 in inverse response
GX = 0.7
TAUX = 0.4
LAMBDA_E = 1.2
RNG = np.random.default_rng(1708)


def G_chi_h(omega, k=0.0, m2=M2, kappa=KAPPA, lam=LAMBDA, tau=TAU):
    """Retarded response of chi to its canonically-normalized conjugate field h.

    Fourier convention: exp(-i omega t).
    Equation: tau chi_ddot + chi_dot + lam (m2 + kappa k^2) chi = lam h.
    """
    return 1.0 / (m2 + kappa*k*k - (tau/lam)*omega*omega - 1j*omega/lam)


def K_x(omega, gx=GX, taux=TAUX):
    """Representative causal branch-odd cross kernel, used only for identification audit."""
    return gx / (1.0 - 1j*omega*taux)


def G_chi_x(omega, k=0.0):
    return G_chi_h(omega, k) * K_x(omega)


def G_gate(omega, k, E):
    """Response to h in a frozen even environment E shifting the quadratic coefficient."""
    return G_chi_h(omega, k, m2=M2 + LAMBDA_E*E)


# 1) Spectral positivity and Kramers-Kronig static-susceptibility check.
omega = np.linspace(0.0, 40.0, 4001)
rows = []
for k in [0.0, 2.0, 5.0]:
    G = G_chi_h(omega, k)
    for om, val in zip(omega, G):
        rows.append({'omega': float(om), 'k': float(k), 'ReG': float(val.real), 'ImG': float(val.imag)})
spec = pd.DataFrame(rows)
spec.to_csv(OUT / 'retarded_spectrum.csv', index=False)

# For exp(-i omega t), Im G^R >= 0 at omega>0 for this passive oscillator convention.
min_im_positive = float(spec[(spec.omega > 0) & (spec.k == 0.0)].ImG.min())

def kk_integrand(w):
    return float(G_chi_h(w, 0.0).imag / w) if w > 0 else 1.0/(LAMBDA*M2*M2)

kk_val, kk_err = quad(kk_integrand, 0.0, np.inf, epsabs=1e-11, epsrel=1e-10, limit=1000)
kk_static = 2.0/np.pi * kk_val
static_exact = float(G_chi_h(0.0, 0.0).real)
kk_relerr = abs(kk_static-static_exact)/static_exact

# 2) Low-frequency inverse-response identification with small deterministic synthetic noise.
# Re G^{-1} = m2 + kappa k^2 - M omega^2; Im G^{-1} = -gamma omega.
omegas = np.linspace(0.04, 1.0, 30)
ks = np.array([0.0, 0.8, 1.6, 2.4])
obs = []
noise_level = 2.0e-3
for k in ks:
    for w in omegas:
        G = G_chi_h(w, k)
        # multiplicative complex perturbation, fixed seed, for an identification stress test
        eps = noise_level*(RNG.normal() + 1j*RNG.normal())
        Gn = G*(1.0 + eps)
        inv = 1.0/Gn
        obs.append({'omega':w,'k':k,'ReInv':inv.real,'ImInv':inv.imag})
obs = pd.DataFrame(obs)
obs.to_csv(OUT/'inverse_response_synthetic.csv', index=False)
A = np.column_stack([np.ones(len(obs)), obs.k.to_numpy()**2, -obs.omega.to_numpy()**2])
coef, *_ = np.linalg.lstsq(A, obs.ReInv.to_numpy(), rcond=None)
m2_hat, kappa_hat, inertia_hat = coef
# Im inv = -gamma omega, no intercept by analyticity at omega=0
gamma_hat = -float(np.dot(obs.omega, obs.ImInv)/np.dot(obs.omega, obs.omega))
lambda_hat = 1.0/gamma_hat
tau_hat = inertia_hat/gamma_hat
c2_hat = kappa_hat/inertia_hat

# 3) Branch-odd cross kernel can be factored out of the phase propagator.
# K_X^R = G_chiX^R / G_chih^R. For the exponential audit, 1/K = 1/g - i omega tauX/g.
wxs = np.linspace(0.0, 1.0, 81)
ratio = G_chi_x(wxs,0.0)/G_chi_h(wxs,0.0)
invK = 1.0/ratio
intercept = float(np.mean(invK.real))
# least-squares imaginary slope through origin
slope_im = float(np.dot(wxs, invK.imag)/np.dot(wxs, wxs))
gx_hat = 1.0/intercept
taux_hat = -slope_im/intercept
pd.DataFrame({'omega':wxs,'ReK':ratio.real,'ImK':ratio.imag,'absK':np.abs(ratio),'phaseK':np.angle(ratio)}).to_csv(OUT/'cross_kernel.csv', index=False)

# 4) Even gate is invisible as a linear odd source at the symmetric parent, but measurable as
# a shift of the inverse odd susceptibility.  lambda_E = d_E Re [G_chih(0,0;E)]^{-1}.
Es = np.linspace(-0.2,0.2,21)
inv_static = np.array([(1.0/G_gate(0.0,0.0,E)).real for E in Es])
p = np.polyfit(Es, inv_static, 1)
lambdaE_hat, m2_gate_hat = map(float,p)
pd.DataFrame({'E':Es,'inverse_static_susceptibility':inv_static}).to_csv(OUT/'even_gate_susceptibility.csv', index=False)

# 5) A compact response dictionary that future microscopic simulations can target.
summary = {
    'benchmark_input': {
        'm2':M2,'kappa':KAPPA,'Lambda':LAMBDA,'tau':TAU,'gamma_inv_Lambda':GAMMA,
        'inertia_tau_over_Lambda':INERTIA,'c2':LAMBDA*KAPPA/TAU,
        'gX':GX,'tauX':TAUX,'lambda_even_gate':LAMBDA_E,
        'synthetic_complex_noise_fraction':noise_level,
    },
    'spectral_checks': {
        'minimum_Im_G_for_positive_omega_k0': min_im_positive,
        'static_susceptibility_exact': static_exact,
        'static_susceptibility_KK': kk_static,
        'KK_quad_error': kk_err,
        'KK_relative_error': kk_relerr,
    },
    'inverse_response_identification': {
        'm2_hat':float(m2_hat),'kappa_hat':float(kappa_hat),'gamma_hat':float(gamma_hat),
        'Lambda_hat':float(lambda_hat),'inertia_hat':float(inertia_hat),'tau_hat':float(tau_hat),
        'c2_hat':float(c2_hat),
    },
    'cross_response_identification': {'gX_hat':float(gx_hat),'tauX_hat':float(taux_hat)},
    'even_gate_identification': {'m2_intercept_hat':float(m2_gate_hat),'lambdaE_hat':float(lambdaE_hat)},
    'structural_selection_rules': {
        'linear_odd_response': 'G_chiX allowed because chi and X are branch-odd',
        'linear_even_to_odd_response_at_symmetric_parent': 'G_chiE = 0 by branch parity',
        'even_gate_measurement': 'lambda_E = d_E Re[(G_chih^R(0,0;E))^{-1}]',
        'odd_kernel_measurement': 'K_X^R = G_chiX^R / G_chih^R',
    }
}
(OUT/'retarded_response_summary.json').write_text(json.dumps(summary,indent=2))

# Figures
fig, ax = plt.subplots(figsize=(7.0,4.3))
for k in [0.0,2.0,5.0]:
    s=spec[spec.k==k]
    ax.plot(s.omega, s.ImG, label=fr'$k={k:g}$')
ax.set_xlim(0,12)
ax.set_xlabel(r'frequency $\omega$ [internal units]')
ax.set_ylabel(r'$\mathrm{Im}\,G^R_{\chi h}$')
ax.set_title('Passive retarded phase response has positive spectral weight')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT/'fig_retarded_spectral_response.pdf')
fig.savefig(OUT/'fig_retarded_spectral_response.png',dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.0,4.3))
ax.plot(wxs, ratio.real, label=r'$\mathrm{Re}\,K_X^R$')
ax.plot(wxs, ratio.imag, label=r'$\mathrm{Im}\,K_X^R$')
ax.set_xlabel(r'frequency $\omega$ [internal units]')
ax.set_ylabel(r'$K_X^R=G^R_{\chi X}/G^R_{\chi h}$')
ax.set_title('Branch-odd cross kernel factors out of the phase propagator')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT/'fig_cross_kernel_identification.pdf')
fig.savefig(OUT/'fig_cross_kernel_identification.png',dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.6,4.2))
ax.plot(Es, inv_static, 'o', label='synthetic response data')
ax.plot(Es, lambdaE_hat*Es+m2_gate_hat, label='linear inverse-susceptibility fit')
ax.set_xlabel(r'branch-even control $\mathcal{E}$ [internal units]')
ax.set_ylabel(r'$[G^R_{\chi h}(0,0;\mathcal{E})]^{-1}$')
ax.set_title('Even gate is measured through the inverse odd susceptibility')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT/'fig_even_gate_kubo.pdf')
fig.savefig(OUT/'fig_even_gate_kubo.png',dpi=180)
plt.close(fig)

print(json.dumps(summary,indent=2))

# Compact three-panel figure for the manuscript.
fig, axs = plt.subplots(1,3,figsize=(12.6,3.7))
# (a) spectral
for k in [0.0,2.0,5.0]:
    s=spec[spec.k==k]
    axs[0].plot(s.omega,s.ImG,label=fr'$k={k:g}$')
axs[0].set_xlim(0,12)
axs[0].set_xlabel(r'$\omega$')
axs[0].set_ylabel(r'$\mathrm{Im}\,G^R_{\chi h}$')
axs[0].set_title('(a) passive spectrum')
axs[0].legend(frameon=False,fontsize=7)
# (b) odd cross kernel
axs[1].plot(wxs,ratio.real,label='Re')
axs[1].plot(wxs,ratio.imag,label='Im')
axs[1].set_xlabel(r'$\omega$')
axs[1].set_ylabel(r'$K_X^R$')
axs[1].set_title('(b) odd cross kernel')
axs[1].legend(frameon=False,fontsize=7)
# (c) even gate
axs[2].plot(Es,inv_static,'o',ms=3,label='response')
axs[2].plot(Es,lambdaE_hat*Es+m2_gate_hat,label='fit')
axs[2].set_xlabel(r'$\mathcal{E}$')
axs[2].set_ylabel(r'$[G^R_{\chi h}(0)]^{-1}$')
axs[2].set_title('(c) even gate')
axs[2].legend(frameon=False,fontsize=7)
fig.tight_layout()
fig.savefig(OUT/'fig_kubo_response_dictionary.pdf')
fig.savefig(OUT/'fig_kubo_response_dictionary.png',dpi=180)
plt.close(fig)
