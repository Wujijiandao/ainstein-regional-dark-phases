from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import newton

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research_v1.9' / 'results' / 'nonlinear_1pi_sk'
OUT.mkdir(parents=True, exist_ok=True)

# Canonical 1PI/Landau normalization of the symmetric sextic used throughout AInstein:
# U = 1/2 r chi^2 + u/4! chi^4 + v/6! chi^6 - h chi.
# W0 = 2 chi^2 - 3 chi^4 + 4/3 chi^6 corresponds exactly to:
R0 = 4.0
U4 = -72.0
U6 = 960.0
KAPPA = 0.02       # internal regression value only
E0 = 1.0           # sets the internal energy-density unit
LAMBDA = 0.5       # same causal benchmark as v1.8
TAU = 0.02
THETA = 0.30       # local-KMS benchmark energy/noise scale; not cosmological inference
RNG = np.random.default_rng(1909)


def Ueff(chi, r=R0, h=0.0, u=U4, v=U6):
    return 0.5*r*chi**2 + (u/24.0)*chi**4 + (v/720.0)*chi**6 - h*chi


def dU(chi, r=R0, h=0.0, u=U4, v=U6):
    return r*chi + (u/6.0)*chi**3 + (v/120.0)*chi**5 - h


def ddU(chi, r=R0, u=U4, v=U6):
    return r + (u/2.0)*chi**2 + (v/24.0)*chi**4


def central_root(h, r=R0, u=U4, v=U6):
    return float(newton(lambda x: dU(x,r,h,u,v), h/r,
                        fprime=lambda x: ddU(x,r,u,v), tol=1e-14, maxiter=100))


def daughter_y_roots(r, u=U4, v=U6):
    disc = u*u - (6.0/5.0)*v*r
    if disc < 0:
        return None
    s = np.sqrt(max(disc, 0.0))
    y_small = 10.0*(-u - s)/v
    y_large = 10.0*(-u + s)/v
    return float(y_small), float(y_large)


# ---- 1) Exact nonlinear/1PI phase-landscape identities ----
r_daughter_spinodal = 5.0*U4**2/(6.0*U6)
r_coex = 5.0*U4**2/(8.0*U6)
chi_coex_sq = -15.0*U4/U6
chi_coex = float(np.sqrt(chi_coex_sq))
parent_to_coex_shift = R0-r_coex

# Exact factorization at symmetric triple coexistence:
# U_coex = (v/720) chi^2 (chi^2-chi0^2)^2.
xs = np.linspace(-1.4*chi_coex, 1.4*chi_coex, 2001)
factorized = (U6/720.0)*xs**2*(xs**2-chi_coex_sq)**2
factorization_maxerr = float(np.max(np.abs(Ueff(xs, r_coex, 0.0)-factorized)))

# Exact coexistence wall constant and 10-90 width in units sqrt(kappa/E0).
C_sigma_sym = (chi_coex_sq**2/4.0)*np.sqrt(U6/360.0)
width_10_90_sym = (
    np.log((0.9**2/(1-0.9**2))/(0.1**2/(1-0.1**2)))
    /(2.0*chi_coex_sq*np.sqrt(U6/360.0))
)

# ---- 2) Higher-order static response / 1PI vertex dictionary ----
# If chi(h)=R1 h + R3 h^3/3! + R5 h^5/5! + ... about the symmetric parent,
# then R1=1/r, R3=-u/r^4, R5=10 u^2/r^7-v/r^6.
R1 = 1.0/R0
R3 = -U4/R0**4
R5 = 10.0*U4**2/R0**7 - U6/R0**6
r_from_R = 1.0/R1
u_from_R = -R3*r_from_R**4
v_from_R = 10.0*u_from_R**2/r_from_R - R5*r_from_R**6

# Synthetic nonlinear response audit over the still-stable central branch.
hs = np.linspace(-0.85,0.85,241)
chis = np.array([central_root(h) for h in hs])
noise_fraction = 2.0e-4
chi_obs = chis*(1.0 + noise_fraction*RNG.normal(size=chis.shape))
A = np.column_stack([chi_obs, chi_obs**3, chi_obs**5])
coef, *_ = np.linalg.lstsq(A, hs, rcond=None)
r_hat = float(coef[0])
u_hat = float(6.0*coef[1])
v_hat = float(120.0*coef[2])
rco_hat = 5.0*u_hat*u_hat/(8.0*v_hat)
rspin_hat = 5.0*u_hat*u_hat/(6.0*v_hat)

pd.DataFrame({'h':hs,'chi_exact':chis,'chi_observed':chi_obs}).to_csv(OUT/'nonlinear_static_response.csv',index=False)

# ---- 3) Even gate becomes a measurable coexistence condition ----
# r(E)=r0+lambda_E E.  Once lambda_E is measured by v1.8 Kubo protocol,
# E_coex=(r_coex-r0)/lambda_E is no longer an arbitrary threshold.
LAMBDA_E = 1.2  # internal response-identification benchmark only
E_coex = (r_coex-R0)/LAMBDA_E
E_daughter_spinodal = (r_daughter_spinodal-R0)/LAMBDA_E

# ---- 4) Thin-wall nucleation barrier near even-gated coexistence ----
# For r = r_coex-delta_r, envelope theorem gives Delta f = E0 * chi0^2 delta_r / 2 + O(delta_r^2).
# sigma = C_sigma_sym E0 xi0, xi0=sqrt(kappa/E0).
# Thus Rc/xi0 ~ 4 C_sigma/(chi0^2 delta_r),
# Delta F*/(E0 xi0^3) ~ [64 pi C_sigma^3/(3 chi0^4)] delta_r^-2.
Rc_coeff = 4.0*C_sigma_sym/chi_coex_sq
Fstar_coeff = 64.0*np.pi*C_sigma_sym**3/(3.0*chi_coex_sq**2)

# Compute exact daughter minima and free-energy advantage to test the near-coexistence asymptotics.
deltas = np.array([0.01,0.02,0.05,0.08,0.12,0.20,0.30])
nuc_rows=[]
for dr in deltas:
    r = r_coex-dr
    roots = daughter_y_roots(r)
    assert roots is not None
    yb, ym = roots
    chi_m=np.sqrt(ym)
    f_d=Ueff(chi_m,r,0.0)
    df=-f_d # parent free energy is 0
    sigma=C_sigma_sym*np.sqrt(KAPPA*E0) # coexistence thin-wall value
    xi0=np.sqrt(KAPPA/E0)
    Rc_exact=2.0*sigma/(E0*df)
    Fstar_exact=16.0*np.pi*sigma**3/(3.0*(E0*df)**2)
    Rc_asym=Rc_coeff*xi0/dr
    Fstar_asym=Fstar_coeff*E0*xi0**3/dr**2
    nuc_rows.append({
        'delta_r':dr,'r':r,'chi_daughter':chi_m,'Deltaf_over_E0':df,
        'Rc_over_xi_exact':Rc_exact/xi0,'Rc_over_xi_asym':Rc_asym/xi0,
        'Fstar_over_E0xi3_exact':Fstar_exact/(E0*xi0**3),
        'Fstar_over_E0xi3_asym':Fstar_asym/(E0*xi0**3),
    })
nuc=pd.DataFrame(nuc_rows)
nuc.to_csv(OUT/'thin_wall_nucleation_audit.csv',index=False)

# ---- 5) Local-KMS/FDT benchmark for the Gaussian SK noise sector ----
# Dividing the v1.7 equation by Lambda gives friction 1/Lambda.  In a local-KMS benchmark:
# <eta eta> = 2 Theta/Lambda delta, hence S_chichi=(2 Theta/Lambda)|G_R|^2
# and Theta = omega S/(2 Im G_R).
def G_R(omega,k=0.0,r=R0,kappa=KAPPA,lam=LAMBDA,tau=TAU):
    return 1.0/(r+kappa*k*k-(tau/lam)*omega*omega-1j*omega/lam)

omegas=np.linspace(0.03,8.0,400)
fdt_rows=[]
for k in [0.0,1.5,3.0]:
    G=G_R(omegas,k)
    S=(2.0*THETA/LAMBDA)*np.abs(G)**2
    theta_rec=omegas*S/(2.0*G.imag)
    for om,gv,sv,tv in zip(omegas,G,S,theta_rec):
        fdt_rows.append({'omega':om,'k':k,'ReG':gv.real,'ImG':gv.imag,'S_chichi':sv,'Theta_recovered':tv})
fdt=pd.DataFrame(fdt_rows)
fdt.to_csv(OUT/'local_kms_fdt.csv',index=False)
theta_maxerr=float(np.max(np.abs(fdt.Theta_recovered-THETA)))

# Dimensionless local-KMS nucleation exponent benchmark (not cosmological prediction).
# theta_bar = Theta/(E0 xi0^3)
xi0=np.sqrt(KAPPA/E0)
theta_bar=THETA/(E0*xi0**3)
nuc['local_KMS_exponent_Fstar_over_Theta']=nuc.Fstar_over_E0xi3_exact/theta_bar
nuc.to_csv(OUT/'thin_wall_nucleation_audit.csv',index=False)

# ---- Figures ----
# Figure 1: 1PI gate landscape / spinodal control plane.
chi_sp=np.linspace(-1.35,1.35,900)
r_sp=-(U4/2.0)*chi_sp**2-(U6/24.0)*chi_sp**4
h_sp=-(U4/3.0)*chi_sp**3-(U6/30.0)*chi_sp**5
mask=r_sp>=-0.3
fig,ax=plt.subplots(figsize=(7.0,4.7))
ax.plot(r_sp[mask],h_sp[mask],label='spinodal envelope')
ax.axvline(r_coex,ls='--',label=fr'symmetric triple coexistence $r={r_coex:.3f}$')
ax.axvline(R0,ls=':',label=fr'baseline $r_0={R0:.1f}$')
ax.scatter([r_coex],[0.0],s=35,zorder=3)
ax.set_xlim(0,5.0)
ax.set_ylim(-1.1,1.1)
ax.set_xlabel(r'even gate / quadratic 1PI vertex $r$')
ax.set_ylabel(r'branch-odd conjugate field $h$')
ax.set_title('Two-control sextic landscape separates gate from daughter direction')
ax.legend(frameon=False,fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_1pi_gate_landscape.pdf')
fig.savefig(OUT/'fig_1pi_gate_landscape.png',dpi=180)
plt.close(fig)

# Figure 2: nonlinear equation-of-state identification.
fig,ax=plt.subplots(figsize=(7.0,4.4))
ax.plot(chis,hs,label='exact response curve')
ax.scatter(chi_obs[::8],hs[::8],s=10,label='perturbed synthetic samples')
chi_grid=np.linspace(chis.min(),chis.max(),500)
h_fit=r_hat*chi_grid+(u_hat/6.0)*chi_grid**3+(v_hat/120.0)*chi_grid**5
ax.plot(chi_grid,h_fit,ls='--',label='recovered 1PI equation of state')
ax.set_xlabel(r'order parameter response $\chi$')
ax.set_ylabel(r'conjugate source $h$')
ax.set_title('Nonlinear response identifies the quartic and sextic 1PI vertices')
ax.legend(frameon=False,fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_nonlinear_1pi_recovery.pdf')
fig.savefig(OUT/'fig_nonlinear_1pi_recovery.png',dpi=180)
plt.close(fig)

# Figure 3: thin-wall nucleation barrier divergence near coexistence.
fig,ax=plt.subplots(figsize=(7.0,4.4))
ax.loglog(nuc.delta_r,nuc.Fstar_over_E0xi3_exact,'o-',label='exact bulk advantage + coexistence wall')
ax.loglog(nuc.delta_r,nuc.Fstar_over_E0xi3_asym,'--',label=r'near-coexistence $\propto\delta r^{-2}$')
ax.set_xlabel(r'gate depth $\delta r=r_{\rm coex}-r$')
ax.set_ylabel(r'$\Delta F_*/(\mathcal{E}_0\xi_0^3)$')
ax.set_title('Even-gated conversion retains a nucleation barrier')
ax.legend(frameon=False,fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_nucleation_barrier.pdf')
fig.savefig(OUT/'fig_nucleation_barrier.png',dpi=180)
plt.close(fig)

# Figure 4: FDT effective-noise recovery.
fig,ax=plt.subplots(figsize=(7.0,4.4))
for k in [0.0,1.5,3.0]:
    s=fdt[fdt.k==k]
    ax.plot(s.omega,s.Theta_recovered,label=fr'$k={k:g}$')
ax.axhline(THETA,ls='--',label=fr'input $\Theta_\chi={THETA:.2f}$')
ax.set_xlabel(r'frequency $\omega$')
ax.set_ylabel(r'$\omega S_{\chi\chi}/(2\,\mathrm{Im}\,G^R_{\chi h})$')
ax.set_title('Local-KMS benchmark makes the Gaussian noise scale measurable')
ax.legend(frameon=False,fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_local_kms_noise.pdf')
fig.savefig(OUT/'fig_local_kms_noise.png',dpi=180)
plt.close(fig)

summary={
    'canonical_1PI_vertices':{'r0':R0,'u4':U4,'u6':U6},
    'exact_phase_landscape':{
        'daughter_spinodal_r':r_daughter_spinodal,
        'symmetric_triple_coexistence_r':r_coex,
        'chi_coexistence_squared':chi_coex_sq,
        'chi_coexistence':chi_coex,
        'baseline_to_coexistence_shift_r0_minus_rcoex':parent_to_coex_shift,
        'factorization_max_abs_error':factorization_maxerr,
        'coexistence_wall_Csigma':C_sigma_sym,
        'coexistence_width_10_90_over_xi0':width_10_90_sym,
        'near_coexistence_Rc_over_xi_coefficient':Rc_coeff,
        'near_coexistence_Fstar_over_E0xi3_coefficient':Fstar_coeff,
    },
    'higher_order_static_response':{
        'R1':R1,'R3':R3,'R5':R5,
        'reconstructed_from_exact_R':{'r':r_from_R,'u4':u_from_R,'u6':v_from_R},
        'synthetic_response_noise_fraction':noise_fraction,
        'fit_from_perturbed_equation_of_state':{'r':r_hat,'u4':u_hat,'u6':v_hat,'r_coex':rco_hat,'r_daughter_spinodal':rspin_hat},
    },
    'even_gate_response_benchmark':{
        'lambda_E':LAMBDA_E,
        'E_coexistence':E_coex,
        'E_daughter_spinodal':E_daughter_spinodal,
        'interpretation':'internal response-identification benchmark only; not a cosmological threshold fit',
    },
    'local_KMS_benchmark':{
        'Theta_input':THETA,'Theta_max_abs_recovery_error':theta_maxerr,
        'xi0_internal':xi0,'Theta_over_E0xi3':theta_bar,
        'interpretation':'optional local-equilibrium/FDT benchmark; AInstein cosmic coarse graining is not assumed thermal',
    },
    'structural_conclusion':[
        'branch-even gate is the quadratic 1PI vertex r(E), not the old Gamma-only selector',
        'branch-odd h selects C versus V after the even landscape has made daughter branches accessible',
        'r,u,v are measurable from nonlinear static response / zero-momentum 1PI vertices',
        'local-KMS noise can be measured from S_chichi and Im G_R; without a KMS/large-deviation closure abundance remains undetermined',
        'thin-wall nucleation gives a rate exponent only after the stochastic/noise sector is specified'
    ]
}
(OUT/'nonlinear_1pi_sk_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

# Compact manuscript figure: landscape, nonlinear 1PI recovery, nucleation barrier.
fig,axs=plt.subplots(1,3,figsize=(12.8,3.8))
axs[0].plot(r_sp[mask],h_sp[mask],label='spinodal')
axs[0].axvline(r_coex,ls='--',label='triple coexistence')
axs[0].axvline(R0,ls=':',label='baseline')
axs[0].set_xlim(0,5.0); axs[0].set_ylim(-1.1,1.1)
axs[0].set_xlabel(r'$r$'); axs[0].set_ylabel(r'$h$')
axs[0].set_title('(a) gate / direction plane')
axs[0].legend(frameon=False,fontsize=7)
axs[1].plot(chis,hs,label='exact')
axs[1].plot(chi_grid,h_fit,'--',label='recovered')
axs[1].set_xlabel(r'$\chi$'); axs[1].set_ylabel(r'$h$')
axs[1].set_title('(b) nonlinear 1PI recovery')
axs[1].legend(frameon=False,fontsize=7)
axs[2].loglog(nuc.delta_r,nuc.Fstar_over_E0xi3_exact,'o-',label='thin-wall')
axs[2].loglog(nuc.delta_r,nuc.Fstar_over_E0xi3_asym,'--',label=r'$\delta r^{-2}$')
axs[2].set_xlabel(r'$\delta r$')
axs[2].set_ylabel(r'$\Delta F_*/(\mathcal{E}_0\xi_0^3)$')
axs[2].set_title('(c) nucleation obstruction')
axs[2].legend(frameon=False,fontsize=7)
fig.tight_layout()
fig.savefig(OUT/'fig_nonlinear_sk_closure.pdf')
fig.savefig(OUT/'fig_nonlinear_sk_closure.png',dpi=180)
plt.close(fig)
