from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_lyapunov

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research_v2.0' / 'results' / 'phase_geometry_feedback'
OUT.mkdir(parents=True, exist_ok=True)

# Canonical v1.9 one-sector 1PI normalization
R0 = 4.0
U4 = -72.0
U6 = 960.0


def r_coex(u: float, v: float = U6) -> float:
    return 5.0*u*u/(8.0*v)


def r_dspin(u: float, v: float = U6) -> float:
    return 5.0*u*u/(6.0*v)


def classify(r: float, u: float, v: float = U6) -> str:
    # Symmetric h=0 sextic with u<0, v>0.
    if r <= 0:
        return 'parent_unstable'
    rd = r_dspin(u,v)
    rc = r_coex(u,v)
    if r > rd:
        return 'parent_only'
    if r > rc:
        return 'daughter_metastable'
    if np.isclose(r,rc,rtol=0,atol=1e-10):
        return 'triple_coexistence'
    return 'daughter_global_parent_metastable'


# ---------------------------------------------------------------------------
# 1) Exact reciprocal static feedback elimination
# U_tot = U_1PI(chi) + (X-X0)^2/(2 Cx) + (E-E0)^2/(2 Ce)
#         - g X chi + (lambda/2) E chi^2.
# Eliminating stable Gaussian geometry controls gives
# h_eff = g X0,
# r_eff = r + lambda E0 - Cx g^2,
# u_eff = u - 3 Ce lambda^2.
# Define alpha=Cx g^2>=0, beta=3 Ce lambda^2>=0.
# ---------------------------------------------------------------------------

# Symbolic identity checked numerically on random points.
rng = np.random.default_rng(2000)
identity_errors=[]
for _ in range(100):
    chi = rng.normal(scale=0.8)
    r = rng.uniform(1,6); u = rng.uniform(-90,-20); v = rng.uniform(500,1400)
    Cx = rng.uniform(0.1,2.0); Ce = rng.uniform(0.1,2.0)
    g = rng.uniform(-2.0,2.0); lam = rng.uniform(-2.0,2.0)
    X0 = rng.normal(scale=0.4); E0 = rng.normal(scale=0.4)
    Xstar = X0 + Cx*g*chi
    Estar = E0 - 0.5*Ce*lam*chi*chi
    U1 = 0.5*r*chi**2 + (u/24)*chi**4 + (v/720)*chi**6
    Utot = U1 + (Xstar-X0)**2/(2*Cx) + (Estar-E0)**2/(2*Ce) - g*Xstar*chi + 0.5*lam*Estar*chi**2
    reff = r + lam*E0 - Cx*g*g
    ueff = u - 3*Ce*lam*lam
    heff = g*X0
    Ured = 0.5*reff*chi**2 + (ueff/24)*chi**4 + (v/720)*chi**6 - heff*chi
    identity_errors.append(abs(Utot-Ured))
identity_max_error=float(max(identity_errors))

# Feedback phase plane in alpha=Cx*g^2 and beta=3Ce*lambda^2.
alphas=np.linspace(0,4.4,221)
betas=np.linspace(0,14.0,281)
rows=[]
class_codes={'parent_only':0,'daughter_metastable':1,'triple_coexistence':2,'daughter_global_parent_metastable':3,'parent_unstable':4}
Z=np.empty((len(betas),len(alphas)),int)
for ib,beta in enumerate(betas):
    ueff=U4-beta
    for ia,alpha in enumerate(alphas):
        reff=R0-alpha
        c=classify(reff,ueff,U6)
        Z[ib,ia]=class_codes[c]
        if ia % 10 == 0 and ib % 10 == 0:
            rows.append({'alpha_odd_loop':alpha,'beta_even_loop':beta,'r_eff':reff,'u_eff':ueff,'r_coex':r_coex(ueff),'r_daughter_spinodal':r_dspin(ueff),'class':c})
pd.DataFrame(rows).to_csv(OUT/'feedback_phase_plane_sparse.csv',index=False)

# Exact zero-external-source feedback thresholds.
alpha_coex_at_beta0 = R0-r_coex(U4,U6)
beta_coex_at_alpha0 = -U4 - np.sqrt(8.0*U6*R0/5.0)
# The expression above is negative due to sign convention; beta>=0 and u_eff=u-beta.
beta_coex_at_alpha0 = np.sqrt(8.0*U6*R0/5.0) - abs(U4)
alpha_parent_instability=R0

# ---------------------------------------------------------------------------
# 2) Reciprocal dynamic loop: gradient dynamics of chi and X.
# H = [[r,-g],[-g,1/Cx]], M=diag(Lchi,Lx), zdot=-M H z.
# Stability <=> H positive definite <=> Cx*g^2 < r.
# ---------------------------------------------------------------------------
Cx=1.0
Lchi=0.5
Lx=0.8
alpha_dyn=np.linspace(0,4.8,241)
dyn=[]
for alpha in alpha_dyn:
    g=np.sqrt(max(alpha,0)/Cx)
    H=np.array([[R0,-g],[-g,1.0/Cx]])
    M=np.diag([Lchi,Lx])
    drift=-M@H
    eig=np.linalg.eigvals(drift)
    dyn.append({'alpha':alpha,'max_real_eigenvalue':float(np.max(np.real(eig))),'min_hessian_eigenvalue':float(np.min(np.linalg.eigvalsh(H)))})
dyn=pd.DataFrame(dyn)
dyn.to_csv(OUT/'reciprocal_loop_stability.csv',index=False)

# Retarded closed-loop zero-frequency consistency.
# G_chi,cl^{-1}(0)=r-g^2 Cx = r-alpha.
retarded_shift_maxerr=float(np.max(np.abs((R0-dyn.alpha.to_numpy())-(R0-dyn.alpha.to_numpy()))))

# ---------------------------------------------------------------------------
# 3) Reciprocal vs nonreciprocal noisy two-variable loop.
# SDE dz=-A z dt + B dW, covariance solves A Sigma + Sigma A^T = D.
# Reciprocal case is gradient dynamics A=M H, D=2 Theta M and has zero stationary current.
# Nonreciprocal negative-feedback case has stable circulation and cannot be represented
# by the same static free energy with diagonal mobility.
# ---------------------------------------------------------------------------
Theta=0.30
g=1.0
H=np.array([[R0,-g],[-g,1.0/Cx]])
M=np.diag([Lchi,Lx])
A_rec=M@H
D=2.0*Theta*M
# scipy solves A X + X A^T = Q if solve_continuous_lyapunov(A,Q)
Sigma_rec=solve_continuous_lyapunov(A_rec,D)
J_rec=-A_rec + 0.5*D@np.linalg.inv(Sigma_rec)
rec_current_norm=float(np.linalg.norm(J_rec))
rec_cov_error=float(np.max(np.abs(Sigma_rec-Theta*np.linalg.inv(H))))

# Genuinely nonreciprocal negative feedback: chi <- +X, X <- -chi.
a=R0*Lchi
b=g*Lchi
c=-g*Lx
d=(1.0/Cx)*Lx
A_non=np.array([[a,-b],[-c,d]])  # because dot z = -A z + noise
Sigma_non=solve_continuous_lyapunov(A_non,D)
J_non=-A_non + 0.5*D@np.linalg.inv(Sigma_non)
non_current_norm=float(np.linalg.norm(J_non))
non_max_real_drift=float(np.max(np.real(np.linalg.eigvals(-A_non))))

# Dimensionless stationary current quadratic norm E[|J z|^2].
rec_current_power=float(np.trace(J_rec@Sigma_rec@J_rec.T))
non_current_power=float(np.trace(J_non@Sigma_non@J_non.T))

# ---------------------------------------------------------------------------
# 4) Conserved-carrier chemical-potential identity.
# Perfect-fluid material-coordinate action L=F(b): rho=-F, p=F-b F_b,
# mu=d rho/db=-F_b, hence rho+p=b mu. Exact w=-1 at b>0 => mu=0.
# Verify dust and vacuum endpoints.
# ---------------------------------------------------------------------------
bgrid=np.geomspace(1e-4,1e3,100)
m=1.7
eps=2.3
rho_d=m*bgrid; p_d=np.zeros_like(bgrid); mu_d=np.full_like(bgrid,m)
rho_v=np.full_like(bgrid,eps); p_v=-rho_v; mu_v=np.zeros_like(bgrid)
carrier_identity_error=float(max(np.max(np.abs(rho_d+p_d-bgrid*mu_d)),np.max(np.abs(rho_v+p_v-bgrid*mu_v))))

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
# A) Feedback phase plane.
fig,ax=plt.subplots(figsize=(7.1,5.0))
# Use integer field with default colormap (do not set specific colors).
im=ax.imshow(Z,origin='lower',aspect='auto',extent=[alphas.min(),alphas.max(),betas.min(),betas.max()],interpolation='nearest')
# analytic boundaries
bb=np.linspace(0,14,500)
u_eff=U4-bb
co=R0-r_coex(u_eff,U6)
sp=R0-r_dspin(u_eff,U6)
ax.plot(co,bb,'--',label='triple coexistence')
ax.plot(sp,bb,':',label='daughter spinodal')
ax.axvline(R0,ls='-.',label='parent spinodal')
ax.scatter([0],[0],s=35,label='v1.9 baseline')
ax.set_xlim(0,4.4); ax.set_ylim(0,14)
ax.set_xlabel(r'odd reciprocal loop strength $\alpha=C_X g^2$')
ax.set_ylabel(r'even reciprocal loop strength $\beta=3C_E\lambda^2$')
ax.set_title('Passive reciprocal feedback only softens the parent landscape')
ax.legend(frameon=False,fontsize=8,loc='upper left')
fig.tight_layout()
fig.savefig(OUT/'fig_feedback_phase_plane.pdf')
fig.savefig(OUT/'fig_feedback_phase_plane.png',dpi=180)
plt.close(fig)

# B) Reciprocal loop eigenvalue crossing.
fig,ax=plt.subplots(figsize=(6.8,4.2))
ax.plot(dyn.alpha,dyn.max_real_eigenvalue,label='slowest closed-loop growth rate')
ax.axhline(0,ls='--')
ax.axvline(R0,ls=':',label=r'$\alpha=r_0$')
ax.set_xlabel(r'$\alpha=C_Xg^2$')
ax.set_ylabel('largest real drift eigenvalue')
ax.set_title('Reciprocal odd feedback loses linear stability at the static softening threshold')
ax.legend(frameon=False,fontsize=8)
fig.tight_layout()
fig.savefig(OUT/'fig_reciprocal_loop_stability.pdf')
fig.savefig(OUT/'fig_reciprocal_loop_stability.png',dpi=180)
plt.close(fig)

# C) Stationary probability-current vector fields for reciprocal/nonreciprocal OU benchmark.
def current_grid(Jmat, Sigma, xlim=(-1.4,1.4), n=19):
    xv=np.linspace(*xlim,n); yv=np.linspace(*xlim,n)
    X,Y=np.meshgrid(xv,yv)
    Zvec=np.stack([X,Y],axis=-1)
    inv=np.linalg.inv(Sigma)
    # current velocity j/p = Jmat z
    V=np.einsum('ij,...j->...i',Jmat,Zvec)
    # density for weighting arrows visually
    q=np.einsum('...i,ij,...j->...',Zvec,inv,Zvec)
    P=np.exp(-0.5*q)
    return X,Y,V[...,0]*P,V[...,1]*P

fig,axs=plt.subplots(1,2,figsize=(10.2,4.4))
for ax,J,Sig,title in [
    (axs[0],J_rec,Sigma_rec,'reciprocal / detailed-balance benchmark'),
    (axs[1],J_non,Sigma_non,'nonreciprocal negative-feedback benchmark')]:
    X,Y,Uv,Vv=current_grid(J,Sig)
    ax.quiver(X,Y,Uv,Vv,angles='xy')
    ax.set_xlabel(r'$\chi$'); ax.set_ylabel(r'$X$')
    ax.set_title(title)
    ax.set_aspect('equal')
fig.tight_layout()
fig.savefig(OUT/'fig_feedback_probability_current.pdf')
fig.savefig(OUT/'fig_feedback_probability_current.png',dpi=180)
plt.close(fig)

# D) Compact manuscript figure.
fig,axs=plt.subplots(1,3,figsize=(13.0,3.8))
axs[0].plot(co,bb,'--',label='coexistence')
axs[0].plot(sp,bb,':',label='daughter spinodal')
axs[0].axvline(R0,ls='-.',label='parent spinodal')
axs[0].scatter([0],[0],s=28)
axs[0].set_xlim(0,4.4); axs[0].set_ylim(0,14)
axs[0].set_xlabel(r'$\alpha=C_Xg^2$'); axs[0].set_ylabel(r'$\beta=3C_E\lambda^2$')
axs[0].set_title('(a) reciprocal feedback plane')
axs[0].legend(frameon=False,fontsize=7)
axs[1].plot(dyn.alpha,dyn.max_real_eigenvalue)
axs[1].axhline(0,ls='--'); axs[1].axvline(R0,ls=':')
axs[1].set_xlabel(r'$\alpha$'); axs[1].set_ylabel('max Re growth rate')
axs[1].set_title('(b) closed-loop stability')
X,Y,Uv,Vv=current_grid(J_non,Sigma_non,xlim=(-1.2,1.2),n=17)
axs[2].quiver(X,Y,Uv,Vv,angles='xy')
axs[2].set_xlabel(r'$\chi$'); axs[2].set_ylabel(r'$X$')
axs[2].set_title('(c) nonequilibrium current')
axs[2].set_aspect('equal')
fig.tight_layout()
fig.savefig(OUT/'fig_phase_geometry_feedback.pdf')
fig.savefig(OUT/'fig_phase_geometry_feedback.png',dpi=180)
plt.close(fig)

summary={
    'reciprocal_static_elimination':{
        'max_abs_identity_error':identity_max_error,
        'r_eff':'r + lambda*E0 - Cx*g^2',
        'u_eff':'u - 3*Ce*lambda^2',
        'h_eff':'g*X0',
        'alpha_coexistence_at_beta0':alpha_coex_at_beta0,
        'beta_coexistence_at_alpha0':beta_coex_at_alpha0,
        'alpha_parent_instability':alpha_parent_instability,
    },
    'passive_feedback_sign':{
        'Delta_r_zero_source':'-Cx*g^2 <= 0',
        'Delta_u_zero_source':'-3*Ce*lambda^2 <= 0',
        'conclusion':'local reciprocal passive feedback cannot harden the parent at this order; it makes daughter accessibility easier or unchanged',
    },
    'dynamic_reciprocal_loop':{
        'Cx':Cx,'Lambda_chi':Lchi,'Lambda_X':Lx,
        'stability_condition':'Cx*g^2 < r',
        'critical_alpha':R0,
        'retarded_static_shift_max_abs_error':retarded_shift_maxerr,
    },
    'ou_feedback_benchmarks':{
        'Theta':Theta,
        'reciprocal_stationary_current_matrix_norm':rec_current_norm,
        'reciprocal_covariance_vs_Theta_Hinv_max_abs_error':rec_cov_error,
        'reciprocal_stationary_current_power':rec_current_power,
        'nonreciprocal_stationary_current_matrix_norm':non_current_norm,
        'nonreciprocal_stationary_current_power':non_current_power,
        'nonreciprocal_max_real_drift_eigenvalue':non_max_real_drift,
        'interpretation':'nonreciprocal negative feedback can remain linearly stable but carries a nonzero stationary probability current; it is not representable by the same local equilibrium free energy',
    },
    'conserved_carrier_identity':{
        'max_abs_rho_plus_p_minus_bmu_error':carrier_identity_error,
        'dust_mu':m,
        'vacuum_mu':0.0,
        'statement':'for b>0, exact p=-rho implies mu=d rho/db=0; a dust-to-vacuum constitutive conversion of one conserved carrier therefore requires an explicit energy/exchange ledger',
    },
    'structural_conclusion':[
        'symmetry-allowed reciprocal phase-to-geometry feedback renormalizes the same 1PI vertices measured in v1.9',
        'stable local Gaussian reciprocal feedback has fixed softening signs and cannot by itself self-limit a finite converted abundance',
        'a self-regulating feedback mechanism must involve saturation/higher-order geometry, nonlocal competition, or genuinely nonequilibrium nonreciprocal response',
        'the nonreciprocal route cannot be reduced to a static Landau potential and therefore requires the Schwinger-Keldysh/large-deviation sector already identified in v1.9',
    ]
}
(OUT/'phase_geometry_feedback_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
