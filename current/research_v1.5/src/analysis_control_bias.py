#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

H_COEX = 0.341688495018
CHI_D = 0.087415890295
CHI_V = 1.035423691165
DELTA_CHI = CHI_V - CHI_D
C_SIGMA = 0.396830116706
XI0_BENCH = 3.12938181925  # h^-1 Mpc
R_BENCH = 15.6             # h^-1 Mpc
OMEGA_M = 0.315
OMEGA_V = 0.685

INTERVAL=(0.75,1.0)

def identities(Gamma, nu):
    qhat = Gamma*(nu-1.0)          # Q/(4 pi G rho_m)
    rhat = Gamma*(nu+3.0)          # -R/(4 pi G rho_m)
    Gamma2 = (rhat-qhat)/4.0
    nu2 = (rhat+3.0*qhat)/(rhat-qhat)
    K_over_rhom = (rhat-3.0*qhat)/12.0
    U_over_rhom = rhat/6.0
    return qhat, rhat, Gamma2, nu2, K_over_rhom, U_over_rhom

def compute_interval(df: pd.DataFrame, a0=0.75, a1=1.0):
    d0=df[np.isclose(df.a,a0)].copy()
    d1=df[np.isclose(df.a,a1)].copy()
    m=d0.merge(d1,on=['seed','final_basin_id'],suffixes=('_0','_1'))
    aratio=(a1/a0)*(m.volume_fraction_1/m.volume_fraction_0)**(1.0/3.0)
    valid=(m.Gamma_0>0)&(m.Gamma_1>0)&(aratio>0)&np.isfinite(aratio)
    m=m.loc[valid].copy(); aratio=aratio.loc[valid]
    m['aD_ratio']=aratio
    m['nu_hist']=np.log(m.Gamma_1/m.Gamma_0)/np.log(aratio)
    m['Gamma_mid']=np.sqrt(m.Gamma_0*m.Gamma_1)
    m['qhat_hist']=m.Gamma_mid*(m.nu_hist-1.0)
    m['rhat_hist']=m.Gamma_mid*(m.nu_hist+3.0)
    return m

def volume_fraction_above(m: pd.DataFrame, threshold: float):
    rows=[]
    for seed,g in m.groupby('seed'):
        rows.append({'seed':int(seed),'threshold':threshold,
                     'absolute_volume_fraction':float(g.loc[g.qhat_hist>threshold,'volume_fraction_1'].sum()),
                     'positive_Gamma_volume_fraction':float(g.volume_fraction_1.sum())})
    return rows

def threshold_for_target(g: pd.DataFrame, target: float):
    s=g[['qhat_hist','volume_fraction_1']].sort_values('qhat_hist',ascending=False)
    cum=s.volume_fraction_1.cumsum().to_numpy()
    if len(cum)==0 or cum[-1] < target:
        return math.nan
    idx=int(np.searchsorted(cum,target,side='left'))
    return float(s.qhat_hist.iloc[idx])

def audit_dataset(path: Path, label: str):
    df=pd.read_csv(path)
    m=compute_interval(df)
    rows=[]
    for seed,g in m.groupby('seed'):
        rows.append({
            'dataset':label,'seed':int(seed),
            'valid_positive_Gamma_volume':float(g.volume_fraction_1.sum()),
            'volume_qhat_gt_0':float(g.loc[g.qhat_hist>0,'volume_fraction_1'].sum()),
            'volume_qhat_gt_0p5':float(g.loc[g.qhat_hist>0.5,'volume_fraction_1'].sum()),
            'volume_qhat_gt_0p743':float(g.loc[g.qhat_hist>0.743,'volume_fraction_1'].sum()),
            'volume_qhat_gt_1':float(g.loc[g.qhat_hist>1.0,'volume_fraction_1'].sum()),
            'threshold_for_1pct_total_volume':threshold_for_target(g,0.01),
            'threshold_for_5pct_total_volume':threshold_for_target(g,0.05),
            'threshold_for_10pct_total_volume':threshold_for_target(g,0.10),
            'threshold_for_20pct_total_volume':threshold_for_target(g,0.20),
            'threshold_for_50pct_total_volume':threshold_for_target(g,0.50),
        })
    return m, pd.DataFrame(rows)

def make_basis_plot(out: Path):
    q=np.linspace(-2.0,2.0,500)
    fig,ax=plt.subplots(figsize=(6.7,5.0))
    ax.plot(q, q+4.0, label=r'$\Gamma=1:\ \hat R-\hat Q=4$')
    ax.plot(q, 3.0*q, linestyle='--', label=r'vacuum-like geometry: $\hat R=3\hat Q$')
    ax.axvline(0.0, linestyle=':', label=r'curvature-only: $\hat Q=0$')
    ax.axvline(1.0, linestyle='-.', label=r'geometry-only acceleration: $\hat Q=1$')
    ax.set_xlim(-1.5,2.0); ax.set_ylim(-2.5,7.0)
    ax.set_xlabel(r'$\hat Q\equiv\mathcal{Q}/(4\pi G\rho_m)$')
    ax.set_ylabel(r'$\hat R\equiv-\langle R\rangle/(4\pi G\rho_m)$')
    ax.set_title('Regional geometry-source basis')
    ax.legend(frameon=False,fontsize=8)
    fig.tight_layout()
    fig.savefig(out/'fig_control_basis.pdf',bbox_inches='tight')
    fig.savefig(out/'fig_control_basis.png',dpi=200,bbox_inches='tight')
    plt.close(fig)

def make_cumulative_plot(m64: pd.DataFrame, m80: pd.DataFrame, out: Path, qcoex: float, qgrow: float):
    fig,ax=plt.subplots(figsize=(7.2,4.8))
    grid=np.linspace(-0.1,1.3,281)
    for seed,g in m64.groupby('seed'):
        y=[g.loc[g.qhat_hist>x,'volume_fraction_1'].sum() for x in grid]
        ax.plot(grid,y,linewidth=1.1,label=f'64^3 seed {seed}')
    # one higher-resolution check
    g=m80
    y=[g.loc[g.qhat_hist>x,'volume_fraction_1'].sum() for x in grid]
    ax.plot(grid,y,linewidth=1.8,label='80^3 seed 11')
    ax.axvline(qcoex,linestyle='--',linewidth=1.2,label='canonical infinite-domain coexistence')
    ax.axvline(qgrow,linestyle=':',linewidth=1.2,label='canonical finite-size benchmark')
    ax.axvline(1.0,linestyle='-.',linewidth=1.0,label='geometry-only acceleration')
    ax.set_xlabel(r'history proxy $\hat Q_{\rm hist}=\Gamma_{\rm mid}(\nu-1)$')
    ax.set_ylabel('final-volume fraction above threshold')
    ax.set_ylim(bottom=0)
    ax.set_title('Curvature-neutral drive audit, final interval')
    ax.legend(frameon=False,fontsize=7.5,ncol=2)
    fig.tight_layout()
    fig.savefig(out/'fig_qdrive_cumulative.pdf',bbox_inches='tight')
    fig.savefig(out/'fig_qdrive_cumulative.png',dpi=200,bbox_inches='tight')
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--n64',type=Path,required=True)
    ap.add_argument('--n80',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)

    # Exact coefficient in finite-size near-coexistence expansion.
    capillary_coeff=2.0*C_SIGMA/DELTA_CHI
    eps_over_rhom=OMEGA_V/OMEGA_M
    qcoex_r0_1=H_COEX*eps_over_rhom
    hgrow=H_COEX+capillary_coeff*XI0_BENCH/R_BENCH
    qgrow_r0_1=hgrow*eps_over_rhom

    m64,s64=audit_dataset(args.n64,'64^3')
    m80,s80=audit_dataset(args.n80,'80^3')
    summary=pd.concat([s64,s80],ignore_index=True)
    summary.to_csv(args.out/'qdrive_volume_audit.csv',index=False)
    m64.to_csv(args.out/'qdrive_interval_n64.csv',index=False)
    m80.to_csv(args.out/'qdrive_interval_n80.csv',index=False)

    # Test algebraic basis on a deterministic grid.
    pts=[]
    for Gamma in [0.1,0.5,1.0,2.0]:
        for nu in [0.5,1.0,2.0,3.0,4.0]:
            qh,rh,G2,n2,K,U=identities(Gamma,nu)
            pts.append({'Gamma':Gamma,'nu':nu,'qhat':qh,'rhat':rh,
                        'Gamma_reconstructed':G2,'nu_reconstructed':n2,
                        'K_over_rhom':K,'U_over_rhom':U})
    pd.DataFrame(pts).to_csv(args.out/'control_basis_identity_grid.csv',index=False)

    outj={
        'exact_identities':{
            'qhat':'Q/(4 pi G rho_m) = Gamma (nu-1)',
            'rhat':'-<R>/(4 pi G rho_m) = Gamma (nu+3)',
            'inverse_Gamma':'(rhat-qhat)/4',
            'inverse_nu':'(rhat+3 qhat)/(rhat-qhat)',
            'K_G_over_rho_m':'(rhat-3 qhat)/12',
            'U_G_over_rho_m':'rhat/6',
            'acceleration':'qhat > 1'
        },
        'leading_bias_basis':{
            'general':'H_ext = b_Q qhat + b_R rhat + higher order/history terms',
            'statement':'Averaged kinematics fixes the basis but not the mixing ratio b_Q/b_R.',
            'curvature_neutral_subclass':'b_R=0 => H_ext proportional to qhat',
            'status':'curvature-neutrality is an additional physical postulate, not a GR theorem'
        },
        'canonical_q_coupling_benchmark':{
            'definition':'H_Q = Q/(4 pi G E0)',
            'r0':'E0/epsilon_V',
            'H_coex':H_COEX,
            'Delta_chi':DELTA_CHI,
            'capillary_coefficient_2Csigma_over_DeltaChi':capillary_coeff,
            'OmegaV_over_OmegaM_reference':eps_over_rhom,
            'qhat_coex_at_r0_1_z0_reference':qcoex_r0_1,
            'benchmark_xi0_hinvMpc':XI0_BENCH,
            'benchmark_R_hinvMpc':R_BENCH,
            'H_grow_benchmark':hgrow,
            'qhat_grow_at_r0_1_z0_reference':qgrow_r0_1,
            'warning':'This is a canonical normalization audit, not a microscopic derivation of the coupling.'
        },
        'volume_audit':summary.to_dict(orient='records'),
        'caveat':'qhat_hist uses finite-difference descendant histories of the operational PM Gamma proxy. It is not a direct relativistic measurement of Buchert Q.'
    }
    (args.out/'control_bias_summary.json').write_text(json.dumps(outj,indent=2),encoding='utf-8')
    make_basis_plot(args.out)
    make_cumulative_plot(m64,m80,args.out,qcoex_r0_1,qgrow_r0_1)
    print(json.dumps(outj['canonical_q_coupling_benchmark'],indent=2))
    print(summary.to_string(index=False))

if __name__=='__main__':
    main()
