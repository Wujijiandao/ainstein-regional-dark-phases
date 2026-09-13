#!/usr/bin/env python3
"""Finite-size selector ceiling audit.

For the v1.2 conditional growth criterion
 Gamma >= 1 + 2 C_sigma (xi0/beta_Gamma)/R,
with positive interface coefficient and positive xi0/beta, the active set is a
subset of {Gamma>=1}. Therefore the measured raw Gamma>=1 filling fraction is
an absolute ceiling for every positive finite-size tuning at fixed dynamics,
coarse graining and regional estimator.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

C_SIGMA=0.396830116706
CURRENT_LAMBDA=3.12938

def audit(df):
    rows=[]
    lambdas=np.r_[0.0,np.geomspace(0.03,30.0,100)]
    curves={}
    for seed,g in df.groupby('seed'):
        q=g[np.isclose(g.a,1.0)].copy()
        vols=q.volume_fraction.to_numpy(float); gamma=q.Gamma.to_numpy(float); R=q.R_eff_hMpc.to_numpy(float)
        phi=[]
        for lam in lambdas:
            active=gamma >= 1.0+2*C_SIGMA*lam/R
            phi.append(float(vols[active].sum()))
        raw=float(vols[gamma>=1.0].sum())
        cur=float(vols[gamma>=1.0+2*C_SIGMA*CURRENT_LAMBDA/R].sum())
        rows.append({'seed':int(seed),'gamma_ge1_ceiling':raw,'current_finite_size_phi':cur,
                     'current_to_ceiling_ratio':cur/raw if raw>0 else np.nan})
        curves[int(seed)]=phi
    return pd.DataFrame(rows),lambdas,curves

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(args.input); tab,lams,curves=audit(d);tab.to_csv(args.out/'selector_ceiling_table.csv',index=False)
    summary={'criterion':'Gamma >= 1 + 2 C_sigma lambda/R, lambda=xi0/beta_Gamma >=0',
             'algebraic_ceiling':'phi_active(lambda) <= phi(Gamma>=1)',
             'seeds':[int(x) for x in tab.seed],
             'gamma_ge1_ceiling_range':[float(tab.gamma_ge1_ceiling.min()),float(tab.gamma_ge1_ceiling.max())],
             'gamma_ge1_ceiling_mean':float(tab.gamma_ge1_ceiling.mean()),
             'current_finite_size_phi_range':[float(tab.current_finite_size_phi.min()),float(tab.current_finite_size_phi.max())],
             'current_finite_size_phi_mean':float(tab.current_finite_size_phi.mean()),
             'current_to_ceiling_ratio_range':[float(tab.current_to_ceiling_ratio.min()),float(tab.current_to_ceiling_ratio.max())]}
    (args.out/'selector_ceiling_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (args.out/'selector_ceiling_summary.txt').write_text('\n'.join(f'{k}={v}' for k,v in summary.items())+'\n',encoding='utf-8')
    fig,ax=plt.subplots(figsize=(6.4,4.2))
    for seed,phi in curves.items(): ax.plot(lams,phi,label=f'seed {seed}')
    ax.axvline(CURRENT_LAMBDA,linewidth=1)
    ax.set_xscale('symlog',linthresh=0.03);ax.set_xlabel(r'$\lambda=\xi_0/\beta_\Gamma\;[h^{-1}{\rm Mpc}]$');ax.set_ylabel('active volume fraction')
    ax.set_title('Finite-size tuning cannot exceed the raw $\\Gamma\\geq1$ ceiling')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(args.out/'fig_selector_ceiling.pdf');fig.savefig(args.out/'fig_selector_ceiling.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
