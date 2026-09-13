#!/usr/bin/env python3
"""Invert a simple Gamma threshold to quantify the shift required for large V filling.

This is not a proposal to change the selector. It is an adversarial diagnostic:
for the measured final PM basin catalogue, find Gamma_* such that the union of
basins with Gamma>=Gamma_* occupies a target volume fraction. It quantifies how
far the matter-curvature-equality landmark Gamma=1 lies from a volume-dominant
classification in the present coarse graining.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TARGETS=(0.10,0.20,0.50,0.70)

def threshold_for_fraction(g,target):
    q=g.sort_values('Gamma',ascending=False)
    c=q.volume_fraction.cumsum().to_numpy(float); gam=q.Gamma.to_numpy(float)
    i=int(np.searchsorted(c,target,side='left'));i=min(i,len(gam)-1)
    return float(gam[i])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(args.input);d=d[np.isclose(d.a,1.0)]
    rows=[]
    fig,ax=plt.subplots(figsize=(6.4,4.2))
    for seed,g in d.groupby('seed'):
        q=g.sort_values('Gamma',ascending=False)
        surv=q.volume_fraction.cumsum().to_numpy(float); gam=q.Gamma.to_numpy(float)
        ax.plot(gam,surv,label=f'seed {int(seed)}')
        row={'seed':int(seed)}
        for t in TARGETS: row[f'Gamma_threshold_for_phi_{int(t*100)}pct']=threshold_for_fraction(g,t)
        rows.append(row)
    tab=pd.DataFrame(rows);tab.to_csv(args.out/'selector_threshold_inversion.csv',index=False)
    summary={'targets':list(TARGETS),'rows':rows,'Gamma_eq_landmark':1.0,
             'Gamma_50pct_range':[float(tab.Gamma_threshold_for_phi_50pct.min()),float(tab.Gamma_threshold_for_phi_50pct.max())],
             'Gamma_70pct_range':[float(tab.Gamma_threshold_for_phi_70pct.min()),float(tab.Gamma_threshold_for_phi_70pct.max())]}
    (args.out/'selector_threshold_inversion.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (args.out/'selector_threshold_inversion.txt').write_text('\n'.join(f'{k}={v}' for k,v in summary.items())+'\n',encoding='utf-8')
    ax.axvline(1.0,linewidth=1);ax.axhline(0.5,linewidth=1);ax.axhline(0.7,linewidth=1)
    ax.set_xlim(min(-0.5,float(d.Gamma.min())),min(2.5,float(d.Gamma.max())))
    ax.set_ylim(0,1);ax.set_xlabel(r'threshold $\Gamma_*$');ax.set_ylabel(r'volume fraction with $\Gamma\geq\Gamma_*$')
    ax.set_title('Threshold inversion of the final nonlinear PM basin catalogue');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(args.out/'fig_selector_threshold_inversion.pdf');fig.savefig(args.out/'fig_selector_threshold_inversion.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
