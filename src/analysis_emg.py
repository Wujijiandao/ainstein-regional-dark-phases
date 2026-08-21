#!/usr/bin/env python3
"""Reproduce the one-parameter SPARC RAR diagnostic used in AInstein/EMG.

The spherical bound-phase formula is algebraically equivalent to the historical
simple-MOND relation; this script treats the fit only as a consistency check.
"""
from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar


def load_rar(path: Path) -> np.ndarray:
    rows=[]
    with path.open('r',encoding='utf-8') as f:
        for line in f:
            if re.match(r'^\s*-?\d', line):
                p=line.split()
                if len(p)>=4:
                    try:
                        rows.append(tuple(map(float,p[:4])))
                    except ValueError:
                        pass
    a=np.asarray(rows,dtype=float)
    if a.ndim != 2 or a.shape[1] != 4:
        raise RuntimeError(f'Could not parse four-column RAR table: {path}')
    return a


def g_emg(gb: np.ndarray, aE: float) -> np.ndarray:
    return 0.5*(gb+np.sqrt(gb*gb+4.0*aE*gb))


def g_rar(gb: np.ndarray, a0: float) -> np.ndarray:
    z=np.sqrt(gb/a0)
    return gb/(1.0-np.exp(-z))


def log_slope(func, gb: np.ndarray, a0: float, eps: float=1e-5) -> np.ndarray:
    return (np.log10(func(gb*10**eps,a0))-np.log10(func(gb/10**eps,a0)))/(2*eps)


def fit_model(x,sx,y,sy,func):
    gb=10.0**x
    def chi2(loga):
        aa=10.0**loga
        ym=np.log10(func(gb,aa))
        sl=log_slope(func,gb,aa)
        var=sy**2+(sl*sx)**2
        return float(np.sum((y-ym)**2/var))
    def rms(loga):
        aa=10.0**loga
        return float(np.sqrt(np.mean((y-np.log10(func(gb,aa)))**2)))
    rw=minimize_scalar(chi2,bounds=(-12,-9),method='bounded')
    ru=minimize_scalar(rms,bounds=(-12,-9),method='bounded')
    return 10**rw.x,chi2(rw.x)/(len(y)-1),rms(rw.x),10**ru.x,rms(ru.x)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data'/'RAR.mrt')
    ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results')
    ap.add_argument('--fig',type=Path,default=Path(__file__).resolve().parents[1]/'figures')
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True); args.fig.mkdir(parents=True,exist_ok=True)
    x,sx,y,sy=load_rar(args.data).T
    models=[('Regional bound-phase spherical relation',g_emg),('Empirical exponential RAR',g_rar)]
    results=[]
    for name,func in models:
        vals=fit_model(x,sx,y,sy,func)
        results.append((name,*vals))
    csv_path=args.out/'rar_fit_summary.csv'
    with csv_path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['model','a_weighted_m_s2','chi2_per_dof','rms_dex_at_weighted','a_unweighted_m_s2','rms_dex_unweighted'])
        for r in results: w.writerow([r[0],*[f'{v:.10g}' for v in r[1:]]])

    best=results[0][1]; gb=10**x
    xx=np.linspace(-13,-8,600); ggrid=10**xx
    fig=plt.figure(figsize=(6.6,6.2)); gs=fig.add_gridspec(2,1,height_ratios=[3.2,1],hspace=0.05)
    ax=fig.add_subplot(gs[0]); ax.scatter(x,y,s=3,alpha=.18,rasterized=True,label='SPARC RAR points')
    ax.plot(xx,np.log10(g_emg(ggrid,best)),lw=2,label='Regional bound-phase relation'); ax.plot(xx,xx,lw=1,ls='--',label='Newtonian')
    ax.set_ylabel(r'$\log_{10} g_{\rm obs}\;[\mathrm{m\,s^{-2}}]$'); ax.set_xlim(-12.2,-8); ax.set_ylim(-12.2,-8); ax.legend(frameon=False,fontsize=8,loc='upper left'); ax.tick_params(labelbottom=False)
    ax2=fig.add_subplot(gs[1],sharex=ax); res=y-np.log10(g_emg(gb,best)); bins=np.linspace(-12.2,-8,18)
    cent=[]; med=[]; lo=[]; hi=[]
    for l,r in zip(bins[:-1],bins[1:]):
        m=(x>=l)&(x<r)
        if m.sum()>5:
            rr=res[m]; cent.append((l+r)/2); med.append(np.median(rr)); lo.append(np.percentile(rr,16)); hi.append(np.percentile(rr,84))
    cent=np.asarray(cent); med=np.asarray(med); lo=np.asarray(lo); hi=np.asarray(hi)
    ax2.axhline(0,lw=1); ax2.fill_between(cent,lo,hi,alpha=.15); ax2.plot(cent,med,marker='o',ms=3,lw=1)
    ax2.set_xlabel(r'$\log_{10} g_{\rm bar}\;[\mathrm{m\,s^{-2}}]$'); ax2.set_ylabel('resid.\n[dex]'); ax2.set_ylim(-.35,.35)
    fig.savefig(args.fig/'fig_rar.pdf',bbox_inches='tight'); plt.close(fig)
    print(f'N = {len(y)}')
    for r in results: print(r)

if __name__=='__main__': main()
