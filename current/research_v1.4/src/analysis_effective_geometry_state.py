#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INTERVALS=[(0.25,0.50),(0.50,0.75),(0.75,1.00)]

def interval_table(df: pd.DataFrame, label: str) -> pd.DataFrame:
    # Current legacy selector defines terminal-active descendants at a=1.
    final=df[np.isclose(df.a,1.0) & (df.instantaneous_active==1)][['seed','final_basin_id','volume_fraction']].copy()
    final=final.rename(columns={'volume_fraction':'final_volume_weight'})
    rows=[]
    for a0,a1 in INTERVALS:
        d0=df[np.isclose(df.a,a0)].merge(final,on=['seed','final_basin_id'])
        d1=df[np.isclose(df.a,a1)].merge(final,on=['seed','final_basin_id'])
        m=d0.merge(d1,on=['seed','final_basin_id','final_volume_weight'],suffixes=('_0','_1'))
        # regional volume scale factor: a_D propto a * (coordinate volume)^(1/3)
        aratio=(a1/a0)*(m.volume_fraction_1/m.volume_fraction_0)**(1/3)
        valid=(m.Gamma_0>0)&(m.Gamma_1>0)&(aratio>0)&np.isfinite(aratio)
        m=m.loc[valid].copy(); aratio=aratio.loc[valid]
        m['aD_ratio']=aratio
        m['nu_gamma']=np.log(m.Gamma_1/m.Gamma_0)/np.log(aratio)
        m['w_gamma_proxy']=-m.nu_gamma/3.0
        # midpoint proxy for the exact Buchert acceleration boundary Gamma*(nu-1)>1
        m['Gamma_mid_geom']=np.sqrt(m.Gamma_0*m.Gamma_1)
        m['acceleration_score_proxy']=m.Gamma_mid_geom*(m.nu_gamma-1.0)-1.0
        m['interval']=f'{a0:g}-{a1:g}'
        m['dataset']=label
        rows.append(m[['dataset','interval','seed','final_basin_id','final_volume_weight','Gamma_0','Gamma_1','aD_ratio','nu_gamma','w_gamma_proxy','Gamma_mid_geom','acceleration_score_proxy']])
    return pd.concat(rows,ignore_index=True)

def weighted_quantile(x,w,q):
    x=np.asarray(x,float); w=np.asarray(w,float); q=np.asarray(q,float)
    o=np.argsort(x); x=x[o]; w=w[o]
    cw=np.cumsum(w); cw=(cw-0.5*w)/w.sum()
    return np.interp(q,cw,x)

def summarize(tab: pd.DataFrame):
    out=[]
    for (dataset,interval),g in tab.groupby(['dataset','interval'],sort=False):
        w=g.final_volume_weight.to_numpy(float)
        nu=g.nu_gamma.to_numpy(float); wp=g.w_gamma_proxy.to_numpy(float); A=g.acceleration_score_proxy.to_numpy(float)
        out.append({
            'dataset':dataset,'interval':interval,'n_terminal_active_basins':int(len(g)),
            'weight_sum':float(w.sum()),
            'nu_volume_weighted_mean':float(np.average(nu,weights=w)),
            'nu_weighted_q10':float(weighted_quantile(nu,w,[.1])[0]),
            'nu_weighted_median':float(weighted_quantile(nu,w,[.5])[0]),
            'nu_weighted_q90':float(weighted_quantile(nu,w,[.9])[0]),
            'w_gamma_volume_weighted_mean':float(np.average(wp,weights=w)),
            'w_gamma_weighted_median':float(weighted_quantile(wp,w,[.5])[0]),
            'fraction_weight_near_vacuum_w_m1_pm0p2':float(np.average(((wp>-1.2)&(wp<-.8)).astype(float),weights=w)),
            'fraction_weight_acceleration_score_positive':float(np.average((A>0).astype(float),weights=w)),
        })
    return pd.DataFrame(out)

def make_plot(summary: pd.DataFrame, out: Path):
    fig,ax=plt.subplots(figsize=(7.2,4.5))
    datasets=list(summary.dataset.unique())
    intervals=list(summary.interval.unique())
    x=np.arange(len(intervals),dtype=float)
    width=0.34
    for j,d in enumerate(datasets):
        s=summary[summary.dataset==d].set_index('interval').loc[intervals]
        ax.bar(x+(j-(len(datasets)-1)/2)*width,s.w_gamma_volume_weighted_mean,width=width,label=d)
    ax.axhline(-1/3,linestyle='--',linewidth=1,label=r'curvature-like $w=-1/3$')
    ax.axhline(-1,linestyle=':',linewidth=1,label=r'vacuum-like $w=-1$')
    ax.set_xticks(x,intervals)
    ax.set_xlabel('scale-factor interval')
    ax.set_ylabel(r'volume-weighted $w_{\Gamma}\equiv-\nu_{\Gamma}/3$')
    ax.set_title('History coordinate of legacy-selector terminal basins')
    ax.legend(frameon=False,fontsize=8)
    fig.tight_layout()
    fig.savefig(out/'fig_effective_geometry_history.pdf',bbox_inches='tight')
    fig.savefig(out/'fig_effective_geometry_history.png',dpi=200,bbox_inches='tight')
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--n64',type=Path,required=True)
    ap.add_argument('--n80',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    d64=pd.read_csv(args.n64); d80=pd.read_csv(args.n80)
    tabs=[interval_table(d64,'64^3, four seeds'),interval_table(d80,'80^3, seed 11')]
    tab=pd.concat(tabs,ignore_index=True)
    summary=summarize(tab)
    tab.to_csv(args.out/'effective_geometry_interval_table.csv',index=False)
    summary.to_csv(args.out/'effective_geometry_summary.csv',index=False)
    rec={'identity':{
        'rho_G':'-(Q+R)/(16 pi G)','p_G':'-(Q-R/3)/(16 pi G)','Gamma':'rho_G/rho_m',
        'nu':'d ln Gamma / d ln a_D','w_G':'-nu/3','vacuum_line':'nu=3','curvature_line':'nu=1',
        'regional_acceleration':'Gamma*(nu-1)>1'},
        'summary':summary.to_dict(orient='records'),
        'caveat':'PM Gamma is an operational Newtonian continuation, not a direct relativistic measurement of Buchert Q and R. w_Gamma is therefore a trajectory-consistency proxy, not a measured physical EOS.'}
    (args.out/'effective_geometry_summary.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
    make_plot(summary,args.out)
    print(summary.to_string(index=False))

if __name__=='__main__': main()
