#!/usr/bin/env python3
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.input)
    seeds=[int(x) for x in sorted(df.seed.unique())]
    rows=[]
    all_traj=[]
    for seed in seeds:
        d=df[df.seed==seed]
        ids=sorted(d.final_basin_id.unique())
        for bid in ids:
            q=d[d.final_basin_id==bid].sort_values('a')
            if len(q)!=4: continue
            gamma=q.Gamma.to_numpy(float); grow=q.Gamma_grow.to_numpy(float); active=q.instantaneous_active.to_numpy(int).astype(bool)
            final_active=bool(active[-1])
            status_monotone=bool(np.all(np.diff(active.astype(int))>=0))
            gamma_monotone=bool(np.all(np.diff(gamma)>=-1e-10))
            fp=float(q.a.to_numpy(float)[np.argmax(active)]) if np.any(active) else np.nan
            final_vol=float(q.volume_fraction.iloc[-1])
            rows.append({'seed':seed,'final_basin_id':bid,'final_active':int(final_active),'first_passage_a':fp,
                         'active_status_monotone':int(status_monotone),'Gamma_monotone':int(gamma_monotone),
                         'final_volume_fraction':final_vol,'final_Gamma':float(gamma[-1]),'final_Gamma_grow':float(grow[-1])})
            if final_active:
                all_traj.append((seed,bid,q.a.to_numpy(float),gamma-grow,final_vol))
    outdf=pd.DataFrame(rows);outdf.to_csv(args.out/'descendant_first_passage_basin_audit.csv',index=False)
    fa=outdf[outdf.final_active==1]
    summary={
        'seeds':seeds,
        'n_basins_total':int(len(outdf)),
        'n_final_active_basins':int(len(fa)),
        'fraction_all_basins_status_monotone':float(outdf.active_status_monotone.mean()),
        'fraction_final_active_status_monotone':float(fa.active_status_monotone.mean()) if len(fa) else None,
        'fraction_final_active_Gamma_monotone':float(fa.Gamma_monotone.mean()) if len(fa) else None,
        'first_passage_count_final_active':{str(a):int(np.count_nonzero(np.isclose(fa.first_passage_a,a,equal_nan=False))) for a in [0.25,0.5,0.75,1.0]},
        'first_passage_final_volume_weight':{}
    }
    denom=fa.final_volume_fraction.sum()
    for a in [0.25,0.5,0.75,1.0]:
        summary['first_passage_final_volume_weight'][str(a)]=float(fa.loc[np.isclose(fa.first_passage_a,a,equal_nan=False),'final_volume_fraction'].sum()/denom) if denom>0 else 0.0
    (args.out/'descendant_first_passage_audit.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    lines=['Descendant first-passage audit']+[f'{k}={v}' for k,v in summary.items()]
    (args.out/'descendant_first_passage_audit.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')

    # Plot Gamma-Gamma_grow trajectories for final-active basins; linewidth encodes final volume weakly.
    fig,ax=plt.subplots(figsize=(6.4,4.2))
    for seed,bid,a,y,v in all_traj:
        ax.plot(a,y,alpha=0.35,linewidth=0.8+10*min(v,0.01))
    ax.axhline(0,linewidth=1)
    ax.set_xlabel('scale factor a');ax.set_ylabel(r'$\Gamma_\Omega-\Gamma_{\rm grow}$')
    ax.set_title('Descendant-anchored trajectories of final-active basins')
    fig.tight_layout();fig.savefig(args.out/'fig_descendant_gamma_first_passage.pdf');fig.savefig(args.out/'fig_descendant_gamma_first_passage.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
