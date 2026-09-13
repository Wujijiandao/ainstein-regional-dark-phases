#!/usr/bin/env python3
"""Matched-estimator control: Zel'dovich vs nonlinear PM at a=1.

Both branches start from the *same* Gaussian delta_L field and are analysed with
identical Eulerian CIC density, 4 Mpc/h velocity smoothing, watershed basin
segmentation, Gamma estimator, and finite-size threshold.  This isolates the
increment due to nonlinear PM evolution from changes in the analysis pipeline.
"""
from __future__ import annotations
import argparse,csv,json,sys,math
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from analysis_nonlinear_pm_basin_lineage import (
    BOX,A_INIT,R_SMOOTH,fourier_grid,gaussian_linear_field,displacement_from_delta,
    grid_positions,pm_step,velocity_divergence_theta,segment,basin_stats,sigma8_grid
)


def analyze(label,pos,p,fg,persistence=0.0):
    rho,rho_sm,theta=velocity_divergence_theta(pos,p,1.0,fg,R_SMOOTH)
    labels,rms=segment(rho_sm-1.0,persistence)
    st=basin_stats(labels,rho,theta,BOX/fg.n)
    counts=st['counts']; act=st['active']; eq=st['gamma']>=1.0
    return {
        'branch':label,'Ngrid':fg.n,'n_basins':st['nlab'],'delta_sm_rms':rms,
        'sigma8_CIC':sigma8_grid(rho-1.0,fg),
        'gamma_ge_1_volume_fraction':float(counts[eq].sum()/labels.size),
        'finite_size_active_volume_fraction':float(counts[act].sum()/labels.size),
        'n_finite_active':int(act.sum()),
        'median_R_active_hMpc':float(np.median(st['reff'][act])) if act.any() else math.nan,
        'median_gamma_active':float(np.median(st['gamma'][act])) if act.any() else math.nan,
    }


def run(seed,n,steps,persistence):
    fg=fourier_grid(n,BOX)
    d0=gaussian_linear_field(seed,fg)
    s=displacement_from_delta(d0,fg)
    q=grid_positions(n,BOX)
    # matched ZA at a=1: x=q+s, p=a^3 H s=s in EdS at a=1.
    za_pos=(q+s)%BOX; za_p=s.copy()
    za=analyze('ZA_matched',za_pos,za_p,fg,persistence)
    # PM from a_init.
    pos=(q+A_INIT*s)%BOX; p=(A_INIT**1.5)*s
    nodes=np.unique(np.round(np.linspace(A_INIT,1.0,steps+1),10)); nodes.sort(); a=A_INIT
    for a2 in nodes[1:]:
        pos,p=pm_step(pos,p,a,float(a2),fg); a=float(a2)
    pm=analyze('PM_nonlinear',pos,p,fg,persistence)
    return [za,pm]


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results/matched_control')
    ap.add_argument('--n',type=int,default=80); ap.add_argument('--steps',type=int,default=38)
    ap.add_argument('--seeds',type=int,nargs='+',default=[11]); ap.add_argument('--persistence',type=float,default=0.0)
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for seed in args.seeds:
        rr=run(seed,args.n,args.steps,args.persistence)
        for r in rr:r['seed']=seed;r['persistence_rms']=args.persistence
        rows.extend(rr)
    with (args.out/'matched_za_pm_control.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    pairs=[]
    for seed in args.seeds:
        z=next(r for r in rows if r['seed']==seed and r['branch']=='ZA_matched')
        p=next(r for r in rows if r['seed']==seed and r['branch']=='PM_nonlinear')
        pairs.append({'seed':seed,'ZA_phi':z['finite_size_active_volume_fraction'],'PM_phi':p['finite_size_active_volume_fraction'],
                      'PM_minus_ZA':p['finite_size_active_volume_fraction']-z['finite_size_active_volume_fraction'],
                      'PM_over_ZA':(p['finite_size_active_volume_fraction']/z['finite_size_active_volume_fraction']) if z['finite_size_active_volume_fraction']>0 else math.nan})
    text='Matched ZA vs PM control\n'+json.dumps({'Ngrid':args.n,'persistence_rms':args.persistence,'pairs':pairs,'rows':rows},indent=2)+'\n'
    (args.out/'matched_za_pm_control_summary.txt').write_text(text,encoding='utf-8')
    print(text)
if __name__=='__main__':main()
