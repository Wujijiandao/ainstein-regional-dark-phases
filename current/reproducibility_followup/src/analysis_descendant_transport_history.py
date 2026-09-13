#!/usr/bin/env python3
"""Descendant-anchored history via transport of the terminal watershed partition.

Compared with inverse-density particle-volume weighting, this construction keeps
exactly the same Eulerian basin estimator as the v1.2/pilot selector. At a=1 the
terminal watershed partition is used verbatim. At earlier snapshots each
particle carries the ID of its a=1 basin backward; the current grid is assigned
by the modal terminal-basin ID of particles in each cell, with empty cells filled
from the nearest occupied cell using a periodic KD-tree. Thus every snapshot is
a full descendant-anchored partition with fixed basin identities and no greedy
snapshot-to-snapshot catalogue matching.

Exploratory obstruction test only; the rasterized transported partition is not
claimed to be a physical material boundary.
"""
from __future__ import annotations
import argparse, csv, json, math, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import analysis_nonlinear_pm_basin_lineage as pm
from analysis_descendant_anchored_history import simulate_trajectory, grid_catalogue_final


def transported_partition(pos:np.ndarray, particle_final_labels:np.ndarray, n:int, box:float, nlab:int):
    dx=box/n
    ijk=np.floor(pos/dx).astype(np.int64)%n
    cell=(ijk[:,0]*n+ijk[:,1])*n+ijk[:,2]
    code=cell.astype(np.int64)*(nlab+1)+particle_final_labels.astype(np.int64)
    u,c=np.unique(code,return_counts=True)
    cells=u//(nlab+1); labs=u%(nlab+1)
    # choose highest-count label per occupied cell, tie -> lower label deterministically
    order=np.lexsort((labs,-c,cells))
    cells_s=cells[order]; labs_s=labs[order]
    first=np.r_[True,cells_s[1:]!=cells_s[:-1]]
    chosen_cells=cells_s[first]; chosen_labs=labs_s[first]
    flat=np.zeros(n**3,dtype=np.int32)
    flat[chosen_cells]=chosen_labs.astype(np.int32)
    empty=np.flatnonzero(flat==0)
    if len(empty):
        occ=chosen_cells
        oz=occ%n; oy=(occ//n)%n; ox=occ//(n*n)
        ecoord=np.column_stack([empty//(n*n),(empty//n)%n,empty%n]).astype(float)
        ocoord=np.column_stack([ox,oy,oz]).astype(float)
        tree=cKDTree(ocoord,boxsize=float(n))
        _,idx=tree.query(ecoord,k=1,workers=-1)
        flat[empty]=chosen_labs[idx]
    return flat.reshape((n,n,n)), float(len(chosen_cells)/(n**3))


def run_seed(seed:int,n:int,steps:int,persistence_values:tuple[float,...],R:float):
    fg,delta0,sigma0,snaps=simulate_trajectory(seed,n,steps)
    final_pos,final_p=snaps[1.0]
    rows=[]; basin_rows=[]; validation={'seed':seed,'Ngrid':n,'R_smooth_hMpc':R,'sigma8_linear_target_a1':sigma0}
    dx=fg.box/fg.n
    for pers in persistence_values:
        final_labels,particle_final,final_st,rms=grid_catalogue_final(final_pos,final_p,fg,pers,R)
        nlab=final_st['nlab']
        ever=np.zeros(nlab,dtype=bool); first=np.full(nlab,np.nan)
        for a in pm.SNAPSHOTS:
            pos,p=snaps[round(a,10)]
            rho,rho_sm,theta=pm.velocity_divergence_theta(pos,p,a,fg,R)
            if abs(a-1.0)<1e-12:
                labels=final_labels
                occupied_frac=1.0
            else:
                labels,occupied_frac=transported_partition(pos,particle_final,n,fg.box,nlab)
            st=pm.basin_stats(labels,rho,theta,dx)
            # st nlab should match terminal; assert fixed identity space
            assert st['nlab']==nlab
            active=st['active']
            newly=active & ~ever; first[newly]=a; ever|=active
            phi=float(st['counts'][active].sum()/labels.size)
            phi_abs=float(st['counts'][ever].sum()/labels.size)
            rows.append({
                'seed':seed,'Ngrid':n,'persistence_rms_final':pers,'R_smooth_hMpc':R,'a':a,
                'n_final_basins':nlab,'directly_occupied_cell_fraction_before_fill':occupied_frac,
                'transport_instantaneous_active_volume_fraction':phi,
                'transport_absorbing_active_volume_fraction':phi_abs,
                'active_basin_count':int(active.sum()),'ever_active_basin_count':int(ever.sum()),
                'median_R_all_hMpc':float(np.median(st['reff'])),
                'median_R_active_hMpc':float(np.median(st['reff'][active])) if np.any(active) else math.nan,
                'median_gamma_all':float(np.median(st['gamma'][np.isfinite(st['gamma'])])) if np.any(np.isfinite(st['gamma'])) else math.nan,
                'median_gamma_active':float(np.median(st['gamma'][active])) if np.any(active) else math.nan,
                'sigma8_nonlinear':pm.sigma8_grid(rho-1.0,fg),
            })
            for i in range(nlab):
                basin_rows.append({
                    'seed':seed,'Ngrid':n,'persistence_rms_final':pers,'R_smooth_hMpc':R,'a':a,
                    'final_basin_id':i+1,'volume_fraction':float(st['counts'][i]/labels.size),
                    'R_eff_hMpc':float(st['reff'][i]),'mean_density_ratio':float(st['mean_rho'][i]),
                    'Hloc_over_H':float(st['hratio'][i]),'Gamma':float(st['gamma'][i]),
                    'Gamma_grow':float(st['gamma_grow'][i]),'instantaneous_active':int(active[i]),
                    'absorbing_active':int(ever[i]),
                })
        validation[f'p{pers:g}_terminal_grid_phi_exact']=float(final_st['counts'][final_st['active']].sum()/final_labels.size)
        validation[f'p{pers:g}_first_passage_counts']={str(a):int(np.count_nonzero(np.isclose(first,a,equal_nan=False))) for a in pm.SNAPSHOTS}
        validation[f'p{pers:g}_ever_active_final_basin_fraction']=float(np.mean(np.isfinite(first)))
    return rows,basin_rows,validation


def write(out:Path,rows,basin_rows,vals):
    out.mkdir(parents=True,exist_ok=True)
    with (out/'descendant_transport_summary_table.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    with (out/'descendant_transport_basin_histories.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(basin_rows[0].keys()));w.writeheader();w.writerows(basin_rows)
    (out/'descendant_transport_validation.json').write_text(json.dumps(vals,indent=2),encoding='utf-8')
    lines=['AInstein descendant-anchored transported terminal-partition history',
           'terminal estimator reproduces the exact a=1 Eulerian watershed catalogue by construction',
           'earlier partition=terminal particle basin IDs rasterized at current positions; periodic nearest-occupied fill',
           'no pairwise snapshot catalogue matching']
    for pers in sorted({r['persistence_rms_final'] for r in rows}):
        sub=[r for r in rows if r['persistence_rms_final']==pers]
        for a in pm.SNAPSHOTS:
            ss=[r for r in sub if abs(r['a']-a)<1e-12]
            for key in ['transport_instantaneous_active_volume_fraction','transport_absorbing_active_volume_fraction']:
                v=np.array([r[key] for r in ss]);lines.append(f'p{pers:g}_a{a:g}_{key}_range={v.min():.8g},{v.max():.8g};mean={v.mean():.8g}')
    (out/'descendant_transport_summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    ap=argparse.ArgumentParser();
    ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results'/'descendant_transport')
    ap.add_argument('--n',type=int,default=64);ap.add_argument('--steps',type=int,default=38)
    ap.add_argument('--seeds',type=int,nargs='+',default=[11,23,47,91]);ap.add_argument('--persistence',type=float,nargs='+',default=[0.0,0.1,0.2])
    ap.add_argument('--R-smooth',type=float,default=pm.R_SMOOTH);args=ap.parse_args()
    rows=[];br=[];vals=[]
    for seed in args.seeds:
        print(f'=== descendant-transport seed={seed} N={args.n} R={args.R_smooth} ===',flush=True)
        r,b,v=run_seed(seed,args.n,args.steps,tuple(args.persistence),args.R_smooth);rows+=r;br+=b;vals.append(v);print(json.dumps(v,indent=2),flush=True)
    write(args.out,rows,br,vals);print((args.out/'descendant_transport_summary.txt').read_text(),flush=True)
if __name__=='__main__':main()
