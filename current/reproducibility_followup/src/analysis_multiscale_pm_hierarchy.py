#!/usr/bin/env python3
"""Multiscale watershed-hierarchy audit on a nonlinear PM snapshot.

This is deliberately a *selector stress test*, not a calibrated cosmological
prediction.  It asks whether the low active fraction of the single-scale v1.2
basin tests is merely an artifact of choosing R_smooth=4 h^-1 Mpc.

A fixed preregistered Gaussian scale ladder is used.  At each scale, watershed
basins are evaluated with the same Gamma and finite-size criterion.  Fine-to-
coarse parent links are assigned by maximum Eulerian volume overlap, following
the standard hierarchical-watershed idea.  We report the union of active nodes
across scales as a deliberately permissive upper envelope, together with a
strict non-overlap maximal-active selection built from coarse-to-fine coverage.
"""
from __future__ import annotations
import argparse,csv,math,json,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from analysis_nonlinear_pm_basin_lineage import (
    BOX,A_INIT,SNAPSHOTS,C_SIGMA,XI0,BETA_GAMMA,
    fourier_grid,gaussian_linear_field,displacement_from_delta,grid_positions,
    pm_step,velocity_divergence_theta,gaussian_smooth,segment,basin_stats,
    sigma8_grid,density_ratio
)

SCALES=(2.0,4.0,8.0,12.0,16.0)
PERSISTENCE=0.0


def evolve_final(seed:int,n:int,steps:int):
    fg=fourier_grid(n,BOX)
    delta0=gaussian_linear_field(seed,fg)
    s=displacement_from_delta(delta0,fg)
    pos=(grid_positions(n,BOX)+A_INIT*s)%BOX
    p=(A_INIT**1.5)*s
    nodes=np.unique(np.round(np.concatenate([np.linspace(A_INIT,1.0,steps+1),[1.0]]),10)); nodes.sort()
    a=A_INIT
    for a2 in nodes[1:]:
        pos,p=pm_step(pos,p,a,float(a2),fg); a=float(a2)
    # use the 4 Mpc velocity smoothing as the operational velocity estimator,
    # then test density-domain scale separately. This avoids changing both
    # state estimator and domain definition at once.
    rho,rho4,theta=velocity_divergence_theta(pos,p,1.0,fg,4.0)
    return fg,rho,theta,pos,{
        'sigma8_linear_field':sigma8_grid(delta0,fg),
        'sigma8_pm_final':sigma8_grid(rho-1.0,fg),
        'mass_closure':float(rho.mean())
    }


def overlap_parent(child:np.ndarray,parent:np.ndarray,nchild:int,nparent:int):
    # volume-overlap dominant parent for every fine-scale child
    packed=child.ravel().astype(np.int64)*(nparent+1)+parent.ravel().astype(np.int64)
    bc=np.bincount(packed)
    best=np.zeros(nchild+1,dtype=np.int32)
    bestn=np.zeros(nchild+1,dtype=np.int64)
    for code in np.nonzero(bc)[0]:
        c=code//(nparent+1); p=code%(nparent+1)
        if c==0 or p==0: continue
        if bc[code]>bestn[c]:
            bestn[c]=bc[code]; best[c]=p
    counts=np.bincount(child.ravel(),minlength=nchild+1)
    frac=np.zeros(nchild+1,float)
    nz=counts>0; frac[nz]=bestn[nz]/counts[nz]
    return best,frac


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results/multiscale')
    ap.add_argument('--n',type=int,default=64)
    ap.add_argument('--steps',type=int,default=38)
    ap.add_argument('--seed',type=int,default=11)
    ap.add_argument('--scales',type=float,nargs='+',default=list(SCALES))
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    scales=tuple(sorted(args.scales))
    fg,rho,theta,pos,val=evolve_final(args.seed,args.n,args.steps)
    dx=BOX/args.n
    cats=[]; rows=[]
    for R in scales:
        dsm=gaussian_smooth(rho-1.0,fg,R)
        labels,rms=segment(dsm,PERSISTENCE)
        st=basin_stats(labels,rho,theta,dx)
        eq=(st['gamma']>=1.0)
        act=st['active']
        phi_eq=float(st['counts'][eq].sum()/labels.size)
        phi_act=float(st['counts'][act].sum()/labels.size)
        cats.append({'R':R,'labels':labels,'stats':st})
        rows.append({'seed':args.seed,'Ngrid':args.n,'smooth_R_hMpc':R,'n_basins':st['nlab'],
                     'delta_sm_rms':rms,'gamma_ge_1_volume_fraction':phi_eq,
                     'finite_size_active_volume_fraction':phi_act,
                     'median_R_basin_hMpc':float(np.median(st['reff'])),
                     'median_R_active_hMpc':float(np.median(st['reff'][act])) if act.any() else math.nan,
                     'median_gamma_active':float(np.median(st['gamma'][act])) if act.any() else math.nan})
    # fine->coarse parent links and hierarchy quality
    links=[]
    for i in range(len(cats)-1):
        fine,coarse=cats[i],cats[i+1]
        parent,frac=overlap_parent(fine['labels'],coarse['labels'],fine['stats']['nlab'],coarse['stats']['nlab'])
        links.append((parent,frac))

    # permissive union upper envelope: any active node at any scale.
    union=np.zeros((args.n,args.n,args.n),dtype=bool)
    raw_union=np.zeros((args.n,args.n,args.n),dtype=bool)
    for c in cats:
        union |= c['stats']['active'][c['labels']-1]
        raw_union |= (c['stats']['gamma']>=1.0)[c['labels']-1]
    phi_union=float(union.mean())
    phi_raw_union=float(raw_union.mean())

    # strict coarse-to-fine maximal non-overlap: once a coarse active node claims a
    # cell, finer active nodes may only add cells not already claimed.
    claimed=np.zeros_like(union)
    additions=[]
    for c in reversed(cats):
        cells=c['stats']['active'][c['labels']-1]
        add=cells & ~claimed
        additions.append((c['R'],float(add.mean())))
        claimed |= cells
    phi_maximal=float(claimed.mean())

    # hierarchy diagnostics: how laminar are adjacent scale partitions?
    linkrows=[]
    for i,(parent,frac) in enumerate(links):
        f=cats[i]
        valid=frac[1:]
        linkrows.append({'fine_R_hMpc':cats[i]['R'],'coarse_R_hMpc':cats[i+1]['R'],
                         'median_dominant_parent_overlap':float(np.median(valid)),
                         'fraction_children_parent_overlap_ge_0p8':float(np.mean(valid>=0.8)),
                         'fraction_children_parent_overlap_ge_0p5':float(np.mean(valid>=0.5))})

    with (args.out/'multiscale_pm_hierarchy_table.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    with (args.out/'multiscale_pm_hierarchy_links.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(linkrows[0].keys()));w.writeheader();w.writerows(linkrows)
    out={
        'seed':args.seed,'Ngrid':args.n,'scales_hMpc':scales,'persistence_rms':PERSISTENCE,
        'validation':val,'single_scale_rows':rows,'adjacent_scale_link_quality':linkrows,
        'permissive_any_scale_gamma_ge1_union_volume_fraction':phi_raw_union,
        'permissive_any_scale_active_union_volume_fraction':phi_union,
        'maximal_nonoverlap_active_union_volume_fraction':phi_maximal,
        'coarse_to_fine_incremental_volume_fractions':additions,
        'warning':'The any-scale union is a preregistered scale-ladder upper envelope, not a calibrated phase abundance. Adjacent watershed partitions are not exactly laminar; overlap diagnostics quantify this.'
    }
    (args.out/'multiscale_pm_hierarchy_summary.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    text=(
        'Nonlinear PM multiscale basin-hierarchy audit\n'
        f'seed={args.seed}\nNgrid={args.n}\nscales_hMpc={scales}\n'
        f'sigma8_pm_final={val["sigma8_pm_final"]:.9f}\n'
        + ''.join(f'R={r["smooth_R_hMpc"]:.1f}: phi_Gamma>=1={r["gamma_ge_1_volume_fraction"]:.6f}, phi_finite={r["finite_size_active_volume_fraction"]:.6f}, nbasin={r["n_basins"]}\n' for r in rows)
        + f'any_scale_gamma_ge1_union_phi={phi_raw_union:.6f}\nany_scale_union_phi={phi_union:.6f}\nmaximal_nonoverlap_union_phi={phi_maximal:.6f}\n'
        + ''.join(f'link {x["fine_R_hMpc"]:.1f}->{x["coarse_R_hMpc"]:.1f}: median_parent_overlap={x["median_dominant_parent_overlap"]:.4f}, frac>=0.8={x["fraction_children_parent_overlap_ge_0p8"]:.4f}\n' for x in linkrows)
        + 'interpretation=If even the preregistered multiscale upper envelope remains small, the obstruction is not a single-smoothing-scale artifact. If it rises sharply, the next task is to replace the permissive scale union by a physically fixed hierarchical selection rule.\n'
    )
    (args.out/'multiscale_pm_hierarchy_summary.txt').write_text(text,encoding='utf-8')
    print(text)

if __name__=='__main__': main()
