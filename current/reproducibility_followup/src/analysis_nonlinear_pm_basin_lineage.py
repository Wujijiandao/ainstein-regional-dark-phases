#!/usr/bin/env python3
"""Exploratory nonlinear PM + persistent basin-lineage stress test.

Purpose
-------
This is a follow-up research prototype to AInstein v1.2.  It asks whether
fully nonlinear collisionless gravitational evolution can materially raise the
finite-size-active regional V-phase filling fraction relative to the frozen
Gaussian-watershed and Zel'dovich-history obstruction tests.

This is NOT a precision LambdaCDM simulation.  It uses an Einstein-de Sitter
particle-mesh (PM) evolution with a Planck-like z=0 linear power-spectrum shape
and sigma8 normalization, then reconstructs Eulerian watershed basins at four
snapshots.  Basin lineage is tracked by particle-ID overlap, so mergers and
history inheritance are explicit rather than imposed on a fixed hierarchy.

Selector
--------
For each Eulerian basin Omega we estimate
    Gamma_Omega = (H_loc/H)^2 * (rho_bg/rho_Omega) - 1,
with H_loc/H = 1 + <div u/(aH)>_Omega/3,
and apply the same conditional finite-size growth threshold used in v1.2,
    Gamma_grow(R) = 1 + 2 C_sigma xi0/(beta_Gamma R).

Two history rules are reported in addition to strict instantaneous activity:
  * lineage50: current basin is active, or >=50% of its particles descend from
    persistent-active parent basins at the preceding snapshot;
  * lineage10: same, with a deliberately permissive 10% inheritance threshold.
The pair is a sensitivity bracket, not a physical law.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from dataclasses import dataclass

import numpy as np
from skimage.morphology import h_minima, local_minima
from skimage.segmentation import watershed

# Frozen v1.2 constants / Planck-like shape parameters.
OMEGA_M_SHAPE = 0.315
H_REDUCED = 0.674
N_S = 0.965
SIGMA8 = 0.811
C_SIGMA = 0.396830116706
XI0 = 3.12938          # h^-1 Mpc, same conditional benchmark as v1.2
BETA_GAMMA = 1.0       # same matched-normalization benchmark as v1.2
BOX = 256.0            # h^-1 Mpc
R_SMOOTH = 4.0         # h^-1 Mpc, same coarse-graining scale as v1.2 3D tests
A_INIT = 0.05
SNAPSHOTS = (0.25, 0.50, 0.75, 1.00)
PERSISTENCE = (0.0, 0.10, 0.20)


def transfer_eh0(k_hmpc: np.ndarray) -> np.ndarray:
    q = np.asarray(k_hmpc, dtype=float) / (OMEGA_M_SHAPE * H_REDUCED)
    L0 = np.log(2.0 * math.e + 1.8 * q)
    C0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
    return L0 / (L0 + C0 * q*q)


def top_hat(x: np.ndarray) -> np.ndarray:
    out = np.ones_like(x, dtype=float)
    m = np.abs(x) > 1e-8
    xm = x[m]
    out[m] = 3.0 * (np.sin(xm) - xm*np.cos(xm)) / xm**3
    return out


@dataclass
class FourierGrid:
    n: int
    box: float
    kx: np.ndarray
    ky: np.ndarray
    kz: np.ndarray
    k2: np.ndarray
    invk2: np.ndarray
    kmag: np.ndarray


def fourier_grid(n: int, box: float) -> FourierGrid:
    dx = box/n
    kx1 = 2.0*math.pi*np.fft.fftfreq(n, d=dx)
    ky1 = 2.0*math.pi*np.fft.fftfreq(n, d=dx)
    kz1 = 2.0*math.pi*np.fft.rfftfreq(n, d=dx)
    kx,ky,kz = np.meshgrid(kx1,ky1,kz1,indexing='ij')
    k2 = kx*kx+ky*ky+kz*kz
    inv = np.zeros_like(k2)
    m=k2>0
    inv[m]=1.0/k2[m]
    return FourierGrid(n,box,kx,ky,kz,k2,inv,np.sqrt(k2))


def gaussian_linear_field(seed: int, fg: FourierGrid) -> np.ndarray:
    """Return delta_L extrapolated to a=1 and normalized to sigma8."""
    rng=np.random.default_rng(seed)
    white=rng.normal(size=(fg.n,fg.n,fg.n))
    W=np.fft.rfftn(white)
    T=transfer_eh0(fg.kmag)
    P=np.zeros_like(fg.kmag)
    m=fg.kmag>0
    P[m]=fg.kmag[m]**N_S*T[m]**2
    d=np.fft.irfftn(W*np.sqrt(P),s=(fg.n,fg.n,fg.n),axes=(0,1,2)).real
    d-=d.mean()
    F=np.fft.rfftn(d)
    sm8=np.fft.irfftn(F*top_hat(fg.kmag*8.0),s=d.shape,axes=(0,1,2)).real
    d*=SIGMA8/sm8.std()
    return d


def gaussian_smooth(arr: np.ndarray, fg: FourierGrid, R: float) -> np.ndarray:
    F=np.fft.rfftn(arr)
    return np.fft.irfftn(F*np.exp(-0.5*(fg.kmag*R)**2),s=arr.shape,axes=(0,1,2)).real


def sigma8_grid(delta: np.ndarray, fg: FourierGrid) -> float:
    sm=np.fft.irfftn(np.fft.rfftn(delta)*top_hat(fg.kmag*8.0),s=delta.shape,axes=(0,1,2)).real
    return float(sm.std())


def displacement_from_delta(delta0: np.ndarray, fg: FourierGrid) -> np.ndarray:
    F=np.fft.rfftn(delta0)
    out=np.empty((delta0.size,3),dtype=np.float64)
    for j,kj in enumerate((fg.kx,fg.ky,fg.kz)):
        sk=np.fft.irfftn(1j*kj*fg.invk2*F,s=delta0.shape,axes=(0,1,2)).real
        out[:,j]=sk.ravel()
    return out


def grid_positions(n:int, box:float) -> np.ndarray:
    dx=box/n
    q=(np.arange(n,dtype=float)+0.5)*dx
    x,y,z=np.meshgrid(q,q,q,indexing='ij')
    return np.column_stack([x.ravel(),y.ravel(),z.ravel()])


def cic_components(pos: np.ndarray, n:int, box:float):
    dx=box/n
    g=pos/dx-0.5  # grid centers at (i+1/2) dx
    i0=np.floor(g).astype(np.int64)
    f=g-i0
    i0%=n
    i1=(i0+1)%n
    # yield flattened index and weight for 8 corners
    for bx in (0,1):
        ix=i1[:,0] if bx else i0[:,0]
        wx=f[:,0] if bx else (1.0-f[:,0])
        for by in (0,1):
            iy=i1[:,1] if by else i0[:,1]
            wy=f[:,1] if by else (1.0-f[:,1])
            for bz in (0,1):
                iz=i1[:,2] if bz else i0[:,2]
                wz=f[:,2] if bz else (1.0-f[:,2])
                idx=(ix*n+iy)*n+iz
                yield idx, wx*wy*wz


def cic_deposit(pos:np.ndarray,n:int,box:float,values:np.ndarray|None=None) -> np.ndarray:
    size=n**3
    out=np.zeros(size,dtype=np.float64)
    if values is None:
        for idx,w in cic_components(pos,n,box):
            out += np.bincount(idx,weights=w,minlength=size)
    else:
        vals=np.asarray(values,dtype=float)
        for idx,w in cic_components(pos,n,box):
            out += np.bincount(idx,weights=w*vals,minlength=size)
    return out.reshape((n,n,n))


def cic_interpolate(grid:np.ndarray,pos:np.ndarray,n:int,box:float) -> np.ndarray:
    flat=grid.ravel()
    out=np.zeros(len(pos),dtype=np.float64)
    for idx,w in cic_components(pos,n,box):
        out += flat[idx]*w
    return out


def density_ratio(pos:np.ndarray,n:int,box:float) -> np.ndarray:
    # one particle per PM cell -> mean deposited mass per cell = 1 exactly
    rho=cic_deposit(pos,n,box)
    return rho/rho.mean()


def force_tilde(delta:np.ndarray,fg:FourierGrid) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Return g_tilde with g(a)=a^-1 g_tilde and H0=1 in EdS."""
    F=np.fft.rfftn(delta)
    coeff=1.5*F*fg.invk2
    forces=[]
    for kj in (fg.kx,fg.ky,fg.kz):
        g=np.fft.irfftn(1j*kj*coeff,s=delta.shape,axes=(0,1,2)).real
        forces.append(g)
    return tuple(forces)


def kick_factor(a1:float,a2:float) -> float:
    # integral da/(a^2 H) in EdS, H=H0 a^-3/2, H0=1
    return 2.0*(math.sqrt(a2)-math.sqrt(a1))


def drift_factor(a1:float,a2:float) -> float:
    # integral da/(a^3 H) in EdS
    return 2.0*(a1**-0.5-a2**-0.5)


def pm_step(pos:np.ndarray,p:np.ndarray,a1:float,a2:float,fg:FourierGrid) -> tuple[np.ndarray,np.ndarray]:
    # kick-drift-kick, split kick integral exactly in half
    amid=((math.sqrt(a1)+math.sqrt(a2))*0.5)**2
    rho=density_ratio(pos,fg.n,fg.box)
    forces=force_tilde(rho-1.0,fg)
    k1=kick_factor(a1,amid)
    for j in range(3):
        p[:,j]+=k1*cic_interpolate(forces[j],pos,fg.n,fg.box)
    pos=(pos+drift_factor(a1,a2)*p)%fg.box
    rho=density_ratio(pos,fg.n,fg.box)
    forces=force_tilde(rho-1.0,fg)
    k2=kick_factor(amid,a2)
    for j in range(3):
        p[:,j]+=k2*cic_interpolate(forces[j],pos,fg.n,fg.box)
    return pos,p


def velocity_divergence_theta(pos:np.ndarray,p:np.ndarray,a:float,fg:FourierGrid,R:float) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Return rho, smoothed rho, and theta=div u/(aH), u=p/a."""
    rho=density_ratio(pos,fg.n,fg.box)
    u=p/a
    mass_sm=gaussian_smooth(rho,fg,R)
    vel=[]
    for j in range(3):
        mom=cic_deposit(pos,fg.n,fg.box,u[:,j])
        mom_sm=gaussian_smooth(mom,fg,R)
        vel.append(mom_sm/np.maximum(mass_sm,1e-8))
    divF=(1j*fg.kx*np.fft.rfftn(vel[0])+
          1j*fg.ky*np.fft.rfftn(vel[1])+
          1j*fg.kz*np.fft.rfftn(vel[2]))
    div=np.fft.irfftn(divF,s=rho.shape,axes=(0,1,2)).real
    H=a**-1.5
    theta=div/(a*H)
    return rho,mass_sm,theta


def segment(delta_sm:np.ndarray,persistence_frac:float) -> tuple[np.ndarray,float]:
    rms=float(delta_sm.std())
    if persistence_frac==0.0:
        minima=local_minima(delta_sm)
    else:
        minima=h_minima(delta_sm,persistence_frac*rms)
    coords=np.argwhere(minima)
    markers=np.zeros(delta_sm.shape,dtype=np.int32)
    for i,c in enumerate(coords,1):
        markers[tuple(c)]=i
    labels=watershed(delta_sm,markers=markers,connectivity=1)
    return labels,rms


def interface_area_density(active_cells:np.ndarray,dx:float,box:float) -> float:
    faces=0
    for axis in range(3):
        faces += int(np.count_nonzero(active_cells != np.roll(active_cells,-1,axis=axis)))
    area=faces*dx*dx
    return float(area/(box**3))


def basin_stats(labels:np.ndarray,rho:np.ndarray,theta:np.ndarray,dx:float):
    lab=labels.ravel()
    nlab=int(lab.max())
    counts=np.bincount(lab,minlength=nlab+1)[1:].astype(float)
    mean_rho=np.bincount(lab,weights=rho.ravel(),minlength=nlab+1)[1:]/counts
    mean_theta=np.bincount(lab,weights=theta.ravel(),minlength=nlab+1)[1:]/counts
    hratio=1.0+mean_theta/3.0
    gamma=np.full(nlab,-np.inf,dtype=float)
    valid=(mean_rho>0)&(hratio>0)
    gamma[valid]=hratio[valid]**2/mean_rho[valid]-1.0
    reff=(3.0*counts*dx**3/(4.0*math.pi))**(1.0/3.0)
    gamma_grow=1.0+2.0*C_SIGMA*XI0/(BETA_GAMMA*reff)
    active=valid & (gamma>=gamma_grow)
    return {
        'nlab':nlab,'counts':counts,'mean_rho':mean_rho,'mean_theta':mean_theta,
        'hratio':hratio,'gamma':gamma,'reff':reff,'gamma_grow':gamma_grow,'active':active
    }


def particle_basin_labels(pos:np.ndarray,labels:np.ndarray,box:float) -> np.ndarray:
    n=labels.shape[0]; dx=box/n
    idx=np.floor(pos/dx).astype(np.int64)%n
    return labels[idx[:,0],idx[:,1],idx[:,2]].astype(np.int32)


def inherited_fraction(prev_particle_labels:np.ndarray,prev_flags:np.ndarray,
                       curr_particle_labels:np.ndarray,ncurr:int) -> np.ndarray:
    parent_active=prev_flags[prev_particle_labels-1]
    totals=np.bincount(curr_particle_labels,minlength=ncurr+1)[1:].astype(float)
    inherited=np.bincount(curr_particle_labels,weights=parent_active.astype(float),minlength=ncurr+1)[1:]
    frac=np.divide(inherited,totals,out=np.zeros_like(inherited),where=totals>0)
    return frac


def merger_metrics(prev_particle_labels:np.ndarray,curr_particle_labels:np.ndarray,nprev:int,ncurr:int):
    # Particle-ID overlap matrix, sparse via packed pair bincount.
    packed=prev_particle_labels.astype(np.int64)*(ncurr+1)+curr_particle_labels.astype(np.int64)
    bc=np.bincount(packed)
    totals=np.bincount(curr_particle_labels,minlength=ncurr+1).astype(float)
    parent_counts=np.zeros(ncurr,dtype=int)
    domfrac=np.zeros(ncurr,dtype=float)
    # iterate only nonzero pair entries
    nz=np.nonzero(bc)[0]
    for code in nz:
        parent=code//(ncurr+1); child=code%(ncurr+1)
        if parent==0 or child==0: continue
        c=bc[code]; tot=totals[child]
        if tot<=0: continue
        frac=c/tot
        if frac>=0.10: parent_counts[child-1]+=1
        if frac>domfrac[child-1]: domfrac[child-1]=frac
    return {
        'fraction_children_mult_parent10':float(np.mean(parent_counts>=2)) if ncurr else 0.0,
        'mean_parent_count10':float(np.mean(parent_counts)) if ncurr else 0.0,
        'median_dominant_parent_fraction':float(np.median(domfrac)) if ncurr else math.nan,
    }


def summarize_snapshot(seed:int,pers:float,a:float,labels:np.ndarray,rho:np.ndarray,theta:np.ndarray,
                       stats:dict,inst_flags:np.ndarray,h50:np.ndarray,h10:np.ndarray,
                       merger:dict|None,fg:FourierGrid,sigma8_now:float,delta_sm_rms:float):
    counts=stats['counts']; ncell=labels.size
    phi_inst=float(counts[inst_flags].sum()/ncell)
    phi50=float(counts[h50].sum()/ncell)
    phi10=float(counts[h10].sum()/ncell)
    active_cells_inst=inst_flags[labels-1]
    active_cells50=h50[labels-1]
    active_cells10=h10[labels-1]
    s_inst=interface_area_density(active_cells_inst,fg.box/fg.n,fg.box)
    s50=interface_area_density(active_cells50,fg.box/fg.n,fg.box)
    s10=interface_area_density(active_cells10,fg.box/fg.n,fg.box)
    def L(phi,s): return 3.0*phi/s if s>0 else math.inf
    act=stats['active']
    row={
        'seed':seed,'persistence_rms':pers,'a':a,'Ngrid':fg.n,
        'n_basins':stats['nlab'],'delta_sm_rms':delta_sm_rms,'sigma8_nonlinear':sigma8_now,
        'instantaneous_active_volume_fraction':phi_inst,
        'lineage50_active_volume_fraction':phi50,
        'lineage10_active_volume_fraction':phi10,
        'instantaneous_interface_area_density_hMpc_inv':s_inst,
        'lineage50_interface_area_density_hMpc_inv':s50,
        'lineage10_interface_area_density_hMpc_inv':s10,
        'instantaneous_hydraulic_scale_hMpc':L(phi_inst,s_inst),
        'lineage50_hydraulic_scale_hMpc':L(phi50,s50),
        'lineage10_hydraulic_scale_hMpc':L(phi10,s10),
        'median_R_all_hMpc':float(np.median(stats['reff'])),
        'median_R_active_hMpc':float(np.median(stats['reff'][act])) if np.any(act) else math.nan,
        'median_gamma_active':float(np.median(stats['gamma'][act])) if np.any(act) else math.nan,
        'mean_density_ratio_all':float(np.average(stats['mean_rho'],weights=counts)),
        'median_density_ratio_basin':float(np.median(stats['mean_rho'])),
        'median_Hratio_basin':float(np.median(stats['hratio'])),
    }
    if merger:
        row.update(merger)
    else:
        row.update({'fraction_children_mult_parent10':math.nan,'mean_parent_count10':math.nan,
                    'median_dominant_parent_fraction':math.nan})
    return row


def run_seed(seed:int,n:int,steps:int,persistence_values:tuple[float,...]):
    fg=fourier_grid(n,BOX); dx=BOX/n
    delta0=gaussian_linear_field(seed,fg)
    sigma0=sigma8_grid(delta0,fg)
    s=displacement_from_delta(delta0,fg)
    q=grid_positions(n,BOX)
    pos=(q+A_INIT*s)%BOX
    # Zel'dovich growing-mode canonical momentum p=a^3 H s=a^(3/2)s for H0=1 EdS.
    p=(A_INIT**1.5)*s

    # initial PM density validation against linear target a_i delta0.
    rho_i=density_ratio(pos,n,BOX)
    sig8_i=sigma8_grid(rho_i-1.0,fg)
    target_i=A_INIT*sigma0

    # fixed linear-a step schedule; force requested snapshots to be exact nodes.
    base=np.linspace(A_INIT,1.0,steps+1)
    nodes=np.unique(np.round(np.concatenate([base,np.array(SNAPSHOTS)]),10))
    nodes.sort()

    state={pers:{'prev_particle_labels':None,'h50':None,'h10':None,'nprev':0} for pers in persistence_values}
    rows=[]
    snapset={round(x,10) for x in SNAPSHOTS}
    a=A_INIT
    for a2 in nodes[1:]:
        pos,p=pm_step(pos,p,a,float(a2),fg)
        a=float(a2)
        if round(a,10) not in snapset:
            continue
        rho,rho_sm,theta=velocity_divergence_theta(pos,p,a,fg,R_SMOOTH)
        delta_sm=rho_sm-1.0
        sig_now=sigma8_grid(rho-1.0,fg)
        for pers in persistence_values:
            labels,rms=segment(delta_sm,pers)
            st=basin_stats(labels,rho,theta,dx)
            curr_particles=particle_basin_labels(pos,labels,BOX)
            info=state[pers]
            current=st['active'].copy()
            merger=None
            if info['prev_particle_labels'] is None:
                h50=current.copy(); h10=current.copy()
            else:
                f50=inherited_fraction(info['prev_particle_labels'],info['h50'],curr_particles,st['nlab'])
                f10=inherited_fraction(info['prev_particle_labels'],info['h10'],curr_particles,st['nlab'])
                h50=current | (f50>=0.50)
                h10=current | (f10>=0.10)
                merger=merger_metrics(info['prev_particle_labels'],curr_particles,info['nprev'],st['nlab'])
            rows.append(summarize_snapshot(seed,pers,a,labels,rho,theta,st,current,h50,h10,merger,fg,sig_now,rms))
            info['prev_particle_labels']=curr_particles
            info['h50']=h50; info['h10']=h10; info['nprev']=st['nlab']
    validation={
        'seed':seed,'Ngrid':n,'sigma8_linear_target_a1':sigma0,
        'initial_PM_sigma8':sig8_i,'initial_linear_sigma8_target':target_i,
        'initial_sigma8_ratio_PM_to_linear':sig8_i/target_i,
        'particle_count':int(len(pos)),'mass_closure_final':float(density_ratio(pos,n,BOX).mean()),
    }
    return rows,validation


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results')
    ap.add_argument('--n',type=int,default=48)
    ap.add_argument('--steps',type=int,default=38)
    ap.add_argument('--seeds',type=int,nargs='+',default=[11,23])
    ap.add_argument('--persistence',type=float,nargs='+',default=list(PERSISTENCE))
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    allrows=[]; validations=[]
    for seed in args.seeds:
        print(f'=== seed {seed}, N={args.n} ===',flush=True)
        rr,v=run_seed(seed,args.n,args.steps,tuple(args.persistence))
        allrows.extend(rr); validations.append(v)
        print(json.dumps(v,indent=2),flush=True)
    table=args.out/'nonlinear_pm_basin_lineage_table.csv'
    with table.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0].keys()))
        w.writeheader(); w.writerows(allrows)
    (args.out/'nonlinear_pm_validation.json').write_text(json.dumps(validations,indent=2),encoding='utf-8')

    final=[r for r in allrows if abs(r['a']-1.0)<1e-10]
    def rng(key):
        x=np.array([r[key] for r in final],dtype=float)
        return float(np.nanmin(x)),float(np.nanmax(x))
    lines=[
        'Exploratory nonlinear PM + persistent basin-lineage stress test',
        f'box_hinvMpc={BOX}',f'Ngrid={args.n}',f'particle_count={args.n**3}',
        f'R_smooth_hinvMpc={R_SMOOTH}',f'a_init={A_INIT}',f'snapshots={SNAPSHOTS}',
        f'seeds={tuple(args.seeds)}',f'persistence_rms={tuple(args.persistence)}',
        f'conditional_xi0_hinvMpc={XI0}',f'conditional_beta_gamma={BETA_GAMMA}',
        'dynamics=Einstein-de Sitter particle-mesh; exploratory, not precision LambdaCDM',
        'lineage=particle-ID overlap between consecutive Eulerian watershed catalogues',
    ]
    for key in ['instantaneous_active_volume_fraction','lineage50_active_volume_fraction','lineage10_active_volume_fraction',
                'instantaneous_interface_area_density_hMpc_inv','lineage50_interface_area_density_hMpc_inv',
                'fraction_children_mult_parent10','median_dominant_parent_fraction']:
        lo,hi=rng(key); lines.append(f'D1_{key}_range={lo:.8g},{hi:.8g}')
    lines.append('interpretation=This pilot tests whether nonlinear gravity plus an evolving basin catalogue and explicit lineage can evade the v1.2 Gaussian/Zeldovich low-filling obstruction. History thresholds are sensitivity brackets, not calibrated phase dynamics.')
    summary='\n'.join(lines)+'\n'
    (args.out/'nonlinear_pm_basin_lineage_summary.txt').write_text(summary,encoding='utf-8')
    print(summary)

if __name__=='__main__':
    main()
