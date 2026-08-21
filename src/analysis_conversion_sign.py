#!/usr/bin/env python3
"""Monotone phase-conversion sign theorem and fixed-scale benchmark.

If a V phase has fixed intrinsic vacuum density epsilon_V and occupies an
Eulerian filling fraction phi_V(a), while parent dust transfers energy into it,
then after subtracting the conserved no-conversion dust baseline C a^-3 the
residual dark-energy density is

  rho_X = 3 epsilon_V a^-3 integral_0^a a'^2 phi_V(a') da',

with pressure p_X = -epsilon_V phi_V(a). Hence

  w_X = - a^3 phi_V(a) / [3 integral_0^a a'^2 phi_V(a') da'].

For nondecreasing phi_V, w_X <= -1.  This is a bookkeeping theorem under the
stated assumptions, not a claim that the physical V phase contains a ghost.
"""
from __future__ import annotations
from pathlib import Path
import csv, math
import numpy as np

from analysis_curvature_selector import eulerian_volume_fractions_curvature, MEDIAN_RL
from analysis_power_clock import growth_D, sigma_R, transfer_bbks, transfer_eh0

ROOT = Path(__file__).resolve().parents[1]


def residual_w_from_history(a: np.ndarray, phi: np.ndarray) -> float:
    a=np.asarray(a,dtype=float); phi=np.asarray(phi,dtype=float)
    if a.ndim!=1 or phi.shape!=a.shape or len(a)<3:
        raise ValueError('a and phi must be matching 1-D arrays')
    if np.any(np.diff(a)<=0) or a[0] < 0 or a[-1] <= 0:
        raise ValueError('a must be strictly increasing and nonnegative')
    I=float(np.trapezoid(a*a*phi,a))
    if I<=0:
        raise ValueError('history must have positive integrated phase volume')
    return -float(a[-1]**3*phi[-1])/(3.0*I)


def benchmark_history(R: float, transfer, a: np.ndarray, D: np.ndarray):
    # External LCDM growth is used only as in the manuscript timing audit.
    S=sigma_R(R,transfer)**2
    phi=np.array([eulerian_volume_fractions_curvature(d*d*S)[1] for d in D])
    return a,phi


def main() -> None:
    out=ROOT/'results'; out.mkdir(exist_ok=True)
    models=[('BBKS-Sugiyama',transfer_bbks),('EH98-zero-baryon',transfer_eh0)]
    radii=[MEDIAN_RL[0],10.0,MEDIAN_RL[1]]
    rows=[]
    # Growth is independent of smoothing scale/transfer choice: compute once.
    a_grid=np.geomspace(0.01,1.0,140)
    D_grid=np.array([growth_D(1/x-1) for x in a_grid])
    for name,T in models:
        for R in radii:
            a,phi=benchmark_history(R,T,a_grid,D_grid)
            w=residual_w_from_history(a,phi)
            rows.append((name,R,phi[-1],w))
    with (out/'conversion_sign_table.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.writer(f); wr.writerow(['transfer_shape','R_L_hinv_Mpc','phiV_z0','w_residual_z0_external_growth'])
        for r in rows: wr.writerow([r[0],f'{r[1]:.9g}',f'{r[2]:.12g}',f'{r[3]:.12g}'])
    # exact synthetic checks: phi=a^m => w=-(m+3)/3 for a_start->0.
    a=np.geomspace(1e-6,1.0,4000)
    w_a2=residual_w_from_history(a,a*a)
    with (out/'conversion_sign_summary.txt').open('w',encoding='utf-8') as f:
        f.write('Monotone V-phase conversion sign audit\n')
        f.write(f'synthetic_phi_equals_a2_w = {w_a2:.12f} (exact -5/3)\n')
        f.write('theorem: nondecreasing phi_V with fixed epsilon_V and negligible interface energy implies w_residual <= -1\n')
        f.write('fixed-scale external-growth benchmark:\n')
        for name,R,ph,w in rows:
            f.write(f'  {name}, R={R:.6f}: phiV0={ph:.9f}, w_res0={w:.9f}\n')
        f.write('interpretation: benchmark values are a sign/timing stress test, not a self-consistent CPL likelihood.\n')
    print((out/'conversion_sign_summary.txt').read_text())

if __name__=='__main__': main()
