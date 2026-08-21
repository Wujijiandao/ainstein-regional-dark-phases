#!/usr/bin/env python3
"""Planar interface audit for the normalized sextic regional phase potential.

This script does not fit an interface tension. It evaluates the dimensionless
Cahn-Hilliard/Ginzburg wall constants at the D-daughter coexistence point.
For a physical phase-field free-energy density
    F = kappa/2 |grad chi|^2 + E0 [W(chi;H)-Wmin],
the planar coexistence wall obeys
    sigma = C_sigma sqrt(kappa E0)
and its 10--90 order-parameter thickness is
    d_10_90 = C_10_90 sqrt(kappa/E0).
"""
from __future__ import annotations
import math
from pathlib import Path
import argparse
from scipy.integrate import quad
from analysis_phase import W0, phase_constants, central_and_selected


def W(x: float, H: float) -> float:
    return W0(x) - H*x


def interface_constants() -> dict[str,float]:
    H, _, _, _, _ = phase_constants()
    (x1,e1),(x2,e2)=central_and_selected(H)
    emin=0.5*(e1+e2)
    if abs(e1-e2)>1e-10:
        raise RuntimeError('coexistence minima are not degenerate')
    I=quad(lambda x: math.sqrt(max(W(x,H)-emin,0.0)),x1,x2,epsabs=1e-13,epsrel=1e-13,limit=300)[0]
    C_sigma=math.sqrt(2.0)*I
    dx=x2-x1
    xa=x1+0.1*dx; xb=x1+0.9*dx
    C_10_90=quad(lambda x: 1.0/math.sqrt(2.0*max(W(x,H)-emin,1e-300)),xa,xb,epsabs=1e-12,epsrel=1e-12,limit=300)[0]
    return {
        'H_coex':H,'chi_D':x1,'chi_daughter':x2,'W_min':emin,
        'C_sigma':C_sigma,'C_10_90':C_10_90,
        'thickness_to_ell_if_E0_eq_epsV':C_10_90/C_sigma,
    }


def required_phase_width(ell_hMpc: float, energy_scale_ratio: float=1.0) -> float:
    """10--90 wall width corresponding to ell=sigma/epsilon_V.

    energy_scale_ratio = E0 / epsilon_V. Since
      ell = C_sigma (E0/epsilon_V) xi0,
      d10-90 = C_10_90 xi0,
    d = (C_10_90/C_sigma) ell / energy_scale_ratio.
    """
    c=interface_constants()
    return c['thickness_to_ell_if_E0_eq_epsV']*ell_hMpc/energy_scale_ratio


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=root/'results')
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    c=interface_constants()
    ell_bench=1.241834074954
    width=required_phase_width(ell_bench,1.0)
    text=(
        'Sextic phase-field planar-interface audit\n'
        f'H_coex = {c["H_coex"]:.12f}\n'
        f'chi_D = {c["chi_D"]:.12f}\n'
        f'chi_daughter = {c["chi_daughter"]:.12f}\n'
        f'C_sigma = {c["C_sigma"]:.12f}\n'
        f'C_10_90 = {c["C_10_90"]:.12f}\n'
        f'd10_90/ell_if_E0_eq_epsilonV = {c["thickness_to_ell_if_E0_eq_epsV"]:.12f}\n'
        f'interface_crossing_ell_upper_benchmark = {ell_bench:.12f} h^-1 Mpc\n'
        f'corresponding_d10_90_if_E0_eq_epsilonV = {width:.12f} h^-1 Mpc\n'
        'interpretation = conditional scale audit only; kappa and E0/epsilonV are not derived\n'
    )
    (args.out/'phase_interface_summary.txt').write_text(text,encoding='utf-8')
    print(text,end='')

if __name__=='__main__':
    main()
