#!/usr/bin/env python3
"""Hamiltonian-to-Landau bias matching audit.

The manuscript uses a dimensionless regional Hamiltonian ratio Gamma.  On the
expanding/V side, Gamma=1 is the proposed matter--geometry equality separatrix.
Instead of introducing an arbitrary function H(Gamma), this audit matches the
bulk free-energy gain of the normalized sextic daughter branch to the
Hamiltonian imbalance,

    Delta w_V(H) = beta * (Gamma - 1),      Gamma >= 1,

where beta is a single dimensionless normalization/coupling.  The phase-field
sharp-interface law then implies the exact finite-size condition

    R_c/xi0 = 2 C_sigma / [beta (Gamma-1)].

Thus Gamma=1 is the infinite-domain coexistence limit, while a finite basin
requires a capillary overshoot above equality.  beta=1 is a matched
normalization benchmark, not a microscopic derivation.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path
from scipy.optimize import root_scalar
from analysis_phase import phase_constants, central_and_selected
from analysis_phase_interface import interface_constants


def void_bulk_gain(H: float) -> float:
    """W_D-W_V for H<0, positive when the V daughter is thermodynamically favored."""
    if H > 0:
        raise ValueError("void_bulk_gain expects H <= 0")
    central, selected = central_and_selected(H)
    return float(central[1] - selected[1])


def H_from_gamma(Gamma: float, beta: float = 1.0) -> float:
    """Implicit matched-bias map on the metastable D--V interval."""
    if beta <= 0:
        raise ValueError("beta must be positive")
    if Gamma < 1:
        raise ValueError("V-side matching requires Gamma >= 1")
    Hcoex, _, _, _, Hspin = phase_constants()
    target = beta*(Gamma-1.0)
    if target == 0:
        return -Hcoex
    max_gain = void_bulk_gain(-Hspin+1e-8)
    if target >= max_gain:
        raise ValueError("target exceeds parent-metastability range")
    f=lambda H: void_bulk_gain(H)-target
    return root_scalar(f, bracket=[-Hspin+1e-8, -Hcoex-1e-12], xtol=1e-13).root


def gamma_for_radius(R_over_xi: float, beta: float = 1.0) -> float:
    """Finite-size growth threshold from capillary balance."""
    if R_over_xi <= 0 or beta <= 0:
        raise ValueError("R_over_xi and beta must be positive")
    Csig=interface_constants()['C_sigma']
    return 1.0 + 2.0*Csig/(beta*R_over_xi)


def radius_over_xi_from_gamma(Gamma: float, beta: float = 1.0) -> float:
    if Gamma <= 1 or beta <= 0:
        return math.inf
    Csig=interface_constants()['C_sigma']
    return 2.0*Csig/(beta*(Gamma-1.0))


def near_coexistence_bias_slope(beta: float = 1.0) -> float:
    """d[-H-Hcoex]/dGamma at coexistence."""
    _,_,_,jump,_=phase_constants()
    return beta/jump


def benchmark() -> dict[str,float]:
    Hcoex,_,_,jump,Hspin=phase_constants()
    Csig=interface_constants()['C_sigma']
    d1090=5.298450046
    C1090=interface_constants()['C_10_90']
    xi=d1090/C1090
    R=15.6
    Rxi=R/xi
    beta=1.0
    Gamma=gamma_for_radius(Rxi,beta)
    H=H_from_gamma(Gamma,beta)
    return {
        'beta_benchmark':beta,
        'H_coexistence_abs':Hcoex,
        'H_parent_spinodal_abs':Hspin,
        'Delta_chi':jump,
        'C_sigma':Csig,
        'near_coexistence_dB_dGamma':near_coexistence_bias_slope(beta),
        'conditional_xi_hMpc':xi,
        'hydraulic_basin_scale_hMpc':R,
        'R_over_xi':Rxi,
        'Gamma_growth_threshold':Gamma,
        'Gamma_overshoot':Gamma-1.0,
        'matched_H_V':H,
        'matched_bias_excess_abs':(-H)-Hcoex,
        'matched_Delta_w':void_bulk_gain(H),
        'Rc_over_xi_backcheck':radius_over_xi_from_gamma(Gamma,beta),
    }


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=root/'results')
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    b=benchmark()
    text='Hamiltonian-to-Landau bias matching audit\n' + '\n'.join(f'{k} = {v:.12g}' for k,v in b.items()) + '\n'
    text += ('interpretation = Delta w = beta (Gamma-1) removes an arbitrary bias function but leaves one normalization/coupling beta; beta=1 is a benchmark, not a microscopic derivation\n')
    (args.out/'bias_matching_summary.txt').write_text(text,encoding='utf-8')
    print(text,end='')

if __name__=='__main__': main()
