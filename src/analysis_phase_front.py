#!/usr/bin/env python3
"""Thin-interface / nucleation audit for the regional sextic phase field.

For Model-A (Allen-Cahn) dynamics
    tau_chi * d_t chi = kappa laplacian chi - E0 d_chi W(chi;H),
the coexistence wall has sigma = C_sigma E0 xi0 with
xi0=sqrt(kappa/E0).  In the thin-interface limit a spherical daughter
bubble in a parent phase obeys
    v_n = (kappa/(tau_chi xi0)) [Delta w/C_sigma - 2/(R/xi0)],
where Delta w = W_parent,min - W_daughter,min > 0.
Thus R_c/xi0 = 2 C_sigma/Delta w and the classical capillary barrier is
Delta F_c/(E0 xi0^3) = 16 pi C_sigma^3/(3 Delta w^2).

This script treats the relation as a deterministic sharp-interface audit.
It does not calculate a thermal/quantum nucleation rate.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path
from scipy.optimize import root_scalar
from analysis_phase import phase_constants, central_and_selected
from analysis_phase_interface import interface_constants


def well_difference(H: float) -> float:
    """Dimensionless bulk free-energy gain W_parent-W_selected (>0 past coexistence)."""
    central, selected = central_and_selected(H)
    return float(central[1] - selected[1])


def critical_radius_over_xi(H: float) -> float:
    c = interface_constants()
    dw = well_difference(H)
    if dw <= 0:
        return math.inf
    return 2.0*c['C_sigma']/dw


def barrier_over_E0_xi3(H: float) -> float:
    c = interface_constants()
    dw = well_difference(H)
    if dw <= 0:
        return math.inf
    return 16.0*math.pi*c['C_sigma']**3/(3.0*dw**2)


def dimensionless_front_speed(H: float, R_over_xi: float) -> float:
    """v_n / [kappa/(tau_chi xi0)] for a spherical daughter bubble."""
    c = interface_constants()
    return well_difference(H)/c['C_sigma'] - 2.0/R_over_xi


def H_for_critical_radius(target_R_over_xi: float) -> float:
    Hcoex, _, _, _, Hspin = phase_constants()
    f=lambda H: critical_radius_over_xi(H)-target_R_over_xi
    return root_scalar(f, bracket=[Hcoex+1e-9, Hspin-1e-7], xtol=1e-13).root


def linear_near_coex_coefficient() -> float:
    Hcoex, xD, xP, jump, _ = phase_constants()
    c=interface_constants()
    # Delta w = jump * Delta H + O(Delta H^2), so Rc/xi ~ 2 Csigma/(jump DeltaH)
    return 2.0*c['C_sigma']/jump


def benchmark() -> dict[str,float]:
    c=interface_constants(); Hcoex,_,_,jump,Hspin=phase_constants()
    # Previous conditional cross-scale benchmark: d_10-90=5.29845 h^-1 Mpc.
    d1090=5.298450046
    xi=d1090/c['C_10_90']
    Rvoid=15.6
    ratio=Rvoid/xi
    Hreq=H_for_critical_radius(ratio)
    dw=well_difference(Hreq)
    return {
        'H_coex':Hcoex,'H_spinodal':Hspin,'jump':jump,
        'C_sigma':c['C_sigma'],'C_10_90':c['C_10_90'],
        'near_coex_Rc_coeff':linear_near_coex_coefficient(),
        'conditional_xi_hMpc':xi,'void_scale_hMpc':Rvoid,'target_R_over_xi':ratio,
        'H_for_target_Rc':Hreq,'Delta_H_for_target_Rc':Hreq-Hcoex,
        'Delta_w_for_target_Rc':dw,'barrier_dimless':barrier_over_E0_xi3(Hreq),
        'speed_at_half_Rc':dimensionless_front_speed(Hreq,0.5*ratio),
        'speed_at_Rc':dimensionless_front_speed(Hreq,ratio),
        'speed_at_twice_Rc':dimensionless_front_speed(Hreq,2.0*ratio),
    }


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,default=root/'results')
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    b=benchmark()
    text='Thin-interface / nucleation audit\n' + '\n'.join(f'{k} = {v:.12g}' for k,v in b.items()) + '\n'
    text += ('interpretation = deterministic thin-interface scale audit only; '
             'tau_chi, kappa/E0, and a microscopic nucleation/noise law remain underived\n')
    (args.out/'phase_front_summary.txt').write_text(text,encoding='utf-8')
    print(text,end='')

if __name__=='__main__': main()
