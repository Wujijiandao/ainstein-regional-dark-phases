#!/usr/bin/env python3
"""Regional expansion-density state coordinate from EdS spherical dynamics.

Define
    Upsilon_D = (H_D/H_bg) * sqrt(rho_bg/rho_D).
For an overdense EdS top-hat the exact identity is Upsilon=cos(theta/2),
and for an underdense top-hat it is Upsilon=cosh(eta/2).

This version also marks *dynamical separatrices* rather than mapping Upsilon to
an arbitrary normalized Landau bias: overdense turnaround (delta_L=1.062),
standard collapse (delta_L=1.686), and mature-void shell crossing
(delta_L=-2.717, rho_void/rho_bg ~= 0.2).
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

DELTA_V_SC = -2.717


def overdense(theta: np.ndarray):
    d=(9/2)*((theta-np.sin(theta))**2)/(1-np.cos(theta))**3
    hratio=(3/2)*(theta-np.sin(theta))*np.sin(theta)/(1-np.cos(theta))**2
    ups=hratio/np.sqrt(d)
    exact=np.cos(theta/2)
    dl=(3/5)*((3/4)*(theta-np.sin(theta)))**(2/3)
    return d,hratio,ups,exact,dl


def underdense(eta: np.ndarray):
    d=(9/2)*((np.sinh(eta)-eta)**2)/(np.cosh(eta)-1)**3
    hratio=(3/2)*(np.sinh(eta)-eta)*np.sinh(eta)/(np.cosh(eta)-1)**2
    ups=hratio/np.sqrt(d)
    exact=np.cosh(eta/2)
    dl=-(3/5)*((3/4)*(np.sinh(eta)-eta))**(2/3)
    return d,hratio,ups,exact,dl


def scalar_void(eta: float):
    d=(9/2)*((math.sinh(eta)-eta)**2)/(math.cosh(eta)-1)**3
    u=math.cosh(eta/2)
    dl=-(3/5)*((3/4)*(math.sinh(eta)-eta))**(2/3)
    return d,u,dl


def curvature_matter_equality_state():
    """EdS spherical void state where local matter and curvature terms are equal."""
    ups=math.sqrt(2.0)
    eta=2.0*math.acosh(ups)
    d,u,dl=scalar_void(eta)
    return eta,d,u,dl


def void_shell_crossing_state(delta_v: float = DELTA_V_SC):
    eta = brentq(lambda e: scalar_void(e)[2] - delta_v, 0.1, 8.0)
    d,u,dl = scalar_void(eta)
    return eta,d,u,dl


def main():
    ap=argparse.ArgumentParser(); root=Path(__file__).resolve().parents[1]
    ap.add_argument('--out',type=Path,default=root/'results'); ap.add_argument('--fig',type=Path,default=root/'figures')
    args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True); args.fig.mkdir(parents=True,exist_ok=True)
    th=np.linspace(.03,2*math.pi-.03,1000); et=np.linspace(.03,4.5,800)
    D_o,H_o,U_o,Ue_o,L_o=overdense(th); D_v,H_v,U_v,Ue_v,L_v=underdense(et)
    maxerr_o=float(np.max(np.abs(U_o-Ue_o))); maxerr_v=float(np.max(np.abs(U_v-Ue_v)))

    # Standard EdS dynamical landmarks.
    theta_ta=math.pi
    delta_ta=(3/5)*((3/4)*(theta_ta-math.sin(theta_ta)))**(2/3)
    eta_sc,Dsc,Usc,DLsc=void_shell_crossing_state()
    eta_keq,Dkeq,Ukeq,DLkeq=curvature_matter_equality_state()

    with (args.out/'regional_state_summary.txt').open('w',encoding='utf-8') as f:
        f.write('Derived state coordinate: Upsilon=(H_D/H_bg)*sqrt(rho_bg/rho_D)\n')
        f.write(f'max_identity_error_overdense = {maxerr_o:.3e}\nmax_identity_error_void = {maxerr_v:.3e}\n')
        f.write('Exact EdS identities: overdense Upsilon=cos(theta/2); void Upsilon=cosh(eta/2).\n')
        f.write('Interpretation: background=1; overdense expansion 0<Upsilon<1; turnaround=0; collapse Upsilon<0; void Upsilon>1.\n\n')
        f.write('Standard dynamical landmarks (not fitted phase thresholds):\n')
        f.write(f'overdense turnaround: theta=pi, Upsilon=0, delta_L={delta_ta:.12f}\n')
        f.write('spherical collapse barrier used in excursion-set prototype: delta_c=1.686\n')
        f.write(f'void shell crossing: delta_L={DLsc:.12f}, eta={eta_sc:.12f}, rho/rho_bg={Dsc:.12f}, delta_NL={Dsc-1:.12f}, Upsilon={Usc:.12f}\n')
        f.write(f'matter-curvature equality: delta_L={DLkeq:.12f}, eta={eta_keq:.12f}, rho/rho_bg={Dkeq:.12f}, delta_NL={Dkeq-1:.12f}, Upsilon={Ukeq:.12f}\n')
        f.write('The void barrier prevents mildly underdense early regions (Upsilon just above 1) from being labelled vacuum phase.\n')

    with (args.out/'regional_state_trajectory.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['branch','parameter','rho_over_bg','Hlocal_over_Hbg','Upsilon','delta_linear'])
        for i in np.linspace(0,len(th)-1,250,dtype=int): w.writerow(['overdense',f'{th[i]:.9g}',f'{D_o[i]:.9g}',f'{H_o[i]:.9g}',f'{U_o[i]:.9g}',f'{L_o[i]:.9g}'])
        for i in np.linspace(0,len(et)-1,180,dtype=int): w.writerow(['void',f'{et[i]:.9g}',f'{D_v[i]:.9g}',f'{H_v[i]:.9g}',f'{U_v[i]:.9g}',f'{L_v[i]:.9g}'])

    fig=plt.figure(figsize=(6.3,4.5)); ax=fig.add_subplot(111)
    # Linear contrast gives a much cleaner common horizontal coordinate for the
    # two standard spherical branches and their finite dynamical landmarks.
    mo=L_o <= 1.6865
    mv=L_v >= -3.05
    ax.plot(L_o[mo],U_o[mo],label='overdense / collapse top-hat')
    ax.plot(L_v[mv],U_v[mv],label='void top-hat')
    ax.axhline(1,lw=1,ls='--'); ax.axhline(0,lw=1,ls=':')
    ax.scatter([DLsc],[Usc],s=28,zorder=5,label='void shell crossing')
    ax.scatter([DLkeq],[Ukeq],s=28,zorder=5,label='matter-curvature equality')
    ax.scatter([delta_ta],[0],s=28,zorder=5,label='turnaround')
    ax.scatter([1.686],[-1],s=28,zorder=5,label='collapse barrier')
    ax.set_xlim(-3.05,1.75); ax.set_ylim(-1.12,3.25)
    ax.set_xlabel(r'Linearly extrapolated density contrast $\delta_L$'); ax.set_ylabel(r'$\Upsilon_{\Omega}$')
    ax.set_title('Regional expansion-density state coordinate'); ax.legend(frameon=False,fontsize=8)
    fig.tight_layout(); fig.savefig(args.fig/'fig_regional_state.pdf',bbox_inches='tight'); plt.close(fig)
    print((args.out/'regional_state_summary.txt').read_text())


if __name__=='__main__': main()
