#!/usr/bin/env python3
"""Regional matter-curvature equality as a parameter-free V-phase separatrix candidate.

This module tests a specific follow-up to the v0.5 shell-crossing obstruction.
For an EdS spherical separate universe, the regional state coordinate obeys

    Upsilon^2 = 1 / Omega_m,local.

An expanding underdensity has Omega_k,local = 1 - Omega_m,local.  The standard
matter-curvature equality Omega_m=Omega_k=1/2 therefore occurs at Upsilon=sqrt(2).
Using the exact EdS spherical void mapping gives a non-fitted linear landmark
(delta_L ~= -0.516745) substantially earlier than mature shell crossing.

We then insert that landmark into the same two-barrier Markov bookkeeping used
elsewhere in the repository.  This is a structural consistency audit, not a
self-consistent cosmological solution: the parent spectrum, basin definition,
phase free energy, interfaces, and backreaction on growth remain open.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from analysis_excursion_phase import DELTA_C, DELTA_NL_C, first_crossing_fractions
from analysis_power_clock import (
    growth_D,
    sigma_R,
    transfer_bbks,
    transfer_eh0,
)
from analysis_regional_state import scalar_void

UPSILON_CURV_EQ = math.sqrt(2.0)
ETA_CURV_EQ = 2.0 * math.acosh(UPSILON_CURV_EQ)
RHO_CURV_EQ, UPSILON_CHECK, DELTA_V_CURV_EQ = scalar_void(ETA_CURV_EQ)

# Observed median Eulerian void radii used in v0.5, mapped by the *mature*
# spherical shell-crossing expansion factor only to define a benchmark
# Lagrangian scale interval.  The new selector itself occurs earlier.
VOID_EXPANSION_SHELL = 1.697
MEDIAN_RE = (15.0, 19.0)  # h^-1 Mpc
MEDIAN_RL = tuple(r / VOID_EXPANSION_SHELL for r in MEDIAN_RE)




def gamma_estimator(delta_mean: float, theta_pec_over_aH: float) -> float:
    """Operational EdS single-stream estimator Gamma=Upsilon^2-1.

    Parameters
    ----------
    delta_mean:
        Regional mean density contrast, rho/rhobar - 1.
    theta_pec_over_aH:
        Regional mean comoving peculiar-velocity divergence divided by a H,
        i.e. <div_x v_p>/(a H).
    """
    if delta_mean <= -1.0:
        raise ValueError("delta_mean must exceed -1")
    hratio = 1.0 + theta_pec_over_aH / 3.0
    ups = hratio / math.sqrt(1.0 + delta_mean)
    return ups * ups - 1.0

def local_omega_m_from_upsilon(upsilon: float) -> float:
    """EdS spherical separate-universe identity Omega_m,local=Upsilon^-2."""
    return 1.0 / (float(upsilon) ** 2)


def eulerian_volume_fractions_curvature(tau: float | np.ndarray):
    """Eulerian C/V/D fractions with the curvature-equality V barrier.

    Collapse uses the usual EdS virial density ratio; the V branch uses the
    exact spherical density ratio at matter-curvature equality rather than the
    mature shell-crossing value 0.2.
    """
    fc, fv, fd = first_crossing_fractions(
        tau, delta_c=DELTA_C, delta_v=DELTA_V_CURV_EQ
    )
    scalar = np.ndim(tau) == 0
    fc = np.atleast_1d(fc)
    fv = np.atleast_1d(fv)
    fd = np.atleast_1d(fd)
    weights = np.vstack([fc / DELTA_NL_C, fv / RHO_CURV_EQ, fd])
    out = weights / np.sum(weights, axis=0)
    if scalar:
        return tuple(float(x) for x in out[:, 0])
    return out[0], out[1], out[2]


def tau_half_volume() -> float:
    return brentq(lambda t: eulerian_volume_fractions_curvature(t)[1] - 0.5,
                  1e-5, 5.0)


TAU_HALF_CURV = tau_half_volume()
SIGMA_HALF_CURV = math.sqrt(TAU_HALF_CURV)


def z_at_half_volume(R_hmpc: float, transfer):
    """Benchmark redshift where the fixed-R prototype crosses phi_V=1/2.

    Uses the Planck-like LambdaCDM growth function from analysis_power_clock.py.
    It is an external consistency benchmark, not a self-consistent prediction
    of the regional-phase theory.
    """
    s2 = sigma_R(R_hmpc, transfer) ** 2
    if s2 < TAU_HALF_CURV:
        return None
    dtarget = math.sqrt(TAU_HALF_CURV / s2)
    if dtarget >= 1.0:
        return 0.0
    return brentq(lambda z: growth_D(z) - dtarget, 0.0, 30.0)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=root / 'results')
    ap.add_argument('--fig', type=Path, default=root / 'figures')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.fig.mkdir(parents=True, exist_ok=True)

    models = [('BBKS-Sugiyama', transfer_bbks), ('EH98-zero-baryon', transfer_eh0)]
    radii = [8.0, MEDIAN_RL[0], 10.0, MEDIAN_RL[1], 12.0, 15.0]

    rows = []
    for name, T in models:
        for R in radii:
            sig = sigma_R(R, T)
            tau0 = sig * sig
            phi_c, phi_v, phi_d = eulerian_volume_fractions_curvature(tau0)
            zhalf = z_at_half_volume(R, T)
            rows.append((name, R, sig, tau0, phi_c, phi_v, phi_d, zhalf))

    with (args.out / 'curvature_selector_table.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow([
            'transfer_shape', 'R_L_hinv_Mpc', 'sigma_R_z0', 'S_R_z0',
            'phiC_z0', 'phiV_z0', 'phiD_z0', 'z_phiV_half_benchmark'
        ])
        for name, R, sig, tau0, pc, pv, pd, zh in rows:
            w.writerow([
                name, f'{R:.9g}', f'{sig:.12g}', f'{tau0:.12g}',
                f'{pc:.12g}', f'{pv:.12g}', f'{pd:.12g}',
                '' if zh is None else f'{zh:.12g}'
            ])

    with (args.out / 'curvature_selector_summary.txt').open('w', encoding='utf-8') as f:
        f.write('Regional matter-curvature equality selector audit\n')
        f.write(f'Upsilon_curvature_equality = {UPSILON_CURV_EQ:.12f}\n')
        f.write(f'eta_curvature_equality = {ETA_CURV_EQ:.12f}\n')
        f.write(f'deltaL_curvature_equality = {DELTA_V_CURV_EQ:.12f}\n')
        f.write(f'rho_over_background_curvature_equality = {RHO_CURV_EQ:.12f}\n')
        f.write(f'Omega_m_local = {local_omega_m_from_upsilon(UPSILON_CURV_EQ):.12f}\n')
        f.write(f'Omega_k_local = {1-local_omega_m_from_upsilon(UPSILON_CURV_EQ):.12f}\n')
        f.write(f'tau_at_50pct_Eulerian_V = {TAU_HALF_CURV:.12f}\n')
        f.write(f'required_sigma_for_half_V = {SIGMA_HALF_CURV:.12f}\n')
        f.write(f'benchmark_median_RL_range_hinvMpc = {MEDIAN_RL[0]:.9f},{MEDIAN_RL[1]:.9f}\n')
        f.write('\nPlanck-normalized transfer-shape benchmarks:\n')
        for name, T in models:
            f.write(f'{name}:\n')
            for R in MEDIAN_RL:
                sig = sigma_R(R, T)
                tau0 = sig * sig
                phi = eulerian_volume_fractions_curvature(tau0)[1]
                zh = z_at_half_volume(R, T)
                f.write(
                    f'  R_L={R:.6f}: sigma={sig:.9f}, S={tau0:.9f}, '
                    f'phiV_z0={phi:.9f}, z(phiV=0.5)=' +
                    ('never by z=0\n' if zh is None else f'{zh:.9f}\n')
                )
        f.write('\nInterpretation:\n')
        f.write('In an EdS spherical separate universe, Upsilon^2=1/Omega_m,local exactly.\n')
        f.write('Matter-curvature equality therefore supplies a non-fitted regional dynamical crossover at Upsilon=sqrt(2).\n')
        f.write('Inserted into the same first-passage bookkeeping, this earlier barrier lies in the variance regime needed for a volume-rich V branch on large-void scales.\n')
        f.write('Curvature domination itself is not acceleration (it is coasting/w=-1/3 in an FRW bookkeeping); the equality is tested only as a phase-selector candidate.\n')
        f.write('The calculation is not a precision cosmology fit and does not derive the vacuum free-energy density.\n')

    # Figure: fixed-scale benchmark histories and present-day scale dependence.
    fig = plt.figure(figsize=(9.4, 4.0))
    ax = fig.add_subplot(1, 2, 1)
    zgrid = np.linspace(0.0, 2.0, 220)
    for R in (8.0, MEDIAN_RL[0], 10.0, MEDIAN_RL[1]):
        ph = []
        for z in zgrid:
            tau = growth_D(z)**2 * sigma_R(R, transfer_bbks)**2
            ph.append(eulerian_volume_fractions_curvature(tau)[1])
        ax.plot(zgrid, ph, label=fr'$R_L={R:.1f}\,h^{{-1}}\,{{\rm Mpc}}$')
    ax.axhline(0.5, lw=0.8, ls=':')
    ax.set_xlim(2.0, 0.0)
    ax.set_ylim(0.0, 0.75)
    ax.set_xlabel('redshift $z$ (external growth benchmark)')
    ax.set_ylabel(r'$\phi_V$')
    ax.set_title('Curvature-equality selector: benchmark histories')
    ax.legend(frameon=False, fontsize=7)

    ax2 = fig.add_subplot(1, 2, 2)
    rgrid = np.linspace(5.0, 16.0, 140)
    for name, T in models:
        vals = [eulerian_volume_fractions_curvature(sigma_R(R, T)**2)[1] for R in rgrid]
        ax2.plot(rgrid, vals, label=name)
    ax2.axhline(0.5, lw=0.8, ls=':')
    ax2.axvspan(MEDIAN_RL[0], MEDIAN_RL[1], alpha=0.12, label='median-void $R_L$ band')
    ax2.set_ylim(0.25, 0.75)
    ax2.set_xlabel(r'$R_L\,[h^{-1}{\rm Mpc}]$')
    ax2.set_ylabel(r'$\phi_V(z=0)$')
    ax2.set_title('Present fixed-scale volume fraction')
    ax2.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(args.fig / 'fig_curvature_selector.pdf', bbox_inches='tight')
    plt.close(fig)

    print((args.out / 'curvature_selector_summary.txt').read_text())


if __name__ == '__main__':
    main()
