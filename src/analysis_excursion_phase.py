#!/usr/bin/env python3
"""Two-barrier first-crossing prototype for regional dark-phase fractions.

This is deliberately a *prototype*, not a precision cosmology calculation.
It uses the standard sharp-k/Markov excursion-set Brownian problem with two
absorbing dynamical barriers:

    +delta_c : mature collapse / halo formation,
    -|delta_v|: mature void shell crossing.

For a Gaussian parent field the variance clock tau = D(z)^2 S_Omega replaces
an arbitrary redshift conversion function.  The cumulative first-hit fractions
are analytic eigenfunction sums.  They are interpreted as Lagrangian phase
fractions only within this idealized model.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

DELTA_C = 1.686
DELTA_V = -2.717
DELTA_NL_C = 18.0 * math.pi**2  # EdS virial overdensity relative to background.
DELTA_NL_V = 0.2                 # Spherical void density at shell crossing.


def first_crossing_fractions(tau: float | np.ndarray,
                             delta_c: float = DELTA_C,
                             delta_v: float = DELTA_V,
                             n_terms: int = 4000):
    """Return (F_C, F_V, F_D) for variance clock tau.

    Brownian walks start at delta=0 between absorbing barriers delta_v<0<delta_c.
    F_C and F_V are cumulative first-passage probabilities to the collapse and
    void barriers. F_D is the survival (unconverted parent) probability.

    The sums are written as asymptotic hitting probabilities minus exponentially
    damped Fourier tails, which is numerically stable as tau -> 0.
    """
    if not (delta_v < 0 < delta_c):
        raise ValueError("Require delta_v < 0 < delta_c")
    arr = np.atleast_1d(np.asarray(tau, dtype=float))
    if np.any(arr < 0):
        raise ValueError("tau must be non-negative")

    A = abs(float(delta_v))
    B = float(delta_c)
    L = A + B
    alpha = A / L
    n = np.arange(1, n_terms + 1, dtype=float)
    sn = np.sin(math.pi * n * alpha)
    alt = (-1.0) ** (n + 1)
    lam = n * n * math.pi**2 / (2.0 * L * L)

    fc = np.empty_like(arr)
    fv = np.empty_like(arr)
    for i, t in enumerate(arr):
        if t == 0.0:
            fc[i] = fv[i] = 0.0
            continue
        e = np.exp(-lam * t)
        # Upper/collapse and lower/void first-hit CDFs.
        fc[i] = A / L - (2.0 / math.pi) * np.sum(alt * sn * e / n)
        fv[i] = B / L - (2.0 / math.pi) * np.sum(sn * e / n)

    # Remove roundoff-level negative values at very small tau.
    fc = np.clip(fc, 0.0, 1.0)
    fv = np.clip(fv, 0.0, 1.0)
    fd = np.clip(1.0 - fc - fv, 0.0, 1.0)
    if np.ndim(tau) == 0:
        return float(fc[0]), float(fv[0]), float(fd[0])
    return fc, fv, fd


def first_crossing_kernels(tau: float | np.ndarray,
                           delta_c: float = DELTA_C,
                           delta_v: float = DELTA_V,
                           n_terms: int = 4000):
    """Return dF_C/dtau and dF_V/dtau for the two-barrier process."""
    arr = np.atleast_1d(np.asarray(tau, dtype=float))
    if np.any(arr <= 0):
        raise ValueError("tau must be positive for first-crossing kernels")
    A = abs(float(delta_v)); B = float(delta_c); L = A + B; alpha = A / L
    n = np.arange(1, n_terms + 1, dtype=float)
    sn = np.sin(math.pi * n * alpha)
    alt = (-1.0) ** (n + 1)
    lam = n * n * math.pi**2 / (2.0 * L * L)
    kc = np.empty_like(arr); kv = np.empty_like(arr)
    for i, t in enumerate(arr):
        e = np.exp(-lam * t)
        pref = math.pi / (L * L)
        kc[i] = pref * np.sum(alt * n * sn * e)
        kv[i] = pref * np.sum(n * sn * e)
    # Tiny ringing at the numerical floor can be clipped; physically the kernels are non-negative.
    kc = np.where(kc > -1e-13, np.maximum(kc, 0.0), kc)
    kv = np.where(kv > -1e-13, np.maximum(kv, 0.0), kv)
    if np.ndim(tau) == 0:
        return float(kc[0]), float(kv[0])
    return kc, kv


def eulerian_volume_fractions(tau: float | np.ndarray,
                              delta_c_nl: float = DELTA_NL_C,
                              delta_v_nl: float = DELTA_NL_V):
    """Illustrative Eulerian volume weighting of the three Lagrangian fractions.

    A phase containing Lagrangian mass fraction F_i at density ratio Delta_i
    occupies volume proportional to F_i/Delta_i.  We use the standard EdS
    spherical values Delta_C=18*pi^2 and Delta_V=0.2, with Delta_D=1.
    """
    fc, fv, fd = first_crossing_fractions(tau)
    scalar = np.ndim(tau) == 0
    fc = np.atleast_1d(fc); fv = np.atleast_1d(fv); fd = np.atleast_1d(fd)
    w = np.vstack([fc / delta_c_nl, fv / delta_v_nl, fd])
    out = w / np.sum(w, axis=0)
    if scalar:
        return tuple(float(x) for x in out[:, 0])
    return out[0], out[1], out[2]



def apparent_w_void(tau: float | np.ndarray, dln_tau_dln_a: float = 2.0,
                      delta_c_nl: float = DELTA_NL_C,
                      delta_v_nl: float = DELTA_NL_V):
    """Apparent independently-conserved EOS of a growing vacuum filling fraction.

    If the intrinsic vacuum branch has fixed energy density epsilon_V and
    rho_V,avg = epsilon_V * phi_V, then an observer who forces this averaged
    component to obey a separate continuity equation infers

        w_app = -1 - (1/3) d ln(phi_V) / d ln(a).

    The default d ln(tau)/d ln(a)=2 is the EdS fixed-S relation tau=D^2 S
    with D proportional to a. This is an illustrative sign/shape prediction,
    not a precision w(z) calculation.
    """
    arr = np.atleast_1d(np.asarray(tau, dtype=float))
    if np.any(arr <= 0):
        raise ValueError("tau must be positive")
    fc, fv, fd = first_crossing_fractions(arr)
    kc, kv = first_crossing_kernels(arr)
    kd = -kc - kv
    wc = fc / delta_c_nl; wv = fv / delta_v_nl; wd = fd
    z = wc + wv + wd
    dz = kc / delta_c_nl + kv / delta_v_nl + kd
    phi = wv / z
    dphi = (kv / delta_v_nl * z - wv * dz) / (z * z)
    out = -1.0 - (dln_tau_dln_a * arr / 3.0) * dphi / phi
    if np.ndim(tau) == 0:
        return float(out[0])
    return out

def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=root / 'results')
    ap.add_argument('--fig', type=Path, default=root / 'figures')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.fig.mkdir(parents=True, exist_ok=True)

    # The variance clock is tau = D^2 S_Omega.  No cosmological power spectrum is specified here; tau is shown directly so
    # that the structural first-passage result is not tied to a fitted redshift.
    taus = np.geomspace(0.03, 20.0, 110)
    fc, fv, fd = first_crossing_fractions(taus)
    vc, vv, vd = eulerian_volume_fractions(taus)

    # Exact asymptotic Brownian hitting probabilities.
    A = abs(DELTA_V); B = DELTA_C; L = A + B
    asym_c = A / L; asym_v = B / L
    tau_half_v = brentq(lambda t: eulerian_volume_fractions(t)[1] - 0.5, 1.0, 6.0)

    with (args.out / 'excursion_phase_summary.txt').open('w', encoding='utf-8') as f:
        f.write('Two-barrier regional phase-fraction prototype\n')
        f.write(f'delta_c = {DELTA_C:.6f}\n')
        f.write(f'delta_v = {DELTA_V:.6f}\n')
        f.write(f'asymptotic_lagrangian_C = {asym_c:.12f}\n')
        f.write(f'asymptotic_lagrangian_V = {asym_v:.12f}\n')
        f.write(f'tau_at_50pct_Eulerian_V = {tau_half_v:.12f}\n')
        f.write(f'illustrative_wV_app_at_tau3_EdS = {apparent_w_void(3.0):.12f}\n')
        f.write('NOTE: tau=D^2 S_Omega is a variance clock, not a fitted redshift.\n')
        f.write('NOTE: Eulerian weighting uses idealized EdS spherical density ratios.\n\n')
        for t in [0.5, 1.0, 2.0, 3.0, 4.0, 8.0, 16.0]:
            c, v, d = first_crossing_fractions(t)
            cV, vV, dV = eulerian_volume_fractions(t)
            wapp = apparent_w_void(t)
            f.write(f'tau={t:4.1f} Lagrangian(C,V,D)=({c:.9f},{v:.9f},{d:.9f}) '
                    f'EulerianVolume(C,V,D)=({cV:.9f},{vV:.9f},{dV:.9f}) wV_app_EdS={wapp:.9f}\n')

    with (args.out / 'excursion_phase_trajectory.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['tau_D2S', 'F_C_lagrangian', 'F_V_lagrangian', 'F_D_parent',
                    'phi_C_eulerian', 'phi_V_eulerian', 'phi_D_eulerian'])
        for row in zip(taus, fc, fv, fd, vc, vv, vd):
            w.writerow([f'{x:.12g}' for x in row])

    fig = plt.figure(figsize=(9.4, 4.0))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(taus, fd, label=r'parent $F_D$')
    ax.plot(taus, fc, label=r'collapse first-hit $F_C$')
    ax.plot(taus, fv, label=r'void first-hit $F_V$')
    ax.set_xscale('log'); ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r'variance clock $\tau=D^2 S_\Omega$')
    ax.set_ylabel('Lagrangian first-hit fraction')
    ax.set_title('Dynamical-barrier phase fractions')
    ax.legend(frameon=False, fontsize=8)

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(taus, vd, label=r'parent $\phi_D$')
    ax2.plot(taus, vc, label=r'bound $\phi_C$')
    ax2.plot(taus, vv, label=r'void $\phi_V$')
    ax2.axhline(0.5, lw=0.8, ls=':')
    ax2.axvline(tau_half_v, lw=0.8, ls='--')
    ax2.set_xscale('log'); ax2.set_ylim(-0.02, 1.02)
    ax2.set_xlabel(r'variance clock $\tau=D^2 S_\Omega$')
    ax2.set_ylabel('illustrative Eulerian volume fraction')
    ax2.set_title('Mass-rich vs volume-rich coexistence')
    ax2.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.fig / 'fig_excursion_phase.pdf', bbox_inches='tight')
    plt.close(fig)

    print((args.out / 'excursion_phase_summary.txt').read_text())


if __name__ == '__main__':
    main()
