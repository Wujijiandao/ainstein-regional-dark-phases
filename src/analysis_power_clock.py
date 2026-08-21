#!/usr/bin/env python3
"""Power-spectrum closure audit for the regional phase first-passage model.

This module does *not* claim a precision Boltzmann calculation.  Its purpose is
structural: connect the variance clock tau=D(z)^2 S(R) to a Planck-normalized
CDM-like linear spectrum and test whether the spherical void-shell-crossing
barrier can make the V phase volume filling on physically large void scales.

Two standard analytic transfer-function shapes are used as a robustness bracket:
  * BBKS with the Sugiyama baryon-corrected shape parameter;
  * the Eisenstein-Hu zero-baryon fitting form.
Both are normalized to sigma_8, so the key conclusion is insensitive to the
absolute primordial amplitude.  A sharper model-independent audit also uses
only sigma_8 and the fact that sigma(R) decreases for R >= 8 h^-1 Mpc in the
standard matter spectrum.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, quad
from scipy.optimize import brentq

from analysis_excursion_phase import eulerian_volume_fractions, first_crossing_fractions, DELTA_NL_C
from analysis_regional_state import scalar_void

# Planck 2018 base-LambdaCDM late-time parameters (rounded to the precision
# relevant for this structural audit).
OMEGA_M = 0.315
OMEGA_B = 0.049
H = 0.674
N_S = 0.965
SIGMA8 = 0.811

# From the v0.4 two-barrier prototype.
TAU_HALF_V = 2.940769839691134
VOID_EXPANSION = 1.697  # spherical shell-crossing R_E/R_L (Jennings et al.)


def window_tophat(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    xs = x[small]
    out[small] = 1.0 - xs * xs / 10.0 + xs**4 / 280.0
    xx = x[~small]
    out[~small] = 3.0 * (np.sin(xx) - xx * np.cos(xx)) / xx**3
    return out


def transfer_bbks(k_hmpc: np.ndarray) -> np.ndarray:
    """BBKS transfer with Sugiyama baryon correction; k in h/Mpc."""
    gamma = OMEGA_M * H * math.exp(-OMEGA_B * (1.0 + math.sqrt(2.0 * H) / OMEGA_M))
    q = np.asarray(k_hmpc, dtype=float) / gamma
    x = 2.34 * q
    first = np.log1p(x) / x
    poly = 1.0 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4
    return first * poly**(-0.25)


def transfer_eh0(k_hmpc: np.ndarray) -> np.ndarray:
    """Eisenstein-Hu zero-baryon transfer-function fit; k in h/Mpc."""
    q = np.asarray(k_hmpc, dtype=float) / (OMEGA_M * H)
    L0 = np.log(2.0 * math.e + 1.8 * q)
    C0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
    return L0 / (L0 + C0 * q*q)


def sigma_ratio(R_hmpc: float, transfer) -> float:
    """Return sigma(R)/sigma8 after normalizing the same P(k) shape at 8 Mpc/h."""
    k = np.logspace(-5, 3, 60000)  # h/Mpc
    T = transfer(k)
    base = k**(N_S + 3.0) * T*T
    lnk = np.log(k)
    I8 = simpson(base * window_tophat(k*8.0)**2, x=lnk)
    IR = simpson(base * window_tophat(k*R_hmpc)**2, x=lnk)
    return math.sqrt(IR / I8)


def sigma_R(R_hmpc: float, transfer) -> float:
    return SIGMA8 * sigma_ratio(R_hmpc, transfer)


def E_lcdm(a: float) -> float:
    return math.sqrt(OMEGA_M / a**3 + (1.0 - OMEGA_M))


def growth_unnormalized(a: float) -> float:
    # Exact growing-mode integral for flat matter+Lambda cosmology.
    val = quad(lambda ap: 1.0 / (ap**3 * E_lcdm(ap)**3), 1e-7, a,
               epsabs=1e-10, epsrel=1e-9, limit=300)[0]
    return 2.5 * OMEGA_M * E_lcdm(a) * val

_G1 = growth_unnormalized(1.0)


def growth_D(z: float) -> float:
    a = 1.0 / (1.0 + float(z))
    return growth_unnormalized(a) / _G1


def tau_R_z(R_hmpc: float, z: float, transfer) -> float:
    s = sigma_R(R_hmpc, transfer)
    return growth_D(z)**2 * s*s


def z_when_tau(R_hmpc: float, tau_target: float, transfer):
    s2 = sigma_R(R_hmpc, transfer)**2
    if s2 < tau_target:
        return None
    Dtarget = math.sqrt(tau_target / s2)
    if Dtarget >= 1.0:
        return 0.0
    return brentq(lambda z: growth_D(z) - Dtarget, 0.0, 30.0)


def R_when_S(tau_target: float, transfer) -> float:
    return brentq(lambda R: sigma_R(R, transfer)**2 - tau_target, 0.1, 20.0)



def spherical_void_state_from_linear(delta_v: float):
    """Return (rho/rhobar, Upsilon) for a negative EdS spherical linear contrast."""
    if not delta_v < 0:
        raise ValueError("delta_v must be negative")
    eta = brentq(lambda e: scalar_void(e)[2] - delta_v, 1e-7, 10.0)
    density_ratio, upsilon, _ = scalar_void(eta)
    return density_ratio, upsilon


def eulerian_void_fraction_selfconsistent(tau: float, delta_v: float, delta_c: float = 1.686) -> float:
    """Eulerian V fraction using the spherical density ratio implied by delta_v."""
    density_v, _ = spherical_void_state_from_linear(delta_v)
    fc, fv, fd = first_crossing_fractions(tau, delta_c=delta_c, delta_v=delta_v)
    weights = np.array([fc / DELTA_NL_C, fv / density_v, fd], dtype=float)
    return float(weights[1] / np.sum(weights))


def required_void_barrier_for_half_volume(tau: float):
    """Inverse target: spherical delta_v required for phi_V=1/2 at a given variance clock."""
    root = brentq(lambda dv: eulerian_void_fraction_selfconsistent(tau, dv) - 0.5, -2.717, -0.05)
    density_v, upsilon = spherical_void_state_from_linear(root)
    return root, density_v, upsilon

def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=root / 'results')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    models = [('BBKS-Sugiyama', transfer_bbks), ('EH98-zero-baryon', transfer_eh0)]
    radii = [1, 2, 3, 5, 8, 10, 15, 20]
    rows = []
    for name, T in models:
        for R in radii:
            sig = sigma_R(R, T)
            tau0 = sig*sig
            phiV = eulerian_volume_fractions(tau0)[1]
            rows.append((name, R, sig, tau0, phiV))

    with (args.out / 'power_clock_table.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['transfer_shape', 'R_L_hinv_Mpc', 'sigma_R_z0', 'S_R_z0',
                    'phiV_spherical_prototype_z0'])
        for row in rows:
            w.writerow([row[0], f'{row[1]:.6g}', f'{row[2]:.12g}', f'{row[3]:.12g}', f'{row[4]:.12g}'])

    # Observed median Eulerian void scale -> spherical shell-crossing Lagrangian scale.
    med_e = (15.0, 19.0)
    med_l = tuple(x / VOID_EXPANSION for x in med_e)

    with (args.out / 'power_clock_summary.txt').open('w', encoding='utf-8') as f:
        f.write('Power-spectrum closure audit for regional V-phase prototype\n')
        f.write(f'Planck-like parameters: Om={OMEGA_M}, Ob={OMEGA_B}, h={H}, ns={N_S}, sigma8={SIGMA8}\n')
        f.write(f'tau_at_50pct_Eulerian_V = {TAU_HALF_V:.12f}\n')
        f.write(f'required_sigma_for_half_V = {math.sqrt(TAU_HALF_V):.12f}\n')
        f.write(f'sigma8^2 = {SIGMA8**2:.12f}\n')
        f.write(f'phiV_at_tau_sigma8sq = {eulerian_volume_fractions(SIGMA8**2)[1]:.12f}\n')
        f.write(f'observed_median_void_RE_range_hinvMpc = {med_e[0]:.3f},{med_e[1]:.3f}\n')
        f.write(f'spherical_RL_range_hinvMpc = {med_l[0]:.6f},{med_l[1]:.6f}\n')
        f.write('\nTransfer-shape sensitivity:\n')
        for name, T in models:
            rhalf = R_when_S(TAU_HALF_V, T)
            f.write(f'{name}: R_L giving S=tau_half today = {rhalf:.9f} h^-1 Mpc\n')
            for R in med_l:
                sig = sigma_R(R, T); tau0 = sig*sig; phi = eulerian_volume_fractions(tau0)[1]
                f.write(f'  R_L={R:.6f}: sigma={sig:.9f}, S={tau0:.9f}, phiV={phi:.12g}, nu_v={2.717/sig:.6f}\n')
            for R in (1.0, 2.0, 3.0):
                zhalf = z_when_tau(R, TAU_HALF_V, T)
                f.write(f'  R_L={R:.1f}: z(phiV=0.5)=' + ('never by z=0\n' if zhalf is None else f'{zhalf:.9f}\n'))
        f.write('\nInverse selector target (self-consistent spherical density ratio):\n')
        for name, T in models:
            f.write(f'{name}:\n')
            for R in med_l:
                tau0 = sigma_R(R, T)**2
                dv_req, rho_req, ups_req = required_void_barrier_for_half_volume(tau0)
                f.write(f'  R_L={R:.6f}: required_deltaL_V={dv_req:.9f}, rho/rhobar={rho_req:.9f}, Upsilon={ups_req:.9f}\n')
        f.write('\nInterpretation:\n')
        f.write('Within the spherical shell-crossing + fixed-scale Markov prototype, mature void scales are far too low-variance to make V volume-dominant by z=0.\n')
        f.write('Therefore the V selector must occur earlier than spherical shell crossing, be genuinely multiscale/topological, or use a parent spectrum/domain rule different from this benchmark.\n')
        f.write('This is an obstruction test, not a precision cosmology likelihood.\n')

    print((args.out / 'power_clock_summary.txt').read_text())


if __name__ == '__main__':
    main()
