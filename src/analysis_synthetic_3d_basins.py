#!/usr/bin/env python3
"""Deterministic 3D synthetic-basin stress test for the regional phase selector.

This is deliberately *not* an N-body calculation.  A periodic Gaussian linear
field with an Eisenstein-Hu-like transfer shape is normalized to sigma8=0.811,
coarse-grained on a few-Mpc scale, segmented into non-spherical watershed
basins, and then tested against the manuscript's conditional finite-size V-phase
criterion.

The purpose is adversarial: quantify whether an instantaneous, single-scale
Gaussian basin catalogue would already make a large V-phase filling fraction.
If it does not, a nonlinear/persistent basin hierarchy and history are genuinely
needed rather than being optional decoration.
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from skimage.morphology import h_minima, local_minima
from skimage.segmentation import watershed

from analysis_power_clock import transfer_eh0, SIGMA8
from analysis_regional_state import scalar_void
from analysis_phase_interface import interface_constants

ROOT = Path(__file__).resolve().parents[1]
BOX = 256.0            # h^-1 Mpc
NGRID = 80
R_SMOOTH = 4.0         # h^-1 Mpc, coarse-graining used only in this stress test
XI0 = 3.12938          # h^-1 Mpc, conditional interface benchmark from manuscript
BETA_GAMMA = 1.0       # matched normalization benchmark, not a fitted cosmological value
SEEDS = (11, 23, 47, 91)
PERSISTENCE = (0.0, 0.05, 0.10, 0.20)  # multiples of rms of smoothed field


def top_hat(x: np.ndarray) -> np.ndarray:
    out = np.ones_like(x)
    m = np.abs(x) > 1e-8
    xm = x[m]
    out[m] = 3.0 * (np.sin(xm) - xm * np.cos(xm)) / xm**3
    return out


def transfer_array(k: np.ndarray) -> np.ndarray:
    out = np.zeros_like(k)
    m = k > 0
    # transfer_eh0 is scalar; the grid is modest and this stress test is run only a few times.
    out[m] = np.array([transfer_eh0(float(x)) for x in k[m]])
    return out


def gaussian_field(seed: int):
    rng = np.random.default_rng(seed)
    dx = BOX / NGRID
    white = rng.normal(size=(NGRID, NGRID, NGRID))
    W = np.fft.fftn(white)
    kval = 2.0 * math.pi * np.fft.fftfreq(NGRID, d=dx)
    kx, ky, kz = np.meshgrid(kval, kval, kval, indexing='ij')
    k = np.sqrt(kx*kx + ky*ky + kz*kz)
    T = transfer_array(k)
    P = np.zeros_like(k)
    m = k > 0
    P[m] = k[m]**0.9665 * T[m]**2
    field = np.fft.ifftn(W * np.sqrt(P)).real
    field -= field.mean()
    F = np.fft.fftn(field)
    sm8 = np.fft.ifftn(F * top_hat(k * 8.0)).real
    field *= SIGMA8 / sm8.std()
    F = np.fft.fftn(field)
    coarse = np.fft.ifftn(F * np.exp(-0.5 * (k * R_SMOOTH)**2)).real
    return field, coarse, dx


def gamma_from_mean_deltaL(deltaL: float) -> float:
    """Spherical-void map used only to turn a negative *regional mean* delta_L into Gamma."""
    if deltaL >= 0:
        return -math.inf
    # grid basins do not approach the numerical singularity in practice; bracket generously.
    eta = brentq(lambda e: scalar_void(e)[2] - deltaL, 1e-7, 12.0)
    _, ups, _ = scalar_void(eta)
    return ups*ups - 1.0


def segment_stats(coarse: np.ndarray, dx: float, persistence_frac: float):
    rms = float(coarse.std())
    if persistence_frac == 0.0:
        minima = local_minima(coarse)
    else:
        minima = h_minima(coarse, persistence_frac * rms)
    coords = np.argwhere(minima)
    markers = np.zeros(coarse.shape, dtype=np.int32)
    for i, c in enumerate(coords, 1):
        markers[tuple(c)] = i
    labels = watershed(coarse, markers=markers, connectivity=1)
    lab = labels.ravel()
    counts = np.bincount(lab)[1:]
    means = np.bincount(lab, weights=coarse.ravel())[1:] / counts
    reff = (3.0 * counts * dx**3 / (4.0 * math.pi))**(1.0/3.0)
    gamma = np.array([gamma_from_mean_deltaL(float(d)) for d in means])
    Csig = interface_constants()['C_sigma']
    gamma_grow = 1.0 + 2.0 * Csig * XI0 / (BETA_GAMMA * reff)
    equality = gamma >= 1.0
    finite = gamma >= gamma_grow
    total = counts.sum()
    def vf(mask): return float(counts[mask].sum() / total)
    return {
        'n_basins': int(len(counts)),
        'coarse_rms': rms,
        'volume_fraction_gamma_ge_1': vf(equality),
        'volume_fraction_finite_size_active': vf(finite),
        'n_gamma_ge_1': int(equality.sum()),
        'n_finite_size_active': int(finite.sum()),
        'median_R_all': float(np.median(reff)),
        'median_R_active': float(np.median(reff[finite])) if finite.any() else math.nan,
        'median_gamma_active': float(np.median(gamma[finite])) if finite.any() else math.nan,
    }


def run_all():
    rows = []
    sigma8_checks = []
    for seed in SEEDS:
        field, coarse, dx = gaussian_field(seed)
        # independent real-space check via Fourier top-hat on the already normalized field
        kval = 2.0 * math.pi * np.fft.fftfreq(NGRID, d=dx)
        kx, ky, kz = np.meshgrid(kval, kval, kval, indexing='ij')
        k = np.sqrt(kx*kx + ky*ky + kz*kz)
        sig8 = float(np.fft.ifftn(np.fft.fftn(field)*top_hat(k*8.0)).real.std())
        sigma8_checks.append(sig8)
        for p in PERSISTENCE:
            s = segment_stats(coarse, dx, p)
            rows.append((seed, p, sig8, s))
    return rows, sigma8_checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=ROOT/'results')
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    rows, checks = run_all()
    csvpath = args.out/'synthetic_3d_basin_table.csv'
    with csvpath.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['seed','persistence_rms','sigma8_check','n_basins','coarse_rms',
                    'volfrac_gamma_ge_1','volfrac_finite_size_active','n_gamma_ge_1',
                    'n_finite_size_active','median_R_all_hMpc','median_R_active_hMpc',
                    'median_gamma_active'])
        for seed,p,s8,s in rows:
            w.writerow([seed,p,f'{s8:.12g}',s['n_basins'],f"{s['coarse_rms']:.12g}",
                        f"{s['volume_fraction_gamma_ge_1']:.12g}",
                        f"{s['volume_fraction_finite_size_active']:.12g}",
                        s['n_gamma_ge_1'],s['n_finite_size_active'],
                        f"{s['median_R_all']:.12g}",f"{s['median_R_active']:.12g}",
                        f"{s['median_gamma_active']:.12g}"])
    vals0 = [s['volume_fraction_finite_size_active'] for seed,p,s8,s in rows if p == 0.0]
    vals_p = [s['volume_fraction_finite_size_active'] for seed,p,s8,s in rows]
    eq0 = [s['volume_fraction_gamma_ge_1'] for seed,p,s8,s in rows if p == 0.0]
    summary = args.out/'synthetic_3d_basin_summary.txt'
    summary.write_text(
        '3D synthetic non-spherical basin stress test\n'
        f'box_hinvMpc={BOX}\nNgrid={NGRID}\nR_smooth_hinvMpc={R_SMOOTH}\n'
        f'seeds={SEEDS}\npersistence_rms={PERSISTENCE}\n'
        f'sigma8_target={SIGMA8:.12f}\nsigma8_check_max_error={max(abs(x-SIGMA8) for x in checks):.3e}\n'
        f'conditional_xi0_hinvMpc={XI0}\nconditional_beta_gamma={BETA_GAMMA}\n'
        f'no_persistence_volfrac_gamma_ge_1_range={min(eq0):.6f},{max(eq0):.6f}\n'
        f'no_persistence_volfrac_finite_size_active_range={min(vals0):.6f},{max(vals0):.6f}\n'
        f'all_segmentations_volfrac_finite_size_active_range={min(vals_p):.6f},{max(vals_p):.6f}\n'
        'interpretation=An instantaneous single-scale Gaussian watershed produces only a small active V fraction in this deliberately conservative toy test. '
        'The result therefore does not validate the cosmological abundance mechanism; instead it shows that nonlinear basin evolution, hierarchy/history, or a different coarse-graining rule is necessary.\n'
        'status=stress-test/obstruction, not an N-body prediction\n', encoding='utf-8')
    print(summary.read_text(encoding='utf-8'))

if __name__ == '__main__':
    main()
