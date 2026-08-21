#!/usr/bin/env python3
"""History-dependent 3D Zel'dovich basin stress test.

This module extends the frozen Gaussian/watershed obstruction test by following
fixed Lagrangian underdensity basins through a Zel'dovich deformation.  It is
still not an N-body simulation.  The purpose is to test whether adding a
minimal, history-bearing gravitational deformation substantially changes the
regional Hamiltonian selector before shell crossing.

Conventions
-----------
The linearly extrapolated density field delta_L(q) is normalized to sigma8 at
D=1.  A smoothed potential Phi obeys nabla^2 Phi = delta_L, and the Zel'dovich
map is x(q,D)=q-D grad Phi.  If lambda_i are eigenvalues of Hess(Phi), the local
Jacobian is J=prod_i(1-D lambda_i).  A fixed Lagrangian basin has Eulerian
volume V_E(D)=int_Omega J d^3q and mass M proportional to its Lagrangian volume
V_L.  In an Einstein-de Sitter timing convention (f=d ln D/d ln a=1),

    H_Omega/H = 1 + (D/3 V_E) dV_E/dD,
    rho_Omega/rho_bg = V_L/V_E,
    Gamma_Omega = (H_Omega/H)^2 (V_E/V_L) - 1.

The V-phase eligibility benchmark is then Gamma >= Gamma_grow(R_eff), with the
same conditional xi0 and beta_Gamma used in the manuscript, and only while the
entire basin remains single stream (D lambda_max < 1 everywhere in the basin).

The test is intentionally adversarial and should not be interpreted as a
self-consistent cosmological abundance calculation.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from skimage.morphology import h_minima, local_minima
from skimage.segmentation import watershed

from analysis_synthetic_3d_basins import (
    BOX, NGRID, R_SMOOTH, XI0, BETA_GAMMA, SEEDS, SIGMA8,
    gaussian_field,
)
from analysis_phase_interface import interface_constants

ROOT = Path(__file__).resolve().parents[1]
PERSISTENCE = (0.0, 0.10, 0.20)
D_SNAPSHOTS = (0.25, 0.50, 0.75, 1.00)


def _fourier_geometry(coarse: np.ndarray, dx: float):
    """Return Hessian invariants and largest eigenvalue of Phi, nabla^2 Phi=delta."""
    n = coarse.shape[0]
    F = np.fft.fftn(coarse)
    kval = 2.0 * math.pi * np.fft.fftfreq(n, d=dx)
    kx, ky, kz = np.meshgrid(kval, kval, kval, indexing="ij")
    k2 = kx*kx + ky*ky + kz*kz
    invk2 = np.zeros_like(k2)
    m = k2 > 0
    invk2[m] = 1.0 / k2[m]

    # Phi_k=-delta_k/k^2; Hess_ij Phi -> k_i k_j delta_k/k^2.
    ks = (kx, ky, kz)
    H = np.empty(coarse.shape + (3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(i, 3):
            hij = np.fft.ifftn(F * (ks[i] * ks[j] * invk2)).real
            H[..., i, j] = hij
            H[..., j, i] = hij

    tr = H[..., 0, 0] + H[..., 1, 1] + H[..., 2, 2]
    trh2 = (
        H[..., 0, 0]**2 + H[..., 1, 1]**2 + H[..., 2, 2]**2
        + 2.0*(H[..., 0, 1]**2 + H[..., 0, 2]**2 + H[..., 1, 2]**2)
    )
    i2 = 0.5 * (tr*tr - trh2)
    # determinant of symmetric 3x3
    a,b,c = H[...,0,0], H[...,0,1], H[...,0,2]
    d,e = H[...,1,1], H[...,1,2]
    f = H[...,2,2]
    i3 = a*(d*f-e*e) - b*(b*f-c*e) + c*(b*e-c*d)
    lam_max = np.linalg.eigvalsh(H.reshape(-1,3,3))[:, -1].reshape(coarse.shape)

    # Internal consistency: trace(H)=delta up to numerical precision.
    trace_error = float(np.max(np.abs(tr - coarse)))
    return tr, i2, i3, lam_max, trace_error


def _segment(coarse: np.ndarray, persistence_frac: float):
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
    return labels, rms


def basin_history(seed: int, persistence_frac: float):
    field, coarse, dx = gaussian_field(seed)
    labels, coarse_rms = _segment(coarse, persistence_frac)
    tr, i2, i3, lam_max, trace_error = _fourier_geometry(coarse, dx)

    lab = labels.ravel()
    nlab = int(lab.max())
    counts = np.bincount(lab, minlength=nlab+1)[1:].astype(float)
    vl = counts * dx**3
    total_box = float(np.prod(coarse.shape) * dx**3)
    Csig = interface_constants()["C_sigma"]

    # Basin-level first-crossing tracker.
    ever_active = np.zeros(nlab, dtype=bool)
    first_D = np.full(nlab, np.nan)
    rows = []
    for D in D_SNAPSHOTS:
        J = 1.0 - D*tr + D*D*i2 - D**3*i3
        dJ = -tr + 2.0*D*i2 - 3.0*D*D*i3
        ve = np.bincount(lab, weights=J.ravel(), minlength=nlab+1)[1:] * dx**3
        dve = np.bincount(lab, weights=dJ.ravel(), minlength=nlab+1)[1:] * dx**3
        # First shell crossing is controlled by the largest deformation eigenvalue.
        local_single = (D * lam_max.ravel()) < 1.0
        single_counts = np.bincount(lab, weights=local_single.astype(float), minlength=nlab+1)[1:]
        fully_single = single_counts == counts

        # Fixed-Lagrangian-basin regional quantities.
        hratio = 1.0 + (D / 3.0) * (dve / ve)
        density_ratio = vl / ve
        ups2 = hratio*hratio * (ve / vl)
        gamma = ups2 - 1.0
        reff = np.full_like(ve, np.nan, dtype=float)
        pos = ve > 0.0
        reff[pos] = (3.0 * ve[pos] / (4.0*math.pi))**(1.0/3.0)
        gamma_grow = 1.0 + 2.0*Csig*XI0/(BETA_GAMMA*reff)
        kinematic = (ve > 0.0) & (hratio > 0.0) & (gamma >= gamma_grow)
        active = fully_single & kinematic

        newly = active & ~ever_active
        first_D[newly] = D
        ever_active |= active

        # Eulerian filling fractions (periodic ZA volume should remain close to box volume).
        valid_ve = np.where(ve > 0.0, ve, 0.0)
        phi_kin = float(valid_ve[kinematic].sum() / valid_ve.sum())
        phi_inst = float(valid_ve[active].sum() / valid_ve.sum())
        phi_abs = float(valid_ve[ever_active].sum() / valid_ve.sum())
        phi_single = float(valid_ve[fully_single].sum() / valid_ve.sum())

        active_gamma = gamma[active]
        active_r = reff[active]
        rows.append({
            "seed": seed,
            "persistence_rms": persistence_frac,
            "D": D,
            "n_basins": nlab,
            "coarse_rms": coarse_rms,
            "trace_error": trace_error,
            "volume_closure": float(valid_ve.sum()/total_box),
            "single_stream_volume_fraction": phi_single,
            "kinematic_upper_bound_volume_fraction": phi_kin,
            "instantaneous_active_volume_fraction": phi_inst,
            "absorbing_active_volume_fraction": phi_abs,
            "n_instantaneous_active": int(active.sum()),
            "n_ever_active": int(ever_active.sum()),
            "median_gamma_active": float(np.median(active_gamma)) if active_gamma.size else math.nan,
            "median_R_active_hMpc": float(np.median(active_r)) if active_r.size else math.nan,
            "median_density_ratio_all": float(np.median(density_ratio)),
            "median_Hratio_all": float(np.median(hratio)),
        })
    return rows, first_D


def run_all():
    all_rows = []
    all_first = []
    for seed in SEEDS:
        for p in PERSISTENCE:
            rows, first = basin_history(seed, p)
            all_rows.extend(rows)
            all_first.append((seed, p, first))
    return all_rows, all_first


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT/"results")
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    rows, firsts = run_all()

    table = args.out/"zeldovich_basin_history_table.csv"
    fields = list(rows[0].keys())
    with table.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    final = [r for r in rows if r["D"] == 1.0]
    kin = np.array([r["kinematic_upper_bound_volume_fraction"] for r in final])
    inst = np.array([r["instantaneous_active_volume_fraction"] for r in final])
    absorb = np.array([r["absorbing_active_volume_fraction"] for r in final])
    single = np.array([r["single_stream_volume_fraction"] for r in final])
    closure = np.array([r["volume_closure"] for r in rows])
    trace = np.array([r["trace_error"] for r in rows])

    # Monotonicity of absorbing history per seed/persistence.
    monotone = True
    for seed in SEEDS:
        for p in PERSISTENCE:
            x = [r["absorbing_active_volume_fraction"] for r in rows
                 if r["seed"] == seed and r["persistence_rms"] == p]
            if any(b < a - 1e-12 for a,b in zip(x,x[1:])):
                monotone = False

    # Compact history figure for manuscript/SI use.
    import matplotlib.pyplot as plt
    Ds = np.array(D_SNAPSHOTS, dtype=float)
    kin_lo=[]; kin_hi=[]; inst_lo=[]; inst_hi=[]; abs_lo=[]; abs_hi=[]
    for D in Ds:
        rr=[r for r in rows if r['D']==float(D)]
        for key,lo,hi in [
            ('kinematic_upper_bound_volume_fraction',kin_lo,kin_hi),
            ('instantaneous_active_volume_fraction',inst_lo,inst_hi),
            ('absorbing_active_volume_fraction',abs_lo,abs_hi)]:
            vals=np.array([r[key] for r in rr])
            lo.append(float(vals.min())); hi.append(float(vals.max()))
    fig,ax=plt.subplots(figsize=(5.4,3.7))
    ax.fill_between(Ds,kin_lo,kin_hi,alpha=0.16,label='kinematic upper envelope')
    ax.plot(Ds,0.5*(np.array(kin_lo)+np.array(kin_hi)),marker='o')
    ax.fill_between(Ds,inst_lo,inst_hi,alpha=0.16,label='strict instantaneous')
    ax.plot(Ds,0.5*(np.array(inst_lo)+np.array(inst_hi)),marker='s')
    ax.fill_between(Ds,abs_lo,abs_hi,alpha=0.16,label='history-dependent')
    ax.plot(Ds,0.5*(np.array(abs_lo)+np.array(abs_hi)),marker='^')
    ax.set_xlabel("Zel'dovich growth factor $D$")
    ax.set_ylabel('finite-size-active Eulerian volume fraction')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False,fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT/'figures/fig_zeldovich_history.pdf')
    plt.close(fig)

    summary = args.out/"zeldovich_basin_history_summary.txt"
    summary.write_text(
        "3D Zel'dovich Lagrangian-basin history stress test\n"
        f"box_hinvMpc={BOX}\nNgrid={NGRID}\nR_smooth_hinvMpc={R_SMOOTH}\n"
        f"seeds={SEEDS}\npersistence_rms={PERSISTENCE}\nD_snapshots={D_SNAPSHOTS}\n"
        f"sigma8_target={SIGMA8}\nconditional_xi0_hinvMpc={XI0}\nconditional_beta_gamma={BETA_GAMMA}\n"
        f"trace_identity_max_error={trace.max():.3e}\n"
        f"ZA_periodic_volume_closure_range={closure.min():.8f},{closure.max():.8f}\n"
        f"D1_single_stream_volume_fraction_range={single.min():.6f},{single.max():.6f}\n"
        f"D1_kinematic_upper_bound_volume_fraction_range={kin.min():.6f},{kin.max():.6f}\n"
        f"D1_instantaneous_active_volume_fraction_range={inst.min():.6f},{inst.max():.6f}\n"
        f"D1_absorbing_active_volume_fraction_range={absorb.min():.6f},{absorb.max():.6f}\n"
        f"absorbing_fraction_monotone={monotone}\n"
        "interpretation=Adding a history-bearing Zel'dovich deformation tests the selector in a minimally nonlinear Lagrangian setting. "
        "The calculation remains pre-shell-crossing and is not an N-body abundance prediction.\n",
        encoding="utf-8",
    )
    print(summary.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
