#!/usr/bin/env python3
"""Interface-coarsening identities and literature benchmark.

This module is deliberately split into (i) exact geometric/effective-EOS
identities and (ii) an external benchmark using published cosmic-web summary
numbers from Cautun et al. (MNRAS 441, 2923, 2014) and the public SDSS DR7
void-catalog summary of Douglass et al. (ApJS 265, 7, 2023).

The benchmark is not a cosmological fit. The NEXUS wall network is only a proxy
for the boundary area of the proposed V-phase union.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Published approximate benchmark values.
PHI_VOID_Z2 = 0.60
PHI_VOID_Z0 = 0.78
A_WALL_Z0_PER_GPC3 = 1.5e8  # (h^-1 Mpc)^2 per (h^-1 Gpc)^3
V_GPC3 = 1.0e9              # (h^-1 Mpc)^3
S_COM_Z0 = A_WALL_Z0_PER_GPC3 / V_GPC3  # (h^-1 Mpc)^-1
VAST_MEDIAN_R_MIN = 15.0     # h^-1 Mpc
VAST_MEDIAN_R_MAX = 19.0     # h^-1 Mpc


def hydraulic_interface_length(phi: float, s_phys: float) -> float:
    """L_Sigma = 3 phi / s_phys; equals sphere radius for equal disjoint spheres."""
    if not (0.0 < phi <= 1.0):
        raise ValueError("phi must lie in (0,1]")
    if s_phys <= 0.0:
        raise ValueError("interface-area density must be positive")
    return 3.0 * phi / s_phys


def apparent_w_phase_sector(p_phi: float, q_L: float, lambda_sigma: float) -> float:
    """Apparent EOS of epsilon*phi + sigma*s with lambda_sigma=(sigma/epsilon)/L.

    p_phi = d ln phi / d ln a
    q_L   = d ln L / d ln a
    lambda_sigma = ell_sigma/L, ell_sigma=sigma/epsilon.
    """
    if lambda_sigma < 0.0:
        raise ValueError("lambda_sigma must be nonnegative for positive tension")
    x = 3.0 * lambda_sigma
    return -1.0 - (p_phi - q_L * x / (1.0 + x)) / 3.0


def crossing_lambda(p_phi: float, q_L: float) -> float:
    """Positive-tension lambda required for w=-1, if q_L>p_phi."""
    if q_L <= p_phi:
        return math.inf
    return p_phi / (3.0 * (q_L - p_phi))


def literature_benchmark() -> dict[str, float]:
    # a(z=2)=1/3, a(z=0)=1. NEXUS reports void volume fraction rising
    # from ~0.60 to 0.78 and total *comoving* wall area declining.
    dln_a = math.log(3.0)
    p_avg = math.log(PHI_VOID_Z0 / PHI_VOID_Z2) / dln_a

    # s_com(z=0) < s_com(z=2). Since L_phys = 3 a phi / s_com,
    # this gives a strict lower bound on the interval-averaged q_L.
    q_lower = math.log(3.0 * PHI_VOID_Z0 / PHI_VOID_Z2) / dln_a

    L0 = hydraulic_interface_length(PHI_VOID_Z0, S_COM_Z0)  # a0=1
    lam_crit_upper = crossing_lambda(p_avg, q_lower)
    ell_crit_upper = lam_crit_upper * L0
    return {
        "p_phi_avg": p_avg,
        "q_L_lower": q_lower,
        "L_sigma_z0_hMpc": L0,
        "lambda_cross_upper": lam_crit_upper,
        "ell_cross_upper_hMpc": ell_crit_upper,
    }


def make_figure(path: Path) -> None:
    b = literature_benchmark()
    lam = np.linspace(0.0, 0.25, 400)
    w = np.array([apparent_w_phase_sector(b["p_phi_avg"], b["q_L_lower"], x) for x in lam])

    fig = plt.figure(figsize=(6.6, 4.2))
    ax = fig.add_subplot(111)
    ax.plot(lam, w, lw=2, label="NEXUS interval lower-bound coarsening")
    ax.axhline(-1.0, ls="--", lw=1)
    ax.axvline(b["lambda_cross_upper"], ls=":", lw=1.5,
               label=fr"crossing at $\ell_\Sigma/L_\Sigma\lesssim {b['lambda_cross_upper']:.3f}$")
    ax.set_xlabel(r"Interface energy length ratio $\ell_\Sigma/L_\Sigma$")
    ax.set_ylabel(r"Apparent $w_{V+\Sigma}$")
    ax.set_ylim(-1.12, -0.72)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Interface coarsening can reverse the conversion-only EOS sign")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    results = root / "results"
    figures = root / "figures"
    results.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)

    b = literature_benchmark()
    lines = [
        "Interface-coarsening benchmark (not a cosmological fit)",
        f"phi_void(z=2) ~= {PHI_VOID_Z2:.3f}",
        f"phi_void(z=0) ~= {PHI_VOID_Z0:.3f}",
        f"wall area density z=0 ~= {S_COM_Z0:.6f} (h^-1 Mpc)^-1",
        f"hydraulic interface length z=0 = {b['L_sigma_z0_hMpc']:.12f} h^-1 Mpc",
        f"VAST median void effective radius = {VAST_MEDIAN_R_MIN:.1f}-{VAST_MEDIAN_R_MAX:.1f} h^-1 Mpc",
        f"interval average p_phi = {b['p_phi_avg']:.12f}",
        f"interval lower bound q_L > {b['q_L_lower']:.12f}",
        f"crossing lambda upper benchmark = {b['lambda_cross_upper']:.12f}",
        f"crossing ell upper benchmark = {b['ell_cross_upper_hMpc']:.12f} h^-1 Mpc",
    ]
    (results / "interface_coarsening_summary.txt").write_text("\n".join(lines) + "\n")

    rows=[]
    for lam in [0.0, 0.04, b["lambda_cross_upper"], 0.10, 0.20]:
        rows.append((lam, apparent_w_phase_sector(b["p_phi_avg"], b["q_L_lower"], lam)))
    np.savetxt(results / "interface_coarsening_table.csv", np.array(rows), delimiter=",",
               header="ell_over_L,w_apparent", comments="")
    make_figure(figures / "fig_interface_coarsening.pdf")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
