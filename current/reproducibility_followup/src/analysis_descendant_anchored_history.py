#!/usr/bin/env python3
"""Descendant-anchored nonlinear basin history for the AInstein follow-up.

This deliberately replaces pairwise snapshot-to-snapshot watershed matching.
A watershed catalogue is defined only at the terminal a=1 snapshot. Particle
membership in each terminal basin is then frozen and traced back through the
same nonlinear PM trajectory.

For each terminal basin and earlier snapshot, Eulerian volume is estimated from
tracked particles using inverse smoothed-density volume weights. Regional mean
density and expansion are then evaluated on that volume-weighted descendant
support. This gives a descendant-anchored history of Gamma, finite-size activity,
first passage, and filling fraction without a greedy catalogue-matching step.

The construction is an exploratory obstruction test, not a precision void
catalogue or a calibrated cosmological phase model.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analysis_nonlinear_pm_basin_lineage as pm


def weighted_bincount(labels: np.ndarray, weights: np.ndarray, nlab: int) -> np.ndarray:
    return np.bincount(labels, weights=weights, minlength=nlab + 1)[1:].astype(float)


def descendant_stats(pos: np.ndarray, p: np.ndarray, a: float, final_particle_labels: np.ndarray,
                     final_nlab: int, fg: pm.FourierGrid, R_smooth: float):
    """Stats for fixed final-basin particle memberships at one earlier snapshot.

    Particle mass is one code unit. A local particle volume proxy is 1/rho_sm(x_i).
    The proxy volumes are globally renormalized so their sum equals the box volume.
    This makes the anchored basin volumes additive by construction while retaining
    the relative inverse-density volume weighting.
    """
    rho, rho_sm, theta = pm.velocity_divergence_theta(pos, p, a, fg, R_smooth)
    rho_p = pm.cic_interpolate(rho_sm, pos, fg.n, fg.box)
    theta_p = pm.cic_interpolate(theta, pos, fg.n, fg.box)
    rho_p = np.maximum(rho_p, 1e-5)

    raw_vol = 1.0 / rho_p
    norm = fg.box**3 / raw_vol.sum()
    vol_w = raw_vol * norm

    counts = np.bincount(final_particle_labels, minlength=final_nlab + 1)[1:].astype(float)
    volumes = weighted_bincount(final_particle_labels, vol_w, final_nlab)
    # mean density relative to background: basin mass / basin volume, with total mass=Np
    mean_rho = (counts / counts.sum()) / (volumes / volumes.sum())
    mean_theta = weighted_bincount(final_particle_labels, vol_w * theta_p, final_nlab) / volumes
    hratio = 1.0 + mean_theta / 3.0
    valid = (volumes > 0) & (mean_rho > 0) & (hratio > 0)
    gamma = np.full(final_nlab, -np.inf, dtype=float)
    gamma[valid] = hratio[valid]**2 / mean_rho[valid] - 1.0
    reff = (3.0 * volumes / (4.0 * math.pi))**(1.0 / 3.0)
    gamma_grow = 1.0 + 2.0 * pm.C_SIGMA * pm.XI0 / (pm.BETA_GAMMA * reff)
    active = valid & (gamma >= gamma_grow)

    return {
        'rho_grid': rho,
        'rho_sm': rho_sm,
        'theta_grid': theta,
        'counts': counts,
        'volumes': volumes,
        'mean_rho': mean_rho,
        'mean_theta': mean_theta,
        'hratio': hratio,
        'gamma': gamma,
        'reff': reff,
        'gamma_grow': gamma_grow,
        'active': active,
        'volume_closure': float(volumes.sum() / fg.box**3),
        'mass_closure': float(counts.sum() / len(pos)),
    }


def simulate_trajectory(seed: int, n: int, steps: int):
    fg = pm.fourier_grid(n, pm.BOX)
    delta0 = pm.gaussian_linear_field(seed, fg)
    sigma0 = pm.sigma8_grid(delta0, fg)
    s = pm.displacement_from_delta(delta0, fg)
    q = pm.grid_positions(n, pm.BOX)
    pos = (q + pm.A_INIT * s) % pm.BOX
    p = (pm.A_INIT ** 1.5) * s

    base = np.linspace(pm.A_INIT, 1.0, steps + 1)
    nodes = np.unique(np.round(np.concatenate([base, np.array(pm.SNAPSHOTS)]), 10))
    nodes.sort()
    snapset = {round(x, 10) for x in pm.SNAPSHOTS}
    snapshots = {}
    a = pm.A_INIT
    for a2 in nodes[1:]:
        pos, p = pm.pm_step(pos, p, a, float(a2), fg)
        a = float(a2)
        if round(a, 10) in snapset:
            snapshots[round(a, 10)] = (pos.copy(), p.copy())
    return fg, delta0, sigma0, snapshots


def grid_catalogue_final(pos: np.ndarray, p: np.ndarray, fg: pm.FourierGrid,
                         persistence: float, R_smooth: float):
    rho, rho_sm, theta = pm.velocity_divergence_theta(pos, p, 1.0, fg, R_smooth)
    labels, rms = pm.segment(rho_sm - 1.0, persistence)
    st = pm.basin_stats(labels, rho, theta, fg.box / fg.n)
    particle_labels = pm.particle_basin_labels(pos, labels, fg.box)
    return labels, particle_labels, st, rms


def overlap_diagnostics(final_labels: np.ndarray, anchored_volumes: np.ndarray, grid_stats: dict, fg: pm.FourierGrid):
    # Compare terminal anchored volume estimates to terminal Eulerian watershed cell volumes.
    cell_vol = (fg.box / fg.n) ** 3
    grid_vol = grid_stats['counts'] * cell_vol
    # both are additive; quantify fractional scatter and rank correlation.
    ratio = np.divide(anchored_volumes, grid_vol, out=np.full_like(grid_vol, np.nan), where=grid_vol > 0)
    finite = np.isfinite(ratio) & (ratio > 0)
    if finite.sum() > 2:
        x = grid_vol[finite]
        y = anchored_volumes[finite]
        corr = float(np.corrcoef(np.log(x), np.log(y))[0, 1])
    else:
        corr = math.nan
    return {
        'terminal_median_anchored_to_grid_volume_ratio': float(np.nanmedian(ratio)),
        'terminal_p16_anchored_to_grid_volume_ratio': float(np.nanpercentile(ratio, 16)),
        'terminal_p84_anchored_to_grid_volume_ratio': float(np.nanpercentile(ratio, 84)),
        'terminal_log_volume_correlation': corr,
    }


def run_seed(seed: int, n: int, steps: int, persistence_values: tuple[float, ...], R_smooth: float):
    fg, delta0, sigma0, snaps = simulate_trajectory(seed, n, steps)
    final_pos, final_p = snaps[1.0]
    rows = []
    basin_rows = []
    validation = {'seed': seed, 'Ngrid': n, 'R_smooth_hMpc': R_smooth,
                  'sigma8_linear_target_a1': sigma0, 'particle_count': int(n**3)}

    for pers in persistence_values:
        final_grid_labels, final_particle_labels, final_grid_stats, rms = grid_catalogue_final(
            final_pos, final_p, fg, pers, R_smooth)
        nlab = final_grid_stats['nlab']
        ever_active = np.zeros(nlab, dtype=bool)
        first_passage = np.full(nlab, np.nan)
        per_a = {}
        terminal_diag = None

        for a in pm.SNAPSHOTS:
            pos, p = snaps[round(a, 10)]
            st = descendant_stats(pos, p, a, final_particle_labels, nlab, fg, R_smooth)
            active = st['active']
            newly = active & ~ever_active
            first_passage[newly] = a
            ever_active |= active
            per_a[a] = st

            phi_inst = float(st['volumes'][active].sum() / fg.box**3)
            phi_absorb = float(st['volumes'][ever_active].sum() / fg.box**3)
            # Geometry-free descendant accounting has no cell boundary; report support volume only.
            rows.append({
                'seed': seed, 'Ngrid': n, 'persistence_rms_final': pers, 'R_smooth_hMpc': R_smooth,
                'a': a, 'n_final_basins': nlab, 'final_delta_sm_rms': rms,
                'anchored_instantaneous_active_volume_fraction': phi_inst,
                'anchored_absorbing_active_volume_fraction': phi_absorb,
                'active_basin_count': int(active.sum()),
                'ever_active_basin_count': int(ever_active.sum()),
                'median_gamma_all': float(np.median(st['gamma'][np.isfinite(st['gamma'])])) if np.any(np.isfinite(st['gamma'])) else math.nan,
                'median_gamma_active': float(np.median(st['gamma'][active])) if np.any(active) else math.nan,
                'median_R_all_hMpc': float(np.median(st['reff'])),
                'median_R_active_hMpc': float(np.median(st['reff'][active])) if np.any(active) else math.nan,
                'volume_closure': st['volume_closure'], 'mass_closure': st['mass_closure'],
            })

            for i in range(nlab):
                basin_rows.append({
                    'seed': seed, 'Ngrid': n, 'persistence_rms_final': pers, 'R_smooth_hMpc': R_smooth,
                    'a': a, 'final_basin_id': i + 1,
                    'mass_fraction': float(st['counts'][i] / st['counts'].sum()),
                    'volume_fraction': float(st['volumes'][i] / fg.box**3),
                    'R_eff_hMpc': float(st['reff'][i]), 'mean_density_ratio': float(st['mean_rho'][i]),
                    'Hloc_over_H': float(st['hratio'][i]), 'Gamma': float(st['gamma'][i]),
                    'Gamma_grow': float(st['gamma_grow'][i]), 'instantaneous_active': int(active[i]),
                    'absorbing_active': int(ever_active[i]),
                })

            if abs(a - 1.0) < 1e-12:
                terminal_diag = overlap_diagnostics(final_grid_labels, st['volumes'], final_grid_stats, fg)
                validation[f'p{pers:g}_terminal_grid_phi'] = float(final_grid_stats['counts'][final_grid_stats['active']].sum() / final_grid_labels.size)
                validation[f'p{pers:g}_terminal_anchored_phi'] = phi_inst
                validation[f'p{pers:g}_terminal_anchor_vs_grid'] = terminal_diag

        # append first-passage information as extra basin rows at a sentinel? keep separate summary later.
        fp = first_passage[np.isfinite(first_passage)]
        validation[f'p{pers:g}_first_passage_counts'] = {
            str(a): int(np.count_nonzero(np.isclose(first_passage, a, equal_nan=False))) for a in pm.SNAPSHOTS
        }
        validation[f'p{pers:g}_ever_active_final_basin_fraction'] = float(np.mean(np.isfinite(first_passage)))
        validation[f'p{pers:g}_median_first_passage_a'] = float(np.median(fp)) if len(fp) else None

    return rows, basin_rows, validation


def write_outputs(out: Path, rows: list[dict], basin_rows: list[dict], validations: list[dict]):
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'descendant_anchored_history_summary_table.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with (out / 'descendant_anchored_basin_histories.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(basin_rows[0].keys())); w.writeheader(); w.writerows(basin_rows)
    (out / 'descendant_anchored_validation.json').write_text(json.dumps(validations, indent=2), encoding='utf-8')

    finals = [r for r in rows if abs(r['a'] - 1.0) < 1e-12]
    lines = [
        'AInstein descendant-anchored basin history stress test',
        'terminal_catalogue=a=1 Eulerian watershed; particle membership frozen and traced backward',
        'regional_volume_estimator=inverse-smoothed-density particle volume weights, globally volume-normalized',
        'history=instantaneous anchored activity plus absorbing first-passage on the same terminal descendants',
    ]
    for pers in sorted({r['persistence_rms_final'] for r in finals}):
        sub = [r for r in finals if r['persistence_rms_final'] == pers]
        for key in ['anchored_instantaneous_active_volume_fraction', 'anchored_absorbing_active_volume_fraction']:
            vals = np.array([r[key] for r in sub])
            lines.append(f'p{pers:g}_D1_{key}_range={vals.min():.8g},{vals.max():.8g};mean={vals.mean():.8g}')
    (out / 'descendant_anchored_history_summary.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1] / 'results' / 'descendant_anchored')
    ap.add_argument('--n', type=int, default=64)
    ap.add_argument('--steps', type=int, default=38)
    ap.add_argument('--seeds', type=int, nargs='+', default=[11, 23, 47, 91])
    ap.add_argument('--persistence', type=float, nargs='+', default=[0.0, 0.10, 0.20])
    ap.add_argument('--R-smooth', type=float, default=pm.R_SMOOTH)
    args = ap.parse_args()
    rows=[]; basin_rows=[]; vals=[]
    for seed in args.seeds:
        print(f'=== descendant-anchored seed={seed} N={args.n} R={args.R_smooth} ===', flush=True)
        r,b,v = run_seed(seed,args.n,args.steps,tuple(args.persistence),args.R_smooth)
        rows.extend(r); basin_rows.extend(b); vals.append(v)
        print(json.dumps(v, indent=2), flush=True)
    write_outputs(args.out, rows, basin_rows, vals)
    print((args.out/'descendant_anchored_history_summary.txt').read_text(), flush=True)

if __name__ == '__main__':
    main()
