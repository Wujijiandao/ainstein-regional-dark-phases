from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research_v1.7' / 'results' / 'causal_kinetics'


def test_persistent_sign_trichotomy_is_exhaustive_and_broad():
    df = pd.read_csv(OUT / 'persistent_sign_trichotomy.csv')
    total = df.persistent_positive_fraction + df.persistent_negative_fraction + df.sign_mixed_fraction
    assert np.allclose(total, 1.0, atol=1e-12)
    assert (df.persistent_positive_fraction.between(0.47, 0.53)).all()
    assert (df.persistent_negative_fraction.between(0.36, 0.40)).all()
    assert (df.persistent_one_sign_fraction > 0.85).all()


def test_passive_memory_sign_support_is_weakly_tau_sensitive():
    df = pd.read_csv(OUT / 'passive_memory_sign_audit.csv')
    assert (df.memory_positive_fraction.between(0.50, 0.58)).all()
    assert (df.memory_negative_fraction.between(0.42, 0.50)).all()


def test_causal_relaxation_has_lyapunov_decay_and_stable_subcausal_modes():
    p = json.loads((OUT / 'causal_kinetics_summary.json').read_text())
    e = p['causal_relaxation_test']
    l = p['linearized_causal_test']
    assert e['max_single_step_extended_energy_increase'] < 1e-12
    assert e['integrated_energy_balance_max_error'] < 1e-7
    assert 0 < l['c_phase_squared'] <= 1
    assert l['overdamped_discriminant_k0'] >= 0
    assert l['green_function_minimum'] > 0
    assert l['max_real_part_of_roots'] < 0


def test_fixed_sextic_threshold_requires_nontrivial_branch_odd_coupling():
    df = pd.read_csv(OUT / 'driver_coupling_inversion.csv')
    co10 = df[(df.threshold == 'coexistence') & np.isclose(df.target_fraction, 0.10)]
    assert np.isfinite(co10.V_required_gX).all()
    assert np.isfinite(co10.C_required_gX).all()
    assert (co10.V_required_gX > 3.0).all()
    assert (co10.C_required_gX > 3.0).all()


def test_spinodal_is_harder_than_coexistence_and_cannot_cover_half_negative_volume():
    df = pd.read_csv(OUT / 'driver_coupling_inversion.csv')
    key = ['dataset','seed','tauN','target_fraction']
    co = df[df.threshold == 'coexistence'].set_index(key)
    sp = df[df.threshold == 'spinodal'].set_index(key)
    common = co.index.intersection(sp.index)
    assert np.all(sp.loc[common,'V_required_gX'].to_numpy() >= co.loc[common,'V_required_gX'].to_numpy())
    halfC = df[(df.threshold == 'coexistence') & np.isclose(df.target_fraction,0.50)]
    assert np.isinf(halfC.C_required_gX).all()
