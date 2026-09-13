import json
from pathlib import Path
import numpy as np
import importlib.util

ROOT=Path(__file__).resolve().parents[2]
MOD=ROOT/'research_v2.0'/'src'/'analysis_phase_geometry_feedback.py'
spec=importlib.util.spec_from_file_location('fb',MOD)
fb=importlib.util.module_from_spec(spec); spec.loader.exec_module(fb)


def test_static_elimination_identity():
    assert fb.identity_max_error < 1e-10


def test_feedback_thresholds():
    assert abs(fb.alpha_coex_at_beta0-0.625) < 1e-12
    assert 6.3 < fb.beta_coex_at_alpha0 < 6.5
    assert abs(fb.alpha_parent_instability-4.0) < 1e-12


def test_reciprocal_detailed_balance():
    assert fb.rec_current_norm < 1e-10
    assert fb.rec_cov_error < 1e-10


def test_nonreciprocal_current_and_stability():
    assert fb.non_current_norm > 1e-3
    assert fb.non_current_power > 1e-5
    assert fb.non_max_real_drift < 0


def test_conserved_carrier_identity():
    assert fb.carrier_identity_error < 1e-12
