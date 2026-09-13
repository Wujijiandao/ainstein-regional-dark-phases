from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'research_v1.8'/'results'/'retarded_response'

import importlib.util
MOD = ROOT/'research_v1.8'/'src'/'analysis_retarded_response.py'
spec = importlib.util.spec_from_file_location('retarded_response_module', MOD)
rr = importlib.util.module_from_spec(spec); spec.loader.exec_module(rr)


def load():
    return json.loads((OUT/'retarded_response_summary.json').read_text())


def test_passive_spectral_weight_and_kramers_kronig():
    p=load()['spectral_checks']
    assert p['minimum_Im_G_for_positive_omega_k0'] > 0
    assert p['KK_relative_error'] < 1e-8


def test_inverse_response_recovers_transport_coefficients():
    p=load(); q=p['inverse_response_identification']; b=p['benchmark_input']
    assert abs(q['m2_hat']/b['m2']-1) < 0.01
    assert abs(q['kappa_hat']/b['kappa']-1) < 0.08
    assert abs(q['Lambda_hat']/b['Lambda']-1) < 0.01
    assert abs(q['tau_hat']/b['tau']-1) < 0.08
    assert abs(q['c2_hat']/b['c2']-1) < 0.12


def test_cross_kernel_ratio_recovers_odd_susceptibility_and_memory():
    p=load(); q=p['cross_response_identification']; b=p['benchmark_input']
    assert abs(q['gX_hat']-b['gX']) < 1e-12
    assert abs(q['tauX_hat']-b['tauX']) < 1e-12


def test_even_gate_is_identified_by_inverse_static_susceptibility_shift():
    p=load(); q=p['even_gate_identification']; b=p['benchmark_input']
    assert abs(q['m2_intercept_hat']-b['m2']) < 1e-12
    assert abs(q['lambdaE_hat']-b['lambda_even_gate']) < 1e-12


def test_selection_rules_are_explicit():
    s=load()['structural_selection_rules']
    assert 'G_chiE = 0' in s['linear_even_to_odd_response_at_symmetric_parent']
    assert 'G_chiX' in s['odd_kernel_measurement']
