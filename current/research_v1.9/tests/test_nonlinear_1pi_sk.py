import importlib.util
from pathlib import Path
import numpy as np

P=Path(__file__).resolve().parents[1]/'src'/'analysis_nonlinear_1pi_sk.py'
spec=importlib.util.spec_from_file_location('a',P)
a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)

def test_exact_1pi_vertices_and_thresholds():
    assert a.R0 == 4.0
    assert a.U4 == -72.0
    assert a.U6 == 960.0
    assert abs(a.r_daughter_spinodal-4.5) < 1e-12
    assert abs(a.r_coex-3.375) < 1e-12
    assert abs(a.chi_coex_sq-1.125) < 1e-12

def test_factorization_and_wall_constants():
    assert a.factorization_maxerr < 1e-12
    assert abs(a.C_sigma_sym-0.5166892426183266) < 1e-12
    assert abs(a.width_10_90_sym-1.645275999066359) < 1e-12

def test_higher_order_response_reconstructs_vertices():
    assert abs(a.r_from_R-a.R0) < 1e-12
    assert abs(a.u_from_R-a.U4) < 1e-12
    assert abs(a.v_from_R-a.U6) < 1e-10
    assert abs(a.r_hat-a.R0) < 5e-3
    assert abs(a.u_hat-a.U4) < 0.2
    assert abs(a.v_hat-a.U6) < 20.0

def test_thin_wall_asymptotics_improve_near_coexistence():
    n=a.nuc.sort_values('delta_r')
    first=n.iloc[0]
    last=n.iloc[-1]
    err_first=abs(first.Fstar_over_E0xi3_exact/first.Fstar_over_E0xi3_asym-1)
    err_last=abs(last.Fstar_over_E0xi3_exact/last.Fstar_over_E0xi3_asym-1)
    assert err_first < err_last
    assert a.Rc_coeff > 0 and a.Fstar_coeff > 0

def test_local_kms_fdt_recovers_noise_scale():
    assert a.theta_maxerr < 1e-12
    assert np.all(a.fdt.S_chichi > 0)
