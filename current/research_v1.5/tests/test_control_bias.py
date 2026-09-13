from pathlib import Path
import importlib.util, math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('cb',ROOT/'src'/'analysis_control_bias.py')
cb=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(cb)

def test_basis_reconstruction():
    for G in [0.1,0.7,1.0,2.3]:
        for nu in [0.2,1.0,2.1,3.0,4.5]:
            q,r,G2,nu2,K,U=cb.identities(G,nu)
            assert abs(G-G2)<1e-12
            assert abs(nu-nu2)<1e-12

def test_curvature_and_vacuum_lines():
    q,r,*_=cb.identities(1.7,1.0)
    assert abs(q)<1e-12
    q,r,*_=cb.identities(1.7,3.0)
    assert abs(r-3*q)<1e-12

def test_acceleration_identity_threshold():
    # qhat is exactly the normalized Buchert kinematical backreaction;
    # qhat=1 is the dust+geometry acceleration boundary.
    q,r,*_=cb.identities(0.8,2.25)
    assert abs(q-1.0)<1e-12

def test_capillary_coefficient():
    expected=2*cb.C_SIGMA/(cb.CHI_V-cb.CHI_D)
    assert abs(expected-0.837187450022718)<1e-12

def test_reference_thresholds():
    eps_over=cb.OMEGA_V/cb.OMEGA_M
    qco=cb.H_COEX*eps_over
    hgrow=cb.H_COEX+(2*cb.C_SIGMA/cb.DELTA_CHI)*cb.XI0_BENCH/cb.R_BENCH
    qgrow=hgrow*eps_over
    assert abs(qco-0.743039)<5e-6
    assert abs(qgrow-1.1086)<5e-4
