import math, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from analysis_leaf_local_field import mn_leaf_moments, q_simple

def test_scaling_solution_state_map():
    # Integrability for Q,R ~ a^n: (n+6)Q+(n+2)R=0.
    # Resulting effective geometric EOS is w=-(n+3)/3 and nu=-3w=n+3.
    for n in [-2.0,0.0,-1.0]:
        R=1.0; Q=-(n+2)/(n+6)*R
        w=(Q-R/3)/(Q+R)
        assert math.isclose(w,-(n+3)/3,rel_tol=1e-12,abs_tol=1e-12)
        nu=-3*w
        assert math.isclose(nu,n+3,rel_tol=1e-12,abs_tol=1e-12)
    # curvature-only and vacuum-like landmarks
    n=-2; Q=0; R=1
    assert math.isclose((Q-R/3)/(Q+R),-1/3)
    n=0; Q=-1/3; R=1
    assert math.isclose((Q-R/3)/(Q+R),-1)

def test_acceleration_boundary_identity():
    # With Gamma=rho_G/rho_m and w=-nu/3, acceleration requires
    # 1+Gamma(1+3w)<0 <=> Gamma(nu-1)>1.
    for G,nu in [(1,3),(0.4,3),(3,1),(2,2)]:
        w=-nu/3
        left=1+G*(1+3*w)
        right=G*(nu-1)-1
        assert math.isclose(right,-left,rel_tol=1e-12,abs_tol=1e-12)

def test_leaf_local_field_identity():
    d=mn_leaf_moments(1.0,0.1)
    gL=d['g_leaf']; aE=0.1
    q=float(q_simple(aE/gL)); alpha=q-1
    assert math.isclose(alpha*(1+alpha),aE/gL,rel_tol=1e-12,abs_tol=1e-12)
    assert 0<d['eta_S']<=1
    assert d['g_leaf_over_g_b_eq']>1

def test_leaf_far_field_spherical_recovery():
    d=mn_leaf_moments(20.0,0.1)
    assert abs(d['eta_S']-1)<1e-7
    assert abs(d['g_leaf_over_g_b_eq']-1)<2e-3

def test_morphon_reconstruction_identity():
    # Algebraic reconstruction of the effective geometry source.
    Gconst=1.0
    for Q,R in [(-0.2,-0.9),(0.1,-0.7),(-1/3,1.0)]:
        rho=-(Q+R)/(16*math.pi*Gconst)
        p=-(Q-R/3)/(16*math.pi*Gconst)
        K=-(Q+R/3)/(16*math.pi*Gconst)
        U=-R/(24*math.pi*Gconst)
        assert math.isclose(K+U,rho,rel_tol=1e-12,abs_tol=1e-12)
        assert math.isclose(K-U,p,rel_tol=1e-12,abs_tol=1e-12)

def test_exclusive_phase_partition_bookkeeping():
    total=set(range(20))
    C={2,3,4,10}
    V0={0,1,2,3,5,6,7,10,11}
    V=V0-C
    D=total-(C|V)
    assert not (C & V) and not (C & D) and not (V & D)
    assert C|V|D == total
    assert 2 not in V and 10 not in V

def test_one_scale_dimensional_closure():
    # If a_E=alpha c^2/L and eps=beta c^4/(G L^2), eliminating L gives eps=(beta/alpha^2)a_E^2/G.
    c=299792458.0; G=6.67430e-11; L=1.3e26; alpha=0.17; beta=0.08
    aE=alpha*c*c/L
    eps=beta*c**4/(G*L**2)
    assert math.isclose(eps,(beta/alpha**2)*aE*aE/G,rel_tol=1e-12)
