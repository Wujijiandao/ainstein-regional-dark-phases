from pathlib import Path
import sys, math, numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from analysis_phase import phase_constants, mn_eta
from analysis_regional_state import overdense, underdense, void_shell_crossing_state
from analysis_excursion_phase import first_crossing_fractions, eulerian_volume_fractions, apparent_w_void
from analysis_power_clock import sigma_R, transfer_bbks, transfer_eh0, R_when_S, TAU_HALF_V, SIGMA8, required_void_barrier_for_half_volume
from analysis_curvature_selector import (UPSILON_CURV_EQ, DELTA_V_CURV_EQ, RHO_CURV_EQ, TAU_HALF_CURV, eulerian_volume_fractions_curvature, gamma_estimator)


def test_phase_constants():
    H,xc,xs,jump,Hspin=phase_constants()
    assert abs(H-0.341688495018)<5e-12
    assert abs(jump-0.948007800871)<5e-12
    assert abs(Hspin-0.928501997833)<5e-12


def test_geometry_factor():
    assert abs(mn_eta(.2,.1)-0.923361539)<5e-9
    assert 0 < mn_eta(1,.3) <= 1


def test_spherical_state_identities():
    th=np.linspace(.08,2*math.pi-.08,200)
    _,_,u,ue,_=overdense(th)
    assert np.max(np.abs(u-ue))<1e-11
    et=np.linspace(.08,3.5,200)
    _,_,u,ue,_=underdense(et)
    assert np.max(np.abs(u-ue))<1e-11


def test_void_shell_crossing_landmark():
    eta,d,u,dl=void_shell_crossing_state()
    assert abs(dl+2.717)<1e-10
    assert 0.19 < d < 0.22
    assert 2.8 < u < 3.1


def test_two_barrier_first_crossing():
    c,v,d=first_crossing_fractions(3.0)
    assert abs(c-0.330308863092)<2e-10
    assert abs(v-0.116288074635)<2e-10
    assert abs(c+v+d-1)<1e-12
    c2,v2,d2=first_crossing_fractions(20.0)
    assert c2>c and v2>v and d2<d


def test_eulerian_void_volume_crosses_half_near_tau3():
    _,v2,_=eulerian_volume_fractions(2.0)
    _,v3,_=eulerian_volume_fractions(3.0)
    assert v2 < 0.5 < v3


def test_growing_void_fraction_is_apparently_phantom():
    assert apparent_w_void(3.0) < -1.0
    assert apparent_w_void(16.0) > apparent_w_void(3.0)
    assert apparent_w_void(16.0) < -1.0


def test_power_clock_shell_crossing_obstruction():
    # The 50% Eulerian-V landmark requires sigma=sqrt(tau_half)>1.7,
    # while Planck-normalized sigma8 is only 0.811.
    assert TAU_HALF_V > 4.0 * SIGMA8**2
    assert eulerian_volume_fractions(SIGMA8**2)[1] < 0.005


def test_power_clock_transfer_shape_robustness():
    r1=R_when_S(TAU_HALF_V, transfer_bbks)
    r2=R_when_S(TAU_HALF_V, transfer_eh0)
    assert 2.0 < r1 < 2.3
    assert 2.3 < r2 < 2.6
    # Typical shell-crossed void Lagrangian scales (~9-11 Mpc/h) are low variance.
    for fn in (transfer_bbks, transfer_eh0):
        assert sigma_R(9.0, fn) < SIGMA8
        assert sigma_R(11.0, fn) < 0.7


def test_inverse_void_selector_target():
    dv,rho,ups=required_void_barrier_for_half_volume(SIGMA8**2)
    assert -0.80 < dv < -0.77
    assert 0.52 < rho < 0.54
    assert 1.60 < ups < 1.64


def test_curvature_matter_equality_landmark():
    assert abs(UPSILON_CURV_EQ-math.sqrt(2.0)) < 1e-14
    assert abs(DELTA_V_CURV_EQ + 0.5167454041616256) < 5e-12
    assert abs(RHO_CURV_EQ - 0.6388164885032406) < 5e-12


def test_curvature_selector_hits_required_variance_regime():
    assert abs(TAU_HALF_CURV - 0.3582167075026434) < 5e-12
    # At the median shell-crossing-mapped void Lagrangian interval, the
    # curvature-equality candidate is already near/above half Eulerian volume.
    for fn in (transfer_bbks, transfer_eh0):
        p_lo=eulerian_volume_fractions_curvature(sigma_R(8.839128, fn)**2)[1]
        p_hi=eulerian_volume_fractions_curvature(sigma_R(11.196229, fn)**2)[1]
        assert 0.60 < p_lo < 0.63
        assert 0.50 < p_hi < 0.55

from analysis_basin_selector import synthetic_hierarchy, maximal_eligible, is_antichain, naive_eligible_volume, antichain_volume


def test_maximal_eligible_basin_rule_is_nonrecursive():
    basins = synthetic_hierarchy()
    selected = maximal_eligible(basins)
    assert [b.name for b in selected] == ["A", "B1"]
    assert is_antichain(selected, basins)
    assert abs(antichain_volume(basins) - 0.75) < 1e-12
    assert naive_eligible_volume(basins) > 1.0  # nested double counting without the rule


def test_operational_gamma_estimator():
    # Background region.
    assert abs(gamma_estimator(0.0, 0.0)) < 1e-14
    # A moderately underdense, faster-expanding region should have Gamma > 0.
    assert gamma_estimator(-0.3, 0.15) > 0.0
    # The estimator is exactly the algebraic Upsilon^2-1 definition.
    delta=-0.2
    q=0.12
    u=(1.0+q/3.0)/math.sqrt(1.0+delta)
    assert abs(gamma_estimator(delta,q)-(u*u-1.0)) < 1e-14


from analysis_conversion_sign import residual_w_from_history

def test_monotone_conversion_phantom_sign_theorem_numeric():
    a=np.geomspace(1e-6,1.0,4000)
    # Exact continuum result for phi=a^2 is w=-(2+3)/3=-5/3.
    w=residual_w_from_history(a,a*a)
    assert abs(w+5.0/3.0) < 2e-5
    # A monotone saturating history remains phantom-or-Lambda.
    phi=1.0-np.exp(-4*a)
    assert residual_w_from_history(a,phi) <= -1.0


from analysis_basin_dynamics import (hierarchy_template, synthetic_history, instantaneous_active, absorbing_active_history, maximal_from_active, union_volume_of_active, switch_count)

def test_history_dependent_basin_activation_is_nonrecursive_and_monotone():
    basins=hierarchy_template()
    states=synthetic_history()
    inst=[instantaneous_active(s,basins) for s in states]
    hist=absorbing_active_history(states,basins)
    phi=[union_volume_of_active(a,basins) for a in hist]
    assert all(b >= a - 1e-14 for a,b in zip(phi,phi[1:]))
    assert switch_count(hist) < switch_count(inst)
    final=maximal_from_active(hist[-1],basins)
    assert [b.name for b in final] == ["A","B"]
    assert abs(phi[-1]-1.0) < 1e-14

def test_maximal_antichain_represents_active_union_under_absorbing_growth():
    basins=hierarchy_template()
    hist=absorbing_active_history(synthetic_history(),basins)
    # Maximal-node volume is the laminar-union volume and never exceeds box volume.
    for active in hist:
        assert 0.0 <= union_volume_of_active(active,basins) <= 1.0 + 1e-14


from analysis_interface_coarsening import (
    hydraulic_interface_length, apparent_w_phase_sector, crossing_lambda,
    literature_benchmark, PHI_VOID_Z0, S_COM_Z0
)

def test_interface_hydraulic_length_cross_catalog_scale():
    b=literature_benchmark()
    assert abs(hydraulic_interface_length(PHI_VOID_Z0,S_COM_Z0)-15.6) < 1e-12
    assert 15.0 <= b["L_sigma_z0_hMpc"] <= 19.0

def test_interface_coarsening_crossing_identity():
    b=literature_benchmark()
    lam=crossing_lambda(b["p_phi_avg"],b["q_L_lower"])
    w=apparent_w_phase_sector(b["p_phi_avg"],b["q_L_lower"],lam)
    assert abs(w+1.0) < 1e-12
    assert apparent_w_phase_sector(b["p_phi_avg"],b["q_L_lower"],0.04) < -1.0
    assert apparent_w_phase_sector(b["p_phi_avg"],b["q_L_lower"],0.12) > -1.0

def test_interface_benchmark_requires_mesoscale_not_horizon_scale():
    b=literature_benchmark()
    assert 0.07 < b["lambda_cross_upper"] < 0.09
    assert 1.0 < b["ell_cross_upper_hMpc"] < 1.5

from analysis_phase_interface import interface_constants, required_phase_width

def test_sextic_phase_field_interface_constants():
    c=interface_constants()
    assert abs(c['C_sigma']-0.396830116706) < 5e-12
    assert abs(c['C_10_90']-1.693129938128) < 5e-12
    assert abs(c['thickness_to_ell_if_E0_eq_epsV']-4.266636696288) < 5e-12

def test_interface_length_maps_to_mpc_wall_width_conditionally():
    d=required_phase_width(1.241834074954,1.0)
    assert 5.29 < d < 5.31


from analysis_phase_front import (
    well_difference, critical_radius_over_xi, barrier_over_E0_xi3,
    dimensionless_front_speed, H_for_critical_radius, linear_near_coex_coefficient, benchmark
)

def test_thin_interface_critical_radius_and_speed_signs():
    b=benchmark()
    H=b['H_for_target_Rc']; rc=b['target_R_over_xi']
    assert abs(critical_radius_over_xi(H)-rc) < 1e-10
    assert dimensionless_front_speed(H,0.5*rc) < 0.0
    assert abs(dimensionless_front_speed(H,rc)) < 1e-10
    assert dimensionless_front_speed(H,2.0*rc) > 0.0

def test_near_coexistence_nucleation_scaling():
    from analysis_phase import phase_constants
    H0,*_=phase_constants()
    coeff=linear_near_coex_coefficient()
    dH=1e-4
    rc=critical_radius_over_xi(H0+dH)
    assert abs(rc*dH/coeff-1.0) < 5e-4
    assert barrier_over_E0_xi3(H0+0.1) > 0.0

def test_conditional_void_radius_requires_moderate_supercoexistence_bias():
    b=benchmark()
    assert 0.16 < b['Delta_H_for_target_Rc'] < 0.18
    assert 4.9 < b['target_R_over_xi'] < 5.1
    assert 0.15 < b['Delta_w_for_target_Rc'] < 0.17


from analysis_bias_matching import (
    H_from_gamma, gamma_for_radius, radius_over_xi_from_gamma,
    near_coexistence_bias_slope, benchmark as bias_matching_benchmark
)

def test_bias_matching_coexistence_and_slope():
    from analysis_phase import phase_constants
    Hc,_,_,jump,_=phase_constants()
    assert abs(H_from_gamma(1.0)+Hc) < 1e-13
    assert abs(near_coexistence_bias_slope(1.0)-1.0/jump) < 1e-13
    eps=1e-5
    B=(-H_from_gamma(1.0+eps))-Hc
    assert abs(B/eps-near_coexistence_bias_slope(1.0)) < 2e-4


def test_bias_matching_finite_size_shift_identity():
    b=bias_matching_benchmark()
    g=b['Gamma_growth_threshold']; r=b['R_over_xi']
    assert abs(radius_over_xi_from_gamma(g,1.0)-r) < 1e-11
    assert abs(gamma_for_radius(r,1.0)-g) < 1e-13
    assert 1.15 < g < 1.17
    assert -0.53 < b['matched_H_V'] < -0.50


from analysis_synthetic_3d_basins import run_all as run_synth3d

def test_synthetic_3d_basin_stress_test_is_not_false_positive():
    rows, checks = run_synth3d()
    assert max(abs(x-0.811) for x in checks) < 1e-12
    no_persist = [s['volume_fraction_finite_size_active'] for seed,p,s8,s in rows if p == 0.0]
    allv = [s['volume_fraction_finite_size_active'] for seed,p,s8,s in rows]
    # Conservative Gaussian single-scale basins should not spuriously yield a volume-dominant V phase.
    assert max(no_persist) < 0.04
    assert max(allv) < 0.04
    assert min(allv) > 0.0
