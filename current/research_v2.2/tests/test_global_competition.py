from pathlib import Path
import importlib.util

ROOT = Path(__file__).resolve().parents[2]
MOD = ROOT/'research_v2.2'/'src'/'analysis_global_competition.py'
spec = importlib.util.spec_from_file_location('gc', MOD)
gc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gc)


def test_gaussian_elimination_identity_and_sign():
    assert gc.identity_max_error < 1e-10
    assert gc.min_G_eigenvalue > -1e-10
    assert gc.max_induced_kernel_eigenvalue < 1e-10
    assert gc.rank1_induced_max_eigenvalue < 1e-10


def test_odd_budget_does_not_fix_total_conversion():
    assert gc.f_total_max_from_m1 - gc.f_total_min_from_m1 > 0.5
    assert gc.simplex_ok
    assert gc.odd_constraint_error < 1e-12


def test_even_budget_fixes_total_conversion_only_if_imposed():
    assert abs(gc.f_total_from_m2 - 0.40) < 1e-12
    assert abs(gc.fD_from_m1m2 + gc.fC_from_m1m2 + gc.fV_from_m1m2 - 1.0) < 1e-12


def test_passive_global_gaussian_is_cooperative_in_frozen_audit():
    assert gc.passive_monotone_all
    assert (gc.summary_by_seed.passive_final_fraction >= gc.summary_by_seed.baseline_fraction).all()


def test_positive_competition_is_distinct_sign():
    assert (gc.summary_by_seed.repulsive_final_fraction <= gc.summary_by_seed.baseline_fraction).all()
