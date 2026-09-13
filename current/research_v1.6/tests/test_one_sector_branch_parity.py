from pathlib import Path
import json
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research_v1.6'/'results'/'branch_parity'

def test_summary_exists_and_ranges():
    df=pd.read_csv(OUT/'branch_odd_volume_audit.csv')
    assert len(df)==5
    for c in ['final_positive_fraction','persistent_all4_positive_fraction','history_integral_positive_fraction']:
        assert ((df[c]>0.45)&(df[c]<0.60)).all()

def test_old_active_nested_in_positive_branch_controls():
    df=pd.read_csv(OUT/'branch_odd_volume_audit.csv')
    assert np.allclose(df['old_active_inside_final_positive'],1.0)
    assert np.allclose(df['old_active_inside_persistent_positive'],1.0)

def test_velocity_gradient_quadratic_source_is_even():
    p=json.loads((OUT/'branch_parity_summary.json').read_text())
    assert p['velocity_gradient_reversal_Q_even_max_abs_error'] < 1e-12

def test_fluid_branch_eos_identities():
    p=json.loads((OUT/'branch_parity_summary.json').read_text())
    assert p['fluid_branch_identity_max_abs_error'] < 1e-14

def test_history_integral_tracks_exact_relative_dilation():
    df=pd.read_csv(OUT/'branch_odd_volume_audit.csv')
    assert (df['history_dilation_corr']>0.90).all()
