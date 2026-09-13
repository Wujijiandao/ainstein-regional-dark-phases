#!/usr/bin/env python3
"""Conditional sparse-V energy contrast bookkeeping audit.

If a V phase carries a fixed positive intrinsic energy density epsilon_V only on
an active volume fraction phi_V and the rest of the volume carries no V energy,
then
    <rho_V> = phi_V epsilon_V,
    sigma(rho_V)/<rho_V> = sqrt((1-phi_V)/phi_V).
Matching a target mean Omega_V therefore requires
    epsilon_V/rho_crit = Omega_V/phi_V.
This is an exact two-state bookkeeping identity, not a relativistic observable
prediction. It quantifies the price of accepting a very sparse V phase rather
than changing the selector.
"""
from pathlib import Path
import argparse,json,math
import pandas as pd

OMEGA_M=0.315
OMEGA_V_TARGET=1.0-OMEGA_M

def metrics(phi):
    return OMEGA_V_TARGET/phi, math.sqrt((1.0-phi)/phi)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--ceiling-table',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(args.ceiling_table)
    rows=[]
    for _,r in d.iterrows():
        for branch,key in [('current_finite','current_finite_size_phi'),('zero_capillary_ceiling','gamma_ge1_ceiling')]:
            phi=float(r[key]); eps,cv=metrics(phi)
            rows.append({'seed':int(r.seed),'branch':branch,'phi_V':phi,'Omega_V_target':OMEGA_V_TARGET,
                         'required_epsilonV_over_rhocrit':eps,'binary_relative_rms':cv})
    out=pd.DataFrame(rows);out.to_csv(args.out/'sparse_v_energy_contrast.csv',index=False)
    summary={'Omega_V_target':OMEGA_V_TARGET,
             'identity_mean':'<rho_V>=phi_V epsilon_V',
             'identity_relative_rms':'sigma/<rho>=sqrt((1-phi)/phi)',
             'current_required_epsilon_range':[float(out[out.branch=='current_finite'].required_epsilonV_over_rhocrit.min()),float(out[out.branch=='current_finite'].required_epsilonV_over_rhocrit.max())],
             'current_relative_rms_range':[float(out[out.branch=='current_finite'].binary_relative_rms.min()),float(out[out.branch=='current_finite'].binary_relative_rms.max())],
             'ceiling_required_epsilon_range':[float(out[out.branch=='zero_capillary_ceiling'].required_epsilonV_over_rhocrit.min()),float(out[out.branch=='zero_capillary_ceiling'].required_epsilonV_over_rhocrit.max())],
             'ceiling_relative_rms_range':[float(out[out.branch=='zero_capillary_ceiling'].binary_relative_rms.min()),float(out[out.branch=='zero_capillary_ceiling'].binary_relative_rms.max())],
             'scope':'conditional two-state bookkeeping only; not a lensing/CMB likelihood or relativistic perturbation calculation'}
    (args.out/'sparse_v_energy_contrast.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (args.out/'sparse_v_energy_contrast.txt').write_text('\n'.join(f'{k}={v}' for k,v in summary.items())+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
