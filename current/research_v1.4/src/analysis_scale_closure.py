#!/usr/bin/env python3
from pathlib import Path
import argparse, json, math
G=6.67430e-11; c=299792458.0; MPC=3.085677581491367e22

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--aE',type=float,default=1.097458984856508e-10); ap.add_argument('--H0',type=float,default=67.4); ap.add_argument('--OmegaV',type=float,default=.685); args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    H=args.H0*1000/MPC
    epsV=args.OmegaV*3*H*H*c*c/(8*math.pi*G)
    scale=args.aE*args.aE/G
    rec={'aE_m_s2':args.aE,'H0_km_s_Mpc':args.H0,'OmegaV_reference':args.OmegaV,
         'aE_over_cH0':args.aE/(c*H),'one_over_2pi':1/(2*math.pi),
         'epsilonV_reference_J_m3':epsV,'aE_squared_over_G_J_m3':scale,
         'epsilonV_over_aE2_over_G':epsV/scale,
         'interpretation':'If a single length L_* controls both a_E=alpha c^2/L_* and epsilon_V=beta c^4/(G L_*^2), then beta/alpha^2 = epsilon_V G/a_E^2. The numerical O(1) value is a dimensional closure target only, not a derivation; the a0-cosmology coincidence is known in MOND literature.'}
    (args.out/'scale_closure.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
    print(json.dumps(rec,indent=2))
if __name__=='__main__':main()
