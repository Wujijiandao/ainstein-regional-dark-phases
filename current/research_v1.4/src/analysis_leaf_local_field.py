#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Units: G M = a = 1 for the Miyamoto-Nagai toy family.
def mn_leaf_moments(R_eq_over_a:float,b_over_a:float=0.1,a:float=1.0):
    b=b_over_a*a; Req=R_eq_over_a*a
    s=math.sqrt(Req*Req+(a+b)**2)
    thmax=math.acos((a+b)/s)
    def raw(th):
        R=s*math.sin(th); u=s*math.cos(th)-a
        z=math.sqrt(max(u*u-b*b,0.0))
        dR=s*math.cos(th)
        dz=-1e100 if z==0 else -(u*s*math.sin(th))/z
        grad=math.sqrt((R/s)**2+(((a+u)/s)*(z/u if u else 0.0))**2)
        gb=grad/s**2
        dA=4*math.pi*R*math.sqrt(dR*dR+dz*dz)
        return gb,dA
    A=quad(lambda t:raw(t)[1],0,thmax,points=[thmax],epsrel=2e-9,limit=300)[0]
    I1=quad(lambda t:raw(t)[0]*raw(t)[1],0,thmax,points=[thmax],epsrel=2e-9,limit=300)[0]
    I2=quad(lambda t:raw(t)[0]**2*raw(t)[1],0,thmax,points=[thmax],epsrel=2e-9,limit=300)[0]
    g_eq=Req/s**3
    eta=I1*I1/(A*I2)
    g_leaf=I2/I1
    return {'R_eq_over_a':R_eq_over_a,'b_over_a':b_over_a,'A':A,'I1':I1,'I2':I2,
            'eta_S':eta,'sqrt_eta_S':math.sqrt(eta),'g_b_eq':g_eq,'g_leaf':g_leaf,
            'g_leaf_over_g_b_eq':g_leaf/g_eq}

def q_simple(x):
    return 0.5*(1+np.sqrt(1+4*x))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True); args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    Rs=np.geomspace(.12,10,90); bs=[.1,.3,1.0]
    rows=[]
    for b in bs:
        for R in Rs: rows.append(mn_leaf_moments(float(R),b))
    df=pd.DataFrame(rows)
    # Illustrative fixed galaxy compactness: a_E a^2/(GM)=0.1. This is not fit and is not used in theorem.
    aE_dimless=0.1
    df['q_leaf']=q_simple(aE_dimless/df.g_leaf)
    df['q_pointwise_equator']=q_simple(aE_dimless/df.g_b_eq)
    df['q_leaf_over_pointwise_equator']=df.q_leaf/df.q_pointwise_equator
    df.to_csv(args.out/'leaf_local_field_miyamoto_nagai.csv',index=False)
    # A small benchmark table at fixed radii.
    bench=df.iloc[0:0].copy()
    b_rows=[]
    for b in bs:
        for R in [.2,.5,1,2,4,8]:
            d=mn_leaf_moments(R,b); d['q_leaf']=float(q_simple(aE_dimless/d['g_leaf'])); d['q_pointwise_equator']=float(q_simple(aE_dimless/d['g_b_eq'])); d['q_leaf_over_pointwise_equator']=d['q_leaf']/d['q_pointwise_equator']; b_rows.append(d)
    pd.DataFrame(b_rows).to_csv(args.out/'leaf_local_field_benchmark.csv',index=False)
    rec={'derived_local_law':{
        'g_leaf':'I2/I1 = <g_b^2>_S/<g_b>_S',
        'alpha_equation':'alpha(1+alpha)=a_E/g_leaf',
        'q_leaf':'1+alpha = [1+sqrt(1+4 a_E/g_leaf)]/2',
        'pointwise_field':'g_tot(x)=q_leaf(S_phi) g_b(x) for x on S_phi',
        'direction':'g_tot is parallel to g_b under leaf locking',
        'spherical_limit':'g_leaf=g_b, recovering the historical simple-MOND algebraic form'},
        'figure_benchmark_aE_a2_over_GM':aE_dimless,
        'status':'Conditional corollary of the finite-amplitude leaf-locked functional; lensing still requires a relativistic completion.'}
    (args.out/'leaf_local_field_identity.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
    fig,ax=plt.subplots(figsize=(6.6,4.4))
    for b in bs:
        d=df[df.b_over_a==b]
        ax.plot(d.R_eq_over_a,d.g_leaf_over_g_b_eq,label=fr'$b/a={b:g}$')
    ax.axhline(1,linestyle='--',linewidth=1)
    ax.set_xscale('log'); ax.set_xlabel(r'equatorial leaf radius $R/a$'); ax.set_ylabel(r'$g_{\rm leaf}/g_{b,\rm eq}$')
    ax.set_title('Leaf-global acceleration moment differs from local equatorial field')
    ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(args.out/'fig_leaf_global_acceleration.pdf',bbox_inches='tight'); fig.savefig(args.out/'fig_leaf_global_acceleration.png',dpi=200,bbox_inches='tight'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.6,4.4))
    for b in bs:
        d=df[df.b_over_a==b]
        ax.plot(d.R_eq_over_a,d.q_leaf_over_pointwise_equator,label=fr'$b/a={b:g}$')
    ax.axhline(1,linestyle='--',linewidth=1)
    ax.set_xscale('log'); ax.set_xlabel(r'equatorial leaf radius $R/a$'); ax.set_ylabel(r'$q_{\rm leaf}/q_{\rm point,eq}$')
    ax.set_title(r'Illustrative local-field distinction, $a_Ea^2/(GM)=0.1$')
    ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(args.out/'fig_leaf_vs_pointwise_amplification.pdf',bbox_inches='tight'); fig.savefig(args.out/'fig_leaf_vs_pointwise_amplification.png',dpi=200,bbox_inches='tight'); plt.close(fig)
    print(pd.DataFrame(b_rows)[['b_over_a','R_eq_over_a','eta_S','g_leaf_over_g_b_eq','q_leaf_over_pointwise_equator']].to_string(index=False))

if __name__=='__main__': main()
