#!/usr/bin/env python3
"""Minimal tri-stable sextic phase model and equipotential geometry factor."""
from __future__ import annotations
import argparse, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import root_scalar


def W0(x): return (4/3)*x**6-3*x**4+2*x**2
def dW0(x): return 8*x**5-12*x**3+4*x
def ddW0(x): return 40*x**4-36*x**2+4

def stationary(H):
    roots=np.roots([8,0,-12,0,4,-H]); out=[]
    for r in roots:
        if abs(r.imag)<1e-9:
            x=float(r.real); out.append((x,W0(x)-H*x,ddW0(x)))
    return sorted(out)

def minima(H): return [(x,e) for x,e,c in stationary(H) if c>0]
def central_and_selected(H):
    ms=minima(H); central=min(ms,key=lambda p:abs(p[0])); selected=max(ms,key=lambda p:p[0]) if H>=0 else min(ms,key=lambda p:p[0]); return central,selected
def energy_diff(H):
    c,s=central_and_selected(H); return s[1]-c[1]
def phase_constants():
    Hcoex=root_scalar(energy_diff,bracket=[.30,.50],xtol=1e-14).root
    c,s=central_and_selected(Hcoex)
    roots=np.roots([40,-36,4]); spin=[]
    for y2 in roots:
        if y2>0:
            for x in [math.sqrt(float(y2)),-math.sqrt(float(y2))]: spin.append((x,dW0(x)))
    Hspin=max(h for x,h in spin if x>0 and h>0)
    return Hcoex,c[0],s[0],s[0]-c[0],Hspin

def mn_eta(R_eq_over_a,b_over_a=.1,a=1.):
    b=b_over_a*a; Req=R_eq_over_a*a; s=math.sqrt(Req**2+(a+b)**2); thmax=math.acos((a+b)/s)
    def raw(th):
        R=s*math.sin(th); u=s*math.cos(th)-a; z=math.sqrt(max(u*u-b*b,0.)); dR=s*math.cos(th); dz=-1e100 if z==0 else -(u*s*math.sin(th))/z
        grad=math.sqrt((R/s)**2+(((a+u)/s)*(z/u if u else 0.))**2); g=grad/s**2; dA=4*math.pi*R*math.sqrt(dR*dR+dz*dz); return g,dA
    A=quad(lambda t:raw(t)[1],0,thmax,points=[thmax],epsrel=2e-9,limit=300)[0]
    I1=quad(lambda t:raw(t)[0]*raw(t)[1],0,thmax,points=[thmax],epsrel=2e-9,limit=300)[0]
    I2=quad(lambda t:raw(t)[0]**2*raw(t)[1],0,thmax,points=[thmax],epsrel=2e-9,limit=300)[0]
    return I1*I1/(A*I2)

def main():
    ap=argparse.ArgumentParser(); root=Path(__file__).resolve().parents[1]
    ap.add_argument('--out',type=Path,default=root/'results'); ap.add_argument('--fig',type=Path,default=root/'figures'); args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True); args.fig.mkdir(parents=True,exist_ok=True)
    Hcoex,xc,xs,jump,Hspin=phase_constants()
    with (args.out/'phase_summary.txt').open('w',encoding='utf-8') as f:
        f.write(f'H_coexistence = {Hcoex:.12f}\nchi_parent_at_coexistence = {xc:.12f}\nchi_selected_at_coexistence = {xs:.12f}\njump = {jump:.12f}\ncentral_spinodal_H = {Hspin:.12f}\n')
        for b in [.1,.3,1.]:
            f.write(f'\nMiyamoto-Nagai b/a={b}\n')
            for r in [.2,.5,1,2,4,8]:
                eta=mn_eta(r,b); f.write(f'R/a={r:4.1f} eta={eta:.9f} sqrt_eta={math.sqrt(eta):.9f}\n')
    xsgrid=np.linspace(-1.6,1.6,1200); fig=plt.figure(figsize=(10,4.2)); ax=fig.add_subplot(1,2,1)
    for H in [-.55,0,.55]: ax.plot(xsgrid,W0(xsgrid)-H*xsgrid,label=fr'$H={H:+.2f}$')
    ax.set_xlabel(r'Regional order parameter $\chi$'); ax.set_ylabel(r'$W(\chi;H)$'); ax.set_ylim(-.45,1.2); ax.legend(frameon=False); ax.set_title('Minimal tri-stable regional potential')
    Hs=np.linspace(-1,1,801); xeq=[]
    for H in Hs: xeq.append(min(minima(float(H)),key=lambda p:p[1])[0])
    ax2=fig.add_subplot(1,2,2); ax2.plot(Hs,xeq); ax2.axvline(Hcoex,ls='--',lw=1); ax2.axvline(-Hcoex,ls='--',lw=1); ax2.set_xlabel('Dynamical bias H'); ax2.set_ylabel(r'Global-minimum $\chi_*$'); ax2.set_title('Discontinuous branch selection'); fig.tight_layout(); fig.savefig(args.fig/'fig_phase_transition.pdf',bbox_inches='tight'); plt.close(fig)
    Rs=np.geomspace(.12,10,70); fig=plt.figure(figsize=(5.5,4.2)); ax=fig.add_subplot(111)
    for b in [.1,.3,1.]: ax.plot(Rs,[math.sqrt(mn_eta(r,b)) for r in Rs],label=fr'$b/a={b:g}$')
    ax.set_xscale('log'); ax.set_ylim(.94,1.002); ax.set_xlabel(r'Equatorial leaf radius $R/a$'); ax.set_ylabel(r'$Q_G/\sqrt{M_bM_A}=\sqrt{\eta_S}$'); ax.set_title('Parameter-free geometric suppression'); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(args.fig/'fig_geometry_factor.pdf',bbox_inches='tight'); plt.close(fig)
    print((args.out/'phase_summary.txt').read_text())
if __name__=='__main__': main()
