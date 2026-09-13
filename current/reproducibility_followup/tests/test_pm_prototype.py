import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from analysis_nonlinear_pm_basin_lineage import (
    BOX,fourier_grid,grid_positions,density_ratio,force_tilde,
    gaussian_linear_field,displacement_from_delta
)

def test_uniform_lattice_mass_closure():
    n=12; q=grid_positions(n,BOX); rho=density_ratio(q,n,BOX)
    assert np.max(np.abs(rho-1.0)) < 1e-12

def test_zero_density_force_zero():
    n=12; fg=fourier_grid(n,BOX); d=np.zeros((n,n,n))
    f=force_tilde(d,fg)
    assert max(np.max(np.abs(x)) for x in f) < 1e-14

def test_displacement_divergence_identity():
    n=16; fg=fourier_grid(n,BOX); d=gaussian_linear_field(11,fg)
    s=displacement_from_delta(d,fg).reshape(n,n,n,3)
    divF=(1j*fg.kx*np.fft.rfftn(s[...,0])+1j*fg.ky*np.fft.rfftn(s[...,1])+1j*fg.kz*np.fft.rfftn(s[...,2]))
    div=np.fft.irfftn(divF,s=(n,n,n),axes=(0,1,2)).real
    # Nyquist modes cannot be represented by the real-valued derivative in the same way;
    # compare after low-pass filtering where the PM force is operationally resolved.
    filt=np.exp(-0.5*(fg.kmag*(BOX/n))**2)
    lhs=np.fft.irfftn(np.fft.rfftn(div)*filt,s=(n,n,n),axes=(0,1,2)).real
    rhs=np.fft.irfftn(np.fft.rfftn(-d)*filt,s=(n,n,n),axes=(0,1,2)).real
    assert np.sqrt(np.mean((lhs-rhs)**2)) < 2e-3
