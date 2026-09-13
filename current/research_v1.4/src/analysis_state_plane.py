from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parents[1] / 'results' / 'state_plane'
out.mkdir(parents=True, exist_ok=True)

gamma = np.logspace(-1.3, 1.0, 500)
nu_acc = 1.0 + 1.0/gamma

fig, ax = plt.subplots(figsize=(7.2, 4.8))
ax.plot(gamma, np.ones_like(gamma), label=r'curvature-like geometry: $\nu=1$')
ax.plot(gamma, 3*np.ones_like(gamma), linestyle='--', label=r'vacuum-like geometry benchmark: $\nu=3$')
ax.plot(gamma, nu_acc, linestyle='-.', label=r'geometry-only acceleration: $\Gamma(\nu-1)=1$')
ax.axvline(1.0, linestyle=':', label=r'matter--geometry amplitude equality: $\Gamma=1$')
ax.set_xscale('log')
ax.set_xlim(gamma.min(), gamma.max())
ax.set_ylim(0, 4.2)
ax.set_xlabel(r'geometry-to-matter amplitude $\Gamma=\rho_G/\rho_m$')
ax.set_ylabel(r'history exponent $\nu=d\ln\Gamma/d\ln a_\Omega$')
ax.set_title('Regional geometry control-state plane')
ax.legend(fontsize=8, loc='upper right')
fig.tight_layout()
fig.savefig(out/'fig_regional_geometry_state_plane.pdf', bbox_inches='tight')
fig.savefig(out/'fig_regional_geometry_state_plane.png', dpi=180, bbox_inches='tight')
print(out)
