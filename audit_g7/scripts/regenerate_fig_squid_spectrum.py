"""
Regenerate Fig 2.12 (asymmetric_squid_spectrum.png) — the asymmetric-SQUID
transmon spectrum vs external flux. The panel (d) title is updated to
remove the generic 'Eq. 1' reference and replace it with a self-contained
formula label.

Numerical parameters (from Tab. 2.4 'squid_design_params' of the thesis):
  E_Sigma = E_J1 + E_J2 = 27 GHz
  E_Delta = E_J1 - E_J2 = 3 GHz  (asymmetry d = E_Delta/E_Sigma = 11.1%)
  E_C     = 200 MHz
"""
import numpy as np
from scipy.sparse import diags
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.dpi': 120,
    'savefig.bbox': 'tight',
})

# Parameters (in GHz)
E_Sigma = 27.0
E_Delta = 3.0
E_C = 0.200

# Effective Josephson energy E_J(Phi)
def E_J(phi_over_phi0):
    """Effective E_J(Phi) for asymmetric SQUID. Phi/Phi_0 in units of one."""
    return np.sqrt(
        E_Sigma**2 * np.cos(np.pi * phi_over_phi0)**2 +
        E_Delta**2 * np.sin(np.pi * phi_over_phi0)**2
    )

# Transmon Hamiltonian in charge basis: H = 4 E_C (n - n_g)^2 - E_J cos(phi)
# Diagonalize in charge basis truncated to |n|<=N
def transmon_levels(EJ, EC=E_C, ng=0.0, N=20, n_levels=5):
    """Return the first n_levels eigenvalues of the transmon Hamiltonian.
    All energies in GHz units.
    """
    dim = 2*N + 1
    n_vals = np.arange(-N, N+1)
    # Diagonal: 4 E_C (n - n_g)^2
    H_diag = 4*EC * (n_vals - ng)**2
    # Off-diagonal cos(phi) couples n -> n+1 and n -> n-1 with -E_J/2
    off = -EJ/2 * np.ones(dim-1)
    H = diags([off, H_diag, off], [-1, 0, 1], format='csr')
    eigs = np.sort(np.real(spla.eigsh(H, k=n_levels, which='SA', tol=1e-10)[0]))
    return eigs

# Compute spectrum vs flux
phi_grid = np.linspace(0, 1, 121)
n_levels = 5
levels = np.zeros((len(phi_grid), n_levels))
for i, p in enumerate(phi_grid):
    EJ = E_J(p)
    levels[i, :] = transmon_levels(EJ, n_levels=n_levels)

# Anchor to E0 = 0
E0 = levels[:, 0:1]
levels_shifted = levels - E0

# Transition frequencies
omega_01 = levels_shifted[:, 1]
omega_12 = levels_shifted[:, 2] - levels_shifted[:, 1]
alpha = omega_12 - omega_01  # GHz, in MHz below

# E_J(Phi) curve
EJ_curve = np.array([E_J(p) for p in phi_grid])

# ─── Build the 2x2 figure ───
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# Panel (a): Energy spectrum
ax = axes[0, 0]
colors = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD']
for i in range(n_levels):
    ax.plot(phi_grid, levels_shifted[:, i],
            color=colors[i], linewidth=2, label=f'$E_{i}$')
ax.axvline(0.5, color='red', linestyle='--', linewidth=1.5,
           alpha=0.7, label='Sweet spot')
ax.set_xlabel(r'Flux $\Phi/\Phi_0$')
ax.set_ylabel('Energy (GHz)')
ax.set_title('Energy spectrum vs flux', fontweight='bold')
ax.legend(loc='upper center', ncol=3, fontsize=10, framealpha=0.95)
ax.grid(alpha=0.3, linestyle=':')

# Panel (b): Transition frequencies
ax = axes[0, 1]
ax.plot(phi_grid, omega_01, color='#1F77B4', linewidth=2.2,
        label=r'$\omega_{01}$')
ax.plot(phi_grid, omega_12, color='#FF7F0E', linewidth=2.2,
        label=r'$\omega_{12}$')
ax.axvline(0.5, color='red', linestyle='--', linewidth=1.5,
           alpha=0.7, label='Sweet spot')
ax.set_xlabel(r'Flux $\Phi/\Phi_0$')
ax.set_ylabel('Transition frequency (GHz)')
ax.set_title('Transition frequencies vs flux', fontweight='bold')
ax.legend(loc='upper center', fontsize=10, framealpha=0.95)
ax.grid(alpha=0.3, linestyle=':')

# Panel (c): Anharmonicity
ax = axes[1, 0]
ax.plot(phi_grid, alpha * 1000, color='#7E1F9E', linewidth=2.2)
ax.axhline(-E_C * 1000, color='gray', linestyle=':', linewidth=1.5,
           label=r'$-E_C$ (approx.)')
ax.axvline(0.5, color='red', linestyle='--', linewidth=1.5,
           alpha=0.7, label='Sweet spot')
ax.set_xlabel(r'Flux $\Phi/\Phi_0$')
ax.set_ylabel(r'Anharmonicity $\alpha$ (MHz)')
ax.set_title('Anharmonicity vs flux', fontweight='bold')
ax.legend(loc='lower center', fontsize=10, framealpha=0.95)
ax.grid(alpha=0.3, linestyle=':')

# Panel (d): E_J(Phi)
ax = axes[1, 1]
ax.plot(phi_grid, EJ_curve, color='#1A8B3B', linewidth=2.2,
        label=r'$E_J(\Phi)$')
ax.axhline(E_Sigma, color='#1F77B4', linestyle=':', linewidth=1.5,
           label=r'$E_\Sigma = E_{J1} + E_{J2}$')
ax.axhline(abs(E_Delta), color='#FF7F0E', linestyle=':', linewidth=1.5,
           label=r'$|E_\Delta| = |E_{J1} - E_{J2}|$')
ax.axvline(0.5, color='red', linestyle='--', linewidth=1.5,
           alpha=0.7, label='Sweet spot')
ax.set_xlabel(r'Flux $\Phi/\Phi_0$')
ax.set_ylabel(r'$E_J(\Phi)$ (GHz)')
ax.set_title(r'Effective Josephson energy $E_J(\Phi)$ '
             r'$= \sqrt{E_\Sigma^2 \cos^2(\pi\Phi/\Phi_0) + E_\Delta^2 \sin^2(\pi\Phi/\Phi_0)}$',
             fontsize=10, fontweight='bold')
ax.legend(loc='center right', fontsize=9, framealpha=0.95)
ax.grid(alpha=0.3, linestyle=':')
ax.set_ylim(0, E_Sigma * 1.05)

plt.tight_layout()
plt.savefig('./asymmetric_squid_spectrum.pdf')
plt.savefig('./asymmetric_squid_spectrum.png', dpi=160)
print("Saved ./asymmetric_squid_spectrum.{pdf,png}")

# Sanity prints
print(f"\nAt sweet spot Phi=Phi_0/2 (index {np.argmin(np.abs(phi_grid-0.5))}):")
i_sw = np.argmin(np.abs(phi_grid - 0.5))
print(f"  E_J(Phi_0/2) = {EJ_curve[i_sw]:.3f} GHz")
print(f"  omega_01 = {omega_01[i_sw]:.4f} GHz")
print(f"  omega_12 = {omega_12[i_sw]:.4f} GHz")
print(f"  alpha    = {alpha[i_sw]*1000:.1f} MHz")
