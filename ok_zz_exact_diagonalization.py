"""
Exact-diagonalization computation of the transverse exchange coupling J and the
longitudinal ZZ interaction for two transmons sharing a cavity bus.

Purpose (wave-physics-cap2): replace the hardcoded, convention-ambiguous values
J = 6.5 MHz and |zeta_zz| = 1.7 MHz that previously permeated Chapter 2 with
numbers derived from first principles, and fix the factor-of-4 convention error
in Eq. (zeta_zz_correct).

Model (three-body Hamiltonian, rotating frame not needed for static spectrum):

    H = sum_i [ omega_i b_i^dag b_i + (alpha_i/2) b_i^dag b_i^dag b_i b_i ]
        + omega_r a^dag a
        + sum_i g_i (b_i^dag a + b_i a^dag)

Transmons are modelled as Duffing oscillators (3 levels each); the cavity is
truncated at Fock 6. Computational eigenstates |n1,n2,0_cav> are identified by
maximum overlap with the bare product states. From the four lowest computational
levels we extract the ZZ level shift

    xi_ZZ = E_11 + E_00 - E_01 - E_10        [standard "ZZ rate", Krantz 2019]

which is what is directly measurable. The sigma_z (x) sigma_z Hamiltonian
coefficient is zeta_zz = xi_ZZ / 4. We report BOTH, explicitly, to close the
convention ambiguity.

Author: thesis revision, wave-physics-cap2.
"""

import numpy as np
import qutip as qt

# ----------------------------------------------------------------------
# Hilbert-space truncation
# ----------------------------------------------------------------------
NT = 3   # transmon levels (0,1,2): |2> is required for the anharmonic ZZ path
NC = 6   # cavity Fock levels


def build_operators():
    """Annihilation operators for (transmon1, transmon2, cavity)."""
    b1 = qt.tensor(qt.destroy(NT), qt.qeye(NT), qt.qeye(NC))
    b2 = qt.tensor(qt.qeye(NT), qt.destroy(NT), qt.qeye(NC))
    a = qt.tensor(qt.qeye(NT), qt.qeye(NT), qt.destroy(NC))
    return b1, b2, a


def hamiltonian(omega_q1, omega_q2, omega_r, g1, g2, alpha1, alpha2):
    """Full three-body Hamiltonian, all arguments in GHz (2*pi factored out:
    energies are frequencies, so results are frequencies too)."""
    b1, b2, a = build_operators()
    H = (omega_q1 * b1.dag() * b1
         + 0.5 * alpha1 * b1.dag() * b1.dag() * b1 * b1
         + omega_q2 * b2.dag() * b2
         + 0.5 * alpha2 * b2.dag() * b2.dag() * b2 * b2
         + omega_r * a.dag() * a
         + g1 * (b1.dag() * a + b1 * a.dag())
         + g2 * (b2.dag() * a + b2 * a.dag()))
    return H


def bare_index(n1, n2, nc):
    """Flat index of bare state |n1,n2,nc> in the tensor ordering."""
    return (n1 * NT + n2) * NC + nc


def dressed_energy(evals, evecs, n1, n2, nc=0):
    """Energy of the dressed eigenstate adiabatically connected to |n1,n2,nc>,
    identified by maximum overlap with the bare state."""
    target = bare_index(n1, n2, nc)
    overlaps = [abs(vec.full().flatten()[target])**2 for vec in evecs]
    k = int(np.argmax(overlaps))
    return float(np.real(evals[k])), overlaps[k]


def compute_zz(omega_q1, omega_q2, omega_r, g1, g2, alpha1, alpha2, label=""):
    H = hamiltonian(omega_q1, omega_q2, omega_r, g1, g2, alpha1, alpha2)
    evals, evecs = H.eigenstates()

    E00, o00 = dressed_energy(evals, evecs, 0, 0)
    E01, o01 = dressed_energy(evals, evecs, 0, 1)
    E10, o10 = dressed_energy(evals, evecs, 1, 0)
    E11, o11 = dressed_energy(evals, evecs, 1, 1)

    # ZZ level shift (standard "ZZ rate"), in GHz -> convert to kHz
    xi_ZZ = (E11 + E00 - E01 - E10)
    zeta_sigma = xi_ZZ / 4.0  # sigma_z (x) sigma_z Hamiltonian coefficient

    # Transverse exchange from the dispersive expression (reference value)
    D1 = omega_q1 - omega_r
    D2 = omega_q2 - omega_r
    J_disp = 0.5 * g1 * g2 * (1.0 / D1 + 1.0 / D2)

    min_overlap = min(o00, o01, o10, o11)

    print(f"\n=== {label} ===")
    print(f"  omega_q1={omega_q1:.4f}  omega_q2={omega_q2:.4f}  omega_r={omega_r:.4f} GHz")
    print(f"  g1={g1*1e3:.1f}  g2={g2*1e3:.1f} MHz   alpha1={alpha1*1e3:.1f}  alpha2={alpha2*1e3:.1f} MHz")
    print(f"  Delta_q = {(omega_q1-omega_q2)*1e3:.1f} MHz   "
          f"Delta_1={D1*1e3:.1f}  Delta_2={D2*1e3:.1f} MHz")
    print(f"  min dressed-state overlap = {min_overlap:.4f}  "
          f"({'OK' if min_overlap > 0.7 else 'WEAK - check truncation/regime'})")
    print(f"  J (dispersive)        = {J_disp*1e3:+.3f} MHz")
    print(f"  xi_ZZ (level shift)   = {xi_ZZ*1e6:+.2f} kHz   [E11+E00-E01-E10]")
    print(f"  zeta_zz (sigma_z conv)= {zeta_sigma*1e6:+.2f} kHz   [= xi_ZZ / 4]")

    return dict(label=label, J_MHz=J_disp*1e3, xi_ZZ_kHz=xi_ZZ*1e6,
                zeta_sigma_kHz=zeta_sigma*1e6, Delta_q_MHz=(omega_q1-omega_q2)*1e3,
                min_overlap=min_overlap)


def perturbative_xi(J_GHz, Delta_q, alpha1, alpha2):
    """Thesis Eq. (zeta_zz_correct), interpreted as the LEVEL SHIFT xi_ZZ:
        xi = 2 J^2 [ 1/(Delta_q - alpha2) - 1/(Delta_q + alpha1) ]
    Returns kHz."""
    val = 2 * J_GHz**2 * (1.0/(Delta_q - alpha2) - 1.0/(Delta_q + alpha1))
    return val * 1e6


if __name__ == "__main__":
    print("="*70)
    print("Exact diagonalization of the two-transmon + cavity system")
    print("QuTiP", qt.__version__)
    print("="*70)

    # ---- Operating point 1: BASELINE (mK validation device) -------------
    base = compute_zz(omega_q1=5.8, omega_q2=5.2, omega_r=6.5,
                      g1=0.080, g2=0.080, alpha1=-0.200, alpha2=-0.200,
                      label="Baseline (mK validation)")
    xi_pert_base = perturbative_xi(base["J_MHz"]/1e3, 0.600, -0.200, -0.200)
    print(f"  perturbative xi (Eq.) = {xi_pert_base:+.2f} kHz   "
          f"(exact/pert = {base['xi_ZZ_kHz']/xi_pert_base:.3f})")

    # ---- Operating point 2: 300 GHz innovation (alpha = -15 GHz) --------
    # Decision 1: E_J/E_C = 50 -> alpha = -15 GHz at 300 GHz.
    # Decision 3: Delta_q = 600 MHz.  Cavity placed 1 GHz above qubits.
    innov = compute_zz(omega_q1=300.3, omega_q2=299.7, omega_r=301.3,
                       g1=0.080, g2=0.080, alpha1=-15.0, alpha2=-15.0,
                       label="Innovation 300 GHz (alpha=-15 GHz)")
    xi_pert_innov = perturbative_xi(innov["J_MHz"]/1e3, 0.600, -15.0, -15.0)
    print(f"  perturbative xi (Eq.) = {xi_pert_innov:+.2f} kHz   "
          f"(exact/pert = {innov['xi_ZZ_kHz']/xi_pert_innov:.3f})")

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for r in (base, innov):
        print(f"{r['label']:<34}  J={r['J_MHz']:+.2f} MHz  "
              f"xi_ZZ={r['xi_ZZ_kHz']:+.1f} kHz  "
              f"|J/xi|={abs(r['J_MHz']*1e3/r['xi_ZZ_kHz']):.1f}")
