"""
qp_kinetics_squid.py
====================

Bilancio cinetico delle quasiparticelle (QP) e impatto su T1 per un transmon
basato su SQUID asimmetrico con elettrodi NbN, operativo a 4 K.

Modello cinetico (stato stazionario):
    dn_qp/dt = G_ph(T) + G_gamma - r n_qp^2 - s n_qp = 0

dove (in unita' di x = n_qp / n_cp, con n_cp = 2 nu0 Delta):
- G_ph(T) generazione termica calibrata da detailed balance per riprodurre
  x_qp_thermal = sqrt(2 pi kBT/Delta) exp(-Delta/kBT)  [Glazman-Catelani 2021]
- G_gamma generazione non-termica (fotoni pair-breaking, raggi cosmici...)
- r       rate di ricombinazione (~1/tau_r0, scaling Kaplan per NbN)
- s       rate di trap (normal-metal trap, vortex trap, gap-asymmetric extraction)

T1 da QP per SQUID asimmetrico [Catelani et al. PRB 84, 064517 (2011)]:
    1/T1_qp = sum_j (omega01/pi) x_qp sqrt(2 Delta_j / hbar omega01)
              * (E_Jj/E_J) cos^2(phi_j/2)
              * exp(-max(0, |Delta1-Delta2| - hbar omega01)/kBT)

Ultimo fattore: soppressione gap-asymmetry [Marchegiani-Amico-Catelani,
PRX Quantum 3, 040338 (2022)]: se |Delta1-Delta2| > hbar omega01, le QP
"fredde" non hanno energia per attraversare la giunzione.

Riferimenti chiave (da bib di tesi se presenti, altrimenti da aggiungere):
- Catelani, Schoelkopf, Devoret, Glazman, PRB 84, 064517 (2011)
- Riwar, Hosseinkhani, Burkhart, Gao, Schoelkopf, Glazman, Catelani,
  PRB 94, 104516 (2016) - normal metal traps
- Glazman & Catelani, SciPost Phys. Lect. Notes 31 (2021) - review
- Marchegiani, Amico, Catelani, PRX Quantum 3, 040338 (2022) - gap asymmetry
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Costanti fisiche (SI)
# ============================================================
hbar = 1.054571817e-34   # J s
kB   = 1.380649e-23      # J/K
e    = 1.602176634e-19   # C
h    = 2 * np.pi * hbar  # J s

# Conversioni utili
meV_J = 1e-3 * e         # 1 meV in Joule
GHz_h = 1e9 * h          # 1 GHz in Joule (energia)

# ============================================================
# 1. Densita' di QP termiche - equilibrium BCS (kBT << Delta)
# ============================================================
def x_qp_thermal(T, Delta):
    """
    Frazione termica di coppie di Cooper rotte (Bogoliubov QP normalizzate).

    Parametri
    ---------
    T     : temperatura del bagno [K]
    Delta : gap superconduttivo [J]

    Ritorna
    -------
    x_qp = n_qp / (2 nu0 Delta) = sqrt(2 pi kBT/Delta) exp(-Delta/kBT)
    """
    arg = Delta / (kB * T)
    return np.sqrt(2 * np.pi / arg) * np.exp(-arg)


# ============================================================
# 2. Bilancio cinetico stazionario
# ============================================================
def x_qp_steady(T, Delta, G_gamma=0.0, s=0.0, r=1.0/100e-9):
    """
    Risolve in stato stazionario: r x^2 + s x - (G_ph + G_gamma) = 0

    Calibrazione: in assenza di G_gamma e s, x_ss -> x_qp_thermal.
    Quindi G_ph e' scelto come G_ph = r * x_th^2 (detailed balance termico).

    Parametri
    ---------
    T       : temperatura [K]
    Delta   : gap [J]
    G_gamma : rate generazione non-termica [unita' di x/s]
    s       : rate di trap [1/s]
    r       : rate di ricombinazione [1/s] (default 1/100 ns ~ NbN; raffinabile)

    Ritorna
    -------
    x_qp stazionario (adimensionale)
    """
    x_th = x_qp_thermal(T, Delta)
    G_ph = r * x_th**2

    a, b, c = r, s, -(G_ph + G_gamma)
    disc = b**2 - 4*a*c
    return (-b + np.sqrt(disc)) / (2*a)


# ============================================================
# 3. SQUID asimmetrico: E_J effettiva e fasi di equilibrio
# ============================================================
def EJ_squid(phi, EJ_sum, d):
    """
    Energia Josephson effettiva dello SQUID asimmetrico.
    phi = pi Phi/Phi_0  (flusso ridotto in radianti)
    d   = (E_J1 - E_J2) / (E_J1 + E_J2)  asimmetria, |d| <= 1
    """
    return EJ_sum * np.sqrt(np.cos(phi)**2 + (d * np.sin(phi))**2)


def phase_split(phi, d):
    """
    Fasi di equilibrio nelle due giunzioni dello SQUID asimmetrico.
    Soddisfano: E_J1 sin(phi1) = E_J2 sin(phi2),  phi1 - phi2 = 2 phi.

    Forma chiusa:
        tan(phi1) = (1-d) sin(2 phi) / [(1+d) - (1-d) cos(2 phi)]
        phi2 = phi1 - 2 phi
    """
    num = (1 - d) * np.sin(2*phi)
    den = (1 + d) - (1 - d) * np.cos(2*phi)
    phi1 = np.arctan2(num, den)
    phi2 = phi1 - 2*phi
    return phi1, phi2


# ============================================================
# 4. Frequenza qubit e T1 da QP per SQUID asimmetrico
# ============================================================
def omega_01(phi, EJ_sum, EC, d):
    """Frequenza 0->1 del transmon: omega_01 ~ (sqrt(8 EJ EC) - EC)/hbar."""
    EJ = EJ_squid(phi, EJ_sum, d)
    return (np.sqrt(8 * EJ * EC) - EC) / hbar


def T1_qp_squid(phi, x_qp, EJ_sum, EC, d, Delta1, Delta2=None, T=4.0):
    """
    T1 da QP per SQUID asimmetrico con possibile gap-asymmetry.

    Parametri
    ---------
    phi    : flusso ridotto pi Phi/Phi_0 [rad]
    x_qp   : densita' di QP (frazione di coppie rotte)
    EJ_sum : E_J1 + E_J2 [J]
    EC     : energia di charging [J]
    d      : asimmetria EJ (geometrica)
    Delta1, Delta2 : gap delle due giunzioni [J]. Se Delta2 None, gap simmetrico.
    T      : temperatura per soppressione gap-asymmetry [K]
    """
    if Delta2 is None:
        Delta2 = Delta1

    EJ = EJ_squid(phi, EJ_sum, d)
    w01 = omega_01(phi, EJ_sum, EC, d)
    phi1, phi2 = phase_split(phi, d)

    EJ1 = 0.5 * (1 + d) * EJ_sum
    EJ2 = 0.5 * (1 - d) * EJ_sum

    # Soppressione gap-asymmetry (Marchegiani-Amico-Catelani 2022)
    gap_gap = abs(Delta1 - Delta2)
    suppress = np.exp(-max(0.0, gap_gap - hbar*w01) / (kB * T))

    G1 = (w01/np.pi) * x_qp * np.sqrt(2*Delta1/(hbar*w01)) \
         * (EJ1/EJ) * np.cos(phi1/2)**2
    G2 = (w01/np.pi) * x_qp * np.sqrt(2*Delta2/(hbar*w01)) \
         * (EJ2/EJ) * np.cos(phi2/2)**2

    return 1.0 / ((G1 + G2) * suppress)


# ============================================================
# 5. Sanity checks
# ============================================================
def run_sanity_checks():
    """Tre controlli prima di accettare i numeri principali."""
    print("=" * 60)
    print("SANITY CHECKS")
    print("=" * 60)

    Delta = 2.8 * meV_J

    # 1. x_qp termico a 4 K NbN -> ordine 1e-4
    x_th_4K = x_qp_thermal(4.0, Delta)
    print(f"[1] x_qp termico a 4 K, Delta=2.8 meV: {x_th_4K:.3e}")
    assert 1e-5 < x_th_4K < 1e-3, "x_qp termico fuori range atteso"

    # 2. x_qp termico a 20 mK -> < 1e-20 (esponenziale collassa)
    x_th_20mK = x_qp_thermal(0.02, Delta)
    print(f"[2] x_qp termico a 20 mK: {x_th_20mK:.3e}")
    assert x_th_20mK < 1e-20, "x_qp termico a mK non collassa"

    # 3. T1 simmetrico Phi=0 con d=0 -> ~25 ns con parametri NbN tipici
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h
    T1_check = T1_qp_squid(0, x_th_4K, EJ_sum, EC, 0, Delta)
    w01_GHz = omega_01(0, EJ_sum, EC, 0) / (2*np.pi*1e9)
    print(f"[3] omega_01/2pi a Phi=0, d=0: {w01_GHz:.2f} GHz")
    print(f"    T1^qp termico (d=0, gap simm.): {T1_check*1e9:.1f} ns")
    assert 5e-9 < T1_check < 200e-9, "T1 termico fuori range atteso"

    # 4. Coerenza dimensionale: confronto formula diretta
    omega = omega_01(0, EJ_sum, EC, 0)
    factor = np.sqrt(2*Delta/(hbar*omega))
    Gamma_ref = (omega/np.pi) * x_th_4K * factor
    T1_ref = 1.0/Gamma_ref
    print(f"[4] T1 formula diretta:           {T1_ref*1e9:.1f} ns  (ref)")
    rel_err = abs(T1_check - T1_ref)/T1_ref
    print(f"    Errore relativo SQUID vs ref: {rel_err:.1%}")
    assert rel_err < 0.05, "Discrepanza tra implementazione SQUID e formula base"

    print("OK - tutti i sanity check passati.\n")
    return True


# ============================================================
# 6. Plot principali
# ============================================================
def make_plots():
    # Parametri di riferimento (NbN ALDmon-like)
    Delta_NbN = 2.8 * meV_J
    T_op = 4.0
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # ---- Pannello A: x_qp(T) con vari livelli di trap ----
    ax = axes[0, 0]
    Ts = np.linspace(0.05, 6, 300)
    x_th = np.array([x_qp_thermal(T, Delta_NbN) for T in Ts])
    ax.semilogy(Ts, x_th, 'k-', lw=2.5, label='termico (no trap)')
    for s_val, ls in zip([1e4, 1e6, 1e8], ['--', '-.', ':']):
        x = np.array([x_qp_steady(T, Delta_NbN, G_gamma=0, s=s_val) for T in Ts])
        ax.semilogy(Ts, x, ls, lw=1.5, label=fr'$s={s_val:.0e}\,\mathrm{{s^{{-1}}}}$')
    ax.axvline(T_op, color='red', alpha=0.3, lw=2)
    ax.text(T_op+0.1, 1e-3, 'T = 4 K', color='red', alpha=0.7)
    ax.set_xlabel('T (K)')
    ax.set_ylabel(r'$x_{\rm qp}$')
    ax.set_title(r'A. $x_{\rm qp}(T)$ con trap di varie efficienze (NbN, $\Delta=2.8$ meV)')
    ax.set_ylim(1e-20, 1e-2)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # ---- Pannello B: omega_01(Phi) per vari d ----
    ax = axes[0, 1]
    phis = np.linspace(0, np.pi, 300)
    for d_val, ls in zip([0.0, 0.3, 0.6, 0.9], ['-', '--', '-.', ':']):
        w = np.array([omega_01(p, EJ_sum, EC, d_val) for p in phis])
        ax.plot(phis/np.pi, w/(2*np.pi*1e9), ls, lw=2, label=f'd = {d_val}')
    ax.axvline(0.0, color='gray', alpha=0.3)
    ax.axvline(0.5, color='gray', alpha=0.3)
    ax.text(0.02, 4.5, 'USS', fontsize=10, color='gray')
    ax.text(0.46, 4.5, 'LSS', fontsize=10, color='gray')
    ax.set_xlabel(r'$\Phi/\Phi_0$')
    ax.set_ylabel(r'$\omega_{01}/2\pi$ (GHz)')
    ax.set_title(r'B. Frequenza qubit vs flusso (SQUID asimmetrico)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # ---- Pannello C: T1^qp(Phi) per vari d, gap simmetrico ----
    ax = axes[1, 0]
    x_qp_4K = x_qp_thermal(T_op, Delta_NbN)
    for d_val, ls in zip([0.0, 0.3, 0.6, 0.9], ['-', '--', '-.', ':']):
        T1 = np.array([T1_qp_squid(p, x_qp_4K, EJ_sum, EC, d_val, Delta_NbN)
                       for p in phis])
        ax.semilogy(phis/np.pi, T1*1e9, ls, lw=2, label=f'd = {d_val}')
    ax.axvline(0.0, color='gray', alpha=0.3)
    ax.axvline(0.5, color='gray', alpha=0.3)
    ax.set_xlabel(r'$\Phi/\Phi_0$')
    ax.set_ylabel(r'$T_1^{\rm qp}$ (ns)')
    ax.set_title(rf'C. $T_1^{{\rm qp}}$ vs flusso, $T=4$ K termico ($x_{{\rm qp}}={x_qp_4K:.1e}$)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # ---- Pannello D: gap asymmetry sweep ----
    ax = axes[1, 1]
    gap_diffs_meV = np.linspace(0, 2.0, 200)
    # Tre punti operativi diversi
    for phi_op, label, ls in zip(
        [0.0, np.pi*0.25, np.pi*0.5],
        [r'USS $\Phi=0$', r'$\Phi=\Phi_0/4$', r'LSS $\Phi=\Phi_0/2$'],
        ['-', '--', '-.']
    ):
        T1_vs_dgap = []
        for dg in gap_diffs_meV:
            Delta2 = Delta_NbN - dg * meV_J
            if Delta2 <= 0:
                T1_vs_dgap.append(np.nan)
                continue
            T1 = T1_qp_squid(phi_op, x_qp_4K, EJ_sum, EC, 0.5,
                             Delta_NbN, Delta2, T=T_op)
            T1_vs_dgap.append(T1)
        ax.semilogy(gap_diffs_meV, np.array(T1_vs_dgap)*1e6, ls, lw=2,
                    label=label)
    ax.set_xlabel(r'$|\Delta_1 - \Delta_2|$ (meV)')
    ax.set_ylabel(r'$T_1^{\rm qp}$ ($\mu$s)')
    ax.set_title(r'D. Soppressione gap-asymmetry ($d=0.5$, $T=4$ K)')
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/claude/qp_kinetics_squid.png', dpi=140, bbox_inches='tight')
    print("Plot salvato: /home/claude/qp_kinetics_squid.png\n")


# ============================================================
# 7. Tabella riepilogo numerico
# ============================================================
def print_summary():
    Delta_NbN = 2.8 * meV_J
    T_op = 4.0
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h

    x_th = x_qp_thermal(T_op, Delta_NbN)

    print("=" * 60)
    print(f"RIEPILOGO @ T = {T_op} K, Delta = 2.8 meV (NbN)")
    print("=" * 60)
    print(f"E_J,sum/h = {EJ_sum/h/1e9:.1f} GHz,  E_C/h = {EC/h/1e9:.2f} GHz")
    print(f"omega_01/2pi (Phi=0, d=0): {omega_01(0,EJ_sum,EC,0)/(2*np.pi*1e9):.2f} GHz")
    print()
    print(f"x_qp termico (limite intrinseco 4 K): {x_th:.2e}")
    print()
    print("--- Scenari engineering ---")
    print(f"{'Scenario':<45} {'x_qp':>10} {'T1 (us)':>10}")
    print("-" * 67)

    scenarios = [
        ("Termico puro, d=0, gap simm.",
         x_th, T1_qp_squid(0, x_th, EJ_sum, EC, 0.0, Delta_NbN)),
        ("Termico, d=0.5, USS",
         x_th, T1_qp_squid(0, x_th, EJ_sum, EC, 0.5, Delta_NbN)),
        ("Termico, d=0.5, LSS (Phi=Phi0/2)",
         x_th, T1_qp_squid(np.pi/2, x_th, EJ_sum, EC, 0.5, Delta_NbN)),
        ("Trap aggressivo s=1e8/s, d=0.5, USS",
         x_qp_steady(T_op, Delta_NbN, s=1e8),
         T1_qp_squid(0, x_qp_steady(T_op, Delta_NbN, s=1e8),
                     EJ_sum, EC, 0.5, Delta_NbN)),
        ("Gap-asym 1 meV (D1=2.8, D2=1.8 meV), d=0.5, USS",
         x_th, T1_qp_squid(0, x_th, EJ_sum, EC, 0.5,
                            2.8*meV_J, 1.8*meV_J, T=T_op)),
        ("Gap-asym 1.5 meV + trap 1e8, d=0.5, USS",
         x_qp_steady(T_op, Delta_NbN, s=1e8),
         T1_qp_squid(0, x_qp_steady(T_op, Delta_NbN, s=1e8),
                     EJ_sum, EC, 0.5,
                     2.8*meV_J, 1.3*meV_J, T=T_op)),
        ("T = 2 K invece di 4 K (no engineering)",
         x_qp_thermal(2.0, Delta_NbN),
         T1_qp_squid(0, x_qp_thermal(2.0, Delta_NbN),
                     EJ_sum, EC, 0.0, Delta_NbN)),
    ]

    for name, xv, T1v in scenarios:
        print(f"{name:<45} {xv:>10.2e} {T1v*1e6:>10.3f}")
    print()


if __name__ == "__main__":
    run_sanity_checks()
    print_summary()
    make_plots()
