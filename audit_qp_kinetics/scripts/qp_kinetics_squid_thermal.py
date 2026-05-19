"""
qp_kinetics_squid_thermal.py
============================

Modello termicamente consistente per T1 da QP in un transmon SQUID asimmetrico
con gap engineering (Delta_1 != Delta_2), operativo a 4 K.

DIFFERENZA RISPETTO A qp_kinetics_squid_4K.py:
nel regime termico (k_BT non trascurabile rispetto a Δ_2), il modello deve
calcolare separatamente i due rate di tunneling:

    Γ_{2→1} ∝ x_qp^th(Δ_2) · e^{-(Δ_1-Δ_2)/k_BT}  (uphill, soppresso)
    Γ_{1→2} ∝ x_qp^th(Δ_1)                        (downhill, non soppresso)

NON è corretto usare x_qp esterno + suppression: doppia conta della fisica.

Il risultato netto, a meno di un prefattore radice, è:
    Γ_tot ≈ const · e^{-Δ_max/k_BT}
ovvero il rate è governato dal gap maggiore. Gap engineering "discendente"
NON aiuta in regime termico, contrariamente al regime non-equilibrio (mK).

Riferimenti:
- Marchegiani, Amico, Catelani, PRX Quantum 3, 040338 (2022), Sezione III
  (regime termico vs non-equilibrio per gap asymmetric junctions)
- Glazman & Catelani, SciPost Lect. Notes 31 (2021)
"""

import numpy as np
import matplotlib.pyplot as plt
from qp_kinetics_squid import (
    hbar, kB, e, h, meV_J, GHz_h,
    x_qp_thermal, x_qp_steady,
    EJ_squid, phase_split, omega_01,
)


# ============================================================
# 1. Δ(T) BCS (riportato qui per autonomia)
# ============================================================
def delta_BCS(T, Tc, Delta_0_meV):
    """Δ(T) ≈ Δ(0) tanh(1.74 √(Tc/T - 1)). Input/output in meV."""
    T_arr = np.atleast_1d(T)
    out = np.zeros_like(T_arr, dtype=float)
    mask = T_arr < Tc
    out[mask] = Delta_0_meV * np.tanh(1.74 * np.sqrt(Tc/T_arr[mask] - 1))
    return out if np.ndim(T) > 0 else out.item()


# ============================================================
# 2. T1 da QP per SQUID asimmetrico, REGIME TERMICO consistente
# ============================================================
def T1_qp_thermal(phi, EJ_sum, EC, d, Delta1_J, Delta2_J, T=4.0):
    """
    T1 da QP per SQUID asimmetrico in regime termico equilibrio.

    Separa i due contributi:
      Γ_{1→2} = (ω_01/π) sqrt(2 Δ_1/ℏω_01) x_qp^th(Δ_1) · weight_1
      Γ_{2→1} = (ω_01/π) sqrt(2 Δ_2/ℏω_01) x_qp^th(Δ_2)
                                 · e^{-max(0, Δ_1-Δ_2 - ℏω_01)/kBT} · weight_2

    weight_j = (E_Jj/E_J) cos^2(φ_j/2)
    """
    EJ = EJ_squid(phi, EJ_sum, d)
    w01 = omega_01(phi, EJ_sum, EC, d)
    phi1, phi2 = phase_split(phi, d)

    EJ1 = 0.5 * (1 + d) * EJ_sum
    EJ2 = 0.5 * (1 - d) * EJ_sum
    w1 = (EJ1/EJ) * np.cos(phi1/2)**2
    w2 = (EJ2/EJ) * np.cos(phi2/2)**2

    # Densità termiche dei due lati
    x1 = x_qp_thermal(T, Delta1_J)
    x2 = x_qp_thermal(T, Delta2_J)

    # Soppressione uphill (solo se Δ_high > Δ_low + ℏω_01)
    if Delta1_J > Delta2_J:
        Delta_low, Delta_high = Delta2_J, Delta1_J
        x_low, x_high = x2, x1
    else:
        Delta_low, Delta_high = Delta1_J, Delta2_J
        x_low, x_high = x1, x2

    suppress = np.exp(-max(0.0, (Delta_high - Delta_low) - hbar*w01)
                       / (kB * T))

    # Rate uphill (lato a gap minore → lato a gap maggiore)
    Gamma_up = (w01/np.pi) * x_low * np.sqrt(2*Delta_low/(hbar*w01)) \
               * suppress

    # Rate downhill (lato a gap maggiore → lato a gap minore): no soppressione
    Gamma_down = (w01/np.pi) * x_high * np.sqrt(2*Delta_high/(hbar*w01))

    # Distribuzione sui due canali geometrici dello SQUID
    # (entrambe le giunzioni contribuiscono con peso (EJj/EJ)cos²(φj/2))
    Gamma_tot = (Gamma_up + Gamma_down) * (w1 + w2)

    return 1.0 / Gamma_tot


# ============================================================
# 3. T1 con trap normal-metal (modifica x_qp sul lato a gap minore)
# ============================================================
def T1_qp_thermal_with_trap(phi, EJ_sum, EC, d, Delta1_J, Delta2_J,
                             T=4.0, s_trap=0.0, r=1e7):
    """
    Come T1_qp_thermal, ma con trap normal-metal che riduce x_qp nel lato
    a gap minore (assunto: il trap è posizionato per catturare QP nel
    materiale a gap minore, che ha più QP termiche).

    s_trap : rate di trap [1/s]
    r      : rate di ricombinazione [1/s]
    """
    EJ = EJ_squid(phi, EJ_sum, d)
    w01 = omega_01(phi, EJ_sum, EC, d)
    phi1, phi2 = phase_split(phi, d)
    EJ1 = 0.5 * (1 + d) * EJ_sum
    EJ2 = 0.5 * (1 - d) * EJ_sum
    w1 = (EJ1/EJ) * np.cos(phi1/2)**2
    w2 = (EJ2/EJ) * np.cos(phi2/2)**2

    # Determina gap_min e gap_max
    if Delta1_J > Delta2_J:
        Delta_low, Delta_high = Delta2_J, Delta1_J
    else:
        Delta_low, Delta_high = Delta1_J, Delta2_J

    # Stato stazionario nel lato a gap minore con trap attivo
    x_low = x_qp_steady(T, Delta_low, G_gamma=0, s=s_trap, r=r)
    # Lato a gap maggiore: x_qp termico (senza trap)
    x_high = x_qp_thermal(T, Delta_high)

    suppress = np.exp(-max(0.0, (Delta_high - Delta_low) - hbar*w01)
                       / (kB * T))

    Gamma_up = (w01/np.pi) * x_low * np.sqrt(2*Delta_low/(hbar*w01)) * suppress
    Gamma_down = (w01/np.pi) * x_high * np.sqrt(2*Delta_high/(hbar*w01))
    Gamma_tot = (Gamma_up + Gamma_down) * (w1 + w2)

    return 1.0/Gamma_tot, x_low, x_high


# ============================================================
# 4. Materiali (gap a T=0 in meV)
# ============================================================
MATERIALS = {
    'NbN bulk':       (16.0, 2.80),
    'NbN thin ~4nm':  (10.0, 1.50),
    'Nb puro':         (9.3, 1.45),
    'NbTiN':         (14.0, 2.30),
    'TaN':            (7.0, 1.00),
}


# ============================================================
# 5. Riepilogo: scenari termici corretti
# ============================================================
def print_summary_correct():
    T_op = 4.0
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h
    d = 0.5
    Tc1, D0_1 = MATERIALS['NbN bulk']
    D1_meV = delta_BCS(T_op, Tc1, D0_1)
    D1_J = D1_meV * meV_J
    x_th_NbN_bulk = x_qp_thermal(T_op, D1_J)
    w01_GHz = omega_01(0, EJ_sum, EC, d)/(2*np.pi*1e9)

    print("=" * 95)
    print(f"MODELLO TERMICO CONSISTENTE @ T = {T_op} K, d = {d}, "
          f"EJsum/h = 25 GHz")
    print(f"ω_01/2π = {w01_GHz:.2f} GHz, Δ_1 NbN bulk = {D1_meV:.3f} meV, "
          f"x_qp^th(NbN bulk) = {x_th_NbN_bulk:.2e}")
    print("=" * 95)

    # --- Solo gap engineering (no trap), modello corretto ---
    print()
    print("--- Solo gap engineering (no trap) ---")
    print(f"{'Materiale 2':<22} {'Δ_2[meV]':>10} {'δΔ[meV]':>10} "
          f"{'x_qp(Δ_2)':>12} {'T1 USS[μs]':>12}")
    print("-" * 75)
    for name, (Tc2, D0_2) in MATERIALS.items():
        D2_meV = delta_BCS(T_op, Tc2, D0_2)
        D2_J = D2_meV * meV_J
        x_low = x_qp_thermal(T_op, D2_J)
        dg = D1_meV - D2_meV
        T1 = T1_qp_thermal(0, EJ_sum, EC, d, D1_J, D2_J, T=T_op)
        print(f"{name:<22} {D2_meV:>10.3f} {dg:>10.3f} "
              f"{x_low:>12.2e} {T1*1e6:>12.3f}")

    # --- Gap engineering + trap nel lato a basso gap ---
    print()
    print("--- Gap engineering + trap normal-metal sul lato a gap minore ---")
    print(f"{'Configurazione':<30} {'s_trap[1/s]':>12} {'x_low_ss':>12} "
          f"{'T1 USS[μs]':>12}")
    print("-" * 75)
    combined = [
        ('NbN thin (no trap)',  'NbN thin ~4nm', 0),
        ('NbN thin + trap 1e5',  'NbN thin ~4nm', 1e5),
        ('NbN thin + trap 1e6',  'NbN thin ~4nm', 1e6),
        ('NbN thin + trap 1e7',  'NbN thin ~4nm', 1e7),
        ('NbN thin + trap 1e8',  'NbN thin ~4nm', 1e8),
        ('Nb + trap 1e7',        'Nb puro',       1e7),
        ('TaN + trap 1e7',       'TaN',           1e7),
        ('NbN bulk simm. + trap 1e8 (no gap eng)', 'NbN bulk', 1e8),
    ]
    for label, mat2, s_val in combined:
        Tc2, D0_2 = MATERIALS[mat2]
        D2_meV = delta_BCS(T_op, Tc2, D0_2)
        D2_J = D2_meV * meV_J
        T1, x_low_ss, _ = T1_qp_thermal_with_trap(
            0, EJ_sum, EC, d, D1_J, D2_J, T=T_op, s_trap=s_val)
        print(f"{label:<30} {s_val:>12.0e} {x_low_ss:>12.2e} "
              f"{T1*1e6:>12.3f}")

    print()
    print("Target: T1 = 49 μs")


# ============================================================
# 6. Plot
# ============================================================
def plot_correct_model():
    T_op = 4.0
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h
    d = 0.5
    Tc1, D0_1 = MATERIALS['NbN bulk']
    D1_J = delta_BCS(T_op, Tc1, D0_1) * meV_J

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # ---- A: T1 vs δΔ, regime termico corretto vs ingenuo ----
    ax = axes[0]
    dg_meV = np.linspace(0, 2.5, 100)
    T1_correct = []
    T1_naive = []
    x_th_NbN = x_qp_thermal(T_op, D1_J)
    for dg in dg_meV:
        D2_J = D1_J - dg * meV_J
        if D2_J <= 0.3 * meV_J:
            T1_correct.append(np.nan)
            T1_naive.append(np.nan)
            continue
        T1c = T1_qp_thermal(0, EJ_sum, EC, d, D1_J, D2_J, T=T_op)
        T1_correct.append(T1c * 1e6)

        # Modello ingenuo: x_qp = x_th(NbN) + suppression manuale
        w01 = omega_01(0, EJ_sum, EC, d)
        sup = np.exp(-max(0.0, dg*meV_J - hbar*w01) / (kB * T_op))
        Gamma_naive = (w01/np.pi) * x_th_NbN * np.sqrt(2*D1_J/(hbar*w01)) * sup
        T1_naive.append(1.0/Gamma_naive * 1e6)
    ax.semilogy(dg_meV, T1_correct, 'C0-', lw=2.5,
                label='Modello termico corretto')
    ax.semilogy(dg_meV, T1_naive, 'C3--', lw=2,
                label='Modello "ingenuo"\n(usa x_qp di NbN + sup)')
    ax.axhline(49, color='red', linestyle=':', lw=1.5, alpha=0.7)
    ax.text(0.05, 60, 'target 49 μs', color='red', fontsize=10)
    ax.axhline(0.02, color='gray', linestyle=':', alpha=0.5)
    ax.text(0.05, 0.025, 'no gap eng. (20 ns)', color='gray', fontsize=9)
    ax.set_xlabel(r'$|\Delta_1 - \Delta_2|$ (meV)')
    ax.set_ylabel(r'$T_1^{\rm qp}$ ($\mu$s)')
    ax.set_title(r'A. Gap engineering puro a 4 K: termico vs ingenuo')
    ax.legend(loc='center right', fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # ---- B: T1 vs s_trap per vari scenari di gap engineering ----
    ax = axes[1]
    s_range = np.logspace(2, 10, 100)
    for mat2, label, color in [
        ('NbN bulk',      'No gap eng. (NbN simm.)', 'k'),
        ('NbN thin ~4nm', 'NbN thin (δΔ=1.3 meV)',  'C0'),
        ('TaN',           'TaN (δΔ=1.9 meV, T/Tc=0.57)', 'C1'),
    ]:
        Tc2, D0_2 = MATERIALS[mat2]
        D2_J = delta_BCS(T_op, Tc2, D0_2) * meV_J
        T1s = []
        for s in s_range:
            T1, _, _ = T1_qp_thermal_with_trap(0, EJ_sum, EC, d, D1_J, D2_J,
                                                T=T_op, s_trap=s)
            T1s.append(T1*1e6)
        ax.loglog(s_range, T1s, color=color, lw=2, label=label)
    ax.axhline(49, color='red', linestyle=':', lw=1.5, alpha=0.7)
    ax.text(1e3, 60, 'target 49 μs', color='red', fontsize=10)
    ax.axvspan(1e4, 1e6, alpha=0.15, color='green',
               label='range realistico\n(Riwar-Catelani 2016)')
    ax.set_xlabel(r'$s_{\rm trap}$ (s$^{-1}$)')
    ax.set_ylabel(r'$T_1^{\rm qp}$ ($\mu$s)')
    ax.set_title(r'B. $T_1$ vs efficienza trap, con gap engineering')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/claude/qp_kinetics_thermal_correct.png', dpi=140,
                bbox_inches='tight')
    print("Plot salvato: /home/claude/qp_kinetics_thermal_correct.png")


if __name__ == "__main__":
    print_summary_correct()
    plot_correct_model()
