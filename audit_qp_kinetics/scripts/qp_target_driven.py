"""
qp_target_driven.py
===================

Analisi inversa: dato il target T1 = 49 us a 4 K, quale combinazione
(δΔ, s_1, s_2) lo realizza?

Modello esteso: bilancio cinetico SEPARATO per i due lati della giunzione
nel SQUID asimmetrico:
    x_qp,j(T) = soluzione di r x^2 + s_j x - G_ph,j(T) = 0

con G_ph,j = r * x_qp^th(Δ_j)^2 (detailed balance su ciascun lato).

Rate di tunneling combinati (Marchegiani regime, valido quando almeno uno
dei due x_qp è abbastanza ridotto da rendere il sistema non-termalizzato):

    Γ_{1→2} = (ω/π) sqrt(2Δ_1/ℏω) · x_qp,1 · (E_J1/E_J) cos²(φ_1/2)
              (downhill, no suppression)
    Γ_{2→1} = (ω/π) sqrt(2Δ_2/ℏω) · x_qp,2 · (E_J2/E_J) cos²(φ_2/2)
              · e^{-max(0, Δ_1-Δ_2 - ℏω)/kBT}
              (uphill, gap-asymmetry suppression)

Riferimenti:
- Marchegiani, Amico, Catelani, PRX Quantum 3, 040338 (2022), Sez. III:
  Marchegiani regime ⇔ x_qp non-equilibrio (trap-limited)
- Riwar, Hosseinkhani, Burkhart, Gao, Schoelkopf, Glazman, Catelani,
  PRB 94, 104516 (2016): rate sperimentali di trap normal-metal in Al,
  range realistico s ~ 10^4 - 10^6 1/s
- Pan, Gershenson et al., Nat. Commun. 13, 7196 (2022): combined gap eng
  + trap in Al transmons
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from qp_kinetics_squid import (
    hbar, kB, e, meV_J, GHz_h,
    x_qp_thermal, x_qp_steady,
    EJ_squid, phase_split, omega_01,
)
from qp_kinetics_squid_thermal import delta_BCS, MATERIALS


# ============================================================
# 1. T1 con trap BILATERALE
# ============================================================
def T1_bilateral(phi, EJ_sum, EC, d, Delta1_J, Delta2_J,
                 T=4.0, s1=0.0, s2=0.0, r=1e7):
    """
    T1^qp con bilancio cinetico separato sui due lati della giunzione.

    Parametri
    ---------
    phi      : flusso ridotto π Φ/Φ_0 [rad]
    EJ_sum   : E_J,Σ [J]
    EC       : E_C [J]
    d        : asimmetria geometrica EJ
    Delta1_J : gap lato 1 a T_op [J] (convenzione: Delta1 ≥ Delta2)
    Delta2_J : gap lato 2 a T_op [J]
    T        : temperatura operativa [K]
    s1, s2   : rate di trap sui due lati [1/s]
    r        : rate di ricombinazione [1/s]

    Ritorna
    -------
    dict con T1, x_qp_1, x_qp_2, Γ_down, Γ_up, suppress
    """
    # x_qp stazionari su ciascun lato (con eventuale trap)
    x1 = x_qp_steady(T, Delta1_J, G_gamma=0, s=s1, r=r)
    x2 = x_qp_steady(T, Delta2_J, G_gamma=0, s=s2, r=r)

    # Geometria SQUID
    EJ = EJ_squid(phi, EJ_sum, d)
    w01 = omega_01(phi, EJ_sum, EC, d)
    phi1, phi2 = phase_split(phi, d)
    EJ1 = 0.5 * (1 + d) * EJ_sum
    EJ2 = 0.5 * (1 - d) * EJ_sum
    w1 = (EJ1/EJ) * np.cos(phi1/2)**2
    w2 = (EJ2/EJ) * np.cos(phi2/2)**2

    # Soppressione Marchegiani per canale uphill
    sup = np.exp(-max(0.0, (Delta1_J - Delta2_J) - hbar*w01) / (kB * T))

    # Rate downhill: dal lato a gap maggiore (1) verso minore (2)
    Gamma_down = (w01/np.pi) * x1 * np.sqrt(2*Delta1_J/(hbar*w01)) * w1

    # Rate uphill: dal lato a gap minore (2) verso maggiore (1), soppresso
    Gamma_up = (w01/np.pi) * x2 * np.sqrt(2*Delta2_J/(hbar*w01)) * w2 * sup

    Gamma_tot = Gamma_down + Gamma_up
    return {
        'T1': 1.0/Gamma_tot,
        'x1': x1, 'x2': x2,
        'Gamma_down': Gamma_down, 'Gamma_up': Gamma_up,
        'suppress': sup, 'w01_GHz': w01/(2*np.pi*1e9)
    }


# ============================================================
# 2. Sanity checks
# ============================================================
def sanity():
    T_op = 4.0
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h
    d = 0.5
    D1_J = delta_BCS(T_op, 16.0, 2.80) * meV_J  # NbN bulk
    D2_J = delta_BCS(T_op,  9.3, 1.45) * meV_J  # Nb puro

    print("=" * 75)
    print("SANITY CHECKS - modello bilaterale a 4 K, NbN/Nb gap engineering")
    print("=" * 75)
    print(f"Δ_1 (NbN) = {D1_J/meV_J:.3f} meV, Δ_2 (Nb) = {D2_J/meV_J:.3f} meV")
    print(f"δΔ = {(D1_J-D2_J)/meV_J:.3f} meV, k_BT = {kB*T_op/e*1e3:.3f} meV")
    print(f"Soppressione Marchegiani: e^(-δΔ/kBT) = "
          f"{np.exp(-(D1_J-D2_J)/(kB*T_op)):.2e}")
    print()

    cases = [
        ('No trap',          0.0,    0.0),
        ('Trap solo lato 1', 1e6,    0.0),
        ('Trap solo lato 2', 0.0,    1e6),
        ('Trap bilaterale 1e5', 1e5, 1e5),
        ('Trap bilaterale 1e6', 1e6, 1e6),
        ('Trap bilaterale 1e7', 1e7, 1e7),
        ('Trap bilaterale 1e8', 1e8, 1e8),
    ]
    print(f"{'Scenario':<22} {'x_qp,1':>10} {'x_qp,2':>10} "
          f"{'Γ_down':>10} {'Γ_up*sup':>10} {'T1[μs]':>10}")
    print("-" * 75)
    for name, s1, s2 in cases:
        res = T1_bilateral(0, EJ_sum, EC, d, D1_J, D2_J, T=T_op,
                            s1=s1, s2=s2)
        print(f"{name:<22} {res['x1']:>10.2e} {res['x2']:>10.2e} "
              f"{res['Gamma_down']:>10.2e} {res['Gamma_up']:>10.2e} "
              f"{res['T1']*1e6:>10.3f}")
    print()


# ============================================================
# 3. Mappa (s_1, s_2) → T1 con curva iso-49 μs
# ============================================================
def plot_engineering_map():
    T_op = 4.0
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h
    d = 0.5

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Caso A: NbN/Nb (gap engineering moderato)
    # Caso B: NbN bulk simm. (no gap engineering)

    cases = [
        (axes[0], 'NbN/Nb (δΔ=1.4 meV)',
         delta_BCS(T_op, 16.0, 2.80) * meV_J,
         delta_BCS(T_op,  9.3, 1.45) * meV_J),
        (axes[1], 'NbN/NbN simm. (no gap eng.)',
         delta_BCS(T_op, 16.0, 2.80) * meV_J,
         delta_BCS(T_op, 16.0, 2.80) * meV_J),
    ]

    s_range = np.logspace(2, 10, 80)
    S1, S2 = np.meshgrid(s_range, s_range)

    for ax, title, D1, D2 in cases:
        T1_map = np.zeros_like(S1)
        for i in range(S1.shape[0]):
            for j in range(S1.shape[1]):
                res = T1_bilateral(0, EJ_sum, EC, d, D1, D2,
                                    T=T_op, s1=S1[i,j], s2=S2[i,j])
                T1_map[i,j] = res['T1'] * 1e6  # μs

        pcm = ax.pcolormesh(S1, S2, T1_map, norm=LogNorm(vmin=1e-3, vmax=1e4),
                            shading='auto', cmap='viridis')
        cbar = plt.colorbar(pcm, ax=ax)
        cbar.set_label(r'$T_1^{\rm qp}$ ($\mu$s)')

        # Iso-49 µs
        cs = ax.contour(S1, S2, T1_map, levels=[49],
                        colors='red', linewidths=2.5)
        ax.clabel(cs, inline=True, fontsize=10, fmt='49 μs')
        # Iso-1 µs e iso-10 µs come riferimento
        cs2 = ax.contour(S1, S2, T1_map, levels=[1, 10, 100],
                         colors='white', linewidths=1, linestyles='--')
        ax.clabel(cs2, inline=True, fontsize=8, fmt='%.0f μs')

        # Box "regime realistico" Riwar-Catelani 2016
        ax.axvspan(1e4, 1e6, alpha=0.0)
        ax.axhspan(1e4, 1e6, alpha=0.0)
        ax.add_patch(plt.Rectangle((1e4, 1e4), 1e6-1e4, 1e6-1e4,
                                    fill=False, edgecolor='cyan', lw=2,
                                    linestyle=':',
                                    label='range dimostrato\nin Al (RC 2016)'))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$s_1$ (trap rate lato a gap maggiore) [s$^{-1}$]')
        ax.set_ylabel(r'$s_2$ (trap rate lato a gap minore) [s$^{-1}$]')
        ax.set_title(title)
        ax.legend(loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.savefig('/home/claude/qp_target_driven_map.png', dpi=140,
                bbox_inches='tight')
    print("Plot salvato: /home/claude/qp_target_driven_map.png")


# ============================================================
# 4. Minimum trap requirement per target 49 µs
# ============================================================
def minimum_trap_requirement():
    """Per ogni scenario di gap engineering, trova s_min (trap bilaterale
    simmetrico s_1=s_2=s) tale che T1 = 49 µs."""
    from scipy.optimize import brentq
    T_op = 4.0
    EJ_sum = 25 * GHz_h
    EC = 0.25 * GHz_h
    d = 0.5
    target_T1 = 49e-6  # secondi

    print("=" * 75)
    print(f"TRAP MINIMO RICHIESTO (s_1 = s_2 = s) per T_1 = {target_T1*1e6} μs")
    print("=" * 75)
    print(f"{'Scenario':<25} {'Δ_1[meV]':>10} {'Δ_2[meV]':>10} "
          f"{'s_min [1/s]':>14} {'note'}")
    print("-" * 90)

    scenarios = [
        ('NbN bulk simm.',   'NbN bulk',       'NbN bulk'),
        ('NbN/NbN thin',     'NbN bulk',       'NbN thin ~4nm'),
        ('NbN/Nb',           'NbN bulk',       'Nb puro'),
        ('NbN/TaN',          'NbN bulk',       'TaN'),
    ]

    def T1_minus_target(log_s, D1, D2):
        s = 10**log_s
        res = T1_bilateral(0, EJ_sum, EC, d, D1, D2, T=T_op, s1=s, s2=s)
        return res['T1'] - target_T1

    for label, mat1, mat2 in scenarios:
        Tc1, D0_1 = MATERIALS[mat1]
        Tc2, D0_2 = MATERIALS[mat2]
        D1_J = delta_BCS(T_op, Tc1, D0_1) * meV_J
        D2_J = delta_BCS(T_op, Tc2, D0_2) * meV_J

        try:
            log_s_min = brentq(T1_minus_target, 2, 12, args=(D1_J, D2_J))
            s_min = 10**log_s_min
            # Classificazione
            if s_min < 1e6:
                note = "✓ realistico (Al lit.)"
            elif s_min < 1e8:
                note = "△ aggressivo, possibile"
            else:
                note = "✗ fuori portata"
        except ValueError:
            s_min = np.inf
            note = "✗ irraggiungibile a qualsiasi s"

        print(f"{label:<25} {D1_J/meV_J:>10.3f} {D2_J/meV_J:>10.3f} "
              f"{s_min:>14.2e}   {note}")
    print()


if __name__ == "__main__":
    sanity()
    minimum_trap_requirement()
    plot_engineering_map()
