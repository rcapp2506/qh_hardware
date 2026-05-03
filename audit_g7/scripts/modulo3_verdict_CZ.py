"""
Modulo 3 — Verdetto sul gate CZ free-evolution
===============================================

Calcolo della fidelity di un gate CZ implementato come free-evolution sotto ζ_zz σ_z⊗σ_z,
con tempo di gate t_CZ = π/(2|ζ_zz|).

Modello di errore: limitazione di T_1 (relaxation) e T_phi (pure dephasing).
Per un gate di durata t_gate, la fidelity è approssimativamente:

   F ≈ 1 - (4/15) (t_gate/T_1) - (4/15) (t_gate/T_phi)        [Abad 2022, Krantz 2019]
                  ↑                          ↑
            relax. errors           dephasing errors

Per F > 99% (soglia surface code) serve t_gate << T_2*.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*78)
print("MODULO 3 — Verdetto: il CZ free-evolution è impraticabile per HEATS-Q")
print("="*78)

# Coherence times: dichiarati nel manoscritto (§3.6, baseline 5-8 GHz)
# T_1 baseline mK: ~ 50 μs (transmon Al, state-of-the-art)
# T_2 (Hahn echo) baseline mK: ~ 30 μs
# T_1 innovation 300 GHz @ 4K: stima ~ 1-10 μs (degradato da quasi-particles thermal)
# T_2 innovation 300 GHz @ 4K: stima ~ 0.5-5 μs

systems = [
    {
        'name': 'Baseline mK (5-8 GHz)',
        'zeta_zz_MHz': 0.124,
        'T1_us': 50,
        'T2_us': 30,
    },
    {
        'name': 'Innovation 300 GHz @ 4 K',
        'zeta_zz_MHz': 0.0009,    # 0.9 kHz
        'T1_us': 5,                # stima ottimistica per sub-THz transmon a 4K
        'T2_us': 1,                # stima
    },
]

def fidelity_decoherence(t_gate_us, T1_us, T2_us):
    """
    Average gate fidelity with depolarizing noise channel from T1, T2.
    Per Abad et al. PRX Quantum 3, 030327 (2022) o Krantz 2019:
        F = 1 - (1/2)·t/T1 - (1/4)·t/T_phi  (semplificato per CZ)
    Una formula più conservativa è:
        F ≈ exp(-t/T_2eff)  con T_2eff^-1 = 1/(2T_1) + 1/T_phi
    """
    Tphi_inv = max(0, 1/T2_us - 1/(2*T1_us))
    if Tphi_inv == 0:
        Tphi_us = float('inf')
    else:
        Tphi_us = 1/Tphi_inv
    # Approx exponential decoherence envelope:
    F = np.exp(-t_gate_us/T1_us) * np.exp(-t_gate_us/Tphi_us)**0.5
    return F, Tphi_us

print(f"\n{'Sistema':<30s} {'ζ_zz [MHz]':>11s} {'t_CZ [μs]':>12s} {'T_1 [μs]':>10s} {'T_2 [μs]':>10s} {'F_CZ':>8s}")
print("-"*90)

results = []
for sys in systems:
    t_CZ = np.pi / (2 * sys['zeta_zz_MHz'])  # in μs (perché 1/MHz = 1 μs)
    F, _ = fidelity_decoherence(t_CZ, sys['T1_us'], sys['T2_us'])
    print(f"{sys['name']:<30s} {sys['zeta_zz_MHz']:>11.4f} {t_CZ:>12.2f} "
          f"{sys['T1_us']:>10.1f} {sys['T2_us']:>10.1f} {F:>8.4f}")
    results.append((sys['name'], t_CZ, F))

print(f"""
─────────────────────────────────────────────────────────────────────────────
Soglia di "usabilità" per QEC: F > 99% (surface code threshold)
                              F > 99.9% (per qubit logici a basso overhead)
─────────────────────────────────────────────────────────────────────────────

VERDETTO:
""")
for name, t_CZ, F in results:
    if F < 0.99:
        if F < 0.5:
            verdict = "❌ IMPRATICABILE (F < 50%)"
        elif F < 0.9:
            verdict = "❌ INADEGUATO (F < 90%)"
        else:
            verdict = "❌ SOTTO SOGLIA QEC (F < 99%)"
    else:
        verdict = "✓ OK"
    print(f"  {name}:   {verdict}    [t_CZ = {t_CZ:.1f} μs, F = {F*100:.1f}%]")

# ────────────────────────────────────────────────────────────────────────────
# CONFRONTO: cosa succederebbe se la formula 2.39 fosse vera?
# ────────────────────────────────────────────────────────────────────────────
print(f"""
──────────────────────────────────────────────────────────────────────────
CONFRONTO: il gate-time DICHIARATO nella tesi originale era
   t_CZ = π/(2|J|) = π/(2·7.03 MHz) = 224 ns      (basato su formula sbagliata)

Con il VERO ζ_zz = 0.124 MHz:
   t_CZ = π/(2·0.124 MHz) = {np.pi/(2*0.124)*1e3:.0f} ns ≈ {np.pi/(2*0.124):.1f} μs
   
   Discrepanza:  {12800/224:.0f}x più lento del previsto.
""")

# ────────────────────────────────────────────────────────────────────────────
# FIGURA: F vs ζ_zz a confronto con T_1
# ────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
zeta_range_MHz = np.logspace(-3, 1.5, 200)  # 0.001 — 30 MHz
t_gate_us = np.pi / (2 * zeta_range_MHz)  # μs (perché 1/MHz = 1 μs)

for sys in systems:
    F = np.array([fidelity_decoherence(t, sys['T1_us'], sys['T2_us'])[0]
                  for t in t_gate_us])
    ax.semilogx(zeta_range_MHz, F*100, label=f"{sys['name']} (T_1={sys['T1_us']} μs)",
                linewidth=2)
    # Marker per il valore vero
    F_actual = fidelity_decoherence(np.pi/(2*sys['zeta_zz_MHz']),
                                     sys['T1_us'], sys['T2_us'])[0]
    ax.scatter([sys['zeta_zz_MHz']], [F_actual*100], s=120, marker='o',
               edgecolor='black', zorder=5, label=f"  → ζ_zz vero = {sys['zeta_zz_MHz']:.4f} MHz")

ax.axhline(99, color='gray', linestyle='--', alpha=0.6, label='Surface-code threshold (99%)')
ax.axvline(0.124, color='C0', linestyle=':', alpha=0.4)
ax.axvline(0.0009, color='C1', linestyle=':', alpha=0.4)
ax.set_xlabel('Cross-Kerr |ζ_zz| (MHz)', fontsize=11)
ax.set_ylabel('CZ gate fidelity F (%)', fontsize=11)
ax.set_title('CZ free-evolution fidelity vs |ζ_zz|\n(t_CZ = π/(2|ζ_zz|), F limited by T₁/T₂)',
             fontsize=11)
ax.set_ylim(0, 100)
ax.legend(loc='lower right', fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('/home/claude/audit_g7/fig_modulo3_F_vs_zeta.pdf', bbox_inches='tight')
print("Figura salvata: fig_modulo3_F_vs_zeta.pdf")
