"""
Modulo 3 — Verdetto sul gate CZ free-evolution
===============================================

Fidelity of a free-evolution CZ driven by the longitudinal cross-Kerr
LEVEL SHIFT xi_ZZ = E11 + E00 - E01 - E10 (ratified convention; the
sigma_z x sigma_z coefficient is zeta_zz = xi_ZZ/4).

For a free-evolution CZ the conditional pi phase accumulates at the rate
set by the level shift, so the gate time is

    t_CZ = pi / |xi_ZZ| .

Error model: T1 (relaxation) and T_phi (pure dephasing) limited,

    F = exp(-t/T1) * exp(-t/T_phi)^(1/2),   1/T_phi = 1/T2 - 1/(2 T1).

Canonical operating points (exact diagonalization, Tab. 2 J_vs_zetazz):
    baseline mK (5.8/5.2 GHz):  |xi_ZZ| = 106 kHz
    300 GHz / 4 K:              |xi_ZZ| =  74 kHz
Both fall well below the surface-code threshold, motivating the
echo cross-resonance gate adopted in the chapter.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*78)
print("MODULO 3 — Verdetto: il CZ free-evolution e' impraticabile per HEATS-Q")
print("="*78)

systems = [
    {
        'name': 'Baseline mK (5.8/5.2 GHz)',
        'xi_ZZ_MHz': 0.106,        # 106 kHz (canonical, exact diag)
        'T1_us': 50,
        'T2_us': 30,
        'color': 'C0',
    },
    {
        'name': 'Elevated-T design (300 GHz, 4 K)',
        'xi_ZZ_MHz': 0.074,        # 74 kHz (canonical, exact diag)
        'T1_us': 15.3,             # total T1 at 300 GHz/4K (canonical)
        'T2_us': 13.1,             # T2 (canonical)
        'color': 'C1',
    },
]

def fidelity_decoherence(t_gate_us, T1_us, T2_us):
    Tphi_inv = max(0.0, 1.0/T2_us - 1.0/(2.0*T1_us))
    Tphi_us = (1.0/Tphi_inv) if Tphi_inv > 0 else float('inf')
    F = np.exp(-t_gate_us/T1_us) * np.exp(-t_gate_us/Tphi_us)**0.5
    return F, Tphi_us

print(f"\n{'System':<34s} {'xi_ZZ [kHz]':>12s} {'t_CZ [us]':>11s} "
      f"{'T1 [us]':>9s} {'T2 [us]':>9s} {'F_CZ':>8s}")
print("-"*88)

results = []
for s in systems:
    t_CZ = np.pi / s['xi_ZZ_MHz']     # us  (t_CZ = pi/|xi_ZZ|, 1/MHz = us)
    F, _ = fidelity_decoherence(t_CZ, s['T1_us'], s['T2_us'])
    print(f"{s['name']:<34s} {s['xi_ZZ_MHz']*1e3:>12.0f} {t_CZ:>11.1f} "
          f"{s['T1_us']:>9.1f} {s['T2_us']:>9.1f} {F*100:>7.1f}%")
    results.append((s['name'], t_CZ, F))

print(f"""
-----------------------------------------------------------------------------
QEC usability threshold: F > 99% (surface-code threshold)
-----------------------------------------------------------------------------
VERDICT:""")
for name, t_CZ, F in results:
    verdict = "OK" if F >= 0.99 else "BELOW QEC THRESHOLD (F < 99%)"
    print(f"  {name}:  {verdict}  [t_CZ = {t_CZ:.1f} us, F = {F*100:.1f}%]")

# ----------------------------------------------------------------------------
# FIGURE: F vs |xi_ZZ|
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
xi_range_MHz = np.logspace(-2, 1.5, 300)        # 10 kHz -- ~30 MHz
t_gate_us = np.pi / xi_range_MHz                # t_CZ = pi/|xi_ZZ|

for s in systems:
    F = np.array([fidelity_decoherence(t, s['T1_us'], s['T2_us'])[0]
                  for t in t_gate_us])
    ax.semilogx(xi_range_MHz*1e3, F*100, color=s['color'],
                label=f"{s['name']} ($T_1={s['T1_us']}$ $\\mu$s)", linewidth=2)
    xi = s['xi_ZZ_MHz']
    F_pt = fidelity_decoherence(np.pi/xi, s['T1_us'], s['T2_us'])[0]
    ax.scatter([xi*1e3], [F_pt*100], s=120, marker='o', color=s['color'],
               edgecolor='black', zorder=5,
               label=f"  $\\to |\\xi_{{ZZ}}|$ = {xi*1e3:.0f} kHz, F = {F_pt*100:.0f}%")
    ax.axvline(xi*1e3, color=s['color'], linestyle=':', alpha=0.4)

ax.axhline(99, color='gray', linestyle='--', alpha=0.6,
           label='Surface-code threshold (99%)')
ax.set_xlabel(r'Cross-Kerr level shift $|\xi_{ZZ}|$ (kHz)', fontsize=11)
ax.set_ylabel(r'CZ gate fidelity $F$ (%)', fontsize=11)
ax.set_title(r'CZ free-evolution fidelity vs $|\xi_{ZZ}|$'
             + '\n' + r'($t_{CZ} = \pi/|\xi_{ZZ}|$, $F$ limited by $T_1/T_2$)',
             fontsize=11)
ax.set_ylim(0, 100)
ax.legend(loc='lower right', fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('./CZ_freevolution_fidelity_vs_zeta.pdf', bbox_inches='tight')
plt.savefig('./CZ_freevolution_fidelity_vs_zeta.png', dpi=160, bbox_inches='tight')
print('Figure saved: CZ_freevolution_fidelity_vs_zeta.{pdf,png}')
