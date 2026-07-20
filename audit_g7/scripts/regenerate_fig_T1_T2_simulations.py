"""
Regenerate thesis Figs 2.16 (T1 relaxation) and 2.17 (T2 Ramsey regimes).

Reconstruction: the original generator was never pushed to qh_hardware
(same provenance gap as the original Fig 2.8 script). The noise
realisation therefore differs from the previously published PNGs; the
fitted T1 values quoted in the thesis caption are updated accordingly.

PHYSICS FIX vs the published Fig 2.17: the left panel claimed a "pure T1
limit" with T2 = T1/2 = 50 ns, but from 1/T2 = 1/(2*T1) + 1/T_phi the
T_phi -> infinity limit gives T2 = 2*T1 = 200 ns. The left panel (and
its envelope) now uses the correct value.

Deterministic: fixed RNG seed (42).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

rng = np.random.default_rng(42)
NOISE_SIGMA = 0.02

# ============================================================
# Figure 1: T1 energy-relaxation measurement simulation
# ============================================================
T1_values = [50, 100, 200]          # ns (illustrative short values)
t = np.linspace(0, 500, 400)        # ns
inv_e = np.exp(-1.0)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fitted = []
for ax, T1 in zip(axes, T1_values):
    P1_theory = np.exp(-t / T1)
    P0_theory = 1.0 - P1_theory
    data = np.clip(P1_theory + rng.normal(0, NOISE_SIGMA, t.size), -0.05, 1.05)

    popt, _ = curve_fit(lambda tt, tau: np.exp(-tt / tau), t, data, p0=[T1])
    T1_fit = popt[0]
    fitted.append(T1_fit)

    ax.plot(t, data, '.', color='blue', markersize=3, alpha=0.55,
            label=r'$|1\rangle$ data')
    ax.plot(t, P1_theory, '-', color='blue', linewidth=2.5,
            label=r'$|1\rangle$ (theory)')
    ax.plot(t, P0_theory, '-', color='red', linewidth=2.5, alpha=0.85,
            label=r'$|0\rangle$ (theory)')
    ax.plot(t, np.exp(-t / T1_fit), 'k--', linewidth=2.5,
            label=f'Fit: T₁ = {T1_fit:.1f} ns')

    # 1/e marker
    ax.plot(T1, inv_e, 'o', color='darkgreen', markersize=11, zorder=5)
    ax.axhline(inv_e, color='green', linestyle=':', linewidth=1.2, alpha=0.7)
    ax.axvline(T1, color='green', linestyle=':', linewidth=1.8, alpha=0.9)
    ax.annotate(f'T₁ = {T1} ns\nP(1) = 1/e', xy=(T1, inv_e),
                xytext=(T1 + 12, inv_e - 0.005), fontsize=11,
                color='darkgreen', fontweight='bold')

    ax.set_title(f'Energy Relaxation\nT₁ = {T1} ns',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Time (ns)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Population', fontsize=13, fontweight='bold')
    ax.set_xlim(0, 500)
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='center right', fontsize=10, framealpha=0.9)

plt.tight_layout()
plt.savefig('./T1_measurement_simulation.png', dpi=300, bbox_inches='tight')
print('Saved T1_measurement_simulation.png; fitted T1 =',
      ', '.join(f'{v:.1f}' for v in fitted), 'ns')

# ============================================================
# Figure 2: Ramsey fringe simulations (T2 regimes), T1 = 100 ns
# ============================================================
T1 = 100.0                # ns
delta_f = 0.008           # GHz -> 8 MHz detuning
tR = np.linspace(0, 300, 360)

scenarios = [
    # (title-line-1, T_phi, box label for T_phi)
    ('Pure T₁ limit (T₂ = 2T₁)', np.inf, 'T_φ → ∞'),
    ('T₁ + weak T_φ', 200.0, 'T_φ = 200 ns'),
    ('T₁ + strong T_φ', 50.0, 'T_φ = 50 ns'),
]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
e_level = 0.5 + 0.5 * inv_e
for ax, (title, T_phi, phi_lbl) in zip(axes, scenarios):
    inv_T2 = 1.0 / (2 * T1) + (0.0 if np.isinf(T_phi) else 1.0 / T_phi)
    T2 = 1.0 / inv_T2

    env_hi = 0.5 + 0.5 * np.exp(-tR / T2)
    env_lo = 0.5 - 0.5 * np.exp(-tR / T2)
    signal = 0.5 + 0.5 * np.exp(-tR / T2) * np.cos(2 * np.pi * delta_f * tR)
    data = signal + rng.normal(0, NOISE_SIGMA, tR.size)

    ax.fill_between(tR, env_lo, env_hi, color='gray', alpha=0.15, zorder=1)
    ax.plot(tR, data, '.', color='blue', markersize=3, alpha=0.45, label='Data')
    ax.plot(tR, signal, '-', color='blue', linewidth=2.2, alpha=0.85,
            label='Ramsey signal')
    ax.plot(tR, env_hi, 'r--', linewidth=2.2,
            label=f'Envelope (T₂ = {T2:.0f} ns)')
    ax.plot(tR, env_lo, 'r--', linewidth=2.2)

    ax.plot(T2, e_level, 'o', color='darkgreen', markersize=11, zorder=5)
    ax.axhline(e_level, color='green', linestyle=':', linewidth=1.2, alpha=0.7)
    ax.axvline(T2, color='green', linestyle=':', linewidth=1.8, alpha=0.9)
    ax.annotate(f'T₂ = {T2:.0f} ns', xy=(T2, e_level),
                xytext=(T2 + 8, e_level + 0.01), fontsize=11,
                color='darkgreen', fontweight='bold')

    box = (f'T₁ = 100 ns\n{phi_lbl}\nT₂ = {T2:.0f} ns')
    ax.text(0.97, 0.97, box, transform=ax.transAxes, fontsize=10.5,
            va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))

    ax.set_title(f'{title}\n1/T₂ = 1/(2T₁) + 1/T_φ',
                 fontsize=13.5, fontweight='bold')
    ax.set_xlabel('Time (ns)', fontsize=13, fontweight='bold')
    ax.set_ylabel(r'P$_1$ (after 2nd $\pi$/2)', fontsize=13, fontweight='bold')
    ax.set_xlim(0, 300)
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower right', fontsize=10, framealpha=0.9)

plt.tight_layout()
plt.savefig('./T2_ramsey_simulation.png', dpi=300, bbox_inches='tight')
print('Saved T2_ramsey_simulation.png; T2 =',
      ', '.join(f"{1.0/(1.0/(2*T1)+(0 if np.isinf(p) else 1.0/p)):.0f}"
                for _, p, _ in scenarios), 'ns')
