"""
qp_analytical_solutions.py
==========================

Verifica numerica delle soluzioni analitiche di Rothwarf-Taylor + trap
derivate analiticamente.

1. Decadimento iperbolico puro (no trap): n(t) = n0/(1+Rn0 t)
2. Decadimento esponenziale (trap dominante): n(t) = n_ss + (n0-n_ss)e^{-st}
3. Soluzione completa Riccati: n(t) con G, s, R
4. Stato stazionario: n^ss = (-s + sqrt(s^2+4RG))/(2R)

Tutti i risultati cross-checkati contro integrazione numerica ODE.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from qp_kinetics_squid import meV_J, kB, x_qp_thermal

# ============================================================
# Parametri tipici NbN @ 4 K
# ============================================================
Delta_NbN = 2.8 * meV_J
T = 4.0
x_eq = x_qp_thermal(T, Delta_NbN)  # ~ 2.7e-4

# In termini di x = n/n_cp (normalizzato a coppie di Cooper):
# G_ph = R * x_eq^2 (detailed balance)
# Useremo unità arbitrarie con R = 1/(100 ns) come riferimento NbN
R = 1.0 / 100e-9      # 1/s
G_th = R * x_eq**2     # da detailed balance termico
print(f"x_eq termico a 4K, NbN: {x_eq:.3e}")
print(f"G_th = R * x_eq² = {G_th:.3e} s⁻¹")
print()


# ============================================================
# 1. Decadimento iperbolico (no trap, no thermal generation)
# ============================================================
def hyperbolic_decay(t, x0, R):
    return x0 / (1 + R*x0*t)

# Verifica numerica
def rhs_norecomb_notrap(x, t):
    return -R * x**2

t = np.linspace(0, 5e-6, 1000)
x0_pulse = 1e-2  # iniezione iniziale "pulse"
x_num = odeint(rhs_norecomb_notrap, x0_pulse, t).flatten()
x_an = hyperbolic_decay(t, x0_pulse, R)
err1 = np.max(np.abs(x_num - x_an)/x_an)
print(f"[Test 1] Decadimento iperbolico:    err max = {err1:.2e}")


# ============================================================
# 2. Decadimento esponenziale (trap dominante)
# ============================================================
def exp_decay_with_trap(t, x0, x_ss, s):
    return x_ss + (x0 - x_ss) * np.exp(-s*t)

s_trap = 1e7  # 1/s
G = G_th
# x_ss approssimato linearmente (regime trap-dominante)
x_ss_lin = G / s_trap

def rhs_lin(x, t):
    return G - s_trap*x  # linearizzato (no R)

x_num2 = odeint(rhs_lin, x0_pulse, t).flatten()
x_an2 = exp_decay_with_trap(t, x0_pulse, x_ss_lin, s_trap)
err2 = np.max(np.abs(x_num2 - x_an2))
print(f"[Test 2] Decadimento esp. (trap):   err max = {err2:.2e}")


# ============================================================
# 3. Soluzione completa (Riccati): G - s n - R n^2 = 0 + ODE
# ============================================================
def steady_state_full(G, s, R):
    """Radice positiva del polinomio quadratico."""
    disc = s**2 + 4*R*G
    return (-s + np.sqrt(disc)) / (2*R)

def riccati_solution(t, x0, G, s, R):
    """Soluzione analitica chiusa di dx/dt = G - sx - Rx^2."""
    x_plus = steady_state_full(G, s, R)
    x_minus = (-s - np.sqrt(s**2 + 4*R*G)) / (2*R)
    # Rate di convergenza esponenziale verso x_plus
    lam = R * (x_plus - x_minus)  # = sqrt(s^2 + 4RG)
    A = (x0 - x_plus) / (x0 - x_minus)
    return (x_plus - x_minus * A * np.exp(-lam*t)) / (1 - A*np.exp(-lam*t))

def rhs_full(x, t):
    return G - s_trap*x - R*x**2

x_num3 = odeint(rhs_full, x0_pulse, t).flatten()
x_an3 = riccati_solution(t, x0_pulse, G, s_trap, R)
err3 = np.max(np.abs(x_num3 - x_an3)/x_num3)
print(f"[Test 3] Riccati completa:          err max = {err3:.2e}")

x_ss_full = steady_state_full(G, s_trap, R)
print(f"        x^ss completo = {x_ss_full:.3e}")
print(f"        x^ss limite trap (G/s) = {G/s_trap:.3e}")
print(f"        x^ss limite recomb (√G/R) = {np.sqrt(G/R):.3e}")


# ============================================================
# 4. Mappa dei regimi nel piano (s, G)
# ============================================================
print()
print("Crossover regimes:")
print("  Trap domina:   s > sqrt(4 R G)  ⇔  G < s²/(4R)")
for s_test in [1e4, 1e6, 1e8]:
    G_cross = s_test**2 / (4*R)
    print(f"    s = {s_test:.0e} s⁻¹  →  G* = {G_cross:.2e} s⁻¹")


# ============================================================
# 5. Plot dei tre regimi
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Pannello A: dinamica temporale - tre scenari
ax = axes[0]
scenarios = [
    ('No trap, no generation\n(iperbolico)',   0,    0,    'C0'),
    ('Trap solo (lineare)',                    0, 1e7,    'C1'),
    ('Full (Riccati): G + s + R',         G_th, 1e7,    'C2'),
    ('Solo R + G (no trap)',              G_th,   0,    'C3'),
]
t_plot = np.logspace(-9, -4, 500)
for label, G_v, s_v, color in scenarios:
    if G_v == 0 and s_v == 0:
        x_t = hyperbolic_decay(t_plot, x0_pulse, R)
    elif R == 0:
        x_t = exp_decay_with_trap(t_plot, x0_pulse, G_v/s_v if s_v else 0, s_v)
    else:
        # ODE integration
        x_t = odeint(lambda x, t: G_v - s_v*x - R*x**2, x0_pulse, t_plot).flatten()
    ax.loglog(t_plot*1e6, x_t, color=color, lw=2, label=label)
ax.axhline(x_eq, color='gray', linestyle=':', alpha=0.6)
ax.text(2e-3, x_eq*1.3, r'$x_{qp}^{th}(4\,\mathrm{K})$', color='gray', fontsize=9)
ax.set_xlabel(r'tempo ($\mu$s)')
ax.set_ylabel(r'$x_{qp}(t)$')
ax.set_title('A. Decadimento dinamico nei tre regimi')
ax.legend(fontsize=8, loc='lower left')
ax.grid(True, which='both', alpha=0.3)

# Pannello B: stato stazionario vs s_trap
ax = axes[1]
s_range = np.logspace(2, 11, 100)
x_ss_array = np.array([steady_state_full(G_th, s, R) for s in s_range])
ax.loglog(s_range, x_ss_array, 'C0-', lw=2.5, label='Soluzione completa')
ax.loglog(s_range, G_th/s_range, 'C1--', lw=1.5, alpha=0.7,
          label=r'Limite trap: $G/s$')
ax.axhline(np.sqrt(G_th/R), color='C3', linestyle='--', lw=1.5, alpha=0.7,
           label=r'Limite recomb: $\sqrt{G/R}$')
ax.axhline(x_eq, color='gray', linestyle=':', alpha=0.6)
ax.text(1e9, x_eq*1.3, r'$x_{qp}^{th}$', color='gray', fontsize=9)
# Crossover
s_crossover = 2*np.sqrt(R*G_th)
ax.axvline(s_crossover, color='black', linestyle=':', alpha=0.5)
ax.text(s_crossover*1.3, 1e-9, 'crossover\n' + fr'$s^*=2\sqrt{{RG}}$',
        fontsize=8)
ax.set_xlabel(r'$s_{\rm trap}$ (s$^{-1}$)')
ax.set_ylabel(r'$x_{qp}^{ss}$')
ax.set_title('B. Stato stazionario vs efficienza trap')
ax.legend(fontsize=9, loc='lower left')
ax.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig('/home/claude/qp_analytical_verification.png', dpi=140,
            bbox_inches='tight')
print("\nPlot salvato: /home/claude/qp_analytical_verification.png")
