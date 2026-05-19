"""
appendix_figures.py
===================

Publication-quality figures for the QP appendix:
  Fig A.1  Pair breaking + BCS DOS + thermal population (3 panels)
  Fig A.2  Kinetic regimes verification (analytical vs numerical)
  Fig A.3  Normal-metal trap mechanism (energy diagram)
  Fig A.4  Trap engineering layout for the asymmetric SQUID transmon

All saved as PDF for LaTeX inclusion (\includegraphics{...}).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle, Ellipse
from matplotlib.lines import Line2D
from scipy.integrate import odeint

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 120,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# ============================================================
# Fig A.1 - Pair breaking, BCS DOS, thermal population
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))

# ---- Panel (a): Pair-breaking schematic (energy diagram) ----
ax = axes[0]
ax.set_xlim(-1, 3.5)
ax.set_ylim(-3.2, 3.2)
ax.set_aspect('equal')

# Condensate (below -Delta)
ax.add_patch(Rectangle((-0.5, -3.0), 3.5, 1.4, facecolor='#B5D4F4',
                        edgecolor='#185FA5', linewidth=0.8))
ax.text(1.25, -2.3, 'Cooper-pair condensate',
        ha='center', va='center', fontsize=9)
ax.text(1.25, -2.7, r'$E \leq -\Delta$', ha='center', va='center',
        fontsize=8, style='italic')

# Gap (forbidden zone, dashed)
ax.add_patch(Rectangle((-0.5, -1.6), 3.5, 3.2, facecolor='#F1EFE8',
                        edgecolor='gray', linewidth=0.5,
                        linestyle='--', alpha=0.4))
ax.text(1.25, 0, r'forbidden gap $2\Delta$',
        ha='center', va='center', fontsize=9, style='italic',
        color='#5F5E5A')

# QP band (above +Delta)
ax.add_patch(Rectangle((-0.5, 1.6), 3.5, 1.4, facecolor='#F0997B',
                        edgecolor='#993C1D', linewidth=0.8))
ax.text(1.25, 2.3, 'Bogoliubov QP band',
        ha='center', va='center', fontsize=9)
ax.text(1.25, 2.7, r'$E \geq +\Delta$', ha='center', va='center',
        fontsize=8, style='italic')

# Cooper pair: two coupled circles
for x_off in [0.4, 1.4]:
    ax.add_patch(Circle((x_off-0.1, -2.0), 0.08, facecolor='#185FA5'))
    ax.add_patch(Circle((x_off+0.1, -2.0), 0.08, facecolor='#185FA5'))
    ax.plot([x_off-0.05, x_off+0.05], [-2.0, -2.0],
            color='#185FA5', linewidth=1.2)

# Broken pair (faded)
ax.add_patch(Circle((2.3, -2.0), 0.08, facecolor='#185FA5', alpha=0.3))
ax.add_patch(Circle((2.5, -2.0), 0.08, facecolor='#185FA5', alpha=0.3))

# Phonon (wavy line)
t = np.linspace(0, 4*np.pi, 100)
ax.plot(3.0 + 0.12*np.sin(t), -0.5 + np.linspace(-1.0, 0.4, 100),
        color='#BA7517', linewidth=1.5)
ax.text(3.3, -0.5, r'phonon' '\n' r'$\hbar\omega \geq 2\Delta$',
        fontsize=8, color='#BA7517')

# Arrows: phonon breaks pair, two QPs go up and down
ax.annotate('', xy=(2.0, 1.7), xytext=(2.7, -0.1),
            arrowprops=dict(arrowstyle='->', color='#BA7517', lw=1.2,
                            linestyle='--'))
ax.annotate('', xy=(2.0, -1.7), xytext=(2.7, -0.1),
            arrowprops=dict(arrowstyle='->', color='#BA7517', lw=1.2,
                            linestyle='--'))

# Generated QPs
ax.add_patch(Circle((1.9, 1.9), 0.12, facecolor='#993C1D',
                     edgecolor='#993C1D'))
ax.text(1.55, 1.9, r'$\gamma_e$', fontsize=9, color='#993C1D', ha='right')

ax.add_patch(Circle((1.9, -1.9), 0.12, facecolor='white',
                     edgecolor='#993C1D', linewidth=1.2))
ax.text(1.55, -1.9, r'$\gamma_h$', fontsize=9, color='#993C1D', ha='right')

# Energy axis
ax.annotate('', xy=(-0.6, 3.0), xytext=(-0.6, -3.0),
            arrowprops=dict(arrowstyle='-|>', color='black', lw=1))
ax.text(-0.75, 3.0, r'$E$', fontsize=10, ha='center')
ax.plot([-0.62, -0.58], [1.6, 1.6], color='black', linewidth=1)
ax.text(-0.65, 1.6, r'$+\Delta$', fontsize=8, ha='right', va='center')
ax.plot([-0.62, -0.58], [-1.6, -1.6], color='black', linewidth=1)
ax.text(-0.65, -1.6, r'$-\Delta$', fontsize=8, ha='right', va='center')
ax.plot([-0.62, -0.58], [0, 0], color='black', linewidth=1)
ax.text(-0.65, 0, r'$E_F$', fontsize=8, ha='right', va='center')

ax.set_xticks([])
ax.set_yticks([])
ax.set_title('(a) pair breaking', loc='left', fontsize=10)
for spine in ax.spines.values():
    spine.set_visible(False)

# ---- Panel (b): BCS DOS ----
ax = axes[1]
E = np.linspace(0, 4, 1000)
Delta = 1.0
# BCS DOS
dos = np.where(E > Delta, E / np.sqrt(np.abs(E**2 - Delta**2 + 1e-9)), 0)
dos = np.clip(dos, 0, 5)
ax.plot(E/Delta, dos, color='#534AB7', linewidth=1.8,
        label=r'$N_S(E) = E/\sqrt{E^2-\Delta^2}$')
ax.fill_between(E[E<Delta]/Delta, 0, 5, color='#F1EFE8', alpha=0.6)
ax.text(0.5, 2.5, 'gap\n(no states)', ha='center', va='center',
        fontsize=9, color='#5F5E5A', style='italic')
ax.axvline(1, color='gray', linestyle=':', linewidth=0.6)
ax.set_xlim(0, 4)
ax.set_ylim(0, 5)
ax.set_xlabel(r'$E / \Delta$')
ax.set_ylabel(r'$N_S(E) / N(0)$')
ax.set_title('(b) BCS density of states', loc='left', fontsize=10)
ax.legend(loc='upper right', frameon=False)
ax.grid(True, alpha=0.3)

# ---- Panel (c): Thermal population N_S * f(E) ----
ax = axes[2]
# at T such that kBT = 0.12 Delta (corresponds to 4K for Delta=2.8meV)
kBT = 0.12 * Delta
f_FD = 1.0 / (np.exp(E/kBT) + 1)
pop = dos * f_FD
ax.plot(E/Delta, pop, color='#0F6E56', linewidth=1.8,
        label=r'$N_S(E)\,f(E)$')
ax.fill_between(E/Delta, 0, pop, color='#5DCAA5', alpha=0.4)
ax.axvline(1, color='gray', linestyle=':', linewidth=0.6)
ax.set_xlim(0, 4)
ax.set_ylim(0, max(pop)*1.3)
ax.set_xlabel(r'$E / \Delta$')
ax.set_ylabel(r'thermal QP density (arb.)')
ax.set_title('(c) population $= N_S \\times f$', loc='left', fontsize=10)
ax.legend(loc='upper right', frameon=False)
ax.text(1.4, max(pop)*0.7,
        r'integrated:' '\n' r'$x_{qp}^{th} \approx \sqrt{2\pi k_BT/\Delta}$' '\n' r'$\times e^{-\Delta/k_BT}$',
        fontsize=8.5, color='#0F6E56', va='center')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/claude/fig_qp_generation_population.pdf')
plt.savefig('/home/claude/fig_qp_generation_population.png', dpi=160)
print("Saved Fig A.1: qp generation + DOS + population")
plt.close()


# ============================================================
# Fig A.2 - Kinetic regimes (already exists, regenerate cleaner)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))

# Parameters
R = 1.0 / 100e-9
x_eq = 2.6e-4
G_th = R * x_eq**2
x0_pulse = 1e-2

# Panel (a): time dynamics
ax = axes[0]
t_plot = np.logspace(-9, -4, 500)
# (i) hyperbolic
x_hyp = x0_pulse / (1 + R*x0_pulse*t_plot)
ax.loglog(t_plot*1e6, x_hyp, color='#534AB7', linewidth=1.8,
          label='no trap, no gen.  $\\sim 1/t$')
# (ii) trap-dominated linear
s_trap = 1e7
x_ss_lin = G_th / s_trap
x_exp = x_ss_lin + (x0_pulse - x_ss_lin) * np.exp(-s_trap*t_plot)
ax.loglog(t_plot*1e6, x_exp, color='#BA7517', linewidth=1.8,
          label='trap only  $\\sim e^{-st}$', linestyle='--')
# (iii) full Riccati
x_full = odeint(lambda x, t: G_th - s_trap*x - R*x**2,
                 x0_pulse, t_plot).flatten()
ax.loglog(t_plot*1e6, x_full, color='#0F6E56', linewidth=1.8,
          label='full: $G - sn - Rn^2$')
# (iv) recomb + gen, no trap
x_rg = odeint(lambda x, t: G_th - R*x**2, x0_pulse, t_plot).flatten()
ax.loglog(t_plot*1e6, x_rg, color='#993C1D', linewidth=1.8,
          label='$R + G$ only', linestyle='-.')

ax.axhline(x_eq, color='gray', linestyle=':', linewidth=0.8)
ax.text(1e-3, x_eq*1.3, r'$x_{qp}^{th}(4\,\mathrm{K})$',
        color='gray', fontsize=9)
ax.set_xlabel(r'time ($\mu$s)')
ax.set_ylabel(r'$x_{qp}(t)$')
ax.set_xlim(1e-3, 100)
ax.set_ylim(1e-10, 0.1)
ax.legend(fontsize=8.5, loc='lower left', framealpha=0.95)
ax.grid(True, which='both', alpha=0.3)
ax.set_title('(a) dynamic decay regimes', loc='left', fontsize=10)

# Panel (b): steady state vs s
ax = axes[1]
s_range = np.logspace(2, 11, 200)

def x_ss(G, s, R):
    return (-s + np.sqrt(s**2 + 4*R*G)) / (2*R)

x_ss_arr = np.array([x_ss(G_th, s, R) for s in s_range])
ax.loglog(s_range, x_ss_arr, color='#185FA5', linewidth=2,
          label=r'full: $\frac{-s + \sqrt{s^2+4RG}}{2R}$')
ax.loglog(s_range, G_th/s_range, color='#BA7517', linestyle='--',
          linewidth=1.5, label=r'trap-limited: $G/s$')
ax.axhline(np.sqrt(G_th/R), color='#993C1D', linestyle='--',
           linewidth=1.5, label=r'recomb-limited: $\sqrt{G/R}$')
ax.axhline(x_eq, color='gray', linestyle=':', linewidth=0.8)
ax.text(1e3, x_eq*1.5, r'$x_{qp}^{th}$', color='gray', fontsize=9)

s_crossover = 2*np.sqrt(R*G_th)
ax.axvline(s_crossover, color='black', linestyle=':', linewidth=0.7,
           alpha=0.6)
ax.annotate(r'$s^* = 2\sqrt{RG}$',
            xy=(s_crossover, 5e-11),
            xytext=(s_crossover*8, 1.5e-10),
            fontsize=9, color='black')

ax.set_xlabel(r'trap rate $s$ (s$^{-1}$)')
ax.set_ylabel(r'$x_{qp}^{ss}$')
ax.legend(fontsize=8.5, loc='lower left', framealpha=0.95)
ax.grid(True, which='both', alpha=0.3)
ax.set_title('(b) steady state across crossover', loc='left', fontsize=10)
ax.set_xlim(1e2, 1e11)
ax.set_ylim(1e-11, 1e-2)

plt.tight_layout()
plt.savefig('/home/claude/fig_kinetic_regimes.pdf')
plt.savefig('/home/claude/fig_kinetic_regimes.png', dpi=160)
print("Saved Fig A.2: kinetic regimes")
plt.close()


# ============================================================
# Fig A.3 - Normal-metal trap mechanism (energy + spatial)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))

# Panel (a): energy diagram - SC gap + NM continuum
ax = axes[0]
ax.set_xlim(0, 4)
ax.set_ylim(-3, 3)

# SC side (left half): condensate + gap + QP band
ax.add_patch(Rectangle((0.2, -3.0), 1.5, 1.5, facecolor='#B5D4F4',
                        edgecolor='#185FA5', linewidth=0.8))
ax.text(0.95, -2.25, 'condensate', ha='center', fontsize=8.5)
ax.add_patch(Rectangle((0.2, -1.0), 1.5, 2.0, facecolor='#F1EFE8',
                        edgecolor='gray', linewidth=0.5,
                        linestyle='--', alpha=0.4))
ax.text(0.95, 0, r'gap $2\Delta_S$', ha='center', fontsize=8.5,
        style='italic', color='#5F5E5A')
ax.add_patch(Rectangle((0.2, 1.0), 1.5, 1.5, facecolor='#F0997B',
                        edgecolor='#993C1D', linewidth=0.8))
ax.text(0.95, 1.75, 'QP band', ha='center', fontsize=8.5)
ax.text(0.95, -2.8, 'Superconductor (NbN)', ha='center', fontsize=9,
        fontweight='bold')

# NM side (right half): continuous DOS, no gap
ax.add_patch(Rectangle((2.3, -3.0), 1.5, 6.0, facecolor='#FAEEDA',
                        edgecolor='#BA7517', linewidth=0.8, alpha=0.6))
ax.text(3.05, 0, 'continuous\nDOS\n(no gap)',
        ha='center', va='center', fontsize=8.5,
        color='#854F0B', style='italic')
ax.text(3.05, -2.8, 'Normal metal (Cu)', ha='center', fontsize=9,
        fontweight='bold')

# Barrier
ax.add_patch(Rectangle((1.7, -3.0), 0.6, 6.0, facecolor='#F1EFE8',
                        edgecolor='gray', linewidth=0.5, hatch='///',
                        alpha=0.3))
ax.text(2.0, 0, 'barrier', ha='center', va='center', fontsize=8,
        rotation=90, color='#5F5E5A')

# Tunneling: QP enters NM at top
ax.add_patch(Circle((1.55, 1.6), 0.08, facecolor='#993C1D'))
ax.annotate('', xy=(2.45, 1.6), xytext=(1.65, 1.6),
            arrowprops=dict(arrowstyle='->', color='#993C1D',
                            lw=1.3, linestyle='--'))
ax.add_patch(Circle((2.5, 1.6), 0.08, facecolor='#993C1D'))
ax.text(2.0, 2.1, 'tunneling', ha='center', fontsize=8, color='#993C1D')

# Relaxation cascade in NM
for y_start, y_end in [(1.6, 1.0), (1.0, 0.5), (0.5, 0.0), (0.0, -0.6)]:
    ax.annotate('', xy=(2.5, y_end), xytext=(2.5, y_start),
                arrowprops=dict(arrowstyle='->', color='#993C1D',
                                lw=0.9, alpha=0.7))
ax.add_patch(Circle((2.5, -0.7), 0.08, facecolor='#993C1D', alpha=0.5))
ax.text(3.0, 0.6, 'inelastic\n$e$-$e$ relax.', fontsize=8,
        color='#993C1D', alpha=0.85)

# Cannot return (X mark)
ax.annotate('', xy=(1.55, -0.7), xytext=(2.45, -0.7),
            arrowprops=dict(arrowstyle='->', color='gray',
                            lw=0.8, linestyle=':', alpha=0.5))
ax.text(2.0, -1.1, 'return blocked\nby gap', ha='center', fontsize=7.5,
        color='#5F5E5A', style='italic')

# Energy axis
ax.annotate('', xy=(0.05, 2.7), xytext=(0.05, -3.0),
            arrowprops=dict(arrowstyle='-|>', color='black', lw=0.9))
ax.text(0.05, 2.9, r'$E$', fontsize=10, ha='center')
ax.set_xticks([])
ax.set_yticks([])
ax.set_ylim(-3, 3.5)
ax.set_title('(a) energy-diagram view of trapping',
             loc='left', fontsize=10, pad=10)
for spine in ax.spines.values():
    spine.set_visible(False)

# Panel (b): spatial view - SC film with NM pad on top
ax = axes[1]
ax.set_xlim(0, 6)
ax.set_ylim(0, 3)
ax.set_aspect('equal')

# SC film
ax.add_patch(Rectangle((0.3, 0.5), 5.4, 0.9, facecolor='#9FE1CB',
                        edgecolor='#0F6E56', linewidth=0.8))
ax.text(1.5, 0.95, 'NbN film (thickness $d_S$)',
        fontsize=9, va='center')

# NM trap on top of part of SC
ax.add_patch(Rectangle((3.4, 1.4), 2.0, 0.6, facecolor='#F5C4B3',
                        edgecolor='#993C1D', linewidth=0.8))
ax.text(4.4, 1.7, 'normal-metal trap', ha='center', fontsize=9,
        va='center')

# Tunnel barrier indication
ax.plot([3.4, 5.4], [1.4, 1.4], color='gray', linewidth=1, alpha=0.6)
ax.text(5.6, 1.4, 'tunnel', fontsize=7.5, color='#5F5E5A',
        va='center')
ax.text(5.6, 1.32, 'contact', fontsize=7.5, color='#5F5E5A',
        va='top')

# QPs diffusing in SC
qp_positions = [(0.7, 0.95), (1.4, 0.85), (2.0, 1.0), (2.6, 0.9),
                (3.2, 0.95), (3.7, 1.1)]
for x, y in qp_positions:
    ax.add_patch(Circle((x, y), 0.07, facecolor='#0F6E56'))

# Diffusion arrows
ax.annotate('', xy=(1.3, 0.95), xytext=(0.8, 0.95),
            arrowprops=dict(arrowstyle='->', color='#0F6E56',
                            lw=0.7, alpha=0.5))
ax.annotate('', xy=(2.4, 0.95), xytext=(2.05, 0.95),
            arrowprops=dict(arrowstyle='->', color='#0F6E56',
                            lw=0.7, alpha=0.5))
ax.annotate('', xy=(3.5, 1.1), xytext=(3.2, 1.0),
            arrowprops=dict(arrowstyle='->', color='#0F6E56',
                            lw=0.7, alpha=0.5))

# Tunneling event into trap
ax.annotate('', xy=(4.0, 1.5), xytext=(3.8, 1.35),
            arrowprops=dict(arrowstyle='->', color='#993C1D',
                            lw=1.2, linestyle='--'))
ax.add_patch(Circle((4.0, 1.6), 0.07, facecolor='#993C1D'))
ax.add_patch(Circle((4.5, 1.7), 0.06, facecolor='#993C1D', alpha=0.7))
ax.add_patch(Circle((4.9, 1.75), 0.06, facecolor='#993C1D', alpha=0.5))

# Length scales
ax.annotate('', xy=(3.4, 0.3), xytext=(0.3, 0.3),
            arrowprops=dict(arrowstyle='<->', color='black', lw=0.8))
ax.text(1.85, 0.18, r'diffusion length $L_D = \sqrt{D/s_{eff}}$',
        ha='center', fontsize=8)

ax.annotate('', xy=(5.4, 2.4), xytext=(3.4, 2.4),
            arrowprops=dict(arrowstyle='<->', color='black', lw=0.8))
ax.text(4.4, 2.55, r'trap area $A_{trap}$',
        ha='center', fontsize=8)

ax.set_xticks([])
ax.set_yticks([])
ax.set_title('(b) spatial view: trap on the qubit pad',
             loc='left', fontsize=10)
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig('/home/claude/fig_trap_mechanism.pdf')
plt.savefig('/home/claude/fig_trap_mechanism.png', dpi=160)
print("Saved Fig A.3: trap mechanism")
plt.close()


# ============================================================
# Fig A.4 - Trap engineering layout for the asymmetric SQUID transmon
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))

# Panel (a): qubit layout schematic
ax = axes[0]
ax.set_xlim(0, 10)
ax.set_ylim(0, 7)
ax.set_aspect('equal')

# Large pads of the transmon (capacitor)
ax.add_patch(Rectangle((0.5, 2.0), 3.0, 3.0, facecolor='#9FE1CB',
                        edgecolor='#0F6E56', linewidth=1.0))
ax.add_patch(Rectangle((6.5, 2.0), 3.0, 3.0, facecolor='#9FE1CB',
                        edgecolor='#0F6E56', linewidth=1.0))
ax.text(2.0, 4.7, 'NbN pad 1', ha='center', fontsize=9)
ax.text(8.0, 4.7, 'NbN pad 2', ha='center', fontsize=9)

# SQUID loop (small) in the middle
loop_xs = [3.5, 6.5, 6.5, 3.5, 3.5]
loop_ys = [3.7, 3.7, 3.3, 3.3, 3.7]
ax.plot(loop_xs, loop_ys, color='#0F6E56', linewidth=1.2)
ax.text(5.0, 3.0, 'asymmetric SQUID loop', ha='center', fontsize=8.5,
        style='italic', color='#5F5E5A')

# The two junctions (asymmetric: bigger one on left)
ax.add_patch(Rectangle((3.45, 3.45), 0.15, 0.30, facecolor='#993C1D',
                        edgecolor='#993C1D'))
ax.text(3.5, 4.0, r'$E_{J1}$', ha='center', fontsize=9,
        color='#993C1D')

ax.add_patch(Rectangle((6.40, 3.55), 0.10, 0.20, facecolor='#993C1D',
                        edgecolor='#993C1D'))
ax.text(6.5, 4.0, r'$E_{J2}$' '\n' r'(smaller)', ha='center',
        fontsize=8.5, color='#993C1D')

# Normal-metal traps on both pads
ax.add_patch(Rectangle((1.0, 2.4), 1.4, 0.5, facecolor='#F5C4B3',
                        edgecolor='#993C1D', linewidth=1.0, alpha=0.85))
ax.text(1.7, 2.65, 'NM trap', ha='center', fontsize=8.5,
        color='#712B13')

ax.add_patch(Rectangle((7.6, 2.4), 1.4, 0.5, facecolor='#F5C4B3',
                        edgecolor='#993C1D', linewidth=1.0, alpha=0.85))
ax.text(8.3, 2.65, 'NM trap', ha='center', fontsize=8.5,
        color='#712B13')

# Distance label
ax.annotate('', xy=(3.5, 2.0), xytext=(2.4, 2.0),
            arrowprops=dict(arrowstyle='<->', color='black', lw=0.7))
ax.text(2.95, 1.8, r'$L$', ha='center', fontsize=9)
ax.text(2.95, 1.45, 'trap-to-junction', ha='center', fontsize=7.5,
        color='#5F5E5A')

# External flux
ax.annotate('', xy=(5.0, 5.6), xytext=(5.0, 6.5),
            arrowprops=dict(arrowstyle='->', color='#185FA5', lw=1.2))
ax.text(5.0, 6.7, r'$\Phi_{ext}$', ha='center', fontsize=10,
        color='#185FA5')

# Capacitor symbols
ax.text(0.2, 3.5, 'C', ha='center', fontsize=11, fontweight='bold',
        color='#0F6E56')
ax.text(9.8, 3.5, 'C', ha='center', fontsize=11, fontweight='bold',
        color='#0F6E56')

# Bottom annotation
ax.text(5.0, 0.5, r'asymmetric SQUID transmon with bilateral normal-metal traps',
        ha='center', fontsize=9, style='italic', color='#5F5E5A')

ax.set_xticks([])
ax.set_yticks([])
ax.set_title('(a) qubit layout', loc='left', fontsize=10)
for spine in ax.spines.values():
    spine.set_visible(False)

# Panel (b): trap rate s vs trap area
ax = axes[1]
# Riwar-Catelani 2016 scaling: s ~ Gamma_tr * V_trap / V_SC
# in the small-trap regime, s saturates at large area to D/L^2 (diffusion-limited)
A_range = np.logspace(-2, 2, 200)  # in μm^2

# Two regimes: linear in A (tunneling-limited) until saturation (diffusion-limited)
A_sat = 5.0  # μm^2
s_max = 1e8  # diffusion-limited ceiling
gamma_per_area = s_max / A_sat
s_vals = s_max * A_range / (A_range + A_sat)

# Multiple curves for different normal metals
for mat, color, sat in [('Cu (clean)', '#185FA5', 1e8),
                          ('Au (clean)', '#BA7517', 5e7),
                          ('Pd (dirty)', '#993C1D', 1e7)]:
    s_vals = sat * A_range / (A_range + 5.0)
    ax.loglog(A_range, s_vals, linewidth=1.8, color=color, label=mat)

# Target band
ax.axhspan(5e6, 5e7, color='#5DCAA5', alpha=0.25,
           label=r'target band $s \sim 10^7$ s$^{-1}$')
ax.axhline(1e7, color='#0F6E56', linestyle='--', linewidth=1.2,
           label=r'target for $T_1 = 49\,\mu$s')

ax.set_xlabel(r'trap area $A_{trap}$ ($\mu$m$^2$)')
ax.set_ylabel(r'effective trap rate $s$ (s$^{-1}$)')
ax.set_title('(b) trap rate scaling with area', loc='left', fontsize=10)
ax.legend(fontsize=8.5, loc='lower right', framealpha=0.95)
ax.grid(True, which='both', alpha=0.3)
ax.set_xlim(1e-2, 1e2)
ax.set_ylim(1e4, 3e8)

plt.tight_layout()
plt.savefig('/home/claude/fig_trap_engineering.pdf')
plt.savefig('/home/claude/fig_trap_engineering.png', dpi=160)
print("Saved Fig A.4: trap engineering layout")
plt.close()

print("\nAll figures generated.")
