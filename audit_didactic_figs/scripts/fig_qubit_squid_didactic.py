"""
fig_qubit_squid_didactic.py
===========================

Generate the didactic figure for §Preliminaries of Cap. 2:
- Panel (a): LC oscillator -> transmon (anharmonic) energy levels
- Panel (b): SQUID flux loop with two junctions, effective E_J(Phi)

Designed for pedagogical clarity. Output: PDF for thesis inclusion.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle, FancyBboxPatch
from matplotlib.lines import Line2D

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


fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))

# ============================================================
# Panel (a): LC oscillator vs transmon energy levels
# ============================================================
ax = axes[0]
ax.set_xlim(0, 10)
ax.set_ylim(-0.5, 7)
ax.set_aspect('auto')

# LC oscillator (left half): equally-spaced potential and levels
phi_lc = np.linspace(-1.5, 1.5, 200)
V_lc = 0.5 * phi_lc**2 + 1.5
ax.plot(phi_lc + 2.0, V_lc, color='#5F5E5A', linewidth=1.2)
# Equally-spaced LC levels
for n in range(5):
    E_n = 0.7 * n + 1.7
    # only draw inside the parabola
    width = 2.0 * np.sqrt(2.0 * (E_n - 1.5)) if E_n > 1.5 else 0
    if width > 0:
        ax.plot([2.0 - width/2, 2.0 + width/2], [E_n, E_n],
                color='#185FA5', linewidth=1.5)
    ax.text(2.0 + width/2 + 0.1, E_n, f'$|{n}\\rangle$',
            fontsize=9, color='#185FA5', va='center')
# LC label
ax.text(2.0, 0.1, 'LC oscillator', ha='center', fontsize=10,
        fontweight='bold')
ax.text(2.0, -0.25, r'$\omega_{n+1,n} = \omega_r$ (equal spacing)',
        ha='center', fontsize=8.5, style='italic', color='#5F5E5A')

# Transmon (right half): cosine potential, anharmonic levels
phi_tr = np.linspace(-np.pi, np.pi, 300)
V_tr = -1.6 * np.cos(phi_tr * 0.95) + 1.5 + 0.18 * phi_tr**2
# offset to right half
ax.plot(phi_tr + 7.5, V_tr, color='#5F5E5A', linewidth=1.2)
# Anharmonic levels: spacing decreases with n
omega_01 = 0.85
alpha = 0.20  # anharmonicity magnitude
E_levels = []
E = 0.5
for n in range(4):
    E += omega_01 - n * alpha
    E_levels.append(E)
for n, E_n in enumerate(E_levels):
    if E_n > V_tr.min() + 0.1:
        # find width at this energy in the cosine potential
        valid = phi_tr[V_tr <= E_n]
        if len(valid) > 0:
            width = valid.max() - valid.min()
            cen = (valid.max() + valid.min()) / 2 + 7.5
            ax.plot([cen - width/2, cen + width/2], [E_n, E_n],
                    color='#993C1D', linewidth=1.5)
            ax.text(cen + width/2 + 0.15, E_n, f'$|{n+1}\\rangle$' if n>0 else r'$|1\rangle$',
                    fontsize=9, color='#993C1D', va='center')
# Ground state
ax.plot([7.5 - 0.5, 7.5 + 0.5], [0.5, 0.5],
        color='#993C1D', linewidth=1.5)
ax.text(7.5 + 0.6, 0.5, r'$|0\rangle$', fontsize=9, color='#993C1D', va='center')

ax.text(7.5, 0.1, 'transmon', ha='center', fontsize=10, fontweight='bold')
ax.text(7.5, -0.25, r'$\alpha = \omega_{21}-\omega_{10} \approx -E_C$',
        ha='center', fontsize=8.5, style='italic', color='#5F5E5A')

# Anharmonicity annotation
ax.annotate('', xy=(8.8, E_levels[1]), xytext=(8.8, E_levels[0]),
            arrowprops=dict(arrowstyle='<->', color='#993C1D', lw=0.8))
ax.annotate('', xy=(8.8, E_levels[0]), xytext=(8.8, 0.5),
            arrowprops=dict(arrowstyle='<->', color='#993C1D', lw=0.8))
ax.text(9.05, (E_levels[0]+0.5)/2, r'$\omega_{10}$', fontsize=8,
        color='#993C1D', va='center')
ax.text(9.05, (E_levels[0]+E_levels[1])/2, r'$\omega_{21}$', fontsize=8,
        color='#993C1D', va='center')

# Separating vertical line
ax.axvline(5.0, color='gray', linestyle=':', linewidth=0.5)

ax.set_xticks([])
ax.set_yticks([])
ax.set_title('(a) anharmonicity makes the transmon a qubit',
             loc='left', fontsize=10)
ax.set_ylabel('energy', fontsize=9)
for spine in ['top', 'right', 'bottom']:
    ax.spines[spine].set_visible(False)


# ============================================================
# Panel (b): SQUID flux loop schematic
# ============================================================
ax = axes[1]
ax.set_xlim(0, 10)
ax.set_ylim(0, 7)
ax.set_aspect('equal')

# SQUID loop: two parallel branches with one junction each
# Outer rectangle (the loop)
loop_left, loop_right = 2.0, 7.0
loop_top, loop_bot = 5.5, 2.5

# Horizontal wires (top and bottom)
ax.plot([loop_left, loop_right], [loop_top, loop_top],
        color='#0F6E56', linewidth=2)
ax.plot([loop_left, loop_right], [loop_bot, loop_bot],
        color='#0F6E56', linewidth=2)
# Vertical wires (sides) — these contain the junctions
# Left side: junction E_J1 in the middle
ax.plot([loop_left, loop_left], [loop_bot, 3.5], color='#0F6E56', linewidth=2)
ax.plot([loop_left, loop_left], [4.5, loop_top], color='#0F6E56', linewidth=2)
# Junction E_J1: "X" symbol or two crosses
jj1_y = 4.0
ax.add_patch(Rectangle((loop_left-0.18, jj1_y-0.2), 0.36, 0.4,
                        facecolor='white', edgecolor='#993C1D', linewidth=1.5))
ax.plot([loop_left-0.12, loop_left+0.12], [jj1_y-0.14, jj1_y+0.14],
        color='#993C1D', linewidth=1.5)
ax.plot([loop_left-0.12, loop_left+0.12], [jj1_y+0.14, jj1_y-0.14],
        color='#993C1D', linewidth=1.5)
ax.text(loop_left - 0.45, jj1_y, r'$E_{J1}$',
        fontsize=10, color='#993C1D', ha='right', va='center')

# Right side: junction E_J2 (smaller, to convey asymmetry)
ax.plot([loop_right, loop_right], [loop_bot, 3.7], color='#0F6E56', linewidth=2)
ax.plot([loop_right, loop_right], [4.3, loop_top], color='#0F6E56', linewidth=2)
jj2_y = 4.0
ax.add_patch(Rectangle((loop_right-0.13, jj2_y-0.15), 0.26, 0.3,
                        facecolor='white', edgecolor='#993C1D', linewidth=1.5))
ax.plot([loop_right-0.09, loop_right+0.09], [jj2_y-0.10, jj2_y+0.10],
        color='#993C1D', linewidth=1.5)
ax.plot([loop_right-0.09, loop_right+0.09], [jj2_y+0.10, jj2_y-0.10],
        color='#993C1D', linewidth=1.5)
ax.text(loop_right + 0.45, jj2_y, r'$E_{J2}$',
        fontsize=10, color='#993C1D', ha='left', va='center')
ax.text(loop_right + 0.45, jj2_y - 0.4, '(smaller)',
        fontsize=7.5, color='#993C1D', ha='left',
        va='center', style='italic')

# Flux Phi through the loop: arrow into the page, symbol
flux_x, flux_y = 4.5, 4.0
ax.add_patch(Circle((flux_x, flux_y), 0.35,
                     facecolor='white', edgecolor='#185FA5', linewidth=1.2))
ax.plot([flux_x-0.18, flux_x+0.18], [flux_y-0.18, flux_y+0.18],
        color='#185FA5', linewidth=1.2)
ax.plot([flux_x-0.18, flux_x+0.18], [flux_y+0.18, flux_y-0.18],
        color='#185FA5', linewidth=1.2)
ax.text(flux_x + 0.55, flux_y, r'$\Phi_{\rm ext}$',
        fontsize=11, color='#185FA5', va='center')

# Asymmetry parameter d annotation
ax.text(4.5, 1.7, r'$d = \dfrac{E_{J1}-E_{J2}}{E_{J1}+E_{J2}}$',
        ha='center', fontsize=10, color='#5F5E5A')

# Effective E_J(Phi) annotation
ax.text(4.5, 6.3, r'$E_J^{\rm eff}(\Phi) = (E_{J1}+E_{J2})|\!\cos(\pi\Phi/\Phi_0)|\sqrt{1+d^2\tan^2(\pi\Phi/\Phi_0)}$',
        ha='center', fontsize=9, color='black')
ax.text(4.5, 0.9,
        r'two flux-insensitive sweet spots per period if $d \neq 0$:',
        ha='center', fontsize=9, style='italic', color='#5F5E5A')
ax.text(4.5, 0.45,
        r'upper at $\Phi=0$, lower at $\Phi=\Phi_0/2$ with $E_J^{\rm eff} = (E_{J1}{+}E_{J2})|d|$',
        ha='center', fontsize=9, style='italic', color='#5F5E5A')

# Connection points to external circuit (capacitor)
ax.plot([loop_left-0.7, loop_left], [loop_top, loop_top],
        color='#0F6E56', linewidth=2)
ax.plot([loop_left-0.7, loop_left], [loop_bot, loop_bot],
        color='#0F6E56', linewidth=2)
ax.plot([loop_right, loop_right+0.7], [loop_top, loop_top],
        color='#0F6E56', linewidth=2)
ax.plot([loop_right, loop_right+0.7], [loop_bot, loop_bot],
        color='#0F6E56', linewidth=2)
# Mini capacitor symbols (parallel plates)
ax.plot([loop_left-0.7, loop_left-0.7], [loop_top-0.15, loop_top+0.15],
        color='#0F6E56', linewidth=2.5)
ax.plot([loop_left-0.85, loop_left-0.85], [loop_top-0.15, loop_top+0.15],
        color='#0F6E56', linewidth=2.5)
ax.plot([loop_right+0.7, loop_right+0.7], [loop_top-0.15, loop_top+0.15],
        color='#0F6E56', linewidth=2.5)
ax.plot([loop_right+0.85, loop_right+0.85], [loop_top-0.15, loop_top+0.15],
        color='#0F6E56', linewidth=2.5)
ax.text(loop_left-1.05, loop_top, '$C$', fontsize=10,
        color='#0F6E56', ha='right', va='center')
ax.text(loop_right+1.05, loop_top, '$C$', fontsize=10,
        color='#0F6E56', ha='left', va='center')

ax.set_xticks([])
ax.set_yticks([])
ax.set_title('(b) asymmetric SQUID with external flux $\\Phi_{\\rm ext}$',
             loc='left', fontsize=10)
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig('/home/claude/fig_qubit_squid_didactic.pdf')
plt.savefig('/home/claude/fig_qubit_squid_didactic.png', dpi=160)
print("Saved fig_qubit_squid_didactic.pdf and .png")
