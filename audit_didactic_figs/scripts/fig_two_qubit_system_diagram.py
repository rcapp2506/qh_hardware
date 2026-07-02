"""
Variant C of Fig 2.2: schematic only, no side table.
Key numbers in floating boxes near each element. Wider aspect ratio.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'figure.dpi': 120,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.08,
})

fig = plt.figure(figsize=(10, 5.5))
ax = fig.add_axes([0.04, 0.06, 0.92, 0.88])
ax.set_xlim(0, 12); ax.set_ylim(0, 6.5); ax.set_aspect('equal'); ax.axis('off')

# Cavity (centered)
cavity_x, cavity_y = 6.0, 3.0
cavity_w, cavity_h = 2.6, 1.6
cav = FancyBboxPatch((cavity_x - cavity_w/2, cavity_y - cavity_h/2),
                      cavity_w, cavity_h, boxstyle="round,pad=0.08",
                      linewidth=1.8, edgecolor='#1F3F6B', facecolor='#DBE7F7')
ax.add_patch(cav)
ax.text(cavity_x, cavity_y, 'Cavity', ha='center', va='center',
        fontsize=14, fontweight='bold', color='#1F3F6B')

# Cavity floating label below
ax.text(cavity_x, cavity_y - 1.20,
        r'$\omega_r/2\pi = 6.5\,$GHz',
        ha='center', va='center', fontsize=10.5, color='#1F3F6B')
ax.text(cavity_x, cavity_y - 1.55,
        r'$\kappa/2\pi = 1.0\,$MHz, $Q = 6500$',
        ha='center', va='center', fontsize=10, color='#1F3F6B')

# Q1
q1_x, q1_y = 2.0, 3.0
q1 = Ellipse((q1_x, q1_y), 1.8, 1.3, linewidth=1.8,
             edgecolor='#923124', facecolor='#F5D9D5')
ax.add_patch(q1)
ax.text(q1_x, q1_y, 'Q1', ha='center', va='center',
        fontsize=15, fontweight='bold', color='#923124')

# Floating box for Q1 parameters (below ellipse)
ax.text(q1_x, q1_y - 1.20,
        r'$\omega_{q1}/2\pi = 5.8\,$GHz',
        ha='center', va='center', fontsize=10.5, color='#923124')
ax.text(q1_x, q1_y - 1.55,
        r'$\Delta_1/2\pi = -0.7\,$GHz',
        ha='center', va='center', fontsize=10, color='#923124')
ax.text(q1_x, q1_y - 1.85,
        r'$\alpha_1/2\pi = -0.20\,$GHz',
        ha='center', va='center', fontsize=10, color='#923124')

# Q2
q2_x, q2_y = 10.0, 3.0
q2 = Ellipse((q2_x, q2_y), 1.8, 1.3, linewidth=1.8,
             edgecolor='#1A5F3F', facecolor='#D6EBDD')
ax.add_patch(q2)
ax.text(q2_x, q2_y, 'Q2', ha='center', va='center',
        fontsize=15, fontweight='bold', color='#1A5F3F')

ax.text(q2_x, q2_y - 1.20,
        r'$\omega_{q2}/2\pi = 5.2\,$GHz',
        ha='center', va='center', fontsize=10.5, color='#1A5F3F')
ax.text(q2_x, q2_y - 1.55,
        r'$\Delta_2/2\pi = -1.3\,$GHz',
        ha='center', va='center', fontsize=10, color='#1A5F3F')
ax.text(q2_x, q2_y - 1.85,
        r'$\alpha_2/2\pi = -0.20\,$GHz',
        ha='center', va='center', fontsize=10, color='#1A5F3F')

# g1 coupling arrow
ax.annotate('', xy=(cavity_x - cavity_w/2, cavity_y),
            xytext=(q1_x + 0.9, q1_y),
            arrowprops=dict(arrowstyle='<->', color='#923124', lw=2.0))
ax.text((q1_x + 0.9 + cavity_x - cavity_w/2)/2, cavity_y + 0.30,
        r'$g_1/2\pi = 80\,$MHz', ha='center', fontsize=10, color='#923124')

# g2 coupling arrow
ax.annotate('', xy=(q2_x - 0.9, q2_y),
            xytext=(cavity_x + cavity_w/2, cavity_y),
            arrowprops=dict(arrowstyle='<->', color='#1A5F3F', lw=2.0))
ax.text((q2_x - 0.9 + cavity_x + cavity_w/2)/2, cavity_y + 0.30,
        r'$g_2/2\pi = 80\,$MHz', ha='center', fontsize=10, color='#1A5F3F')

# Effective J coupling (curved dashed arrow above)
arc = FancyArrowPatch(posA=(q1_x, q1_y + 0.75), posB=(q2_x, q2_y + 0.75),
                      connectionstyle="arc3,rad=-0.4",
                      arrowstyle='<->', color='#5E1F9E',
                      linewidth=1.8, linestyle='--')
ax.add_patch(arc)
ax.text(cavity_x, 4.75, 'effective transverse exchange',
        ha='center', va='center', fontsize=10.5, style='italic', color='#5E1F9E')
ax.text(cavity_x, 4.40, r'$|J|/2\pi \simeq 7.0\,$MHz',
        ha='center', va='center', fontsize=11, color='#5E1F9E')

# Cross-Kerr label
ax.text(cavity_x, 4.10,
        r'(residual cross-Kerr $|\xi_{ZZ}|/2\pi \simeq 106\,$kHz)',
        ha='center', va='center', fontsize=9, style='italic', color='#5E1F9E')

# Readout port
ax.annotate('', xy=(cavity_x, 0.55), xytext=(cavity_x, cavity_y - cavity_h/2 - 0.05),
            arrowprops=dict(arrowstyle='->', color='#5F5E5A', lw=1.2))
ax.text(cavity_x, 0.25, 'dispersive readout', ha='center', fontsize=9.5,
        color='#5F5E5A', style='italic')

plt.savefig('/tmp/fig1_system_schematic_clean.pdf')
plt.savefig('/tmp/fig1_system_schematic_clean.png', dpi=160)
print("Saved PDF and PNG")
