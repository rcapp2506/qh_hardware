"""
Regenerate fig2_energy_levels.png (Fig 2.6 of Cap.2 of the thesis)
with a zoomed y-axis so that the three energy levels at 5.2 / 5.8 /
6.5 GHz are visually distinguishable. The previous version had y in
[0, 7] GHz with all the relevant content compressed in the upper
1.5 GHz, leaving 5 GHz of empty white space.

Numerical values are taken from the system parameters table
(quantum_hw.tex Tab 2.2):
  Qubit 2: omega_q2/2pi = 5.2 GHz
  Qubit 1: omega_q1/2pi = 5.8 GHz
  Cavity:  omega_r/2pi  = 6.5 GHz
Detunings:
  Delta_1 = omega_q1 - omega_r = -0.7 GHz
  Delta_2 = omega_q2 - omega_r = -1.3 GHz
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.dpi': 120,
    'savefig.bbox': 'tight',
})

omega_q2 = 5.2  # GHz
omega_q1 = 5.8  # GHz
omega_r  = 6.5  # GHz
Delta_1 = omega_q1 - omega_r  # -0.7 GHz
Delta_2 = omega_q2 - omega_r  # -1.3 GHz

fig, ax = plt.subplots(figsize=(8, 5))

# Energy levels
levels = [
    (omega_q2, 'Qubit 2  ' + r'$|1\rangle$  ' + f'({omega_q2} GHz)', '#1A8B3B'),
    (omega_q1, 'Qubit 1  ' + r'$|1\rangle$  ' + f'({omega_q1} GHz)', '#B22222'),
    (omega_r,  'Cavity  '  + r'$|1\rangle$  ' + f'({omega_r} GHz)',  '#1F4E8C'),
]

for energy, label, color in levels:
    ax.hlines(energy, 0.10, 0.90, colors=color, linewidth=5, alpha=0.9)
    ax.text(0.92, energy, label, va='center', fontsize=12,
            color=color, fontweight='bold')

# Detuning arrows (between cavity and each qubit, drawn at x=0.30 and x=0.50 for clarity)
# Delta_1: between Cavity and Qubit 1
ax.annotate('', xy=(0.35, omega_r), xytext=(0.35, omega_q1),
            arrowprops=dict(arrowstyle='<->', color='#FF7700',
                            lw=2.0, ls='--'))
ax.text(0.31, (omega_r + omega_q1)/2,
        r'$\Delta_1 = ' + f'{Delta_1:+.1f}' + r'$ GHz',
        va='center', ha='right', fontsize=11,
        fontweight='bold', color='#FF7700')

# Delta_2: between Cavity and Qubit 2
ax.annotate('', xy=(0.55, omega_r), xytext=(0.55, omega_q2),
            arrowprops=dict(arrowstyle='<->', color='#7E1F9E',
                            lw=2.0, ls='--'))
ax.text(0.59, (omega_r + omega_q2)/2,
        r'$\Delta_2 = ' + f'{Delta_2:+.1f}' + r'$ GHz',
        va='center', ha='left', fontsize=11,
        fontweight='bold', color='#7E1F9E')

# Bare-qubit transition lines on the left (ground state to excited)
# Showing the |g> -> |1> excitation as colored arrows
ax.annotate('', xy=(0.15, omega_q1), xytext=(0.15, 4.9),
            arrowprops=dict(arrowstyle='->', color='#B22222',
                            lw=1.5, alpha=0.5))
ax.annotate('', xy=(0.15, omega_q2), xytext=(0.15, 4.9),
            arrowprops=dict(arrowstyle='->', color='#1A8B3B',
                            lw=1.5, alpha=0.5))

ax.text(0.05, 4.92, r'$|g\rangle$', va='top', ha='left',
        fontsize=12, color='#555555', fontweight='bold')
ax.hlines(4.92, 0.10, 0.90, colors='#555555', linewidth=3, alpha=0.7,
          linestyles='solid')

ax.set_xlim(0, 1.85)
ax.set_ylim(4.85, 6.7)
ax.set_ylabel('Energy (GHz)', fontsize=12, fontweight='bold')
ax.set_title('Energy-level structure of the two-transmon cavity-QED system',
             fontweight='bold')
ax.set_xticks([])
ax.grid(True, alpha=0.3, axis='y', linestyle=':')

plt.tight_layout()
plt.savefig('./fig2_energy_levels.pdf')
plt.savefig('./fig2_energy_levels.png', dpi=160)
print("Saved fig2_energy_levels.pdf and .png")
