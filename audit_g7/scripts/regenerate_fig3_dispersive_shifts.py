"""
Regenerate fig3_dispersive_shifts.png with a zoomed y-axis on
panel (b) 'Two-Qubit State Readout', so that the four cavity-frequency
values 6.4973, 6.4986, 6.5014, 6.5027 GHz are visually distinguishable
(they were not in the previous version, where the axis ran 0--6.5 GHz
and the four bars looked identical).

Numerical values: omega_r = 6.5 GHz, chi_1 = -2.03 MHz, chi_2 = -0.66 MHz.
Two-qubit readout: omega_c = omega_r + chi_1 s_z1 + chi_2 s_z2  with s_z = +/-1.
  |00>: omega_r - chi_1 - chi_2 = 6.5 + 2.69e-3 = 6.50269 GHz
  |01>: omega_r - chi_1 + chi_2 = 6.5 + 1.37e-3 = 6.50137 GHz
  |10>: omega_r + chi_1 - chi_2 = 6.5 - 1.37e-3 = 6.49863 GHz
  |11>: omega_r + chi_1 + chi_2 = 6.5 - 2.69e-3 = 6.49731 GHz

(Convention: s_z = +1 for |0>, s_z = -1 for |1>, matching the
labels of the original figure.)
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

omega_r = 6.500  # GHz
chi1 = -2.03e-3  # GHz
chi2 = -0.66e-3  # GHz

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5))

# === Panel (a): single-qubit dispersive shifts ===
# Bare cavity (dashed line)
ax_a.axhline(omega_r, color='#3D7AB8', linestyle='--', linewidth=1.6,
             label='Bare cavity')

# Four bars: Q1 in |0> / |1>, Q2 in |0> / |1>
# When qubit is in |0>: cavity sees +chi (positive shift)
# When qubit is in |1>: cavity sees -chi (negative shift)
# But chi1, chi2 are negative -> |0> -> shift negative, |1> -> shift positive
freqs_a = [omega_r + chi1, omega_r + chi2, omega_r - chi1, omega_r - chi2]
labels_a = [r'Q1', r'Q2', r'Q1', r'Q2']
colors_a = ['#E84A40', '#3CB371', '#E84A40', '#3CB371']
alphas_a = [0.45, 0.45, 0.95, 0.95]

x_pos = [0.7, 1.3, 2.7, 3.3]
group_centers = [1.0, 3.0]

for x, f, col, alph in zip(x_pos, freqs_a, colors_a, alphas_a):
    ax_a.bar(x, f - 6.492, bottom=6.492, color=col, alpha=alph,
             edgecolor='black', linewidth=1.4, width=0.5)
    shift_MHz = (f - omega_r) * 1000
    ax_a.text(x, f + 0.0003, f'{shift_MHz:+.1f} MHz',
              ha='center', va='bottom', fontsize=10, fontweight='bold',
              color=col)

ax_a.set_xticks(group_centers)
ax_a.set_xticklabels([r'$|0\rangle$', r'$|1\rangle$'], fontsize=14)
ax_a.set_xlabel('Qubit State')
ax_a.set_ylabel('Cavity Frequency (GHz)')
ax_a.set_title('(a) Single-Qubit Dispersive Shifts', fontweight='bold')
ax_a.set_ylim(6.492, 6.508)
ax_a.set_xlim(0.2, 3.8)
ax_a.grid(axis='y', alpha=0.3, linestyle=':')

# Legend with proxy patches
from matplotlib.patches import Patch
proxy = [plt.Line2D([0], [0], color='#3D7AB8', linestyle='--', linewidth=1.6),
         Patch(facecolor='#E84A40', alpha=0.7, edgecolor='black'),
         Patch(facecolor='#3CB371', alpha=0.7, edgecolor='black')]
ax_a.legend(proxy, ['Bare cavity', 'Qubit 1', 'Qubit 2'],
            loc='upper left', fontsize=10, framealpha=0.95)


# === Panel (b): two-qubit state readout, with zoomed y-axis ===
# Convention used in the original figure: |00> is the highest frequency
# (both qubits in ground state -> chi-shift positive in both)
states = [r'$|00\rangle$', r'$|01\rangle$', r'$|10\rangle$', r'$|11\rangle$']
freqs_b = [omega_r - chi1 - chi2,   # |00>: 6.50269
           omega_r - chi1 + chi2,   # |01>: 6.50137
           omega_r + chi1 - chi2,   # |10>: 6.49863
           omega_r + chi1 + chi2]   # |11>: 6.49731
shifts_MHz = [(f - omega_r) * 1000 for f in freqs_b]
colors_b = ['#3D7AB8', '#3CB371', '#E84A40', '#9D55D9']

# Bare cavity reference line (dashed grey)
ax_b.axhline(omega_r, color='gray', linestyle='--', linewidth=1.6, alpha=0.6,
             label=f'Bare cavity ({omega_r:.3f} GHz)', zorder=1)

# Bars: plot as deviations from a baseline at omega_r - 0.004
y_baseline = 6.4960
for i, (f, col) in enumerate(zip(freqs_b, colors_b)):
    ax_b.bar(i, f - y_baseline, bottom=y_baseline,
             color=col, alpha=0.85, edgecolor='black', linewidth=1.6,
             width=0.65, zorder=2)

# Annotations
for i, (f, s) in enumerate(zip(freqs_b, shifts_MHz)):
    ax_b.text(i, f + 0.00018, f'{f:.4f} GHz\n({s:+.2f} MHz)',
              ha='center', va='bottom', fontsize=10, fontweight='bold')

ax_b.set_xticks(range(4))
ax_b.set_xticklabels(states, fontsize=14)
ax_b.set_xlabel('Two-Qubit State')
ax_b.set_ylabel('Cavity Frequency (GHz)')
ax_b.set_title('(b) Two-Qubit State Readout', fontweight='bold')
ax_b.set_ylim(6.4960, 6.5045)   # zoomed
ax_b.grid(axis='y', alpha=0.3, linestyle=':')

# Minimum separation annotation
min_sep_MHz = min(abs(shifts_MHz[i+1] - shifts_MHz[i]) for i in range(3))
ax_b.text(0.5, 0.04,
          f'Minimum separation: {min_sep_MHz:.2f} MHz $>$ $\\kappa$ (1 MHz) $\\checkmark$',
          transform=ax_b.transAxes, ha='center', fontsize=10.5,
          bbox=dict(boxstyle='round,pad=0.4',
                    facecolor='#FFF8C5', edgecolor='#9B7B0A',
                    linewidth=1.0))

ax_b.legend(loc='upper right', fontsize=9.5, framealpha=0.95)

plt.tight_layout()
plt.savefig('./fig3_dispersive_shifts.png', dpi=160)
plt.savefig('./fig3_dispersive_shifts.pdf')
print("Done: ./fig3_dispersive_shifts.{png,pdf}")
