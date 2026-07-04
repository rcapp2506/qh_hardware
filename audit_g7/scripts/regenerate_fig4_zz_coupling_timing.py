"""
Regenerate Fig 2.8 (zz_coupling_timing) with a clearer y-axis on
panel (b): linear scale 0-350 ns instead of log scale 10-500.
This preserves panel (a) unchanged.

The original script for this figure was not pushed to qh_hardware, so
this is a reconstruction based on the values rendered in the current
PNG (system parameters from Tab. 2.2 of the thesis).
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'figure.dpi': 120,
    'savefig.bbox': 'tight',
})

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5))

# === Panel (a): eigenstate energy shifts under H_J = J(s+^1 s-^2 + h.c.) ===
# At qubit-qubit degeneracy the single-excitation manifold hybridises into
# (|01> +/- |10>)/sqrt(2), split by +/-|J|; |00> and |11> are unshifted at
# this order. Canonical baseline value |J|/2pi = 7.0 MHz (canonical_params).
J_MHz = 7.0
J_GHz = J_MHz * 1e-3

states = [r'$|00\rangle$',
          r'$\frac{|01\rangle+|10\rangle}{\sqrt{2}}$',
          r'$\frac{|01\rangle-|10\rangle}{\sqrt{2}}$',
          r'$|11\rangle$']
shifts = [0.0, +J_GHz, -J_GHz, 0.0]
colors = ['#1F77B4', '#2CA02C', '#D62728', '#9467BD']

# Background shading for +/- regions
ax_a.axhspan(0, 0.02, facecolor='#E6F2EA', alpha=0.6, zorder=0)
ax_a.axhspan(-0.02, 0, facecolor='#FBE9E7', alpha=0.6, zorder=0)

bars = ax_a.bar(range(4), shifts, color=colors, edgecolor='black', linewidth=0.8, width=0.6)
for i, (rect, sh) in enumerate(zip(bars, shifts)):
    lbl = '0' if sh == 0 else (f'+{J_MHz} MHz' if sh > 0 else f'\u2212{J_MHz} MHz')
    txt_y = sh + (0.001 if sh >= 0 else -0.001)
    ax_a.text(rect.get_x() + rect.get_width()/2, txt_y, lbl,
              ha='center', va='bottom' if sh >= 0 else 'top',
              fontsize=10, fontweight='bold')

ax_a.set_xticks(range(4))
ax_a.set_xticklabels(states, fontsize=12)
ax_a.set_xlabel('Two-qubit eigenstate (qubit\u2013qubit degeneracy)')
ax_a.set_ylabel(r'Energy shift $\pm|J|$ [GHz]')
ax_a.set_title(r'(a) Cavity-mediated transverse exchange $|J|/2\pi = 7.0$ MHz',
               fontsize=12, fontweight='bold')
ax_a.set_ylim(-0.018, 0.018)
ax_a.axhline(0, color='black', linewidth=0.6)
ax_a.grid(axis='y', alpha=0.3, linestyle=':')


# === Panel (b): Gate timing comparison — linear y-axis 0-350 ns ===
operations = ['Single qubit', 'Hadamard', 'CNOT (echo-CR)', 'Readout']
times_ns = [20, 20, 750, 300]
bar_colors = ['#3CB371', '#3CB371', '#FF8C00', '#1F77B4']

bars_b = ax_b.bar(range(4), times_ns, color=bar_colors, edgecolor='black', linewidth=0.8, width=0.6)
for rect, t in zip(bars_b, times_ns):
    ax_b.text(rect.get_x() + rect.get_width()/2, t + 6,
              f'{t} ns',
              ha='center', va='bottom', fontsize=10, fontweight='bold')

# Fast-gate threshold line at 200 ns
ax_b.axhline(200, color='red', linestyle='--', linewidth=1.2,
             label='Fast-gate threshold (200 ns)')

ax_b.set_xticks(range(4))
ax_b.set_xticklabels(operations, fontsize=10)
ax_b.set_xlabel('Operation')
ax_b.set_ylabel('Gate time [ns]')
ax_b.set_title('(b) Gate timing comparison', fontsize=12, fontweight='bold')
ax_b.set_ylim(0, 350)
ax_b.grid(axis='y', alpha=0.3, linestyle=':')
ax_b.legend(loc='upper left', framealpha=0.9, fontsize=9)

plt.tight_layout()
out_pdf = './fig4_zz_coupling_timing_corrected.pdf'
out_png = './fig4_zz_coupling_timing_corrected.png'
plt.savefig(out_pdf)
plt.savefig(out_png, dpi=160)
print(f"Saved: {out_pdf}, {out_png}")
