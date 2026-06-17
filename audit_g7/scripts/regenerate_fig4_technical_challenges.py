"""
Regenerate fig4_technical_challenges.png (Fig 2.20 in current numbering)
with two readability fixes:

1. Remove the inset 'KEY CHALLENGE / MITIGATION' yellow box that was
   overlapping the leftmost bar of panel (a), and the inset
   'RECOMMENDATION: Sapphire or high-rho Si' green box in panel (b)
   that was overlapping the 'Standard Si' bar. The information they
   carried is in the caption, not in the figure.
2. Simplify the title hierarchy (was 4 nested levels) to a single
   suptitle plus per-panel titles.
3. Drop the 'Innovation' label terminology from the figure (residual
   from the pre-HEATS-Q cleanup).

Numerical content unchanged:
  Panel (a): Dielectric 52x, Radiative 2675x, Quasiparticle 1x
  Panel (b): Standard Si 10 us, High-rho Si 200 us, Sapphire 100 us,
             Target 500 us; minimum threshold 30 us
  Panel (c): Current 97.90%, Better coherence 98.20%, Optimized control 98.25%,
             Thermal ceiling 98.59% (= 100% - eps_thermal at fixed 300 GHz / 4 K);
             threshold 98%. No temperature reduction (4 K is the design floor).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'figure.dpi': 120,
    'savefig.bbox': 'tight',
})

fig = plt.figure(figsize=(13, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.25,
              height_ratios=[1, 1])

# ============================================================
# Panel (a) — Loss mechanism scaling (full width on top)
# ============================================================
ax_a = fig.add_subplot(gs[0, :])

mechanisms = ['Dielectric\n($\\propto \\omega$)',
              'Radiative\n($\\propto \\omega^2$)',
              'Quasiparticle\n(const)']
baseline = [1.0, 1.0, 1.0]
elevated = [52.0, 2675.0, 1.0]

x_pos = np.arange(len(mechanisms))
width = 0.36

bars1 = ax_a.bar(x_pos - width/2, baseline, width,
                  color='#5BAA66', edgecolor='black', linewidth=1.2,
                  label='Baseline (5.8 GHz, 20 mK)')
bars2 = ax_a.bar(x_pos + width/2, elevated, width,
                  color='#D86060', edgecolor='black', linewidth=1.2,
                  label='Elevated-T design (300 GHz, 4 K)')

# Annotations above each bar
for i, (b1, b2, val_b, val_e) in enumerate(zip(bars1, bars2, baseline, elevated)):
    ax_a.text(b1.get_x() + b1.get_width()/2, val_b * 1.5,
              f'{val_b:.0f}$\\times$',
              ha='center', va='bottom', fontsize=10, color='#3A6B43')
    ax_a.text(b2.get_x() + b2.get_width()/2, val_e * 1.5,
              f'{val_e:.0f}$\\times$',
              ha='center', va='bottom', fontsize=10, fontweight='bold',
              color='#8C2828')

ax_a.set_yscale('log')
ax_a.set_xticks(x_pos)
ax_a.set_xticklabels(mechanisms)
ax_a.set_ylabel('Relative loss rate')
ax_a.set_title('(a) Loss-mechanism scaling: baseline vs.\\ elevated-T operating point',
               fontweight='bold')
ax_a.set_ylim(0.5, 8000)
ax_a.legend(loc='upper right', fontsize=10, framealpha=0.95)
ax_a.grid(axis='y', alpha=0.3, linestyle=':')


# ============================================================
# Panel (b) — Substrate comparison
# ============================================================
ax_b = fig.add_subplot(gs[1, 0])

substrates = ['Standard\nSi', 'High-$\\rho$\nSi', 'Sapphire', 'Target']
t1_us = [10, 200, 100, 500]
tan_delta_text = [r'$\tan\delta = 10^{-6}$',
                  r'$\tan\delta = 5\!\times\!10^{-8}$',
                  r'$\tan\delta = 10^{-7}$',
                  r'$\tan\delta = 10^{-8}$']
colors_b = ['#D85959', '#E59A40', '#5BAA66', '#3A6B43']

# Green shaded region above 30 us
ax_b.axhspan(30, 600, color='#D6EBDD', alpha=0.5, zorder=0)

bars = ax_b.bar(range(4), t1_us, color=colors_b, edgecolor='black',
                linewidth=1.4, width=0.65)

# Threshold line
ax_b.axhline(30, color='#D85959', linestyle='--', linewidth=1.8,
             label=r'30 $\mu$s minimum')

# Numerical labels above each bar
for i, (bar, t1, tan_d) in enumerate(zip(bars, t1_us, tan_delta_text)):
    ax_b.text(bar.get_x() + bar.get_width()/2, t1 + 20,
              f'{t1} $\\mu$s\n{tan_d}',
              ha='center', va='bottom', fontsize=9, fontweight='bold')

ax_b.set_xticks(range(4))
ax_b.set_xticklabels(substrates)
ax_b.set_ylabel(r'Projected $T_1$ ($\mu$s)')
ax_b.set_title('(b) Substrate comparison at 300 GHz', fontweight='bold')
ax_b.set_ylim(0, 600)
ax_b.legend(loc='upper left', fontsize=10, framealpha=0.95)
ax_b.grid(axis='y', alpha=0.3, linestyle=':')


# ============================================================
# Panel (c) — Improvement pathways
# ============================================================
ax_c = fig.add_subplot(gs[1, 1])

pathways = ['Current', 'Better\ncoherence', 'Optimized\ncontrol', 'Thermal\nceiling']
fidelities = [97.90, 98.20, 98.25, 98.59]
deltas = ['—', '+0.30%', '+0.35%', 'ceiling']
colors_c = ['#E59A40', '#92C394', '#5BAA66', '#3A6B43']

ax_c.axhspan(98, 100, color='#D6EBDD', alpha=0.5, zorder=0)

bars = ax_c.bar(range(len(fidelities)), fidelities, color=colors_c,
                edgecolor='black', linewidth=1.4, width=0.6)

ax_c.axhline(97.90, color='#D85959', linewidth=1.6, label='Current: 97.90%')
ax_c.axhline(98.0, color='#3A6B43', linestyle='--', linewidth=1.6,
             label='Threshold: 98%')

for i, (bar, f, d) in enumerate(zip(bars, fidelities, deltas)):
    ax_c.text(bar.get_x() + bar.get_width()/2, f + 0.03,
              f'{f:.2f}%\n({d})',
              ha='center', va='bottom', fontsize=8.5, fontweight='bold')

ax_c.set_xticks(range(len(fidelities)))
ax_c.set_xticklabels(pathways, fontsize=9)
ax_c.set_ylabel('CNOT fidelity (%)')
ax_c.set_title('(c) Improvement pathways at fixed 4 K', fontweight='bold')
ax_c.set_ylim(97.0, 99.0)
ax_c.legend(loc='upper left', fontsize=9, framealpha=0.95)
ax_c.grid(axis='y', alpha=0.3, linestyle=':')


plt.savefig('./fig4_technical_challenges.png', dpi=160)
plt.savefig('./fig4_technical_challenges.pdf')
print("Saved ./fig4_technical_challenges.{png,pdf}")
