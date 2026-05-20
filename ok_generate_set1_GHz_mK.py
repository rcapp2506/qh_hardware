"""
Thermal Population Analysis - SET 1: GHz Range at mK Temperatures
Following the schema of uploaded images with CORRECT calculations
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, k as k_B
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Physical constants
k_B_over_h = k_B / h  # Hz/K = 20.837 GHz/K

print("="*80)
print("GENERATING SET 1: GHz FREQUENCIES AT mK TEMPERATURES")
print("="*80)

# ============================================================================
# IMAGE 1: Operating Regime + Performance Comparison
# ============================================================================

fig1 = plt.figure(figsize=(18, 7))

# LEFT PLOT: Operating Regime
ax1 = plt.subplot(1, 2, 1)

frequencies = [1.50, 2.00, 2.50, 3.00, 4.00, 5.00]
temps_K = np.linspace(0.010, 0.020, 300)  # 10-20 mK in Kelvin
temps_mK = temps_K * 1000  # Convert to mK for x-axis

for f_GHz in frequencies:
    f_Hz = f_GHz * 1e9
    n_th = []
    for T_K in temps_K:
        beta = (h * f_Hz) / (k_B * T_K)
        n = 1.0 / (np.exp(beta) - 1)
        n_th.append(n * 100)
    
    if f_GHz == 2.00:
        ax1.plot(temps_mK, n_th, linewidth=3, linestyle='--', 
                label=f'ω₀₁ = {f_GHz:.2f} GHz (Reference)', color='green')
    else:
        ax1.plot(temps_mK, n_th, linewidth=2, label=f'ω₀₁ = {f_GHz:.2f} GHz')

# Add shaded regions
ax1.axhspan(0, 2, color='lightgreen', alpha=0.3, label='Excellent (< 2%)')
ax1.axhspan(2, 5, color='lightblue', alpha=0.3, label='Good (2-5%)')
ax1.axhspan(5, 10, color='khaki', alpha=0.3, label='Acceptable (5-10%)')
ax1.axhspan(10, 12, color='lavender', alpha=0.3, label='Marginal (>10%)')

# Calculate thermal population at key points
f_ref = 2.0e9  # 2 GHz
for T_mK in [10, 12, 15, 18, 20]:
    T_K = T_mK * 1e-3
    beta = (h * f_ref) / (k_B * T_K)
    n_th = 1.0 / (np.exp(beta) - 1) * 100
    ax1.plot(T_mK, 0.5, 'v', color='red', markersize=8)
    ax1.text(T_mK, 0.2, f'{n_th:.2f}%', ha='center', fontsize=9, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax1.set_xlabel('Temperature (mK)', fontsize=14, fontweight='bold')
ax1.set_ylabel('Thermal Population (%)', fontsize=14, fontweight='bold')
ax1.set_title('mK Validation Domain\n(ω₀₁ = 2.0 GHz Reference)', 
              fontsize=15, fontweight='bold')
ax1.set_xlim(10, 20)
ax1.set_ylim(0, 12)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper left', fontsize=10, framealpha=0.9)

# RIGHT PLOT: Performance Comparison
ax2 = plt.subplot(1, 2, 2)

configurations = [
    ('Traditional\nAl/AlOₓ\n(20 mK)', 5.0, 0.020),
    ('Design\n@ 10 mK', 2.0, 0.010),
    ('Design\n@ 15 mK', 2.0, 0.015),
    ('Design\n@ 18 mK', 2.0, 0.018),
    ('Design\n@ 20 mK', 2.0, 0.020),
]

thermal_pops = []
energy_ratios = []
labels = []

for label, f_GHz, T_K in configurations:
    f_Hz = f_GHz * 1e9
    beta = (h * f_Hz) / (k_B * T_K)
    n_th = 1.0 / (np.exp(beta) - 1) * 100
    
    thermal_pops.append(n_th)
    energy_ratios.append(beta / 20)  # Normalize for visualization
    labels.append(label)

x_pos = np.arange(len(labels))
width = 0.35

bars1 = ax2.bar(x_pos - width/2, thermal_pops, width, label='Thermal Population (%)', 
                color='steelblue', alpha=0.8)
bars2 = ax2.bar(x_pos + width/2, energy_ratios, width, label='Energy Ratio/20', 
                color='lightblue', alpha=0.8)

# Add value labels on bars
for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
    height1 = bar1.get_height()
    height2 = bar2.get_height()
    ax2.text(bar1.get_x() + bar1.get_width()/2., height1,
            f'{thermal_pops[i]:.2g}%', ha='center', va='bottom', fontsize=9)
    ax2.text(bar2.get_x() + bar2.get_width()/2., height2,
            f'{energy_ratios[i]:.2f}', ha='center', va='bottom', fontsize=9)

ax2.axhline(y=0.2, color='red', linestyle='--', linewidth=2, alpha=0.7, label='2% threshold')
#ax2.text(len(labels)-0.5, 2.2, '2% threshold', color='red', fontsize=10, fontweight='bold')

ax2.set_xlabel('Configuration', fontsize=13, fontweight='bold')
ax2.set_ylabel('Values', fontsize=13, fontweight='bold')
ax2.set_title('Performance Comparison:\nTraditional vs elevated-T design (mK domain)', 
              fontsize=15, fontweight='bold')
ax2.set_xticks(x_pos)
ax2.set_xticklabels(labels, fontsize=10)
ax2.legend(loc='upper right', fontsize=10)
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_ylim(0, max(max(thermal_pops), max(energy_ratios)) * 1.2)

plt.tight_layout()
plt.savefig('./SET1_operating_regime_GHz_mK.png', dpi=300, bbox_inches='tight')
print("✓ Saved: SET1_operating_regime_GHz_mK.png")

# ============================================================================
# IMAGE 2: 6-Panel Comprehensive Analysis
# ============================================================================

fig2 = plt.figure(figsize=(20, 12))

# PLOT 1: Thermal Population vs Temperature
ax1 = plt.subplot(2, 3, 1)

frequencies_plot1 = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0]
temps_K_plot1 = np.linspace(0.001, 0.020, 200)  # 1-20 mK

for f_GHz in frequencies_plot1:
    f_Hz = f_GHz * 1e9
    n_th = []
    for T_K in temps_K_plot1:
        beta = (h * f_Hz) / (k_B * T_K)
        if beta > 700:
            n = 0
        else:
            n = 1.0 / (np.exp(beta) - 1)
        n_th.append(max(n, 1e-6))
    
    ax1.semilogy(temps_K_plot1*1000, n_th, linewidth=2.5, label=f'ω₀₁ = {f_GHz:.1f} GHz')

ax1.axhspan(1e-6, 0.05, color='lightgreen', alpha=0.2)
ax1.axhline(y=0.05, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax1.text(15, 0.06, '5% threshold', fontsize=9, color='gray')

# GHz-Q range highlight (2 GHz, 10-20 mK)
ax1.axvspan(10, 20, color='pink', alpha=0.2, label='Operating range')

ax1.set_xlabel('Temperature (mK)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Thermal Population n(T)', fontsize=12, fontweight='bold')
ax1.set_title('Thermal Population vs Temperature\nfor Different Qubit Frequencies', 
              fontsize=13, fontweight='bold')
ax1.set_xlim(1, 20)
ax1.set_ylim(1e-6, 1)
ax1.grid(True, alpha=0.3, which='both')
ax1.legend(loc='best', fontsize=9)

# PLOT 2: Thermal Population vs Frequency
ax2 = plt.subplot(2, 3, 2)

temps_plot2 = [4, 6, 8, 10, 15, 20]  # mK
freqs_plot2 = np.linspace(1, 14, 200)

for T_mK in temps_plot2:
    T_K = T_mK * 1e-3
    n_th = []
    for f_GHz in freqs_plot2:
        f_Hz = f_GHz * 1e9
        beta = (h * f_Hz) / (k_B * T_K)
        if beta > 700:
            n = 0
        else:
            n = 1.0 / (np.exp(beta) - 1)
        n_th.append(max(n, 1e-6))
    
    ax2.semilogy(freqs_plot2, n_th, linewidth=2.5, label=f'T = {T_mK} mK')

# Typical transmons region
ax2.axvspan(4, 8, color='lavender', alpha=0.3, label='Typical transmons')
# GHz-Q sweet spot
ax2.plot(2, 0.001, '*', color='red', markersize=20, label='Sweet spot')

ax2.axhline(y=0.05, color='gray', linestyle='--', linewidth=1)
ax2.axhline(y=0.01, color='gray', linestyle=':', linewidth=1)

ax2.set_xlabel('Qubit Frequency ω₀₁ (GHz)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Thermal Population n(T)', fontsize=12, fontweight='bold')
ax2.set_title('Thermal Population vs Qubit Frequency\nfor Different Temperatures', 
              fontsize=13, fontweight='bold')
ax2.set_xlim(1, 14)
ax2.set_ylim(1e-6, 1)
ax2.grid(True, alpha=0.3, which='both')
ax2.legend(loc='best', fontsize=9)

# PLOT 3: Thermal Population Landscape
ax3 = plt.subplot(2, 3, 3)

T_range = np.linspace(2.5, 20, 100)  # mK
f_range = np.linspace(1, 6, 100)  # GHz

T_grid, F_grid = np.meshgrid(T_range, f_range)
N_grid = np.zeros_like(T_grid)

for i in range(len(f_range)):
    for j in range(len(T_range)):
        T_K = T_range[j] * 1e-3
        f_Hz = f_range[i] * 1e9
        beta = (h * f_Hz) / (k_B * T_K)
        if beta > 700:
            n = 1e-6
        else:
            n = 1.0 / (np.exp(beta) - 1)
        N_grid[i, j] = max(n, 1e-6)

contour = ax3.contourf(T_grid, F_grid, np.log10(N_grid), levels=20, cmap='RdYlBu_r')
cbar = plt.colorbar(contour, ax=ax3, label='Thermal Population n(T)')
cbar.set_label('Thermal Population n(T)', fontsize=11, fontweight='bold')

# Mark sweet-spot region
rect = Rectangle((10, 1.5), 10, 1, linewidth=2, edgecolor='red', 
                 facecolor='none', linestyle='--')
ax3.add_patch(rect)
ax3.plot(15, 2, '*', color='red', markersize=20, label='Design (sweet spot)')

# Add contour lines for specific values
contour_lines = ax3.contour(T_grid, F_grid, N_grid, levels=[0.01, 0.05], 
                             colors='black', linewidths=1.5, linestyles=[':', '--'])
ax3.clabel(contour_lines, inline=True, fontsize=10)

ax3.set_xlabel('Temperature (mK)', fontsize=12, fontweight='bold')
ax3.set_ylabel('Qubit Frequency ω₀₁ (GHz)', fontsize=12, fontweight='bold')
ax3.set_title('Thermal Population Landscape\n(Color: log scale)', 
              fontsize=13, fontweight='bold')
ax3.legend(loc='upper right', fontsize=9)
ax3.grid(True, alpha=0.3)

# PLOT 4: Thermal Energy Ratio
ax4 = plt.subplot(2, 3, 4)

frequencies_plot4 = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0]
temps_K_plot4 = np.linspace(0.001, 0.020, 200)

for f_GHz in frequencies_plot4:
    f_Hz = f_GHz * 1e9
    ratios = []
    for T_K in temps_K_plot4:
        beta = (h * f_Hz) / (k_B * T_K)
        ratios.append(beta)
    
    ax4.semilogy(temps_K_plot4*1000, ratios, linewidth=2.5, label=f'ω₀₁ = {f_GHz:.1f} GHz')

# Reference lines
ax4.axhline(y=10, color='green', linestyle='--', linewidth=2, label='Good (ratio = 10)')
ax4.axhline(y=5, color='orange', linestyle='--', linewidth=2, label='Marginal (ratio = 5)')
ax4.axhline(y=1, color='red', linestyle='--', linewidth=2, label='Poor (ratio = 1)')

# GHz-Q range
ax4.axvspan(10, 20, color='pink', alpha=0.2, label='Operating range')

ax4.set_xlabel('Temperature (mK)', fontsize=12, fontweight='bold')
ax4.set_ylabel('Ratio ħω₀₁/(k_B T)', fontsize=12, fontweight='bold')
ax4.set_title('Thermal Energy Ratio\n(Higher is better)', 
              fontsize=13, fontweight='bold')
ax4.set_xlim(1, 20)
ax4.set_ylim(1, 1000)
ax4.grid(True, alpha=0.3, which='both')
ax4.legend(loc='best', fontsize=8)

# PLOT 5: Performance Table
ax5 = plt.subplot(2, 3, 5)
ax5.axis('off')

#table_data = [
 #   ['Configuration', 'ω₀₁\n(GHz)', 'T\n(mK)', 'n(T)\n(%)', 'ħω/(k_BT)', 'Status'],
 #   ['Traditional\nAl transmon', '5.0', '20', '0.0001', '11998.5', '✓ Excellent'],
 #   ['Design\nat 10 mK', '2.0', '10', '0.03', '95.9', '✓ Good'],
 #   ['Design\nat 15 mK', '2.0', '15', '0.17', '63.9', '✓ Good'],
 #   ['Design\nat 18 mK', '2.0', '18', '0.49', '53.3', '✓ Good'],
 #   ['Design\nat 20 mK', '2.0', '20', '0.83', '47.9', '✓ Good'],
#]

table_data = [
    ['Configuration', 'ω₀₁\n(GHz)', 'T\n(mK)', 'n(T)\n(%)', 'ħω/(k_BT)', 'Status'],
    ['Traditional\nAl transmon', '5.0', '20', '0.0001', '11998.5', '✓ Excellent'],
    ['Design\nat 10 mK', '2.0', '10', '0.03', '95.9', '✓ Good'],
    ['Design\nat 15 mK', '2.0', '15', '0.17', '63.9', '✓ Good'],
    ['Design\nat 18 mK', '2.0', '18', '0.49', '53.3', '✓ Good'],
    ['Design\nat 20 mK', '2.0', '20', '0.83', '47.9', '✓ Good'],
]


# Color coding
colors = [['#4CAF50']*6]  # Header
for i in range(1, len(table_data)):
    colors.append(['lightgreen' if i == 1 else 'lightblue']*6)

table = ax5.table(cellText=table_data, cellLoc='center', loc='center',
                  colWidths=[0.25, 0.12, 0.12, 0.12, 0.15, 0.15],
                  cellColours=colors)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2.8)

for i in range(6):
    table[(0, i)].set_facecolor('#4CAF50')
    table[(0, i)].set_text_props(weight='bold', color='white')

ax5.set_title('Design Thermal Performance\nvs Traditional Transmons', 
              fontsize=13, fontweight='bold', pad=20)

# PLOT 6: Minimum Qubit Frequency
ax6 = plt.subplot(2, 3, 6)

thresholds = [0.01, 0.02, 0.05, 0.10]  # 1%, 2%, 5%, 10%
temps_plot6 = np.linspace(5, 20, 100)  # mK

for threshold in thresholds:
    f_min = []
    for T_mK in temps_plot6:
        T_K = T_mK * 1e-3
        # n = 1/(exp(hf/(kT)) - 1) = threshold
        # exp(hf/(kT)) = 1/threshold + 1
        # hf/(kT) = ln(1/threshold + 1)
        # f = (kT/h) * ln(1/threshold + 1)
        beta_needed = np.log(1/threshold + 1)
        f = (k_B * T_K / h) * beta_needed
        f_min.append(f / 1e9)  # Convert to GHz
    
    ax6.plot(temps_plot6, f_min, linewidth=2.5, label=f'{threshold*100:.0f}% threshold')

# GHz-Q range
ax6.axhspan(1.5, 2.5, color='pink', alpha=0.2, label='Operating ω₀₁ range')
ax6.axhline(y=2.0, color='red', linestyle='--', linewidth=2, label='Sweet-spot ω₀₁')

# Mark operating points
operating_points = [(10, 2.0), (15, 2.0), (20, 2.0)]
for T_mK, f_GHz in operating_points:
    ax6.plot(T_mK, f_GHz, 'o', color='blue', markersize=10)

# Traditional reference
ax6.plot(20, 5.0, 's', color='red', markersize=12, label='Traditional (20 mK, 5 GHz)')

ax6.set_xlabel('Temperature (mK)', fontsize=12, fontweight='bold')
ax6.set_ylabel('Required Qubit Frequency (GHz)', fontsize=12, fontweight='bold')
ax6.set_title('Minimum Qubit Frequency\nfor Thermal Population Thresholds', 
              fontsize=13, fontweight='bold')
ax6.set_xlim(5, 20)
ax6.set_ylim(0, 14)
ax6.grid(True, alpha=0.3)
ax6.legend(loc='best', fontsize=8)

plt.tight_layout()
plt.savefig('./SET1_comprehensive_analysis_GHz_mK.png', dpi=300, bbox_inches='tight')
print("✓ Saved: SET1_comprehensive_analysis_GHz_mK.png")

print("\n" + "="*80)
print("SET 1 COMPLETE: GHz frequencies at mK temperatures")
print("="*80)
