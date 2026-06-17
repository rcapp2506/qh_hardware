"""
Thermal Population Analysis - SET 2: 300 GHz at 4-10 K Temperatures
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
print("GENERATING SET 2: 300 GHz AT K TEMPERATURES")
print("="*80)

# ============================================================================
# IMAGE 1: Operating Regime + Performance Comparison
# ============================================================================

fig1 = plt.figure(figsize=(18, 7))

# LEFT PLOT: Operating Regime
ax1 = plt.subplot(1, 2, 1)

frequencies = [250, 270, 290, 300, 310, 330, 350]
temps_K = np.linspace(4, 10, 300)

for f_GHz in frequencies:
    f_Hz = f_GHz * 1e9
    n_th = []
    for T_K in temps_K:
        beta = (h * f_Hz) / (k_B * T_K)
        n = 1.0 / (np.exp(beta) - 1)
        n_th.append(n * 100)
    
    if f_GHz == 300:
        ax1.plot(temps_K, n_th, linewidth=3, linestyle='--', 
                label=f'ω₀₁ = {f_GHz} GHz (Target)', color='green')
    else:
        ax1.plot(temps_K, n_th, linewidth=2, label=f'ω₀₁ = {f_GHz} GHz')

# Add shaded regions
ax1.axhspan(0, 2, color='lightgreen', alpha=0.3, label='Excellent (< 2%)')
ax1.axhspan(2, 5, color='lightblue', alpha=0.3, label='Good (2-5%)')
ax1.axhspan(5, 10, color='khaki', alpha=0.3, label='Acceptable (5-10%)')
ax1.axhspan(10, 35, color='lavender', alpha=0.3, label='Marginal (>10%)')

# Calculate thermal population at key points
f_ref = 300e9  # 300 GHz
for T_K in [4, 5, 6, 7, 8]:
    beta = (h * f_ref) / (k_B * T_K)
    n_th = 1.0 / (np.exp(beta) - 1) * 100
    ax1.plot(T_K, 1, 'v', color='red', markersize=8)
    ax1.text(T_K, 0.5, f'{n_th:.2f}%', ha='center', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax1.set_xlabel('Temperature (K)', fontsize=14, fontweight='bold')
ax1.set_ylabel('Thermal Population (%)', fontsize=14, fontweight='bold')
ax1.set_title('300 GHz Qubits Operating Regime\n(ω₀₁ = 300 GHz at 4K Sweet Spot)', 
              fontsize=15, fontweight='bold')
ax1.set_xlim(4, 10)
ax1.set_ylim(0, 35)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper left', fontsize=10, framealpha=0.9)

# RIGHT PLOT: Performance Comparison
ax2 = plt.subplot(1, 2, 2)

configurations = [
    ('Traditional\nAl/AlOₓ\n(20 mK)', 5.8, 0.020),
    ('Design\n@ 4 K', 300, 4.0),
    ('Design\n@ 6 K', 300, 6.0),
    ('Design\n@ 8 K', 300, 8.0),
    ('Design\n@ 10 K', 300, 10.0),
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
    if i == 0:
        ax2.text(bar1.get_x() + bar1.get_width()/2., height1,
                f'{thermal_pops[i]:.5f}%', ha='center', va='bottom', fontsize=9)
    else:
        ax2.text(bar1.get_x() + bar1.get_width()/2., height1,
                f'{thermal_pops[i]:.2f}%', ha='center', va='bottom', fontsize=9)
    ax2.text(bar2.get_x() + bar2.get_width()/2., height2,
            f'{energy_ratios[i]:.1f}', ha='center', va='bottom', fontsize=9)

ax2.axhline(y=2, color='red', linestyle='--', linewidth=2, alpha=0.7, label='2% threshold')
ax2.text(len(labels)-0.5, 2.5, '2% threshold', color='red', fontsize=10, fontweight='bold')

ax2.set_xlabel('Configuration', fontsize=13, fontweight='bold')
ax2.set_ylabel('Values', fontsize=13, fontweight='bold')
ax2.set_title('Performance Comparison:\nTraditional vs elevated-T design (K domain)', 
              fontsize=15, fontweight='bold')
ax2.set_xticks(x_pos)
ax2.set_xticklabels(labels, fontsize=10)
ax2.legend(loc='upper left', fontsize=10)
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_ylim(0, max(max(thermal_pops), max(energy_ratios)) * 1.2)

plt.tight_layout()
plt.savefig('./SET2_operating_regime_300GHz_K.png', dpi=300, bbox_inches='tight')
print("✓ Saved: SET2_operating_regime_300GHz_K.png")

# ============================================================================
# IMAGE 2: 6-Panel Comprehensive Analysis
# ============================================================================

fig2 = plt.figure(figsize=(20, 12))

# PLOT 1: Thermal Population vs Temperature
ax1 = plt.subplot(2, 3, 1)

frequencies_plot1 = [100, 200, 250, 300, 350, 400]
temps_K_plot1 = np.linspace(2.5, 20, 200)

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
    
    ax1.semilogy(temps_K_plot1, n_th, linewidth=2.5, label=f'ω₀₁ = {f_GHz} GHz')

ax1.axhspan(1e-6, 0.05, color='lightgreen', alpha=0.2)
ax1.axhline(y=0.05, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax1.text(15, 0.06, '5% threshold', fontsize=9, color='gray')

# 300GHz-Q range highlight (300 GHz, 4-10 K)
ax1.axvspan(4, 10, color='pink', alpha=0.2, label='Operating range')

ax1.set_xlabel('Temperature (K)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Thermal Population n(T)', fontsize=12, fontweight='bold')
ax1.set_title('Thermal Population vs Temperature\nfor Different Qubit Frequencies', 
              fontsize=13, fontweight='bold')
ax1.set_xlim(2.5, 20)
ax1.set_ylim(1e-6, 1)
ax1.grid(True, alpha=0.3, which='both')
ax1.legend(loc='best', fontsize=9)

# PLOT 2: Thermal Population vs Frequency
ax2 = plt.subplot(2, 3, 2)

temps_plot2 = [4, 6, 8, 10, 15, 20]  # K
freqs_plot2 = np.linspace(100, 500, 200)

for T_K in temps_plot2:
    n_th = []
    for f_GHz in freqs_plot2:
        f_Hz = f_GHz * 1e9
        beta = (h * f_Hz) / (k_B * T_K)
        if beta > 700:
            n = 0
        else:
            n = 1.0 / (np.exp(beta) - 1)
        n_th.append(max(n, 1e-6))
    
    ax2.semilogy(freqs_plot2, n_th, linewidth=2.5, label=f'T = {T_K} K')

# Typical transmons region
ax2.axvspan(4, 8, color='lavender', alpha=0.3, label='Typical transmons')
# 300GHz-Q sweet spot
ax2.plot(300, 0.03, '*', color='red', markersize=20, label='Design (sweet spot)')

ax2.axhline(y=0.05, color='gray', linestyle='--', linewidth=1)
ax2.axhline(y=0.01, color='gray', linestyle=':', linewidth=1)

ax2.set_xlabel('Qubit Frequency ω₀₁ (GHz)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Thermal Population n(T)', fontsize=12, fontweight='bold')
ax2.set_title('Thermal Population vs Qubit Frequency\nfor Different Temperatures', 
              fontsize=13, fontweight='bold')
ax2.set_xlim(100, 500)
ax2.set_ylim(1e-6, 1)
ax2.grid(True, alpha=0.3, which='both')
ax2.legend(loc='best', fontsize=9)

# PLOT 3: Thermal Population Landscape
ax3 = plt.subplot(2, 3, 3)

T_range = np.linspace(2.5, 20, 100)  # K
f_range = np.linspace(100, 500, 100)  # GHz

T_grid, F_grid = np.meshgrid(T_range, f_range)
N_grid = np.zeros_like(T_grid)

for i in range(len(f_range)):
    for j in range(len(T_range)):
        T_K = T_range[j]
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

# Mark 300GHz-Q region
rect = Rectangle((4, 250), 6, 100, linewidth=2, edgecolor='red', 
                 facecolor='none', linestyle='--')
ax3.add_patch(rect)
ax3.plot(7, 300, '*', color='red', markersize=20, label='Design (sweet spot)')

# Add contour lines for specific values
contour_lines = ax3.contour(T_grid, F_grid, N_grid, levels=[0.01, 0.05], 
                             colors='black', linewidths=1.5, linestyles=[':', '--'])
ax3.clabel(contour_lines, inline=True, fontsize=10)

ax3.set_xlabel('Temperature (K)', fontsize=12, fontweight='bold')
ax3.set_ylabel('Qubit Frequency ω₀₁ (GHz)', fontsize=12, fontweight='bold')
ax3.set_title('Thermal Population Landscape\n(Color: log scale)', 
              fontsize=13, fontweight='bold')
ax3.legend(loc='upper right', fontsize=9)
ax3.grid(True, alpha=0.3)

# PLOT 4: Thermal Energy Ratio
ax4 = plt.subplot(2, 3, 4)

frequencies_plot4 = [100, 200, 250, 300, 350, 400]
temps_K_plot4 = np.linspace(2.5, 20, 200)

for f_GHz in frequencies_plot4:
    f_Hz = f_GHz * 1e9
    ratios = []
    for T_K in temps_K_plot4:
        beta = (h * f_Hz) / (k_B * T_K)
        ratios.append(beta)
    
    ax4.semilogy(temps_K_plot4, ratios, linewidth=2.5, label=f'ω₀₁ = {f_GHz} GHz')

# Reference lines
ax4.axhline(y=10, color='green', linestyle='--', linewidth=2, label='Good (ratio = 10)')
ax4.axhline(y=5, color='orange', linestyle='--', linewidth=2, label='Marginal (ratio = 5)')
ax4.axhline(y=1, color='red', linestyle='--', linewidth=2, label='Poor (ratio = 1)')

# 300GHz-Q range
ax4.axvspan(4, 10, color='pink', alpha=0.2, label='Operating range')

ax4.set_xlabel('Temperature (K)', fontsize=12, fontweight='bold')
ax4.set_ylabel('Ratio ħω₀₁/(k_B T)', fontsize=12, fontweight='bold')
ax4.set_title('Thermal Energy Ratio\n(Higher is better)', 
              fontsize=13, fontweight='bold')
ax4.set_xlim(2.5, 20)
ax4.set_ylim(1, 1000)
ax4.grid(True, alpha=0.3, which='both')
ax4.legend(loc='best', fontsize=8)

# PLOT 5: Performance Table
ax5 = plt.subplot(2, 3, 5)
ax5.axis('off')

table_data = [
    ['Configuration', 'ω₀₁\n(GHz)', 'T\n(K)', 'n(T)\n(%)', 'ħω/(k_BT)', 'Status'],
    ['Traditional\nAl transmon', '5.8', '0.020', '0.0001', '13.9', '✓ Excellent'],
    ['Design\nat 4 K', '300', '4', '2.81', '3.60', '✓ Good'],
    ['Design\nat 6 K', '300', '6', '9.98', '2.40', '△ Marginal'],
    ['Design\nat 8 K', '300', '8', '19.81', '1.80', '✗ Poor'],
    ['Design\nat 10 K', '300', '10', '31.06', '1.44', '✗ Poor'],
]

# Color coding
colors = [['#4CAF50']*6]  # Header
colors.append(['lightgreen']*6)  # Traditional
colors.append(['lightblue']*6)   # 4K
colors.append(['khaki']*6)       # 6K
colors.append(['#FFB6C1']*6)     # 8K
colors.append(['#FFB6C1']*6)     # 10K

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
temps_plot6 = np.linspace(2.5, 20, 100)  # K

for threshold in thresholds:
    f_min = []
    for T_K in temps_plot6:
        beta_needed = np.log(1/threshold + 1)
        f = (k_B * T_K / h) * beta_needed
        f_min.append(f / 1e9)  # Convert to GHz
    
    ax6.plot(temps_plot6, f_min, linewidth=2.5, label=f'{threshold*100:.0f}% threshold')

# 300GHz-Q range
ax6.axhspan(250, 350, color='pink', alpha=0.2, label='Operating ω₀₁ range')
ax6.axhline(y=300, color='red', linestyle='--', linewidth=2, label='Sweet-spot ω₀₁')

# Mark operating points
operating_points = [(4, 300), (6, 300), (8, 300), (10, 300)]
for T_K, f_GHz in operating_points:
    ax6.plot(T_K, f_GHz, 'o', color='blue', markersize=10)

# Traditional reference
ax6.plot(0.020*1000, 5.0, 's', color='green', markersize=12, label='Traditional (20 mK, 5 GHz)')

ax6.set_xlabel('Temperature (K)', fontsize=12, fontweight='bold')
ax6.set_ylabel('Required Qubit Frequency (GHz)', fontsize=12, fontweight='bold')
ax6.set_title('Minimum Qubit Frequency\nfor Thermal Population Thresholds', 
              fontsize=13, fontweight='bold')
ax6.set_xlim(2.5, 20)
ax6.set_ylim(0, 500)
ax6.grid(True, alpha=0.3)
ax6.legend(loc='best', fontsize=8)

plt.tight_layout()
plt.savefig('./SET2_comprehensive_analysis_300GHz_K.png', dpi=300, bbox_inches='tight')
print("✓ Saved: SET2_comprehensive_analysis_300GHz_K.png")

print("\n" + "="*80)
print("SET 2 COMPLETE: 300 GHz at K temperatures")
print("="*80)
