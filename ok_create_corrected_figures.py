"""
CORRECTED FIGURES - All values verified
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.constants import h, k as k_B, hbar, e

k_B_freq = k_B / h
Phi_0 = h / (2*e)

# Set style
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'font.family': 'sans-serif',
})

def thermal_occupation(freq, T):
    if T < 0.001:
        return 0.0
    beta = freq / (k_B_freq * T)
    try:
        return 1.0 / (np.exp(beta) - 1)
    except:
        return 0.0

# Load corrected results
I_c = 1.510e-6  # A (CORRECTED)
T1_total = 15.3e-6  # s (CORRECTED)
T2_total = 13.2e-6  # s (CORRECTED)
F_single = 0.99813  # CORRECTED
F_CNOT = 0.9595  # CORRECTED

print("Creating corrected figures with verified values...")
print(f"I_c = {I_c*1e6:.3f} µA")
print(f"T1 = {T1_total*1e6:.1f} µs")
print(f"T2 = {T2_total*1e6:.1f} µs")
print(f"F_CNOT = {F_CNOT*100:.2f}%")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: System Parameters (CORRECTED)
# ═══════════════════════════════════════════════════════════════════════════════

fig1 = plt.figure(figsize=(16, 10))
gs1 = GridSpec(2, 3, figure=fig1, hspace=0.3, wspace=0.3)

# Panel A: Quantum ratio vs temperature
ax1 = fig1.add_subplot(gs1[0, :2])
temps = np.linspace(0.02, 10, 200)
freqs = [5.8, 50, 100, 200, 300, 500]
colors = ['#C62828', '#EF5350', '#FF9800', '#FF6F00', '#2E7D32', '#1B5E20']

for freq, color in zip(freqs, colors):
    ratios = [freq * 1e9 / (k_B_freq * T) for T in temps]
    lw = 3 if freq in [5.8, 300] else 2
    ls = '--' if freq == 5.8 else '-'
    ax1.semilogy(temps, ratios, linewidth=lw, linestyle=ls,
                 color=color, label=f'{freq:.0f} GHz', alpha=0.8)

ax1.axhline(5, color='green', linestyle='--', linewidth=2, alpha=0.5)
ax1.axhline(3, color='orange', linestyle='--', linewidth=2, alpha=0.5)
ax1.plot(0.02, 5.8e9 / (k_B_freq * 0.02), 'o', markersize=12,
         color='#2E7D32', markeredgecolor='black', markeredgewidth=2)
ax1.plot(4, 300e9 / (k_B_freq * 4), '*', markersize=18,
         color='#FF6F00', markeredgecolor='black', markeredgewidth=2)

ax1.set_xlabel('Temperature (K)', fontweight='bold')
ax1.set_ylabel('Quantum Ratio ℏω/(kT)', fontweight='bold')
ax1.set_title('(a) Quantum Regime: Frequency Requirements', fontweight='bold', pad=10)
ax1.legend(loc='upper right', ncol=2)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 10)
ax1.set_ylim(0.3, 1000)

# Panel B: Thermal population  
ax2 = fig1.add_subplot(gs1[0, 2])
temps_th = np.array([2, 3, 4, 5, 6, 8, 10])
n_th_vals = [thermal_occupation(300e9, T) * 100 for T in temps_th]

bars = ax2.bar(range(len(temps_th)), n_th_vals, color='#FF6F00',
               edgecolor='black', linewidth=2, alpha=0.7)
ax2.axhline(5, color='green', linestyle='--', linewidth=2, alpha=0.7)
ax2.set_ylabel('Thermal Population (%)', fontweight='bold')
ax2.set_xlabel('Temperature (K)', fontweight='bold')
ax2.set_title('(b) Thermal Population\nat 300 GHz', fontweight='bold', pad=10)
ax2.set_xticks(range(len(temps_th)))
ax2.set_xticklabels([f'{T:.0f}' for T in temps_th])
ax2.grid(True, alpha=0.3, axis='y')

# Panel C: Junction size vs current density (CORRECTED)
ax3 = fig1.add_subplot(gs1[1, 0])

J_c_values = np.array([100, 500, 1000, 2000, 5000])
I_c_correct = 1.510e-6  # CORRECTED VALUE

diameters_corrected = []
for J_c in J_c_values:
    A = I_c_correct / (J_c * 1e4)
    d = np.sqrt(4 * A / np.pi) * 1e9
    diameters_corrected.append(d)

colors_jc = []
for d in diameters_corrected:
    if d < 200:
        colors_jc.append('#C62828')
    elif d < 500:
        colors_jc.append('#FF9800')
    elif d < 1000:
        colors_jc.append('#FFC107')
    else:
        colors_jc.append('#EF5350')

bars = ax3.barh(range(len(J_c_values)), diameters_corrected, color=colors_jc,
                 edgecolor='black', linewidth=2, alpha=0.7)

ax3.axvline(200, color='orange', linestyle='--', linewidth=2, alpha=0.7,
            label='200 nm (target)')
ax3.axvline(500, color='green', linestyle='--', linewidth=2, alpha=0.7,
            label='500 nm (easier)')

ax3.set_xlabel('Junction Diameter (nm)', fontweight='bold')
ax3.set_ylabel('Current Density (A/cm²)', fontweight='bold')
ax3.set_title(f'(c) Junction Size (I_c={I_c_correct*1e6:.2f} µA)', fontweight='bold', pad=10)
ax3.set_yticks(range(len(J_c_values)))
ax3.set_yticklabels([f'{J:.0f}' for J in J_c_values])
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3, axis='x')

for i, (d, J_c) in enumerate(zip(diameters_corrected, J_c_values)):
    ax3.text(d + 50, i, f'{d:.0f} nm', va='center', fontsize=9, fontweight='bold')

# Panel D: Dispersive parameters
ax4 = fig1.add_subplot(gs1[1, 1])

params_names = ['|χ₁|', '|χ₂|', 'J_zz', 'Δ₁/g', 'Δ₂/g']
params_values = [5.36, 1.70, 9.38, 40, 80]  # CORRECTED chi values
colors_params = ['#42A5F5', '#29B6F6', '#26C6DA', '#66BB6A', '#9CCC65']

bars = ax4.bar(range(len(params_names)), params_values, color=colors_params,
                edgecolor='black', linewidth=2, alpha=0.7)

ax4.set_ylabel('Value (MHz or ratio)', fontweight='bold')
ax4.set_title('(d) Dispersive Parameters', fontweight='bold', pad=10)
ax4.set_xticks(range(len(params_names)))
ax4.set_xticklabels(params_names, fontweight='bold')
ax4.grid(True, alpha=0.3, axis='y')
ax4.set_yscale('log')

for i, (val, name) in enumerate(zip(params_values, params_names)):
    if 'Δ' in name:
        label = f'{val:.0f}'
    else:
        label = f'{val:.1f} MHz'
    ax4.text(i, val * 1.3, label, ha='center', fontsize=9, fontweight='bold')

# Panel E: Energy levels
ax5 = fig1.add_subplot(gs1[1, 2])

omega_q = 300e9
alpha = -15e9
E0, E1 = 0, 1
E2 = 2 + alpha/(omega_q * 2*np.pi)
E3 = 3 + 3*alpha/(omega_q * 2*np.pi)

levels = [E0, E1, E2, E3]
labels = ['|0⟩', '|1⟩', '|2⟩', '|3⟩']

for i, (E, label) in enumerate(zip(levels, labels)):
    ax5.plot([0, 1], [E, E], 'k-', linewidth=3)
    ax5.text(1.1, E, label, fontsize=14, fontweight='bold', va='center')
    pop = np.exp(-i * omega_q*2*np.pi / (k_B_freq * 4)) * 100
    ax5.text(-0.3, E, f'{pop:.1f}%', fontsize=8, ha='right', va='center')

ax5.annotate('', xy=(0.5, E1), xytext=(0.5, E0),
             arrowprops=dict(arrowstyle='<->', lw=2, color='#2E7D32'))
ax5.text(0.6, (E0+E1)/2, 'ω_q', fontsize=11, fontweight='bold', color='#2E7D32')

ax5.set_xlim(-0.4, 1.5)
ax5.set_ylim(-0.2, 3.5)
ax5.set_ylabel('Energy (ℏω_q units)', fontweight='bold')
ax5.set_title('(e) Transmon Energy Levels', fontweight='bold', pad=10)
ax5.set_xticks([])
ax5.grid(True, alpha=0.3, axis='y')

plt.suptitle('300 GHz / 4K System Parameters', fontsize=16, fontweight='bold', y=0.98)
plt.savefig('./CORRECTED_fig1_system_parameters.png',
            dpi=300, bbox_inches='tight', facecolor='white')
print("✓ Saved: CORRECTED_fig1_system_parameters.png")
plt.close()

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Decoherence (CORRECTED)
# ═══════════════════════════════════════════════════════════════════════════════

fig2 = plt.figure(figsize=(16, 12))
gs2 = GridSpec(3, 3, figure=fig2, hspace=0.3, wspace=0.3)

# Panel A: T1 vs frequency
ax1 = fig2.add_subplot(gs2[0, :])

freqs_scan = np.linspace(5, 500, 100)
T1_diel_5 = 100
T1_diel_scan = T1_diel_5 * (5 / freqs_scan)
T1_qp_scan = np.ones_like(freqs_scan) * 200
T1_rad_unshielded_scan = T1_diel_5 * (5 / freqs_scan)**2
T1_rad_shielded_scan = T1_rad_unshielded_scan * 1000

ax1.plot(freqs_scan, T1_diel_scan, '-', linewidth=3, label='Dielectric (∝ 1/ω)', color='#EF5350')
ax1.plot(freqs_scan, T1_qp_scan, '-', linewidth=3, label='Quasiparticle (const)', color='#4CAF50')
ax1.plot(freqs_scan, T1_rad_unshielded_scan, '--', linewidth=2, label='Radiative unshielded', color='#C62828', alpha=0.7)
ax1.plot(freqs_scan, T1_rad_shielded_scan, '-', linewidth=3, label='Radiative 3D cavity', color='#26C6DA')

ax1.plot(5.8, 100, 'o', markersize=12, color='#2E7D32', markeredgecolor='black', markeredgewidth=2)
ax1.plot(300, T1_total*1e6, '*', markersize=18, color='#FF6F00', markeredgecolor='black', markeredgewidth=2)

ax1.set_xlabel('Frequency (GHz)', fontweight='bold', fontsize=13)
ax1.set_ylabel('T₁ (µs)', fontweight='bold', fontsize=13)
ax1.set_title('(a) T₁ Decoherence Mechanisms vs Frequency', fontweight='bold', pad=10)
ax1.legend(fontsize=10, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')
ax1.set_ylim(0.1, 1000)

ax1.fill_between([0, 500], 15, 1000, alpha=0.1, color='green')

# Panel B: T1 contributions (CORRECTED)
ax2 = fig2.add_subplot(gs2[1, 0])

mechanisms = ['Dielectric', 'Quasiparticle', 'Radiative\n(3D)', 'TLS', 'Purcell']
T1_values = [53.1, 48.6, 373.8, 50.0, 320.0]  # CORRECTED
colors_mech = ['#EF5350', '#4CAF50', '#26C6DA', '#AB47BC', '#FF9800']

contributions = [1/T1 for T1 in T1_values]
total_rate = sum(contributions)
percentages = [c/total_rate * 100 for c in contributions]

bars = ax2.barh(range(len(mechanisms)), percentages, color=colors_mech,
                 edgecolor='black', linewidth=2, alpha=0.7)

ax2.set_xlabel('Contribution to 1/T₁ (%)', fontweight='bold')
ax2.set_title('(b) T₁ Breakdown', fontweight='bold', pad=10)
ax2.set_yticks(range(len(mechanisms)))
ax2.set_yticklabels(mechanisms, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='x')

for i, (bar, T1, pct) in enumerate(zip(bars, T1_values, percentages)):
    ax2.text(pct + 1, i, f'{T1:.0f} µs\n({pct:.1f}%)',
             va='center', fontsize=9, fontweight='bold')

ax2.text(0.5, 0.95, f'Total T₁ = {T1_total*1e6:.1f} µs',
         transform=ax2.transAxes, fontsize=11, fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8), ha='center')

# Panel C: Gate fidelities (CORRECTED)
ax3 = fig2.add_subplot(gs2[1, 1])

gates = ['Single-qubit\n(20 ns)', 'CNOT\n(168 ns)']
fidelities = [F_single*100, F_CNOT*100]  # CORRECTED
colors_gates = ['#42A5F5', '#FF6F00']

bars = ax3.bar(range(len(gates)), fidelities, color=colors_gates,
                edgecolor='black', linewidth=2, alpha=0.7)

ax3.axhline(98, color='red', linestyle='--', linewidth=3, alpha=0.7, label='98% threshold')
ax3.fill_between([-0.5, 1.5], 98, 100, alpha=0.15, color='green')
ax3.fill_between([-0.5, 1.5], 95, 98, alpha=0.15, color='yellow')

ax3.set_ylabel('Fidelity (%)', fontweight='bold')
ax3.set_title('(c) Gate Fidelities', fontweight='bold', pad=10)
ax3.set_xticks(range(len(gates)))
ax3.set_xticklabels(gates, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3, axis='y')
ax3.set_ylim(94, 100)

for i, (F, gate) in enumerate(zip(fidelities, gates)):
    status = '✓' if F >= 98 else '⚠️'
    #ax3.text(i, F + 0.3, f'{F:.2f}%\n{status}', ha='center',
    ax3.text(i, F - 1.0, f'{F:.2f}%\n{status}', ha='center',
             fontsize=10, fontweight='bold')

# Panel D: CNOT error budget (CORRECTED)
ax4 = fig2.add_subplot(gs2[1, 2])

error_sources = ['Decoherence', 'Thermal', 'Leakage', 'Control']
error_values = [2.545, 1.405, 0.000, 0.100]  # CORRECTED
colors_err = ['#EF5350', '#FFA726', '#AB47BC', '#42A5F5']

bars = ax4.barh(range(len(error_sources)), error_values, color=colors_err,
                 edgecolor='black', linewidth=2, alpha=0.7)

ax4.set_xlabel('Error (%)', fontweight='bold')
ax4.set_title('(d) CNOT Error Budget', fontweight='bold', pad=10)
ax4.set_yticks(range(len(error_sources)))
ax4.set_yticklabels(error_sources, fontweight='bold')
ax4.grid(True, alpha=0.3, axis='x')

total_error = sum(error_values)
for i, val in enumerate(error_values):
    ax4.text(val + 0.1, i, f'{val:.2f}%', va='center', fontsize=10, fontweight='bold')

ax4.axvline(2.0, color='red', linestyle='--', linewidth=3, alpha=0.7)
ax4.text(2.0, 3.5, '2% target', fontsize=9, fontweight='bold')

ax4.text(0.98, 0.05, f'Total: {total_error:.2f}%\nF={100-total_error:.2f}%',
         transform=ax4.transAxes, fontsize=11, fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='orange', alpha=0.8),
         ha='right', va='bottom')

# Panel E: Fidelity vs Temperature (CORRECTED)
ax5 = fig2.add_subplot(gs2[2, :])

temps_fid = np.array([2, 3, 4, 5, 6, 8, 10])
T2_corrected = 13.2e-6
t_CNOT = 167.6e-9
alpha_hz = 15e9
J_zz = 9.38e6

fidelities_temp = []
for T_scan in temps_fid:
    n_th_scan = thermal_occupation(300e9, T_scan)
    eps_decoh = t_CNOT / (T2_corrected / 2)
    eps_thermal = n_th_scan * 0.5
    eps_leakage = (J_zz / abs(alpha_hz))**2
    eps_control = 0.001
    F_scan = (1 - (eps_decoh + eps_thermal + eps_leakage + eps_control)) * 100
    fidelities_temp.append(F_scan)

ax5.plot(temps_fid, fidelities_temp, 'o-', linewidth=3, markersize=10,
         color='#FF6F00', markeredgecolor='black', markeredgewidth=2)

ax5.axhline(98, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Threshold')
ax5.fill_between([0, 12], 98, 100, alpha=0.15, color='green')
ax5.fill_between([0, 12], 95, 98, alpha=0.15, color='yellow')

ax5.plot(4, fidelities_temp[2], '*', markersize=20, color='#2E7D32',
         markeredgecolor='black', markeredgewidth=2)

ax5.set_xlabel('Temperature (K)', fontweight='bold')
ax5.set_ylabel('CNOT Fidelity (%)', fontweight='bold')
ax5.set_title('(e) Fidelity vs Temperature', fontweight='bold', pad=10)
ax5.legend(fontsize=10)
ax5.grid(True, alpha=0.3)
ax5.set_xlim(1, 11)
ax5.set_ylim(80, 100)

plt.suptitle('Decoherence and Gate Fidelity Analysis', fontsize=16, fontweight='bold', y=0.995)
plt.savefig('./CORRECTED_fig2_decoherence_fidelity.png',
            dpi=300, bbox_inches='tight', facecolor='white')
print("✓ Saved: CORRECTED_fig2_decoherence_fidelity.png")
plt.close()

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Junction Arrays (NEW)
# ═══════════════════════════════════════════════════════════════════════════════

fig3 = plt.figure(figsize=(14, 10))
gs3 = GridSpec(2, 2, figure=fig3, hspace=0.3, wspace=0.3)

# Panel A: Array configuration comparison
ax1 = fig3.add_subplot(gs3[0, :])

array_configs = ['1×1\n(single)', '3×3', '5×5', '7×7', '10×10']
N_junctions = [1, 9, 25, 49, 100]
I_c_total = 1.51  # µA

d_single_array = []
for N in N_junctions:
    I_single = I_c_total / N
    A_single = I_single*1e-6 / (1000 * 1e4)  # J_c = 1000 A/cm²
    d = np.sqrt(4 * A_single / np.pi) * 1e9
    d_single_array.append(d)

colors_array = []
for d in d_single_array:
    if d < 50:
        colors_array.append('#C62828')
    elif d < 100:
        colors_array.append('#FF9800')
    elif d < 200:
        colors_array.append('#4CAF50')
    else:
        colors_array.append('#FFC107')

bars = ax1.bar(range(len(array_configs)), d_single_array, color=colors_array,
                edgecolor='black', linewidth=2, alpha=0.7)

ax1.axhline(50, color='red', linestyle='--', linewidth=2, alpha=0.7, label='50 nm (limit)')
ax1.axhline(100, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='100 nm (challenging)')
ax1.axhline(200, color='green', linestyle='--', linewidth=2, alpha=0.7, label='200 nm (feasible)')

ax1.set_ylabel('Individual Junction Diameter (nm)', fontweight='bold')
ax1.set_xlabel('Array Configuration', fontweight='bold')
ax1.set_title('(a) Junction Array Solution: Individual Junction Size vs Array Size',
              fontweight='bold', pad=10)
ax1.set_xticks(range(len(array_configs)))
ax1.set_xticklabels(array_configs, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3, axis='y')
ax1.set_ylim(0, 500)

for i, (d, config) in enumerate(zip(d_single_array, array_configs)):
    ax1.text(i, d + 20, f'{d:.0f} nm', ha='center', fontsize=10, fontweight='bold')

# Panel B: Array advantages
ax2 = fig3.add_subplot(gs3[1, 0])

advantages = ['Fab Yield', 'Tunability', 'Redundancy', 'Uniformity']
scores = [9, 8, 9, 7]  # Out of 10

bars = ax2.barh(range(len(advantages)), scores, color='#4CAF50',
                 edgecolor='black', linewidth=2, alpha=0.7)

ax2.set_xlabel('Advantage Score (0-10)', fontweight='bold')
ax2.set_title('(b) Junction Array Advantages', fontweight='bold', pad=10)
ax2.set_yticks(range(len(advantages)))
ax2.set_yticklabels(advantages, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='x')
ax2.set_xlim(0, 10)

for i, score in enumerate(scores):
    ax2.text(score + 0.3, i, f'{score}/10', va='center', fontsize=10, fontweight='bold')

# Panel C: E_J/E_C trade-off
ax3 = fig3.add_subplot(gs3[1, 1])

J_c_scan = np.array([500, 1000, 2000, 3000, 5000])
d_feasible = 200  # nm target

E_J_E_C_ratios = []
for J_c in J_c_scan:
    I_c_max = (J_c * 1e4) * (np.pi * (d_feasible*1e-9/2)**2)
    E_J_max = (Phi_0 / (2*np.pi)) * I_c_max
    E_C_needed = 9.94e-24  # Fixed by 300 GHz requirement
    ratio = E_J_max / E_C_needed
    E_J_E_C_ratios.append(ratio)

colors_ratio = []
for ratio in E_J_E_C_ratios:
    if ratio < 10:
        colors_ratio.append('#C62828')
    elif ratio < 30:
        colors_ratio.append('#FF9800')
    else:
        colors_ratio.append('#4CAF50')

bars = ax3.bar(range(len(J_c_scan)), E_J_E_C_ratios, color=colors_ratio,
                edgecolor='black', linewidth=2, alpha=0.7)

ax3.axhline(10, color='red', linestyle='--', linewidth=2, alpha=0.7, label='E_J/E_C = 10 (min)')
ax3.axhline(30, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='E_J/E_C = 30 (good)')
ax3.axhline(50, color='green', linestyle='--', linewidth=2, alpha=0.7, label='E_J/E_C = 50 (target)')

ax3.set_ylabel('Achievable E_J/E_C', fontweight='bold')
ax3.set_xlabel('Current Density (A/cm²)', fontweight='bold')
ax3.set_title(f'(c) E_J/E_C vs J_c (d={d_feasible} nm)', fontweight='bold', pad=10)
ax3.set_xticks(range(len(J_c_scan)))
ax3.set_xticklabels([f'{J:.0f}' for J in J_c_scan])
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3, axis='y')
ax3.set_ylim(0, 60)

for i, ratio in enumerate(E_J_E_C_ratios):
    ax3.text(i, ratio + 2, f'{ratio:.1f}', ha='center', fontsize=10, fontweight='bold')

plt.suptitle('Junction Array Solution for 300 GHz Operation', fontsize=16, fontweight='bold', y=0.98)
plt.savefig('./CORRECTED_fig3_junction_arrays.png',
            dpi=300, bbox_inches='tight', facecolor='white')
print("✓ Saved: CORRECTED_fig3_junction_arrays.png")
plt.close()

print("\n✓ All corrected figures generated!")
print(f"\nKEY CORRECTED VALUES IN FIGURES:")
print(f"  I_c = {I_c*1e6:.3f} µA")
print(f"  Junction diameter: 196-1387 nm (depending on J_c)")
print(f"  T₁ = {T1_total*1e6:.1f} µs")
print(f"  T₂ = {T2_total*1e6:.1f} µs")
print(f"  F_CNOT = {F_CNOT*100:.2f}%")
print(f"  Gap to threshold: {F_CNOT*100 - 98:.2f}%")
