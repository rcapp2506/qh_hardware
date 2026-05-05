"""
Two-Qubit Cavity QED System - Visualization and Summary
Complete analysis with plots and LaTeX documentation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import patches
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("GENERATING COMPREHENSIVE VISUALIZATIONS")
print("="*80)

# System parameters
omega_r = 2 * np.pi * 6.5  # GHz
omega_q1 = 2 * np.pi * 5.8
omega_q2 = 2 * np.pi * 5.2
#g1 = 2 * np.pi * 0.05
#g2 = 2 * np.pi * 0.05
g1 = 2 * np.pi * 0.08  # GHz
g2 = 2 * np.pi * 0.08  # GHz

#################################### aggiunto #################################
# alpha in pulsazione angolare per coerenza con g, Delta (tutti in 2π·GHz)
alpha_1 = 2 * np.pi * (-0.200)   # Qubit 1 anharmonicity (-200 MHz)
alpha_2 = 2 * np.pi * (-0.200)   # Qubit 2 anharmonicity (-200 MHz)
######################################

Delta1 = omega_q1 - omega_r
Delta2 = omega_q2 - omega_r
#chi1 = g1**2 / Delta1
#chi2 = g2**2 / Delta2
# Transverse exchange J (Eq. eq:J_dispersive in the thesis, post-G7 review).
# The pre-G7 version used g1*g2*(1/D1+1/D2) (factor 2 missing); the correct
# second-order Schrieffer-Wolff result of Majer-Gambetta has /2:
J_xchange = g1 * g2 * (1/Delta1 + 1/Delta2) / 2
# Canonical thesis values after full anharmonic correction
# (Sec. cavity_mediated_couplings, Eq. zeta_zz_correct):
J_thesis_MHz = 6.5
zeta_thesis_MHz = 1.7
# Use thesis-canonical for downstream displays
zeta12 = -J_thesis_MHz * 1e-3 * 2 * np.pi  # GHz, kept name for back-compat
kappa = 2 * np.pi * 0.001  # Cavity decay rate

################################
def dispersive_shift(g, Delta, alpha):
    """Calculate dispersive shift χ = (g²/Δ) × α/(Δ+α)"""
    return (g**2 / Delta) * (alpha / (Delta + alpha))

chi1 = dispersive_shift(g1, Delta1, alpha_1)
chi2 = dispersive_shift(g2, Delta2, alpha_2)
################################


# Calculate CNOT gate time via cross-Kerr
t_cnot_crosskerr = np.pi / (4 * abs(zeta12))

# Create comprehensive figure
fig = plt.figure(figsize=(18, 12))

# ============================================================================
# Plot 1: Energy Level Diagram
# ============================================================================
ax1 = plt.subplot(2, 3, 1)

# Draw energy levels
cavity_levels = np.arange(5)
qubit1_levels = np.array([0, 1])
qubit2_levels = np.array([0, 1])

# Cavity
for n in cavity_levels:
    energy = omega_r/(2*np.pi) * (n + 0.5)
    ax1.plot([0, 0.8], [energy, energy], 'b-', linewidth=2)
    ax1.text(0.9, energy, f'|{n}⟩', fontsize=10, va='center')

# Qubit 1
for n in qubit1_levels:
    energy = omega_q1/(2*np.pi) * (n + 0.5)
    ax1.plot([2, 2.8], [energy, energy], 'r-', linewidth=2)
    ax1.text(2.9, energy, f'|{n}⟩', fontsize=10, va='center')

# Qubit 2
for n in qubit2_levels:
    energy = omega_q2/(2*np.pi) * (n + 0.5)
    ax1.plot([4, 4.8], [energy, energy], 'g-', linewidth=2)
    ax1.text(4.9, energy, f'|{n}⟩', fontsize=10, va='center')

# Draw coupling arrows
ax1.annotate('', xy=(2.1, omega_q1/(2*np.pi)*0.5), xytext=(0.7, omega_r/(2*np.pi)*0.5),
            arrowprops=dict(arrowstyle='<->', color='purple', lw=2))
ax1.text(1.4, (omega_q1/(2*np.pi) + omega_r/(2*np.pi))/4, 'g₁', fontsize=11, 
         ha='center', color='purple', fontweight='bold')

ax1.annotate('', xy=(4.1, omega_q2/(2*np.pi)*0.5), xytext=(0.7, omega_r/(2*np.pi)*1.5),
            arrowprops=dict(arrowstyle='<->', color='orange', lw=2))
ax1.text(2.4, (omega_q2/(2*np.pi) + omega_r/(2*np.pi)*1.5)/2, 'g₂', fontsize=11,
         ha='center', color='orange', fontweight='bold')

ax1.set_ylim(0, 7)
ax1.set_xlim(-0.5, 5.5)
ax1.set_ylabel('Energy (GHz)', fontsize=12, fontweight='bold')
ax1.set_title('Energy Level Diagram\nJaynes-Cummings Coupling', fontsize=13, fontweight='bold')
ax1.set_xticks([0.4, 2.4, 4.4])
ax1.set_xticklabels(['Cavity', 'Qubit 1', 'Qubit 2'], fontsize=11)
ax1.grid(True, alpha=0.3)

# ============================================================================
# Plot 2: Dispersive Shift Visualization
# ============================================================================
ax2 = plt.subplot(2, 3, 2)

# Cavity frequency for different qubit states
omega_00 = omega_r - chi1 - chi2
omega_01 = omega_r - chi1 + chi2
omega_10 = omega_r + chi1 - chi2
omega_11 = omega_r + chi1 + chi2

states = ['|00⟩', '|01⟩', '|10⟩', '|11⟩']
frequencies = [omega_00, omega_01, omega_10, omega_11]
colors = ['blue', 'green', 'red', 'purple']

y_pos = np.arange(len(states))
for i, (state, freq, color) in enumerate(zip(states, frequencies, colors)):
    ax2.barh(i, (freq - omega_r)/(2*np.pi)*1000, left=omega_r/(2*np.pi),
             color=color, alpha=0.6, edgecolor='black', linewidth=2)
    ax2.text(omega_r/(2*np.pi) + (freq - omega_r)/(2*np.pi)*1000/2, i,
             f'{freq/(2*np.pi):.3f} GHz', ha='center', va='center',
             fontsize=9, fontweight='bold')

ax2.axvline(omega_r/(2*np.pi), color='black', linestyle='--', linewidth=2,
            label='ω_r (bare)')
ax2.set_yticks(y_pos)
ax2.set_yticklabels(states, fontsize=11)
ax2.set_xlabel('Cavity Frequency (GHz)', fontsize=12, fontweight='bold')
ax2.set_title('Dispersive Frequency Shifts\nχᵢ-dependent', fontsize=13, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3, axis='x')

# ============================================================================
# Plot 3: CNOT Gate Circuit and Truth Table
# ============================================================================
ax3 = plt.subplot(2, 3, 3)
ax3.axis('off')

# Draw circuit
y_q1 = 0.7
y_q2 = 0.3

# Qubit lines
ax3.plot([0.1, 0.9], [y_q1, y_q1], 'k-', linewidth=2)
ax3.plot([0.1, 0.9], [y_q2, y_q2], 'k-', linewidth=2)

# Labels
ax3.text(0.05, y_q1, 'q₁ (control)', ha='right', va='center', fontsize=11, fontweight='bold')
ax3.text(0.05, y_q2, 'q₂ (target)', ha='right', va='center', fontsize=11, fontweight='bold')

# Control dot
ax3.plot(0.5, y_q1, 'ko', markersize=15)

# Target (⊕)
circle = plt.Circle((0.5, y_q2), 0.04, fill=False, edgecolor='black', linewidth=2)
ax3.add_patch(circle)
ax3.plot([0.5, 0.5], [y_q2-0.04, y_q2+0.04], 'k-', linewidth=2)
ax3.plot([0.5-0.04, 0.5+0.04], [y_q2, y_q2], 'k-', linewidth=2)

# Vertical connection
ax3.plot([0.5, 0.5], [y_q2+0.04, y_q1], 'k-', linewidth=2)

# Truth table
table_data = [
    ['Input', 'Output'],
    ['|00⟩', '|00⟩'],
    ['|01⟩', '|01⟩'],
    ['|10⟩', '|11⟩ ✓'],
    ['|11⟩', '|10⟩ ✓'],
]

table_x = 0.15
table_y = 0.05
cell_height = 0.03
cell_width = 0.15

for i, row in enumerate(table_data):
    y = table_y - i * cell_height
    for j, cell in enumerate(row):
        x = table_x + j * cell_width
        if i == 0:
            ax3.text(x, y, cell, ha='left', va='center', fontsize=10, fontweight='bold')
        else:
            ax3.text(x, y, cell, ha='left', va='center', fontsize=9,
                    color='red' if '✓' in cell else 'black')

ax3.set_xlim(0, 1)
ax3.set_ylim(0, 1)
ax3.set_title('CNOT Gate\nCircuit and Truth Table', fontsize=13, fontweight='bold',
              y=0.95)

# ============================================================================
# Plot 4: Gate Time Comparison
# ============================================================================
ax4 = plt.subplot(2, 3, 4)

gate_types = ['X/Y/Z\n(single)', 'Hadamard\n(single)', 'CNOT\n(CR)', 'CNOT\n(Dispersive)']
gate_times = [20, 20, 150, t_cnot_crosskerr]  # ns
colors_gates = ['lightblue', 'lightgreen', 'orange', 'red']

bars = ax4.bar(gate_types, gate_times, color=colors_gates, edgecolor='black', linewidth=2)

for bar, time in zip(bars, gate_times):
    height = bar.get_height()
    if time > 1000:
        label = f'{time/1000:.1f} μs'
    else:
        label = f'{time:.0f} ns'
    ax4.text(bar.get_x() + bar.get_width()/2., height,
            label, ha='center', va='bottom', fontsize=10, fontweight='bold')

ax4.set_ylabel('Gate Time', fontsize=12, fontweight='bold')
ax4.set_yscale('log')
ax4.set_title('Quantum Gate Times\n(Logarithmic Scale)', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3, axis='y')
ax4.tick_params(axis='x', labelsize=9)

# Add annotations
ax4.axhline(100, color='green', linestyle='--', alpha=0.5, linewidth=2)
ax4.text(3.5, 120, 'Fast gates\n(< 100 ns)', fontsize=9, color='green',
         ha='right', fontweight='bold')

ax4.axhline(1000, color='orange', linestyle='--', alpha=0.5, linewidth=2)
ax4.text(3.5, 1200, 'Moderate\n(~1 μs)', fontsize=9, color='orange',
         ha='right', fontweight='bold')

# ============================================================================
# Plot 5: Readout Frequency Spectrum
# ============================================================================
ax5 = plt.subplot(2, 3, 5)

# Frequency axis
freq_range = np.linspace(omega_r/(2*np.pi) - 0.01, omega_r/(2*np.pi) + 0.01, 1000)

# Lorentzian peaks for each state
linewidth = kappa / (2*np.pi)

def lorentzian(x, x0, gamma):
    return gamma / ((x - x0)**2 + gamma**2)

# Plot each state's peak
omega_states = [omega_00, omega_01, omega_10, omega_11]
labels_states = ['|00⟩', '|01⟩', '|10⟩', '|11⟩']
colors_states = ['blue', 'green', 'red', 'purple']

for omega, label, color in zip(omega_states, labels_states, colors_states):
    spectrum = lorentzian(freq_range, omega/(2*np.pi), linewidth)
    ax5.plot(freq_range, spectrum/spectrum.max(), linewidth=2, label=label, color=color)
    ax5.axvline(omega/(2*np.pi), color=color, linestyle='--', alpha=0.5)

ax5.axvline(omega_r/(2*np.pi), color='black', linestyle=':', linewidth=2,
            label='ω_r (bare)', alpha=0.7)

ax5.set_xlabel('Frequency (GHz)', fontsize=12, fontweight='bold')
ax5.set_ylabel('Readout Signal (a.u.)', fontsize=12, fontweight='bold')
ax5.set_title('Dispersive Readout Spectrum\nState-Dependent Cavity Response',
              fontsize=13, fontweight='bold')
ax5.legend(fontsize=9, loc='upper right')
ax5.grid(True, alpha=0.3)

# ============================================================================
# Plot 6: Bell State Preparation Sequence
# ============================================================================
ax6 = plt.subplot(2, 3, 6)

# Time evolution of populations during Bell state prep
# Conceptual plot showing |00⟩ → (|00⟩ + |11⟩)/√2

time_steps = np.linspace(0, 1, 100)

# After H on q1: equal superposition
p00_t = 0.5 * np.ones_like(time_steps)
p10_t = 0.5 * np.ones_like(time_steps)
p01_t = np.zeros_like(time_steps)
p11_t = np.zeros_like(time_steps)

# After CNOT (at t > 0.5): |00⟩ → |00⟩, |10⟩ → |11⟩
transition_idx = int(len(time_steps) * 0.5)
p00_t[transition_idx:] = 0.5
p10_t[transition_idx:] = 0
p01_t[transition_idx:] = 0
p11_t[transition_idx:] = 0.5

ax6.fill_between(time_steps, 0, p00_t, label='|00⟩', alpha=0.7, color='blue')
ax6.fill_between(time_steps, p00_t, p00_t + p01_t, label='|01⟩', alpha=0.7, color='green')
ax6.fill_between(time_steps, p00_t + p01_t, p00_t + p01_t + p10_t, label='|10⟩', alpha=0.7, color='red')
ax6.fill_between(time_steps, p00_t + p01_t + p10_t, p00_t + p01_t + p10_t + p11_t,
                label='|11⟩', alpha=0.7, color='purple')

# Mark gates
ax6.axvline(0.25, color='black', linestyle='--', linewidth=2, alpha=0.5)
ax6.text(0.25, 1.05, 'H₁', ha='center', fontsize=12, fontweight='bold')

ax6.axvline(0.5, color='black', linestyle='--', linewidth=2, alpha=0.5)
ax6.text(0.5, 1.05, 'CNOT', ha='center', fontsize=12, fontweight='bold')

ax6.set_xlabel('Normalized Time', fontsize=12, fontweight='bold')
ax6.set_ylabel('Population', fontsize=12, fontweight='bold')
ax6.set_title('Bell State Preparation\n|ψ⟩ = (|00⟩ + |11⟩)/√2', fontsize=13, fontweight='bold')
ax6.legend(fontsize=10, loc='right')
ax6.set_ylim(0, 1.1)
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('./two_qubit_cavity_qed_analysis.png',
            dpi=300, bbox_inches='tight')
print("\n✓ Comprehensive analysis figure saved!")

# ============================================================================
# Create Circuit Diagram Figure
# ============================================================================
fig2, axes = plt.subplots(1, 2, figsize=(10, 5))

# Left: Detailed circuit diagram
ax = axes[0]
ax.set_xlim(0, 10)
ax.set_ylim(0, 4)
ax.axis('off')

# Title
ax.text(5, 3.8, 'Two-Qubit Cavity QED System', ha='center', fontsize=14, fontweight='bold')

# Draw cavity
cavity_box = FancyBboxPatch((3.5, 1.5), 3, 1.5, boxstyle="round,pad=0.1",
                           facecolor='lightblue', edgecolor='blue', linewidth=2)
ax.add_patch(cavity_box)
ax.text(5, 2.25, 'Cavity\nω_r', ha='center', va='center', fontsize=12, fontweight='bold')

# Draw qubits
qubit1_box = FancyBboxPatch((0.5, 2.5), 1.5, 0.7, boxstyle="round,pad=0.05",
                           facecolor='lightcoral', edgecolor='red', linewidth=2)
ax.add_patch(qubit1_box)
ax.text(1.25, 2.85, 'Qubit 1\nω_q1', ha='center', va='center', fontsize=10, fontweight='bold')

qubit2_box = FancyBboxPatch((0.5, 0.8), 1.5, 0.7, boxstyle="round,pad=0.05",
                           facecolor='lightgreen', edgecolor='green', linewidth=2)
ax.add_patch(qubit2_box)
ax.text(1.25, 1.15, 'Qubit 2\nω_q2', ha='center', va='center', fontsize=10, fontweight='bold')

# Draw couplings
arrow1 = FancyArrowPatch((2, 2.85), (3.5, 2.5), arrowstyle='<->', mutation_scale=20,
                        color='red', linewidth=2)
ax.add_patch(arrow1)
ax.text(2.75, 2.9, 'g₁', fontsize=11, fontweight='bold', color='red')

arrow2 = FancyArrowPatch((2, 1.15), (3.5, 1.7), arrowstyle='<->', mutation_scale=20,
                        color='green', linewidth=2)
ax.add_patch(arrow2)
ax.text(2.75, 1.2, 'g₂', fontsize=11, fontweight='bold', color='green')

# Virtual coupling
arrow3 = FancyArrowPatch((0.75, 2.5), (0.75, 1.5), arrowstyle='<->', mutation_scale=15,
                        color='purple', linewidth=2, linestyle='--')
ax.add_patch(arrow3)
ax.text(0.2, 2, 'J\n(transv.\nexch.)', fontsize=9, fontweight='bold', color='purple', ha='center')

# Add readout
ax.text(7.5, 2.25, '→ Readout', fontsize=11, fontweight='bold', color='blue')

# Right: Parameter summary
ax = axes[1]
ax.axis('off')

summary_text = f"""
SYSTEM PARAMETERS

Cavity:
  ω_r/2π = {omega_r/(2*np.pi):.2f} GHz
  κ/2π = {kappa/(2*np.pi)*1000:.1f} kHz
  Q = {omega_r/kappa:.0f}

Qubit 1:
  ω_q1/2π = {omega_q1/(2*np.pi):.2f} GHz
  g₁/2π = {g1/(2*np.pi)*1000:.1f} MHz
  Δ₁/2π = {Delta1/(2*np.pi):.2f} GHz

Qubit 2:
  ω_q2/2π = {omega_q2/(2*np.pi):.2f} GHz
  g₂/2π = {g2/(2*np.pi)*1000:.1f} MHz
  Δ₂/2π = {Delta2/(2*np.pi):.2f} GHz

Dispersive Parameters:
  χ₁/2π = {chi1/(2*np.pi)*1000:.2f} MHz
  χ₂/2π = {chi2/(2*np.pi)*1000:.2f} MHz
  J/2π  = {J_thesis_MHz:.2f} MHz   (transv. exch.)
  |ζ_zz|/2π = {zeta_thesis_MHz:.2f} MHz   (long. cross-Kerr)

Regime Check:
  |Δ₁|/g₁ = {abs(Delta1)/g1:.1f} ✓
  |Δ₂|/g₂ = {abs(Delta2)/g2:.1f} ✓
  χ₁/κ = {abs(chi1)/kappa:.1f} ✓
"""

ax.text(0.1, 0.95, summary_text, fontsize=10, verticalalignment='top',
        family='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('./two_qubit_system_diagram.png',
            dpi=300, bbox_inches='tight')
print("✓ System diagram figure saved!")

print("\n" + "="*80)
print("ALL VISUALIZATIONS COMPLETE")
print("="*80)
print("\nGenerated figures:")
print("  1. two_qubit_cavity_qed_analysis.png (6-panel comprehensive)")
print("  2. two_qubit_system_diagram.png (system schematic)")
