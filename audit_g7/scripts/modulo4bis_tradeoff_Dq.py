"""
Modulo 4-bis — Trade-off di Δ_q (qubit-qubit detuning) per HEATS-Q baseline
============================================================================

Variabile di design: Δ_q = ω_q1 - ω_q2  (range: 100 MHz → 1500 MHz)
Parametri fissi: ω_r=6.5 GHz, g=80 MHz, α=-200 MHz, T_1=50 μs, T_2=30 μs

Quantità tracciate al variare di Δ_q:
  1. ζ_zz statico (cross-Kerr "always-on") — più piccolo è meglio (non genera errori parassiti)
  2. Ω_ZX rate del CR gate                  — più grande è meglio (gate più veloce)
  3. t_CR = π/(4|Ω_ZX|)                     — più piccolo è meglio
  4. Crosstalk single-qubit drive            — più piccolo è meglio
  5. Spectral selectivity (Δ_q vs σ_drive)   — più alto è meglio
  6. F_CR atteso da T_1/T_2                  — più alto è meglio

Generiamo grafico riassuntivo + tabella di valori chiave.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parametri fissi
omega_r = 6.5    # GHz
g = 0.080        # 80 MHz
alpha = -0.200   # -200 MHz
T1 = 50          # μs
T2 = 30          # μs
Omega_drive = 0.080  # 80 MHz (limite "safe" = α/2.5, regola Sheldon 2016)

# Centro le frequenze attorno a ω̄ = 5.5 GHz (intermediate)
omega_bar = 5.5

print("="*78)
print("TRADE-OFF Δ_q: scan da 100 MHz a 1500 MHz")
print("="*78)
print(f"Parametri fissi: g = {g*1e3:.0f} MHz, α = {alpha*1e3:.0f} MHz, ω_r = {omega_r} GHz")
print(f"Centro frequenze qubit: ω̄ = {omega_bar} GHz, Ω_drive = {Omega_drive*1e3:.0f} MHz")
print()

Delta_q_arr = np.linspace(0.100, 1.500, 30)  # GHz

results = []
for Dq in Delta_q_arr:
    om_q1 = omega_bar + Dq/2
    om_q2 = omega_bar - Dq/2
    D1 = om_q1 - omega_r
    D2 = om_q2 - omega_r
    
    # 1. J transverse exchange
    J = (g**2 / 2) * (1/D1 + 1/D2)
    
    # 2. ζ_zz cross-Kerr statico (Blais 2021)
    zeta_zz = 2 * J**2 * (1/(Dq - alpha) - 1/(Dq + alpha))
    
    # 3. χ_i dispersive shifts (transmon-corrected)
    chi_1 = (g**2 / D1) * (alpha / (D1 + alpha))
    chi_2 = (g**2 / D2) * (alpha / (D2 + alpha))
    
    # 4. Ω_ZX rate del CR (Sheldon 2016)
    Omega_ZX = -J * Omega_drive * alpha / (Dq * (Dq + alpha))
    
    # 5. t_CR = π/(4|Ω_ZX|) in μs (zeta in GHz → 1/GHz = ns; *1e3 → μs)
    if abs(Omega_ZX) > 1e-9:
        t_CR_us = np.pi / (4 * abs(Omega_ZX) * 1e3)  # GHz → ns conversion
    else:
        t_CR_us = float('inf')
    
    # 6. Single-qubit cross-talk: rapporto Ω/Δ_q per drive selectivity
    # Quando Δ_q < ~5*σ_drive_pulse, il drive on Q1 eccita anche Q2.
    # Per pulse di t_pi=20 ns, σ ≈ 1/t_pi = 50 MHz. Selectivity = Δ_q/σ
    sigma_pulse = 50  # MHz, larghezza spettrale per π-pulse 20 ns
    selectivity = Dq * 1e3 / sigma_pulse  # adimensionale
    
    # 7. Crosstalk: ampiezza del drive su Q1 sentita da Q2 (Lorentziana)
    # Reduction factor = (σ/Δ_q)^2
    crosstalk_pct = 100 * (sigma_pulse / (Dq * 1e3))**2
    
    # 8. Fidelity del CR gate (limitata da T_1/T_2)
    Tphi_inv = max(0, 1/T2 - 1/(2*T1))
    Tphi = 1/Tphi_inv if Tphi_inv > 0 else float('inf')
    F_CR = np.exp(-t_CR_us/T1) * np.exp(-t_CR_us/Tphi)**0.5 if t_CR_us != float('inf') else 0
    
    # 9. ZZ-induced gate error (always-on ζ_zz applied for t_CR generates spurious phase)
    zz_error = (zeta_zz * t_CR_us)**2 if t_CR_us != float('inf') else float('inf')  # rad²
    
    results.append({
        'Dq_MHz': Dq*1e3,
        'J_MHz': J*1e3,
        'zeta_zz_kHz': zeta_zz*1e6,
        'chi_1_MHz': chi_1*1e3,
        'chi_2_MHz': chi_2*1e3,
        'Omega_ZX_MHz': Omega_ZX*1e3,
        't_CR_ns': t_CR_us*1e3,
        'selectivity': selectivity,
        'crosstalk_pct': crosstalk_pct,
        'F_CR': F_CR,
        'zz_err_rad2': zz_error,
    })

# Stampa tabella riassuntiva (solo punti chiave)
print(f"{'Δ_q (MHz)':>10s} {'J (MHz)':>9s} {'ζ_zz (kHz)':>11s} {'Ω_ZX (MHz)':>11s} "
      f"{'t_CR (ns)':>10s} {'selectivity':>12s} {'crosstalk':>10s} {'F_CR':>7s}")
print("-"*88)
for r in results[::3]:  # ogni 3 valori
    print(f"{r['Dq_MHz']:>10.0f} {abs(r['J_MHz']):>9.2f} {abs(r['zeta_zz_kHz']):>11.2f} "
          f"{abs(r['Omega_ZX_MHz']):>11.3f} {r['t_CR_ns']:>10.0f} "
          f"{r['selectivity']:>12.1f} {r['crosstalk_pct']:>9.2f}% {r['F_CR']:>7.4f}")

# ─── Trova il punto ottimo (massimo F_CR con vincoli) ───
print(f"\n{'─'*78}")
print("OTTIMIZZAZIONE: trovo Δ_q ottimo con vincoli pratici")
print(f"{'─'*78}")
print("Vincoli:")
print("  - selectivity ≥ 4 (cioè Δ_q ≥ 4×σ_pulse = 200 MHz)")
print("  - crosstalk ≤ 5% (cioè Δ_q ≥ 224 MHz)")
print("  - leakage ζ_zz·T₁ ≤ 0.01 rad (limite single-shot drift)")
print("  - F_CR ≥ 0.99 (soglia surface code)")

best = None
for r in results:
    if (r['selectivity'] >= 4 and r['crosstalk_pct'] <= 5
        and abs(r['zeta_zz_kHz'])*1e-3*T1 <= 0.01*1e3   # in kHz·μs = mrad
        and r['F_CR'] >= 0.99):
        if best is None or r['t_CR_ns'] < best['t_CR_ns']:
            best = r

if best:
    print(f"\n  → Δ_q ottimale = {best['Dq_MHz']:.0f} MHz")
    print(f"     • J = {best['J_MHz']:.2f} MHz  (vs originale -7.03 MHz @ Δ_q=600)")
    print(f"     • ζ_zz statico = {best['zeta_zz_kHz']:.1f} kHz  (vs -124 kHz originale)")
    print(f"     • χ_1 = {best['chi_1_MHz']:.2f} MHz, χ_2 = {best['chi_2_MHz']:.2f} MHz")
    print(f"     • t_CR = {best['t_CR_ns']:.0f} ns  ← gate time praticabile!")
    print(f"     • F_CR ≈ {best['F_CR']*100:.2f}%")
else:
    print("\n  ✗ Nessun valore di Δ_q soddisfa tutti i vincoli con questi parametri.")
    print("    Serve aumentare g, ridurre |α|, o migliorare T_1.")

# ─── PLOT: 4 panel con trade-off ───
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
Dqs = np.array([r['Dq_MHz'] for r in results])

# Panel A: rate del gate vs Δ_q
ax = axes[0, 0]
ax.semilogy(Dqs, [r['t_CR_ns'] for r in results], 'b-', linewidth=2, label='t_CR cross-resonance')
ax.axhline(200, color='green', linestyle='--', alpha=0.6, label='target < 200 ns')
ax.axhline(500, color='orange', linestyle=':', alpha=0.6, label='marginale 500 ns')
ax.set_xlabel('Δ_q (MHz)'); ax.set_ylabel('t_CR (ns)')
ax.set_title('A. Gate time vs qubit-qubit detuning')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9); ax.set_ylim(50, 1e5)

# Panel B: cross-Kerr statico parassita
ax = axes[0, 1]
ax.semilogy(Dqs, [abs(r['zeta_zz_kHz']) for r in results], 'r-', linewidth=2)
ax.set_xlabel('Δ_q (MHz)'); ax.set_ylabel('|ζ_zz| (kHz)')
ax.set_title('B. Cross-Kerr statico parassita\n(piccolo è meglio)')
ax.grid(True, alpha=0.3)

# Panel C: spectral selectivity & crosstalk
ax = axes[1, 0]
ax2 = ax.twinx()
ax.plot(Dqs, [r['selectivity'] for r in results], 'g-', linewidth=2, label='selectivity')
ax2.plot(Dqs, [r['crosstalk_pct'] for r in results], 'orange', linewidth=2, label='crosstalk %')
ax.axhline(4, color='gray', linestyle='--', alpha=0.5)
ax2.axhline(5, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('Δ_q (MHz)'); ax.set_ylabel('Selectivity', color='g')
ax2.set_ylabel('Crosstalk (%)', color='orange')
ax.set_title('C. Single-qubit drive selectivity & crosstalk')
ax.grid(True, alpha=0.3); ax.set_ylim(0, 30); ax2.set_ylim(0, 30)

# Panel D: fidelity totale
ax = axes[1, 1]
ax.plot(Dqs, [r['F_CR']*100 for r in results], 'b-', linewidth=2, label='F_CR (T₁/T₂ limit)')
ax.axhline(99, color='gray', linestyle='--', alpha=0.6, label='soglia 99%')
ax.axhline(99.9, color='gray', linestyle=':', alpha=0.6, label='soglia 99.9%')
if best:
    ax.scatter([best['Dq_MHz']], [best['F_CR']*100], s=120, c='red', zorder=5,
               label=f'Δ_q* = {best["Dq_MHz"]:.0f} MHz')
ax.set_xlabel('Δ_q (MHz)'); ax.set_ylabel('F_CR (%)')
ax.set_title('D. Fidelity attesa')
ax.set_ylim(95, 100); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

plt.suptitle('Trade-off di Δ_q per HEATS-Q baseline (g=80 MHz, α=-200 MHz)',
             fontsize=12, y=1.00)
plt.tight_layout()
plt.savefig('/home/claude/audit_g7/fig_modulo4bis_tradeoff_Dq.pdf', bbox_inches='tight')
print(f"\nFigura salvata: fig_modulo4bis_tradeoff_Dq.pdf")
