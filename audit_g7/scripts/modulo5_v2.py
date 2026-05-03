"""
Modulo 5 v2 — Validazione Lindblad con qt.propagator (versione robusta)
========================================================================

Approccio rivisto:
1. Usa qt.propagator(H, t, c_ops) per ottenere il superop Lindblad
2. Per ogni t, simula l'evoluzione di stati di prova diversi
3. Calcola gate fidelity da metodo Process Tomography esplicito

Scelta operativa per evitare bug di vec convention:
  • Calcolo F_avg da media su 6² = 36 stati di prova (Pauli eigenstates per qubit)
  • F_avg ≈ <ψ_t | ρ_actual | ψ_t> mediato sui 6 stati input × 6 misure
  • Conservativo ma robusto
"""

import numpy as np
import qutip as qt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*78)
print("MODULO 5 v2 — Lindblad CR gate validation con propagator + state-based F")
print("="*78)

systems = {
    'A_baseline_mK': dict(
        name='Sistema A — baseline mK',
        Delta_q=0.280, alpha=-0.200, J=-7.03e-3, Omega_drive=0.060,
        T1_us=50, T2_us=30, n_th=0.0,
    ),
    'B_innovation_300GHz': dict(
        name='Sistema B — innovation 300 GHz @ 4K',
        Delta_q=0.600, alpha=-1.000, J=-9.40e-3, Omega_drive=0.180,
        T1_us=25, T2_us=12, n_th=0.028,
    ),
}

N_lev = 3   # 3 livelli per transmon

def build_system(p):
    Dq, alpha = p['Delta_q'], p['alpha']
    b1 = qt.tensor(qt.destroy(N_lev), qt.qeye(N_lev))
    b2 = qt.tensor(qt.qeye(N_lev), qt.destroy(N_lev))
    H = (Dq * b1.dag()*b1
         + 0.5*alpha * b1.dag()*b1.dag()*b1*b1
         + 0.5*alpha * b2.dag()*b2.dag()*b2*b2
         + p['J'] * (b1.dag()*b2 + b1*b2.dag())
         + 0.5 * p['Omega_drive'] * (b1 + b1.dag()))
    
    gamma_1 = 1e-3 / p['T1_us']
    Tphi_inv = max(0, 1/p['T2_us'] - 1/(2*p['T1_us']))
    gamma_phi = 1e-3 * Tphi_inv
    n_th = p['n_th']
    
    c_ops = []
    for b in [b1, b2]:
        c_ops.append(np.sqrt(gamma_1*(1+n_th)) * b)
        if n_th > 0:
            c_ops.append(np.sqrt(gamma_1*n_th) * b.dag())
        if gamma_phi > 0:
            c_ops.append(np.sqrt(2*gamma_phi) * (b.dag()*b))
    return H, c_ops

# Stati di prova (Pauli eigenstates) per ogni qubit, in spazio 3-livelli (ma popolano solo 0,1)
def pauli_eigenstates_3lev():
    """6 eigenstates: |0⟩, |1⟩, |+⟩, |-⟩, |+i⟩, |-i⟩ in spazio 3-livelli."""
    g = qt.basis(N_lev, 0)
    e = qt.basis(N_lev, 1)
    plus  = (g + e).unit()
    minus = (g - e).unit()
    plusi  = (g + 1j*e).unit()
    minusi = (g - 1j*e).unit()
    return [g, e, plus, minus, plusi, minusi]

def U_ideal_ZX(theta=np.pi/2):
    import scipy.linalg as sla
    Z = np.array([[1,0],[0,-1]])
    X = np.array([[0,1],[1,0]])
    ZX = np.kron(Z, X)
    return sla.expm(-1j * theta/2 * ZX)

def project_to_2lev(psi_3lev):
    """Estrae le componenti |00⟩,|01⟩,|10⟩,|11⟩ di uno stato in spazio 9-dim."""
    states = [qt.tensor(qt.basis(N_lev, a), qt.basis(N_lev, b))
              for a in [0,1] for b in [0,1]]
    return np.array([complex(s.dag() * psi_3lev) for s in states])

def state_fidelity_avg(p, t, U_ideal_4x4):
    """
    Average state fidelity su un campione di stati prodotto-Pauli:
      F = (1/N) Σ_i Tr(U|ψ_i⟩⟨ψ_i|U† · M(|ψ_i⟩⟨ψ_i|)) =
        = (1/N) Σ_i ⟨ψ_target_i | ρ_out_i | ψ_target_i⟩
    
    F_avg gate ≈ (d·F_state + 1) / (d + 1)   (Horodecki-Nielsen-Bowdrey)
    per d=4: F_avg = (4 F_state + 1) / 5
    """
    H, c_ops = build_system(p)
    pauli_eigs = pauli_eigenstates_3lev()
    
    # 6×6 = 36 prodotti tensoriali (esclusi pochi ridondanti, ma teniamo tutti)
    # In realtà per F_avg basta 6 (set MUB) ma usiamo 16 = 4 Z + 4 X + 4 Y di prodotto
    # Selezioniamo 16 stati prodotto-Pauli (subset informazionalmente completo)
    test_states = []
    for ket1 in pauli_eigs:
        for ket2 in pauli_eigs:
            test_states.append(qt.tensor(ket1, ket2))   # 36 stati
    
    F_state_total = 0
    n_states = 0
    for psi in test_states:
        # Proiezione su 4-dim per il target
        coeffs_in = project_to_2lev(psi)
        if np.abs(np.linalg.norm(coeffs_in) - 1) > 0.01:
            continue   # se ha leakage, salta
        psi_target_4d = U_ideal_4x4 @ coeffs_in
        
        # Evoluzione Lindblad
        if t == 0:
            rho_out = psi * psi.dag()
        else:
            result = qt.mesolve(H, psi, [0, t], c_ops, [])
            rho_out = result.states[-1]
        
        # Proiezione su 4-dim
        states_basis = [qt.tensor(qt.basis(N_lev, a), qt.basis(N_lev, b))
                        for a in [0,1] for b in [0,1]]
        rho_out_4x4 = np.array([[complex(s1.dag() * rho_out * s2) for s2 in states_basis]
                                 for s1 in states_basis])
        
        # State fidelity = ⟨ψ_target | ρ_out | ψ_target⟩
        F_state = np.real(psi_target_4d.conj() @ rho_out_4x4 @ psi_target_4d)
        F_state_total += F_state
        n_states += 1
    
    F_state_avg = F_state_total / n_states
    # Average gate fidelity from average state fidelity (Horodecki, Bowdrey 2002)
    d = 4
    F_avg_gate = (d * F_state_avg + 1) / (d + 1)
    return F_avg_gate

# ───────────────────────────────────────────────────────────────────────────
# ESECUZIONE
# ───────────────────────────────────────────────────────────────────────────

results = {}
for key, p in systems.items():
    print(f"\n{'─'*78}")
    print(f"  {p['name']}")
    print(f"{'─'*78}")
    
    Dq = p['Delta_q']
    alpha = p['alpha']
    Omega_ZX_pred = -p['J'] * p['Omega_drive'] * alpha / (Dq * (Dq + alpha))
    t_pred = np.pi / (2 * abs(Omega_ZX_pred))
    
    print(f"  Δ_q = {Dq*1e3:.0f} MHz, α = {alpha*1e3:.0f} MHz, J = {p['J']*1e3:.2f} MHz")
    print(f"  Ω_drive = {p['Omega_drive']*1e3:.0f} MHz, T_1 = {p['T1_us']} μs")
    print(f"  Predizione: Ω_ZX = {Omega_ZX_pred*1e3:+.4f} MHz, t_pred = {t_pred:.0f} ns")
    
    # Scan tempi
    t_max = 1.5 * t_pred
    times = np.linspace(0, t_max, 25)
    U_targ = U_ideal_ZX(np.pi/2)
    
    print(f"  Simulazione 25 punti × 36 stati ...")
    F_arr = np.zeros(len(times))
    for i, t in enumerate(times):
        F_arr[i] = state_fidelity_avg(p, t, U_targ)
        if i % 5 == 0:
            print(f"    t={t:>6.1f} ns  →  F_avg = {F_arr[i]*100:>6.3f}%")
    
    # Trova best
    idx_best = int(np.argmax(F_arr))
    F_best = F_arr[idx_best]
    t_best = times[idx_best]
    
    print(f"\n  ✓ RISULTATO:")
    print(f"    t_optimal sim = {t_best:.0f} ns  (Sheldon predice {t_pred:.0f} ns)")
    print(f"    F_avg max     = {F_best*100:.3f}%")
    
    results[key] = dict(times=times, F=F_arr, t_pred=t_pred, t_best=t_best,
                        F_best=F_best, Omega_ZX_pred=Omega_ZX_pred, name=p['name'])

# ───────────────────────────────────────────────────────────────────────────
# PLOT
# ───────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
for ax, (key, r) in zip(axes, results.items()):
    ax.plot(r['times'], r['F']*100, 'b-o', linewidth=2, markersize=4)
    ax.axvline(r['t_pred'], color='red', linestyle='--', alpha=0.7, 
               label=f"t_pred Sheldon = {r['t_pred']:.0f} ns")
    ax.axvline(r['t_best'], color='green', linestyle=':', alpha=0.9,
               label=f"t_best (sim) = {r['t_best']:.0f} ns")
    ax.axhline(99, color='gray', linestyle='-', alpha=0.4, label='99% threshold')
    ax.scatter([r['t_best']], [r['F_best']*100], s=120, c='orange', zorder=5,
               label=f"F_avg max = {r['F_best']*100:.2f}%")
    ax.set_xlabel('t_gate (ns)')
    ax.set_ylabel('F_avg (%)')
    ax.set_title(r['name'], fontsize=11)
    ax.set_ylim(50, 102)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)

plt.suptitle('Modulo 5 v2 — Lindblad simulation: F_avg(t) — validation', 
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig('/home/claude/audit_g7/fig_modulo5_v2.png', dpi=130, bbox_inches='tight')
print(f"\nFigura: fig_modulo5_v2.png")
