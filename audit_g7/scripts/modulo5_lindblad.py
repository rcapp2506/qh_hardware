"""
Modulo 5 — Validazione Lindblad del cross-resonance gate per HEATS-Q
=====================================================================

Simula via mesolve di QuTiP 5.2 il gate CR per:
  Sistema A: baseline mK   (Δ_q=280 MHz, T_1=50 μs, T_2=30 μs, n̄_th=0)
  Sistema B: innovation 300 GHz @ 4K (Δ_q=600 MHz, T_1=25 μs, T_2=12 μs, n̄_th=0.028)

Modello: 2 transmon a 3 livelli (Duffing) coupled da J transverse exchange
(post-SW della cavità, vedi Modulo 1) + drive cross-resonance sul controllo.

Frame: rotating frame del target Q2.

Output:
- Ω_ZX dinamico (estratto da fit dei livelli)
- t_gate ottimo (massima F_avg)
- F_avg al sweet-spot
- Confronto con predizioni analitiche dei moduli precedenti
"""

import numpy as np
import qutip as qt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("="*78)
print("MODULO 5 — Validazione Lindblad del cross-resonance gate")
print("="*78)

# ───────────────────────────────────────────────────────────────────────────
# DEFINIZIONE DEI DUE SISTEMI
# ───────────────────────────────────────────────────────────────────────────

systems = {
    'A_baseline_mK': {
        'name': 'Sistema A — baseline mK',
        'Delta_q':     0.280,    # GHz qubit-qubit detuning (re-design)
        'alpha':      -0.200,    # GHz anharmonicity transmon
        'J':          -7.03e-3,  # GHz transverse exchange (dal Modulo 1)
        'Omega_drive': 0.060,    # 60 MHz, < min(α/3, Δq/3) = min(67, 93) = 67
        'T1_us':       50,
        'T2_us':       30,
        'n_th':        0.0,
        'cap_t_ns':    1500,
    },
    'B_innovation_300GHz': {
        'name': 'Sistema B — innovation 300 GHz @ 4K',
        'Delta_q':     0.600,    # 600 MHz (re-design da 20 GHz)
        'alpha':      -1.000,    # GHz
        'J':          -9.40e-3,  # GHz
        'Omega_drive': 0.180,    # 180 MHz, < min(α/3, Δq/3) = min(333, 200) = 200
        'T1_us':       25,
        'T2_us':       12,
        'n_th':        0.028,
        'cap_t_ns':    600,
    }
}

N_lev = 3   # 3 livelli per ciascun transmon (cattura leakage |2⟩)

def build_system(p):
    """Costruisce H, c_ops nel rotating frame del target."""
    Dq = p['Delta_q']
    alpha = p['alpha']
    
    b1 = qt.tensor(qt.destroy(N_lev), qt.qeye(N_lev))
    b2 = qt.tensor(qt.qeye(N_lev), qt.destroy(N_lev))
    
    # Hamiltoniano nel frame rotante del target (ω_drive = ω_q2)
    # Q1 ha frequenza relativa Δ_q, Q2 ha frequenza 0
    H_static = (Dq * b1.dag()*b1
                + 0.5*alpha * b1.dag()*b1.dag()*b1*b1
                + 0.5*alpha * b2.dag()*b2.dag()*b2*b2
                + p['J'] * (b1.dag()*b2 + b1*b2.dag()))
    
    # Drive sul controllo (Q1) alla frequenza del target
    # Nel frame rotante: drive è statico, ampiezza Ω
    H_drive = 0.5 * p['Omega_drive'] * (b1 + b1.dag())
    
    H = H_static + H_drive
    
    # Collapse operators (rates in GHz units, perché H è in GHz)
    gamma_1 = 1e-3 / p['T1_us']   # 1/T1[μs] in GHz
    Tphi_inv = max(0, 1/p['T2_us'] - 1/(2*p['T1_us']))
    gamma_phi = 1e-3 * Tphi_inv   # in GHz
    n_th = p['n_th']
    
    c_ops = []
    for b in [b1, b2]:
        # Relaxation: γ_1·(1+n̄)
        c_ops.append(np.sqrt(gamma_1 * (1 + n_th)) * b)
        # Thermal excitation: γ_1·n̄
        if n_th > 0:
            c_ops.append(np.sqrt(gamma_1 * n_th) * b.dag())
        # Pure dephasing: convertito in collapse via 2γ_φ b†b
        if gamma_phi > 0:
            c_ops.append(np.sqrt(2 * gamma_phi) * (b.dag() * b))
    
    return H, c_ops, b1, b2

# ───────────────────────────────────────────────────────────────────────────
# FIDELITY UTILITIES
# ───────────────────────────────────────────────────────────────────────────

def projector_to_computational():
    """Proiettore sul sottospazio computazionale 4D (livelli 0,1 dei due qubit)."""
    P1 = qt.basis(N_lev,0)*qt.basis(N_lev,0).dag() + qt.basis(N_lev,1)*qt.basis(N_lev,1).dag()
    return qt.tensor(P1, P1)

def truncate_to_4x4(rho_full):
    """Estrae blocco 4x4 (computazionale) di rho 9x9."""
    states = [qt.tensor(qt.basis(N_lev,i), qt.basis(N_lev,j)) for i in [0,1] for j in [0,1]]
    M = np.array([[complex(s1.dag() * rho_full * s2) for s2 in states] for s1 in states])
    return M

def U_ideal_ZX(theta=np.pi/2):
    """Ideal ZX(θ) = exp(-i θ/2 σ_z ⊗ σ_x). Notation: |gg⟩,|ge⟩,|eg⟩,|ee⟩.
    qutip basis: σ_z|0⟩=+|0⟩, |0⟩=|g⟩."""
    import scipy.linalg as sla
    Z = np.array([[1,0],[0,-1]])
    X = np.array([[0,1],[1,0]])
    ZX = np.kron(Z, X)
    return sla.expm(-1j * theta/2 * ZX)

def average_gate_fidelity_2q(M_actual, U_ideal):
    """
    Average gate fidelity per un canale 2-qubit.
    F_avg = (Tr(U_ideal† · Λ(U_ideal) ) + d) / (d² + d)  ma è equivalente a:
    F_avg = (Tr(M_actual · M_ideal†) + d) / (d² + d)  per quasi-unitary
    
    Più semplice: per stati puri input, F_avg ≈ ⟨ψ_out|ψ_target⟩|² mediato.
    
    Usiamo la formula da Bowdrey 2002 / Nielsen 2002:
    F_avg = (1/d²) Σ_i Tr(U_ideal · σ_i · U_ideal† · M(σ_i)) / d + 1/(d+1)
    
    Per noi è più diretto fare tomografia su 4 stati input:
    state |i⟩ → ρ_target_i = U|i⟩⟨i|U†, ρ_actual_i = M(|i⟩⟨i|)
    F_state = ⟨ψ_target|ρ_actual|ψ_target⟩
    F_avg ≈ media (più correzioni)
    """
    # Costruisce process matrix da U_actual (assumendo channel è quasi-unitary)
    d = 4
    U_actual = M_actual   # se input era |ψ⟩⟨ψ| e output è ρ, M è il propagatore
    # Average gate fidelity (Horodecki 1999):
    # F_avg = (|Tr(U_ideal† · U_actual)|² + d) / (d² + d)
    overlap = np.abs(np.trace(U_ideal.conj().T @ U_actual))**2
    F_avg = (overlap + d) / (d*(d + 1))
    return F_avg

# ───────────────────────────────────────────────────────────────────────────
# SIMULAZIONE: estrai propagatore 4x4 (con dissipazione) e fidelity
# ───────────────────────────────────────────────────────────────────────────

def simulate_gate(p, n_times=80, max_t=None):
    """
    Simula il sistema, propaga la matrice densità per ognuno dei 16 elementi
    di base del sottospazio computazionale, ricostruisce il superoperatore
    Λ(t), e calcola F_avg(t) vs ZX(π/2).
    """
    H, c_ops, b1, b2 = build_system(p)
    Dq = p['Delta_q']
    alpha = p['alpha']
    
    # Predizione analitica del Ω_ZX (Sheldon)
    Omega_ZX_pred = -p['J'] * p['Omega_drive'] * alpha / (Dq * (Dq + alpha))
    t_pred = np.pi / (2 * abs(Omega_ZX_pred))   # ns
    
    if max_t is None:
        max_t = min(2 * t_pred, p['cap_t_ns'])
    times = np.linspace(0, max_t, n_times)
    
    # Target unitary (4x4 = 2 qubits)
    U_targ_np = U_ideal_ZX(np.pi/2)
    # Convert to qutip Qobj with dims for 2-qubit operator
    U_targ = qt.Qobj(U_targ_np, dims=[[2,2],[2,2]])
    
    # Costruisco i 16 elementi base |i⟩⟨j| del sottospazio computazionale (in 9x9 ambient)
    states_basis = [qt.tensor(qt.basis(N_lev, a), qt.basis(N_lev, b)) for a in [0,1] for b in [0,1]]
    
    F_vs_t = np.zeros(n_times)
    L = qt.liouvillian(H, c_ops)
    
    for k, t in enumerate(times):
        if t == 0:
            # F_avg(identity, U_ZX) per d=4 (cfr Nielsen 2002): (|Tr(U)|²+d)/(d²+d)
            tr = np.trace(U_targ_np)
            F_vs_t[k] = (np.abs(tr)**2 + 4) / (16 + 4)
            continue
        super_op = (L * t).expm()
        # Ricostruzione del 4x4 superoperatore (Choi-like) 
        # via propagazione degli E_ij = |i⟩⟨j|
        chi = np.zeros((4, 4, 4, 4), dtype=complex)
        for i, si in enumerate(states_basis):
            for j, sj in enumerate(states_basis):
                rho_in = si * sj.dag()
                rho_out = qt.vector_to_operator(super_op * qt.operator_to_vector(rho_in))
                M_out = truncate_to_4x4(rho_out)
                chi[i, j, :, :] = M_out
        
        # Ricostruisco superoperatore 16x16 (in colonna-stack basis)
        # super_M[(k,l),(i,j)] = ⟨k|M(|i⟩⟨j|)|l⟩ = chi[i,j,k,l]
        super_M = np.zeros((16, 16), dtype=complex)
        for i in range(4):
            for j in range(4):
                for k_ in range(4):
                    for l_ in range(4):
                        super_M[k_*4+l_, i*4+j] = chi[i,j,k_,l_]
        
        # Convert to Qobj as supermatrix (column-stacking convention)
        super_qobj = qt.Qobj(super_M, dims=[[[2,2],[2,2]], [[2,2],[2,2]]])
        
        # Compute average gate fidelity
        try:
            F = qt.average_gate_fidelity(super_qobj, U_targ)
        except Exception as e:
            F = np.nan
        F_vs_t[k] = float(np.real(F))
    
    return times, F_vs_t, t_pred, Omega_ZX_pred

# ───────────────────────────────────────────────────────────────────────────
# ESECUZIONE PER ENTRAMBI I SISTEMI
# ───────────────────────────────────────────────────────────────────────────

results = {}
for key, p in systems.items():
    print(f"\n{'─'*78}")
    print(f"  {p['name']}")
    print(f"{'─'*78}")
    print(f"  Δ_q = {p['Delta_q']*1e3:.0f} MHz  (qubit-qubit detuning)")
    print(f"  α = {p['alpha']*1e3:.0f} MHz  (anarmonicità transmon)")
    print(f"  J = {p['J']*1e3:.2f} MHz  (transverse exchange)")
    print(f"  Ω_drive = {p['Omega_drive']*1e3:.0f} MHz")
    print(f"  T_1 = {p['T1_us']} μs, T_2 = {p['T2_us']} μs, n̄_th = {p['n_th']:.4f}")
    print(f"\n  Simulazione in corso (8x8 superoperator × {80} tempi)...")
    
    times, F_vs_t, t_pred, Omega_ZX_pred = simulate_gate(p, n_times=80)
    
    # Best fidelity
    idx_best = int(np.argmax(F_vs_t))
    F_best = F_vs_t[idx_best]
    t_best = times[idx_best]
    
    print(f"\n  RISULTATI:")
    print(f"  ────────────────────────────────────────────────────")
    print(f"   Predizione analitica Sheldon 2016:")
    print(f"     Ω_ZX = {Omega_ZX_pred*1e3:+.4f} MHz")
    print(f"     t_pred (ZX_π/2) = {t_pred:.0f} ns")
    print(f"")
    print(f"   Risultato simulazione Lindblad QuTiP:")
    print(f"     t_optimal = {t_best:.0f} ns")
    print(f"     F_avg(t*) = {F_best*100:.3f}%   ← validazione del numero target")
    print(f"     Discrepanza analytical/numerical t: {(t_best-t_pred)/t_pred*100:+.1f}%")
    
    results[key] = {
        'times': times, 'F_vs_t': F_vs_t,
        't_pred': t_pred, 'Omega_ZX_pred': Omega_ZX_pred,
        't_best': t_best, 'F_best': F_best,
        'name': p['name']
    }

# ───────────────────────────────────────────────────────────────────────────
# PLOT
# ───────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
for ax, (key, r) in zip(axes, results.items()):
    ax.plot(r['times'], r['F_vs_t']*100, 'b-', linewidth=2)
    ax.axvline(r['t_pred'], color='red', linestyle='--', alpha=0.7, 
               label=f"t_pred Sheldon = {r['t_pred']:.0f} ns")
    ax.axvline(r['t_best'], color='green', linestyle=':', alpha=0.9,
               label=f"t_best (sim.) = {r['t_best']:.0f} ns")
    ax.axhline(99, color='gray', linestyle='-', alpha=0.4)
    ax.scatter([r['t_best']], [r['F_best']*100], s=100, c='orange', zorder=5,
               label=f"F_max = {r['F_best']*100:.2f}%")
    ax.set_xlabel('t_gate (ns)')
    ax.set_ylabel('F_avg (%)')
    ax.set_title(r['name'], fontsize=10)
    ax.set_ylim(85, 100)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)

plt.suptitle('Modulo 5 — Lindblad simulation: F_avg(t) per CR gate', fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig('/home/claude/audit_g7/fig_modulo5_lindblad_F_vs_t.png', dpi=130, bbox_inches='tight')
print(f"\nFigura salvata: fig_modulo5_lindblad_F_vs_t.png")
