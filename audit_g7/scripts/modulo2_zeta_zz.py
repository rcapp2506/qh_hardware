"""
Modulo 2 — Cross-Kerr ζ_zz dal transmon (anarmonicità)
=======================================================

Approccio "double Schrieffer-Wolff":
  PRIMA SW: elimina i gradi di libertà cavità → H_eff (a) sui due transmon
            con J exchange (Modulo 1) e auto-Kerr (anarmonicità) ereditati.
  SECONDA SW: nel sottospazio computazionale dei due qubit, elimina J usando
              il qubit-qubit detuning Δ_q = ω_1 − ω_2 come energy gap.
              Genera ζ_zz al secondo ordine in J/Δ_q, modulato dall'anarmonicità.

Risultato canonico (Blais 2021 RMP §IV.B, anche Yan 2018, Krinner 2020):
    
    ζ_zz = 2 J² [ 1/(Δ_q - α_2) - 1/(Δ_q + α_1) ]

Verifichiamolo numericamente con qutip su transmon a 3 livelli.
"""

import numpy as np
import qutip as qt

print("="*78)
print("MODULO 2 — Cross-Kerr ζ_zz vero (transmon con anarmonicità)")
print("="*78)

# Parametri HEATS-Q baseline (GHz)
omega_r  = 6.5
omega_q1 = 5.8
omega_q2 = 5.2
g1 = 0.080
g2 = 0.080
alpha_1 = -0.200
alpha_2 = -0.200

Delta_1 = omega_q1 - omega_r
Delta_2 = omega_q2 - omega_r
Delta_q = omega_q1 - omega_q2  # qubit-qubit detuning

# J trasversale dal Modulo 1
J = (g1 * g2 / 2) * (1/Delta_1 + 1/Delta_2)

# Formula canonica per ζ_zz (Blais 2021)
zeta_blais = 2 * J**2 * (1/(Delta_q - alpha_2) - 1/(Delta_q + alpha_1))

print(f"""
Parametri:
  ω_r = {omega_r} GHz, ω_q1 = {omega_q1}, ω_q2 = {omega_q2}, Δ_q = {Delta_q*1e3:+.0f} MHz
  g_1 = g_2 = {g1*1e3:.0f} MHz
  α_1 = α_2 = {alpha_1*1e3:.0f} MHz   (anarmonicità transmon)

PREDIZIONI ANALITICHE:
  J trasversale (dal Modulo 1)                       = {J*1e3:+.4f} MHz
  ζ_zz canonica (Blais 2021):
     2J² [1/(Δ_q − α_2) − 1/(Δ_q + α_1)]             = {zeta_blais*1e3:+.4f} MHz
""")

# ───────────────────────────────────────────────────────────────────────────
# VERIFICA NUMERICA: diagonalizzazione esatta di transmon (3 lev) + cavità
# ───────────────────────────────────────────────────────────────────────────
N_cav = 6
N_q   = 4   # 4 livelli per ciascun transmon (per catturare leakage virtuale)

a   = qt.tensor(qt.destroy(N_cav), qt.qeye(N_q),    qt.qeye(N_q))
b1  = qt.tensor(qt.qeye(N_cav),    qt.destroy(N_q), qt.qeye(N_q))
b2  = qt.tensor(qt.qeye(N_cav),    qt.qeye(N_q),    qt.destroy(N_q))

# Hamiltoniano transmon Duffing: H_t = ω b†b + (α/2) b†b†bb
H_cav = omega_r * a.dag() * a
H_q1  = omega_q1 * b1.dag()*b1 + 0.5*alpha_1 * b1.dag()*b1.dag()*b1*b1
H_q2  = omega_q2 * b2.dag()*b2 + 0.5*alpha_2 * b2.dag()*b2.dag()*b2*b2

# Coupling charge-cavity: g(a + a†)(b + b†)
H_int = g1 * (a + a.dag()) * (b1 + b1.dag()) \
      + g2 * (a + a.dag()) * (b2 + b2.dag())

H = H_cav + H_q1 + H_q2 + H_int
evals, evecs = H.eigenstates()

def find_state(n_c, q1_lev, q2_lev):
    psi = qt.tensor(qt.basis(N_cav, n_c), qt.basis(N_q, q1_lev), qt.basis(N_q, q2_lev))
    overlaps = np.array([abs(complex(psi.dag() * v)) for v in evecs])
    idx = int(np.argmax(overlaps))
    return evals[idx], overlaps[idx]

E00, ov00 = find_state(0, 0, 0)
E01, ov01 = find_state(0, 0, 1)
E10, ov10 = find_state(0, 1, 0)
E11, ov11 = find_state(0, 1, 1)

# Definizione standard: ζ_zz = E11 - E10 - E01 + E00
zeta_zz_exact = E11 - E10 - E01 + E00

print(f"""
DIAGONALIZZAZIONE ESATTA (transmon 4 livelli + cavità Fock-6):
  E_|00⟩ = {E00:>10.6f} GHz  (overlap {ov00:.4f})
  E_|01⟩ = {E01:>10.6f} GHz  (overlap {ov01:.4f})
  E_|10⟩ = {E10:>10.6f} GHz  (overlap {ov10:.4f})
  E_|11⟩ = {E11:>10.6f} GHz  (overlap {ov11:.4f})

  ζ_zz_esatto = E11 - E10 - E01 + E00 = {zeta_zz_exact*1e3:+.4f} MHz

CONFRONTO:
  Predizione Blais 2021                             = {zeta_blais*1e3:+.4f} MHz
  Estratto da diagonalizzazione esatta              = {zeta_zz_exact*1e3:+.4f} MHz
  Rapporto esatto/Blais                             = {zeta_zz_exact/zeta_blais:.4f}
""")

# ───────────────────────────────────────────────────────────────────────────
# TEST DI ROBUSTEZZA: scaling con anarmonicità
# ───────────────────────────────────────────────────────────────────────────
print("="*78)
print("TEST DI ROBUSTEZZA: scaling con anarmonicità α")
print("="*78)
print(f"\n{'α (MHz)':>10s} {'ζ_zz esatto (MHz)':>22s} {'Blais 2021 (MHz)':>20s} {'rapporto':>10s}")
for alpha_test_MHz in [-1000, -500, -300, -200, -150, -100, -50]:
    a_t = alpha_test_MHz * 1e-3
    H_q1t = omega_q1 * b1.dag()*b1 + 0.5*a_t * b1.dag()*b1.dag()*b1*b1
    H_q2t = omega_q2 * b2.dag()*b2 + 0.5*a_t * b2.dag()*b2.dag()*b2*b2
    H_t = H_cav + H_q1t + H_q2t + H_int
    ev, evec = H_t.eigenstates()
    
    def fs(n_c, q1, q2):
        psi = qt.tensor(qt.basis(N_cav, n_c), qt.basis(N_q, q1), qt.basis(N_q, q2))
        return ev[int(np.argmax([abs(complex(psi.dag()*v)) for v in evec]))]
    
    zz_e = fs(0,1,1) - fs(0,1,0) - fs(0,0,1) + fs(0,0,0)
    zz_b = 2 * J**2 * (1/(Delta_q - a_t) - 1/(Delta_q + a_t))
    print(f"{alpha_test_MHz:>10d} {zz_e*1e3:>22.4f} {zz_b*1e3:>20.4f} {zz_e/zz_b if zz_b!=0 else float('inf'):>10.4f}")

# ───────────────────────────────────────────────────────────────────────────
# CASO INNOVATION 300 GHz (per coerenza con il capitolo)
# ───────────────────────────────────────────────────────────────────────────
print("\n" + "="*78)
print("CASO INNOVATION (300 GHz transmon @ 4 K)")
print("="*78)
omega_r_inn  = 320.0
omega_q1_inn = 300.0
omega_q2_inn = 280.0
g_inn        = 0.500   # 500 MHz
alpha_inn    = -1.0    # 1 GHz anharmonicity (citato nel testo)

D1_inn = omega_q1_inn - omega_r_inn
D2_inn = omega_q2_inn - omega_r_inn
Dq_inn = omega_q1_inn - omega_q2_inn
J_inn  = (g_inn**2 / 2) * (1/D1_inn + 1/D2_inn)
zeta_inn = 2 * J_inn**2 * (1/(Dq_inn - alpha_inn) - 1/(Dq_inn + alpha_inn))
t_CZ_inn = np.pi / (2 * abs(zeta_inn) * 1e9) * 1e9  # in ns (zeta in GHz)

print(f"""
  ω_r = {omega_r_inn} GHz, ω_q1 = {omega_q1_inn}, ω_q2 = {omega_q2_inn}, Δ_q = {Dq_inn*1e3:+.0f} MHz
  g = {g_inn*1e3:.0f} MHz, α = {alpha_inn*1e3:.0f} MHz

  J transverse exchange    = {J_inn*1e3:+.4f} MHz
  ζ_zz cross-Kerr (Blais)  = {zeta_inn*1e3:+.4f} MHz
  
  t_CZ free-evolution      = π/(2|ζ_zz|) = {t_CZ_inn:.0f} ns
""")
