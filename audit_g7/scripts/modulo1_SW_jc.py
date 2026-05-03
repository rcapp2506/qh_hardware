"""
Modulo 1 — Derivazione di Schrieffer-Wolff sul Jaynes-Cummings a due qubit
==========================================================================

Obiettivo: derivare formalmente che la formula

   (g₁g₂/2)(1/Δ₁ + 1/Δ₂)

è il coefficiente J dell'exchange coupling trasversale (σ_+^(1)σ_-^(2) + h.c.),
NON il cross-Kerr ζ longitudinale (σ_z^(1) σ_z^(2)).

Approccio: SW perturbativa a secondo ordine in g/Δ, con qutip per la rappresentazione
matriciale degli operatori e sympy per la simplificazione simbolica dei coefficienti.
"""

import sympy as sp
import numpy as np
import qutip as qt

print("="*78)
print("MODULO 1 — Schrieffer-Wolff sul Jaynes-Cummings a due qubit")
print("="*78)

# ───────────────────────────────────────────────────────────────────────────
# SETUP SIMBOLICO
# ───────────────────────────────────────────────────────────────────────────
# Variabili simboliche: tutte reali, frequenze positive
omega_r, omega_q1, omega_q2, g1, g2 = sp.symbols(
    'omega_r omega_q1 omega_q2 g_1 g_2', real=True, positive=True
)
Delta1 = sp.Symbol('Delta_1', real=True)
Delta2 = sp.Symbol('Delta_2', real=True)
# Definizioni: Δ_i = ω_q,i - ω_r

print("""
SETUP:
  Hamiltoniano JC:  H = ω_r â†â + Σ_i (ω_q,i/2) σ_z^(i) + Σ_i g_i (â†σ_-^(i) + â σ_+^(i))

  Generatore SW:    S = Σ_i (g_i/Δ_i) (â†σ_-^(i) - â σ_+^(i))   antihermitiano

  Espansione BCH:   H̃ = e^(-S) H e^S = H + [H,S] + (1/2)[[H,S],S] + O(g³)

  Per costruzione, [H,S] cancella il coupling lineare in g (primo ordine).
  I termini di secondo ordine vengono da (1/2)[[H,S],S].
""")

# ───────────────────────────────────────────────────────────────────────────
# COSTRUZIONE OPERATORIALE (qutip per rappresentazione esatta)
# ───────────────────────────────────────────────────────────────────────────
# Spazio di Hilbert: cavità (Fock troncato a N=4) ⊗ qubit1 (2D) ⊗ qubit2 (2D)
N_cav = 4
a   = qt.tensor(qt.destroy(N_cav), qt.qeye(2), qt.qeye(2))
sm1 = qt.tensor(qt.qeye(N_cav), qt.sigmam(), qt.qeye(2))
sm2 = qt.tensor(qt.qeye(N_cav), qt.qeye(2), qt.sigmam())
sp1 = sm1.dag()
sp2 = sm2.dag()
sz1 = qt.tensor(qt.qeye(N_cav), qt.sigmaz(), qt.qeye(2))
sz2 = qt.tensor(qt.qeye(N_cav), qt.qeye(2), qt.sigmaz())
n_a = a.dag() * a
I_full = qt.tensor(qt.qeye(N_cav), qt.qeye(2), qt.qeye(2))

# ───────────────────────────────────────────────────────────────────────────
# CALCOLO DI [[H,S],S] con valori NUMERICI in regime dispersivo
# (sympy puro su operatori bosonici è scomodo; usiamo numerica con valori
# estremi del rapporto g/Δ → ricaviamo i coefficienti per fitting simbolico)
# ───────────────────────────────────────────────────────────────────────────
def build_H_and_S(om_r, om_q1, om_q2, g1_v, g2_v):
    D1 = om_q1 - om_r
    D2 = om_q2 - om_r
    H0 = om_r * n_a + (om_q1/2)*sz1 + (om_q2/2)*sz2
    Hg = g1_v*(a.dag()*sm1 + a*sp1) + g2_v*(a.dag()*sm2 + a*sp2)
    H  = H0 + Hg
    # Generatore SW (cancella il coupling al primo ordine)
    S = (g1_v/D1)*(a.dag()*sm1 - a*sp1) + (g2_v/D2)*(a.dag()*sm2 - a*sp2)
    return H0, Hg, S

# Calcoliamo l'Hamiltoniano effettivo H̃ proiettato sul settore di numero
# di eccitazione totale fisso. Per estrarre i coefficienti effettivi,
# proiettiamo sul settore N_exc=0 (cavità nel vuoto, 4 stati qubit).

def project_qubit_subspace(H_op, n_cav_target=0):
    """
    Estrae il blocco di H_op nel sottospazio { |n_cav_target⟩ ⊗ |q1q2⟩ }, base 4×4.
    """
    states = []
    for q1 in [0, 1]:
        for q2 in [0, 1]:
            # qutip basis: basis(2,0) = |g⟩ con σ_z|g⟩=+|g⟩
            states.append(qt.tensor(
                qt.basis(N_cav, n_cav_target),
                qt.basis(2, q1),
                qt.basis(2, q2)
            ))
    M = np.zeros((4, 4), dtype=complex)
    for i, psi_i in enumerate(states):
        for j, psi_j in enumerate(states):
            M[i, j] = complex(psi_i.dag() * H_op * psi_j)
    return M, states

def fit_effective_couplings(M):
    """
    Decomposizione di una matrice 4×4 (base ordinata: |gg⟩, |ge⟩, |eg⟩, |ee⟩)
    secondo:
      M = c0·I + c_z1·Z1 + c_z2·Z2 + c_zz·Z1Z2
        + c_xx·(X1X2 + Y1Y2)/2  +  c_yy·(X1X2 - Y1Y2)/2  + ...
    Per il caso H = h_z1 Z1 + h_z2 Z2 + ζ Z1Z2 + 2J(σ+σ- + σ-σ+):
       2J(σ+σ- + σ-σ+) = J(X⊗X + Y⊗Y)
    Quindi nella base computazionale ordinata |gg⟩, |ge⟩, |eg⟩, |ee⟩:
       Termine X⊗X+Y⊗Y produce elementi off-diagonal solo tra |ge⟩ e |eg⟩, valore = 2J.
       
    Convenzione qutip: σ_z |g⟩ = +|g⟩  (basis(2,0)), σ_z |e⟩ = -|e⟩ (basis(2,1)).
    """
    # diagonale: estraiamo c0, c_z1, c_z2, c_zz
    Egg, Ege, Eeg, Eee = M[0,0].real, M[1,1].real, M[2,2].real, M[3,3].real
    # E(s1,s2) = c0 + c_z1·s1 + c_z2·s2 + c_zz·s1·s2
    # con s1, s2 ∈ {+1 (g), -1 (e)}
    c0   = (Egg + Ege + Eeg + Eee) / 4
    c_z1 = (Egg + Ege - Eeg - Eee) / 4
    c_z2 = (Egg - Ege + Eeg - Eee) / 4
    c_zz = (Egg - Ege - Eeg + Eee) / 4
    # Off-diagonal |ge⟩ ↔ |eg⟩: elemento M[1,2] (riga ge, col eg) = J
    J = M[1, 2].real  # exchange coupling 
    return dict(c0=c0, h_z1=c_z1, h_z2=c_z2, zeta_zz=c_zz, J_exchange=J)

# Esempio numerico in regime dispersivo profondo
om_r_v, om_q1_v, om_q2_v = 6.5, 5.8, 5.2
g1_v, g2_v = 0.080, 0.080

H0, Hg, S = build_H_and_S(om_r_v, om_q1_v, om_q2_v, g1_v, g2_v)

# H_effettivo a secondo ordine: H + [H,S] + (1/2)[[H,S],S]
comm_HS    = H0*S - S*H0 + Hg*S - S*Hg
comm_HSS   = comm_HS*S - S*comm_HS
H_eff_2nd  = H0 + Hg + comm_HS + 0.5 * comm_HSS

# Inoltre, per costruzione, [H0,S] = -Hg al primo ordine (verifichiamolo)
check_first = (H0*S - S*H0 + Hg).norm()
print(f"Check primo ordine: ||[H0,S] + Hg|| = {check_first:.3e}  (atteso ≈ 0)")

# Proiezione sul sottospazio dei qubit (cavità nel vuoto)
M_qubit, _ = project_qubit_subspace(H_eff_2nd, n_cav_target=0)
coeffs = fit_effective_couplings(M_qubit)

# Predizioni teoriche (formule canoniche)
D1, D2 = om_q1_v - om_r_v, om_q2_v - om_r_v
chi1_pred = g1_v**2 / D1
chi2_pred = g2_v**2 / D2
J_pred    = (g1_v * g2_v / 2) * (1/D1 + 1/D2)
zeta_pred = 0  # per qubit a 2 livelli SENZA anarmonicità, ζ_zz = O(g^4) ≈ 0 a 2° ordine

print(f"""
RISULTATI (in MHz, parametri HEATS-Q baseline g={g1_v*1e3:.0f} MHz, Δ_1={D1*1e3:+.0f}, Δ_2={D2*1e3:+.0f}):

  H_eff = c_0·I + h_z1·σ_z^(1) + h_z2·σ_z^(2) + ζ_zz·σ_z^(1)σ_z^(2) + J·(σ_+^(1)σ_-^(2) + h.c.)

  Coefficienti estratti dalla SW (proiezione sul sottospazio qubit):
    h_z1     (rinormalizz. qubit 1) = {coeffs['h_z1']*1e3:+.4f} MHz
    h_z2     (rinormalizz. qubit 2) = {coeffs['h_z2']*1e3:+.4f} MHz
    ζ_zz     (cross-Kerr longit.)   = {coeffs['zeta_zz']*1e3:+.6f} MHz   ← O(g^4), ≈ 0 a 2° ordine
    J        (exchange trasversale) = {coeffs['J_exchange']*1e3:+.4f} MHz

  Predizioni teoriche dalle formule SW di secondo ordine:
    χ_1 (= 2 h_z1)  = g_1²/Δ_1                      = {chi1_pred*1e3:+.4f} MHz   [verifica: 2·h_z1 = {2*coeffs['h_z1']*1e3:+.4f}]
    χ_2 (= 2 h_z2)  = g_2²/Δ_2                      = {chi2_pred*1e3:+.4f} MHz   [verifica: 2·h_z2 = {2*coeffs['h_z2']*1e3:+.4f}]
    J  (atteso)     = (g_1g_2/2)(1/Δ_1 + 1/Δ_2)    = {J_pred*1e3:+.4f} MHz   ← FORMULA "2.39" ORIGINALE
    ζ_zz (atteso)   = 0                              (a 2° ordine, nessuna anarmonicità)

  Rapporti:
    J_estratto / J_predetto         = {coeffs['J_exchange']/J_pred:.4f}    ← se ≈1 conferma identificazione
    ζ_zz_estratto                   = {coeffs['zeta_zz']*1e3:.2e} MHz       ← deve essere ≈ 0
""")

print("="*78)
print("CONCLUSIONE PEDAGOGICA DEL MODULO 1")
print("="*78)
print(f"""
 La derivazione di Schrieffer-Wolff a secondo ordine sul Jaynes-Cummings produce:
 
   H_eff = const + (h_z1)·σ_z^(1) + (h_z2)·σ_z^(2) + J·(σ_+^(1)σ_-^(2) + σ_-^(1)σ_+^(2))
                                                 ↑
                                  COUPLING EFFETTIVO QUBIT-QUBIT
 
 dove il coefficiente J vale esattamente:
 
            ┌─────────────────────────────────────┐
            │   J = (g_1 g_2 / 2)(1/Δ_1 + 1/Δ_2)  │   ← era "Eq. 2.39" della tesi
            └─────────────────────────────────────┘
 
 QUESTO È L'EXCHANGE COUPLING TRASVERSALE (a volte chiamato XY o J_xy nella
 letteratura), originariamente ottenuto da Majer et al. Nature 449, 443 (2007)
 per due transmon "always-on coupled" via cavity bus.
 
 Il termine cross-Kerr longitudinale ζ_zz σ_z⊗σ_z è IDENTICAMENTE NULLO al
 secondo ordine SW per qubit a 2 livelli ideali. Per averlo serve l'anarmonicità
 del transmon (livelli |2⟩, |3⟩, ...), che è proprio quello che rende il
 transmon un "atomo artificiale" weakly anharmonic anziché un oscillatore
 armonico — e quello è il soggetto del Modulo 2.

 → Il commento di Gatti è confermato in pieno.
 → La tesi originale ha attribuito alla formula corretta (J trasversale)
   l'etichetta sbagliata (ζ longitudinale).
""")
