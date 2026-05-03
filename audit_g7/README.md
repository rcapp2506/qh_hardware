# G7 audit scripts — closure of Gatti review comment on Eq. 2.39

This folder contains the numerical audit scripts that support the G7 revision
of `chapters/quantum_hw.tex` in the
[PhDThesis](https://github.com/rcapp2506/PhDThesis) repository.

## Background

The external reviewer (C. Gatti) noted that Eq. 2.39 of the original
manuscript,
$$\zeta = \frac{g_1 g_2}{2}\left(\frac{1}{\Delta_1}+\frac{1}{\Delta_2}\right),$$
was labelled as the coefficient of a longitudinal $\sigma_z\otimes\sigma_z$
coupling, but is in fact the **transverse exchange coupling** $J$ of
Majer *et al.* (Nature 449, 443, 2007). The correction propagates: the
free-evolution CZ gate based on this coupling is unfeasible, and the
manuscript was revised to use the **cross-resonance gate with echo
correction** (Sheldon *et al.*, PRA 93, 060302, 2016) on a redesigned
qubit-qubit detuning $\Delta_q$.

## Scripts (in `scripts/`)

| File | Module | Output |
|---|---|---|
| `modulo1_SW_jc.py` | Schrieffer-Wolff verification on Jaynes-Cummings | confirms `(g₁g₂/2)(1/Δ₁+1/Δ₂)` is transverse J, ζ_zz = 0 at 2nd order |
| `modulo2_zeta_zz.py` | True cross-Kerr from transmon anharmonicity | validates Blais 2021 Eq. 130 vs exact diagonalization |
| `modulo3_verdict_CZ.py` | Verdict on free-evolution CZ gate | t_CZ = 12.7 µs (mK), 1.7 ms (300 GHz) — unfeasible |
| `modulo4_cross_resonance.py` | Cross-resonance derivation | Sheldon 2016 ZX rate formula |
| `modulo4bis_tradeoff_Dq.py` | Trade-off analysis of Δ_q | sweet-spot Δ_q ≈ 280 MHz (baseline mK) |
| `modulo4ter_300GHz.py` | Three strategies for 300 GHz | RIP gate, re-tune, open challenge |
| `modulo4q_sweetspot_4K.py` | Sweet-spot map for F > 99% at 4 K | T₁ ≥ 16 µs at 300 GHz, redesigned Δ_q = 600 MHz |
| `modulo5_lindblad.py` | Lindblad master equation, naïve attempt | (process tomography setup, see v3 for the working version) |
| `modulo5_v2.py` | Lindblad with state-fidelity averaging | F_max ≈ 78% for bare CR |
| `modulo5_v3_echoCR.py` | **Lindblad CR vs Echo CR comparison** | F_bare = 83%, F_echo = 91% (square pulse) |

## Figures (in `figures/`)

The four PNG figures used in the revised `chapters/quantum_hw.tex`:

| File | Used in section |
|---|---|
| `g7_fig1_CZ_freevolution_unfeasible.png` | §3.7 (cavity-mediated couplings) |
| `g7_fig2_tradeoff_Dq.png` | §3.8 (gate strategy) |
| `g7_fig3_sweetspot_4K.png` | §3.10 (technological prerequisites) |
| `g7_fig4_lindblad_echoCR.png` | §3.8 (gate strategy) |

## Reproducibility

Dependencies: Python ≥ 3.10, `numpy`, `scipy`, `qutip ≥ 5.2`, `sympy`,
`matplotlib`. To reproduce the figures:

```bash
cd scripts/
python3 modulo3_verdict_CZ.py
python3 modulo4bis_tradeoff_Dq.py
python3 modulo4q_sweetspot_4K.py
python3 modulo5_v3_echoCR.py    # ~10 min for full Lindblad scan
```

## Key physical results

- The redesigned operating point **Δ_q = 280 MHz** (baseline mK) /
  **Δ_q = 600 MHz** (innovation 300 GHz) lies in the *stretched-CR sweet
  spot*, far enough from the straddling resonance Δ_q = |α| to avoid the
  divergence of static cross-Kerr ζ_zz.
- Echo CR (square drive) achieves F_avg ≈ 91% in the Lindblad simulation;
  closing the remaining gap to F > 99% requires Gaussian-Flat-Gaussian
  pulse shaping with DRAG correction (Motzoi 2009, Gambetta 2011), as
  routinely demonstrated in transmon processors (Sheldon 2016).
- The **claim of high-temperature universal computing at 4 K is preserved**.
  The required T₁ at 300 GHz, 4 K is ≥ 16 µs, corresponding to a dielectric
  loss tangent tan δ ≤ 3 × 10⁻⁸ — at the optimistic end of state-of-the-art
  microwave dielectrics but consistent with extrapolations from
  Romanenko 2020 and Read 2023.
