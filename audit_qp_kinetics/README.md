# QP kinetics audit — quasiparticle generation, recombination, and trapping at 4 K

This folder contains the numerical audit scripts that support the
quasiparticle appendix `chapters/appendix_quasiparticles.tex` of the
[PhDThesis](https://github.com/rcapp2506/PhDThesis) repository.

## Background

The choice of an elevated operating temperature (\$T = 4\,\$K, accessible
with a two-stage pulse-tube cryocooler without sub-Kelvin refrigeration)
imposes an intrinsic thermal floor on the quasiparticle density,
$x_{qp}^{th}(\Delta_{NbN},T) \simeq 2.7\times 10^{-4}$, which translates
into a quasiparticle-limited relaxation time $T_1^{qp,\text{intrinsic}}
\simeq 20\,\text{ns}$. Reaching the operational target $T_1 = 49\,\mu$s
required by the algorithms of Chapter 3 calls for a reduction of $x_{qp}$
by approximately three orders of magnitude below the thermal floor. The
scripts in this folder examine two engineering routes — gap engineering
and trap engineering — and conclude that only the latter is effective in
the thermal regime.

## Scripts (in `scripts/`)

| File | Module | Output |
|---|---|---|
| `qp_kinetics_squid.py` | Base model: BCS gap, thermal $x_{qp}$, SQUID asymmetric $T_1$ | three-panel diagnostic plot (`qp_kinetics_squid.png`) |
| `qp_kinetics_squid_thermal.py` | Thermal-consistent version with $\Delta(T)$ BCS, material scan (NbN, NbN thin, Nb, NbTiN, TaN) | shows the cancellation $e^{-\Delta_2/k_BT}\cdot e^{-(\Delta_1-\Delta_2)/k_BT}=e^{-\Delta_1/k_BT}$ |
| `qp_target_driven.py` | Inverse analysis: bilateral trap rates $(s_1,s_2)$ → $T_1$ map with iso-49 μs contour | confirms target requires $s\sim 10^7\,\text{s}^{-1}$ |
| `qp_analytical_solutions.py` | Closed-form solutions in three regimes (hyperbolic, exponential, Riccati) verified against numerical ODE | error max $10^{-2}$ across crossover |
| `appendix_figures.py` | Publication-quality figures for the appendix (4 PDFs) | regenerates all figures under `figures/` |

## Figures (in `figures/`)

The four figures used in the appendix:

| File | Used in appendix section |
|---|---|
| `fig_qp_generation_population.pdf` | §A.2 Generation and thermal population (Fig. 1) |
| `fig_kinetic_regimes.pdf` | §A.3 Kinetic equations (Fig. 2) |
| `fig_trap_mechanism.pdf` | §A.6 Trap engineering (Fig. 3) |
| `fig_trap_engineering.pdf` | §A.6.1 Implementation strategy (Fig. 4) |

PDFs are intended for inclusion in the thesis via
`\includegraphics{chapters/qh_figures/...}` and are also mirrored in the
PhDThesis repo under `chapters/qh_figures/`. PNG copies are provided here
for quick visual inspection without a PDF viewer.

## Reproducibility

The kinetic scripts rely only on `numpy`, `scipy`, and `matplotlib`. To
regenerate all figures:

```bash
cd scripts/
python qp_analytical_solutions.py   # verification plot
python appendix_figures.py          # the four appendix figures
```

Parameters used throughout:

| Quantity | Value | Source |
|---|---|---|
| $\Delta_{\rm NbN}$ (bulk) | 2.8 meV | Chockalingam *et al.* PRB 77, 214503 (2008) |
| $T_c$ (bulk NbN) | 16 K | same |
| Recombination time $\tau_R^{eq}$ | 100 ns (scaled from Kaplan 1976 for NbN) | Kaplan *et al.* PRB 14, 4854 (1976) |
| $T_{\rm op}$ | 4 K | cryocooler floor |
| $\omega_{01}/2\pi$ | 5 GHz | nominal transmon |
| $E_{J,\Sigma}/h$, $E_C/h$ | 25 GHz, 0.25 GHz | Cap. 2 baseline |
| Geometric asymmetry $d$ | 0.5 | Cap. 2 baseline |

## References used in the appendix

1. A. Rothwarf and B. N. Taylor, *Measurement of recombination lifetimes
   in superconductors*, **Phys. Rev. Lett. 19, 27 (1967)** — kinetic
   framework with phonon bottleneck.
2. S. B. Kaplan *et al.*, *Quasiparticle and phonon lifetimes in
   superconductors*, **Phys. Rev. B 14, 4854 (1976)** — microscopic
   calculation of $\tau_R$.
3. G. Catelani, R. J. Schoelkopf, M. H. Devoret, L. I. Glazman,
   *Relaxation and frequency shifts induced by quasiparticles in
   superconducting qubits*, **Phys. Rev. B 84, 064517 (2011)** — $T_1^{qp}$
   formula used in §A.4.
4. R.-P. Riwar *et al.*, *Normal-metal quasiparticle traps for
   superconducting qubits*, **Phys. Rev. B 94, 104516 (2016)** — trap
   model and demonstrated rates in Al.
5. A. Hosseinkhani *et al.*, *Optimal configurations for normal-metal
   traps in transmon qubits*, **Phys. Rev. Applied 8, 064028 (2017)** —
   geometric optimisation of traps.
6. G. Marchegiani, L. Amico, G. Catelani, *Quasiparticles in
   superconducting qubits with asymmetric junctions*, **PRX Quantum 3,
   040338 (2022)** — gap-asymmetry suppression, non-equilibrium regime.
7. P. Kamenov *et al.*, **Phys. Rev. Applied 21, 054030 (2024)
   [arXiv:2309.02655]** — experimental gap engineering for parity
   preservation.
8. M. McEwen *et al.*, **Phys. Rev. Lett. 133, 240601 (2024)** — gap
   engineering against cosmic-ray bursts.

## Status

Material exploratory: not part of the delivered thesis. To integrate the
appendix into the manuscript:

1. add to `main.tex`, inside the `\appendix` block:
   ```latex
   \include{chapters/appendix_quasiparticles}
   ```
2. ensure `chapters/qh_figures/` contains the four PDF figures (mirrored
   here under `figures/`);
3. add a back-pointer from Chapter 3 (Algorithms) connecting the
   algorithm-level fidelity requirements to the appendix; a draft of this
   wording is provided alongside the audit material.
