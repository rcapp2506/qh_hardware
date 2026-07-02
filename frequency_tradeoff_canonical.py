#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Canonical frequency trade-off for Sec. 2.9 (Frequency Trade-off at Elevated
Temperature) of the PhD thesis, Ch. 2.

Model (dielectric + thermal only, by design of the section):
    Gamma_diel(f)   = (omega/Q_eff) * (1 + n_th(f,T))        [dielectric decay]
    T1_diel(f)      = 1 / Gamma_diel(f)                       -> MONOTONE in f
    eps(f)          = t_CNOT * Gamma_diel(f) + 0.5 * n_th(f)  [gate error]
                      \_ decoherence term _/   \_ thermal _/
The interior optimum lives in eps(f), not in T1_diel(f):
    d/df [f (1+n)] = (1+n)(1 - x n) > 0 for all x = hf/kT > 0,
since x n = x/(e^x - 1) < 1.  Hence T1_diel has no interior maximum.

All inputs derive from canonical_params.json (SSoT):
    - Q_eff from T1_dielectric = 53.1 us at 300 GHz / 4 K,
    - t_CNOT = epsilon_decoh * T2_echo = 152.5 ns,
    - tan-delta factor-3 uncertainty band (G9) -> Q_eff/3 lower edge.

Outputs:
    - frequency_tradeoff_analysis.png  (Fig. 2.19, same filename as before)
    - LaTeX rows for tab:frequency_comparison (printed to stdout)
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ----------------------------------------------------------------- constants
h, kB = 6.62607015e-34, 1.380649e-23

ROOT = Path(__file__).resolve().parent
par = json.load(open(ROOT / "canonical_params.json"))

T = 4.0                                            # K, design temperature
f_op = 300.0                                       # GHz, operating point

def n_th(f_GHz, T=T):
    with np.errstate(over="ignore"):
        return 1.0 / np.expm1(h * f_GHz * 1e9 / (kB * T))

# --- Q_eff from the canonical dielectric budget entry (53.1 us @ 300 GHz/4 K)
T1_diel_op = par["decoherence_breakdown_innovation_us"]["T1_dielectric"] * 1e-6
Q_eff = T1_diel_op * 2 * np.pi * f_op * 1e9 * (1 + n_th(f_op))
# --- t_CNOT derived exactly as ratified: epsilon_decoh * T2_echo
t_cnot = (par["gate_fidelity"]["epsilon_decoh_percent"] / 100.0
          * par["decoherence_breakdown_innovation_us"]["T2_echo_us"] * 1e-6)

print(f"Q_eff  = {Q_eff:.3e}   (from T1_diel = {T1_diel_op*1e6:.1f} us @ 300 GHz)")
print(f"t_CNOT = {t_cnot*1e9:.1f} ns")

def gamma_diel(f_GHz, Q=Q_eff, T=T):
    return 2 * np.pi * f_GHz * 1e9 / Q * (1 + n_th(f_GHz, T))

def eps_total(f_GHz, Q=Q_eff, T=T):
    return t_cnot * gamma_diel(f_GHz, Q, T) + 0.5 * n_th(f_GHz, T)

# ------------------------------------------------------------ sanity asserts
fs = np.linspace(2, 1500, 60000)
T1 = 1.0 / gamma_diel(fs)
assert np.all(np.diff(T1) < 0), "T1_diel must be strictly decreasing"
x = 1.26
assert abs(x / np.expm1(x) - 0.499) < 5e-3          # Eq. 2.146 has no solution
eps = eps_total(fs)
i_opt = int(np.argmin(eps))
f_opt, eps_opt = fs[i_opt], eps[i_opt]
assert 400 < f_opt < 700, f_opt
eps_300 = eps_total(f_op)
print(f"f_opt  = {f_opt:.0f} GHz, eps_min = {eps_opt*100:.2f} %  (F = {100-eps_opt*100:.2f} %)")
print(f"eps(300 GHz) = {eps_300*100:.2f} %  (F = {100-eps_300*100:.2f} %)  "
      f"[full budget: eps_decoh+eps_th = 0.58+1.41 = 1.99 %]")

# tan-delta factor-3 uncertainty band (G9): Q in [Q_eff/3, Q_eff]
Q_lo = Q_eff / 3.0
eps_lo = eps_total(fs, Q=Q_lo)
f_opt_lo = fs[int(np.argmin(eps_lo))]

# --------------------------------------------------------------- the figure
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12,
                     "axes.labelsize": 11.5, "savefig.bbox": "tight"})
fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))
(ax_a, ax_b), (ax_c, ax_d) = axs
C300, COPT = "#FF6F00", "#6A1B9A"

# (a) thermal occupation
ax_a.semilogy(fs, n_th(fs), color="#C62828", lw=2.5)
ax_a.axhline(1, color="gray", ls=":", lw=1)
ax_a.axvspan(2, fs[np.searchsorted(-(n_th(fs) - 1), 0)], color="#FBE9E7",
             alpha=0.6, label=r"classical ($n_{\rm th}>1$)")
ax_a.axvline(f_op, color=C300, ls="--", lw=1.5)
ax_a.annotate(r"$\bar n_{\rm th}(300\,{\rm GHz})=0.028$",
              xy=(f_op, n_th(f_op)), xytext=(420, 0.15),
              arrowprops=dict(arrowstyle="->", color=C300), color=C300)
ax_a.set_xlabel("Frequency (GHz)"); ax_a.set_ylabel(r"$\bar n_{\rm th}$")
ax_a.set_title("(a) Thermal occupation at 4 K")
ax_a.set_xlim(0, 1000); ax_a.set_ylim(1e-5, 50)
ax_a.legend(loc="upper right", fontsize=9); ax_a.grid(alpha=0.3)

# (b) monotone dielectric-limited T1 (4 K and 20 mK)
ax_b.loglog(fs, T1 * 1e6, color="#EF5350", lw=2.5, label=r"4 K: $Q/[\omega(1+\bar n_{\rm th})]$")
T1_mK = Q_eff / (2 * np.pi * fs * 1e9 * (1 + n_th(fs, 0.020))) * 1e6
ax_b.loglog(fs, T1_mK, color="#1565C0", lw=2, ls="--", label=r"20 mK: $Q/\omega$")
ax_b.plot(f_op, T1_diel_op * 1e6, "*", ms=16, color=C300, mec="k",
          label=r"budget value 53.1 $\mu$s")
ax_b.set_xlabel("Frequency (GHz)"); ax_b.set_ylabel(r"$T_1^{\rm diel}$ ($\mu$s)")
ax_b.set_title(r"(b) $T_1^{\rm diel}(f)$ is monotone (plateau at low $f$, no maximum)")
ax_b.set_xlim(2, 1500); ax_b.legend(fontsize=9); ax_b.grid(alpha=0.3, which="both")

# (c) error budget vs frequency: the interior optimum
eps_dec = t_cnot * gamma_diel(fs)
eps_th = 0.5 * n_th(fs)
ax_c.semilogy(fs, eps_dec * 100, color="#1565C0", lw=2,
              label=r"decoherence  $t_{\rm CNOT}\,\omega(1+\bar n_{\rm th})/Q_{\rm eff}$")
ax_c.semilogy(fs, eps_th * 100, color="#C62828", lw=2,
              label=r"thermal  $\bar n_{\rm th}/2$")
ax_c.semilogy(fs, eps * 100, color="k", lw=3, label=r"total $\varepsilon(f)$")
ax_c.plot(f_opt, eps_opt * 100, "o", ms=10, color=COPT, mec="k")
ax_c.plot(f_op, eps_300 * 100, "*", ms=16, color=C300, mec="k")
ax_c.annotate(rf"formal optimum: {f_opt:.0f} GHz, $\varepsilon={eps_opt*100:.2f}\%$",
              xy=(f_opt, eps_opt * 100), xytext=(560, 6), color=COPT,
              arrowprops=dict(arrowstyle="->", color=COPT))
ax_c.annotate(rf"design point: 300 GHz, $\varepsilon={eps_300*100:.2f}\%$",
              xy=(f_op, eps_300 * 100), xytext=(90, 22), color=C300,
              arrowprops=dict(arrowstyle="->", color=C300))
ax_c.set_xlabel("Frequency (GHz)"); ax_c.set_ylabel(r"gate error $\varepsilon$ (%)")
ax_c.set_title(r"(c) Fidelity trade-off: $\varepsilon(f)=t_{\rm CNOT}\Gamma_{\rm diel}+\bar n_{\rm th}/2$")
ax_c.set_xlim(0, 1200); ax_c.set_ylim(0.2, 300)
ax_c.legend(fontsize=9, loc="upper right"); ax_c.grid(alpha=0.3)

# (d) sensitivity to the tan-delta uncertainty (G9 factor-3 band)
ax_d.semilogy(fs, eps * 100, color="k", lw=2.5,
              label=rf"$Q_{{\rm eff}}=1.0\times10^{{8}}$ ($\tan\delta=10^{{-7}}$)")
ax_d.semilogy(fs, eps_lo * 100, color="k", lw=1.5, ls="--",
              label=rf"$Q_{{\rm eff}}/3$ ($\tan\delta=3\times10^{{-7}}$)")
ax_d.fill_between(fs, eps * 100, eps_lo * 100, color="gray", alpha=0.25)
for ff, ee, col, lbl in [(f_opt, eps_opt, COPT, None), (f_opt_lo, eps_lo.min(), COPT, None)]:
    ax_d.plot(ff, ee * 100, "o", ms=8, color=col, mec="k")
ax_d.plot([f_op, f_op], [eps_300 * 100, eps_total(f_op, Q=Q_lo) * 100],
          color=C300, lw=3, solid_capstyle="round",
          label="300 GHz design point (band)")
ax_d.set_xlabel("Frequency (GHz)"); ax_d.set_ylabel(r"gate error $\varepsilon$ (%)")
ax_d.set_title(r"(d) Sensitivity to the sub-THz $\tan\delta$ extrapolation (factor 3)")
ax_d.set_xlim(0, 1200); ax_d.set_ylim(0.2, 300)
ax_d.legend(fontsize=9, loc="upper right"); ax_d.grid(alpha=0.3)

fig.tight_layout()
out = ROOT / "frequency_tradeoff_analysis.png"
fig.savefig(out, dpi=170)
print(f"saved {out}")

# -------------------------------------------------- LaTeX rows for Tab. 2.25
rows_f = [50, 100, 200, 300, 400, round(f_opt), 800]
print("\n% --- tab:frequency_comparison rows (canonical model) ---")
for f in rows_f:
    n = n_th(f)
    t1 = 1.0 / gamma_diel(f) * 1e6
    ed, et = t_cnot * gamma_diel(f) * 100, 0.5 * n * 100
    e = ed + et
    mark = r" (\emph{optimum})" if f == round(f_opt) else (r" (design)" if f == 300 else "")
    print(f"{f} & {n:.3g} & {t1:.1f} & {ed:.2f} & {et:.2f} & {e:.2f} & {100-e:.1f}{mark} \\\\")
