# Didactic figures for the thesis

This folder collects the source scripts for didactic-purpose figures
introduced in the thesis chapters during the post-defense W&C style
cleanup.

## Figures

| Script | Output | Used in |
|---|---|---|
| `fig_qubit_squid_didactic.py` | `fig_qubit_squid_didactic.pdf` | PhDThesis Cap. 2 Preliminaries (Fig. 2.1) |

## Reproducibility

```bash
cd scripts/
python fig_qubit_squid_didactic.py
```

Output goes to current working directory; copy the PDF to
`PhDThesis/chapters/qh_figures/` for thesis inclusion.

Dependencies: matplotlib, numpy. No special skill libraries required.
