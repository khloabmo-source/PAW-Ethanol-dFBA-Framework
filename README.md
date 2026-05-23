# Figure-generation scripts and data

This folder contains the code and data used to produce **every figure and the
graphical abstract** in the manuscript. All figures are generated programmatically
with **matplotlib** from numerical simulation outputs.

## Environment

```
Python 3.10+
numpy, pandas, matplotlib
(for the full dFBA pipeline: cobra==0.29.1, optlang; see the main repository)
```

Install plotting dependencies:
```
pip install numpy pandas matplotlib
```

## Which script produces which figure

| Script | Produces | Reads from `data/` |
|---|---|---|
| `make_results_figures.py` | **Figure 1** (cross-chassis dose-response), **Figure 2** (controller decomposition), **Figure 3** (carbon redistribution & yield), **Figure 4** (Monte Carlo & safety), **Figure 5** (concordance & benchmarking) | `sweep_v4_combined.csv`, `scenarios.csv`, `mc_full.csv`, `mc_KO11.csv`, `mc_LY01.csv` |
| `adversarial_simulation.py` | Runs the AS-1/AS-2/AS-3 surrogate simulations → writes `adversarial_sweeps.csv`, `adversarial_summary.csv` | (self-contained) |
| `redraw_figure6.py` | **Figure 6** (adversarial sensitivity, overlap-corrected) | `adversarial_sweeps.csv` |
| `controller_benchmark.py` | Runs single-enzyme vs redox-partitioning benchmark → writes `controller_benchmark_sweeps.csv`, `controller_benchmark_summary.csv` | (self-contained) |
| `make_figure7.py` | **Figure 7** (controller benchmark, 3 panels) | `controller_benchmark_sweeps.csv` |
| `make_graphical_abstract_bej.py` | **Graphical abstract** (horizontal, 2.5:1, for Biochemical Engineering Journal) | (self-contained) |

## Reproduction order

```bash
# 1. Figures 1–5 (from pre-computed pipeline outputs in data/)
python make_results_figures.py

# 2. Figure 6 (adversarial) — run the simulation first, then plot
python adversarial_simulation.py        # writes data/adversarial_*.csv
python redraw_figure6.py                 # writes Figure_6_*.png

# 3. Figure 7 (controller benchmark) — run the benchmark, then plot
python controller_benchmark.py           # writes data/controller_benchmark_*.csv
python make_figure7.py                   # writes Figure_7_*.png

# 4. Graphical abstract
python make_graphical_abstract_bej.py    # writes Graphical_Abstract_BEJ.png
```

Output PNGs are written at **300 dpi on white backgrounds**, matching Elsevier's
artwork requirements.

## Notes on the data files

- `sweep_v4_combined.csv` — dose-response sweeps for all three chassis on both GSMMs (Figures 1, 3).
- `scenarios.csv` — the five-scenario controller decomposition at the operating point (Figure 2, Table 3).
- `mc_full.csv`, `mc_KO11.csv`, `mc_LY01.csv` — Monte Carlo samples per chassis (Figure 4).
- `adversarial_sweeps.csv` / `_summary.csv` — AS-1/AS-2/AS-3 outputs (Figure 6, §3.9).
- `controller_benchmark_sweeps.csv` / `_summary.csv` — single-enzyme vs redox-partitioning outputs (Figure 7, §3.10).

The **full** genome-scale dFBA pipeline (COBRApy + iJO1366/iML1515) that produced
the primary `data/` files is provided as a separate Jupyter notebook in the main
project repository; these scripts reproduce the figures from its saved outputs.

## Provenance statement (for the AI-use declaration)

All plots in this folder are produced by deterministic matplotlib code operating on
CSV simulation outputs. 
