# A cross-model hybrid dynamic flux balance framework for plasma-activated water exposure in *Escherichia coli* ethanol fermentation

Companion code and data for the manuscript:

> *A cross-model hybrid dynamic flux balance framework as a hypothesis-generating tool for plasma-activated water exposure in Escherichia coli ethanol fermentation.*

This repository contains everything required to reproduce the computational results, **except the two genome-scale metabolic reconstructions** (`iJO1366`, `iML1515`), which are loaded automatically from the BiGG / COBRApy model registry at run time (see *Genome-scale models* below).

---

## What this framework does

A hybrid **dynamic flux balance analysis (dFBA)** pipeline couples three layers:

1. **An analytical PAW-chemistry layer** — first-order decay of the four dominant reactive species (H₂O₂, NO₂⁻, NO₃⁻, and peroxynitrite ONOO⁻), aggregated into a single dimensionless stress index scaled by a dose parameter `α_PAW`.
2. **A network-level redox-partitioning controller** — a phenomenological constraint on the ethanol/by-product exchange-flux ratio that represents stress-driven NADH/NAD⁺ rebalancing, plus a peroxynitrite-dependent inhibition of the ethanol-forming enzyme (AdhE / `ALCD2x`).
3. **Two genome-scale reconstructions** (`iJO1366`, `iML1515`) across three chassis proxies (WT, KO11, LY160), re-solved at each time step over a 10 h anaerobic batch.

The framework predicts a **hormetic ethanol-enhancement window** at `α_PAW ≈ 0.60` in WT/KO11 and no response in the acetate-disabled LY160 proxy, and includes Monte-Carlo, Morris, and adversarial-sensitivity diagnostics that identify the load-bearing assumptions.

---

## Repository structure

```
PAW_Ethanol_dFBA/
├── README.md                     ← this file
├── notebook/
│   └── PAW_Ethanol_dFBA.ipynb    ← main reproducible notebook (Colab-ready)
├── code/
│   ├── dfba_framework.py         ← standalone dFBA framework (models, stress, controller, kd)
│   └── adhe_inhibition.py        ← AdhE peroxynitrite-inhibition module (adversarial scenario AS-4)
├── data/                         ← all precomputed result tables (see DATA_DICTIONARY.md)
│   ├── sweep_v4_combined.csv
│   ├── mc_full.csv / mc_KO11.csv / mc_LY160.csv
│   ├── adversarial_sweeps.csv / adversarial_summary.csv
│   ├── controller_benchmark_sweeps.csv / controller_benchmark_summary.csv
│   ├── adhe_AS4_summary.csv
│   ├── adhe_iJO1366.json / adhe_iML1515.json
│   ├── adhe_sustained_iJO1366.json / adhe_sustained_iML1515.json
│   └── NEW_NUMBERS.json

    
```

---

## Requirements

- Python ≥ 3.9
- [`cobra`](https://opencobra.github.io/cobrapy/) == 0.29.1
- `numpy`, `pandas`, `scipy`, `matplotlib`, `openpyxl`, `tqdm`

Install:

```bash
pip install cobra==0.29.1 numpy pandas scipy matplotlib openpyxl tqdm
```

---

## Genome-scale models (NOT included in this package)

The two reconstructions are **deliberately excluded** from this package. They are loaded directly from the COBRApy / BiGG registry by name:

```python
from cobra.io import load_model
m_J = load_model("iJO1366")   # Orth et al. 2011
m_M = load_model("iML1515")   # Monk et al. 2017
```

`load_model()` downloads and caches them automatically on first call — no local file is needed. If you prefer local SBML copies, download them from the BiGG Models database (http://bigg.ucsd.edu/models/iJO1366 and http://bigg.ucsd.edu/models/iML1515) and load with `cobra.io.read_sbml_model(path)`.

Both models are used with their **native default constraints**, including their respective non-growth-associated ATP-maintenance fluxes (ATPM = 3.15 for iJO1366, 6.86 for iML1515 mmol gDW⁻¹ h⁻¹). Anaerobic batch conditions are set by fixing the O₂ exchange lower bound to 0 and the glucose exchange lower bound to −10 mmol gDW⁻¹ h⁻¹.

---

## How to reproduce

**Option A — notebook (recommended).** Open `notebook/PAW_Ethanol_dFBA.ipynb` in Jupyter or Google Colab and run all cells. The first cell installs dependencies; the models are fetched automatically. The notebook regenerates every sweep, the Monte-Carlo ensemble, the adversarial scenarios (AS-1…AS-4), the controller benchmark, and all manuscript figures. (The notebook was developed in Colab; if running locally, remove or skip the `from google.colab import drive` cell and adjust output paths.)

**Option B — standalone modules.** `code/dfba_framework.py` exposes the core API:

```python
from dfba_framework import load_model_cfg, make_configured, run_dfba_cfg, P

m = make_configured(load_model_cfg("iJO1366"), chassis="WT")
result = run_dfba_cfg(m, alpha=0.60, scenario="nominal")
```

`code/adhe_inhibition.py` adds the peroxynitrite-dependent AdhE cap used in adversarial scenario AS-4 (transient vs. sustained peroxynitrite decay).

---

## Key results (from `data/NEW_NUMBERS.json`)

- Baseline (α = 0) ethanol titer: WT/iJO1366 18.0 mM, WT/iML1515 11.0 mM.
- Hormetic peak at α* = 0.60: 21.05 mM (iJO1366, +16.9 %) and 14.0 mM (iML1515, +26.9 %).
- LY160 proxy: unresponsive (internal negative control).
- Monte Carlo (1200 samples): productive window stable, CV ≈ 6 %, 100 % stress-threshold compliant.
- Adversarial analysis: the hormesis is carried by the **redox-partitioning constraint** (not single-enzyme forcing, not the PFL-inhibition arm); peroxynitrite weighting is the critical structural assumption.
- AS-4 (AdhE inhibition): the window is robust to even severe (90 %) AdhE inhibition under **transient** peroxynitrite, but collapses under **sustained** peroxynitrite — a falsifiable prediction.

---

## License & citation

Please cite the manuscript above when using this code or data. The genome-scale models retain their original BiGG licenses.
