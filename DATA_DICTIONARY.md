# Data dictionary

All result tables in `data/`. Ethanol/by-product titers are in mM; `α_PAW` (alpha) is the dimensionless PAW dose (α = 0.60 corresponds to 3.0 mM H₂O₂, 1.5 mM NO₂⁻, 0.6 mM NO₃⁻, 0.3 mM ONOO⁻).

## sweep_v4_combined.csv
Full dose-response sweep across chassis × controller-variant × reconstruction.
| column | meaning |
|---|---|
| `chassis` | WT, KO11, or LY160 |
| `variant` | controller variant label |
| `gsmm` | `iJO1366` or `iML1515` |
| `alpha` | PAW dose α_PAW |
| `X` | final biomass (gDW L⁻¹) |
| `glc` | residual glucose (mM) |
| `EtOH`, `For`, `Ac` | ethanol, formate, acetate titers (mM) |
| `yield` | ethanol yield |
| `peak_S` | peak aggregate stress index |

## mc_full.csv / mc_KO11.csv / mc_LY160.csv
Monte-Carlo ensemble (1200 samples) for WT, KO11, LY160 respectively. Columns include the sampled parameters (`r_max`, `kd_max`, `beta_ATPM`, `Ka_r`, `Ki_r`, `PFL_inh_K`, `gamma_NADH`) plus outcomes (`EtOH`, `peak_S`, `X_final`, `yield`) at the sampled `alpha`.

## adversarial_sweeps.csv / adversarial_summary.csv
Adversarial scenarios AS-1 (×10 peroxynitrite weight), AS-2 (no high-stress damping), AS-3 (PFL inhibition removed).
`adversarial_summary.csv` columns: `scenario`, `baseline_EtOH`, `peak_EtOH`, `peak_alpha`, `fold_change`.

## controller_benchmark_sweeps.csv / controller_benchmark_summary.csv
Controller A (single-enzyme activation, forcing `ALCD2x`/`ACALD`) vs. Controller B (redox-partitioning, this work) under matched conditions.
`controller_benchmark_summary.csv` reports fold-change and titers per reconstruction; Controller A collapses ethanol (18.0 → 7.85 mM at α = 0.60), Controller B redistributes carbon (18.0 → 21.05 mM).

## adhe_AS4_summary.csv
Adversarial scenario AS-4 — peroxynitrite-dependent AdhE inhibition.
| column | meaning |
|---|---|
| `GSMM` | reconstruction |
| `AdhE_inhibition` | none / moderate / strong / severe |
| `I_AdhE_at_op` | residual ethanol-enzyme capacity at the operating point (1 = none) |
| `fold_fastdecay`, `astar_fastdecay` | hormetic fold-change and optimal α under nominal (transient) peroxynitrite decay |
| `fold_sustained`, `astar_sustained` | the same under ten-fold slower (sustained) peroxynitrite decay |

## adhe_iJO1366.json / adhe_iML1515.json (and `*_sustained_*`)
Full AdhE-inhibition dose-response curves. Each file: `{"alphas": [...], "rows": [{inhibition, K_AdhE, GSMM, baseline, peak, fold, alpha_star, I_op, curve:[EtOH per alpha]}, ...]}`. The `_sustained_` files use the ten-fold slower peroxynitrite decay.

## NEW_NUMBERS.json
Consolidated key numbers reported in the manuscript (baseline titers, fold-changes, peak α, cross-GSMM ratios, Monte-Carlo statistics, adversarial and controller-benchmark outcomes).
