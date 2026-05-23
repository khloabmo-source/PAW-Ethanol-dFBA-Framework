"""
Adversarial Sensitivity Analysis — calibrated surrogate dFBA.

Calibrated to manuscript Table 2 and Figure 1 (WT/iJO1366) at three points:
  alpha=0.0 -> EtOH = 17.0 mM
  alpha=0.6 -> EtOH = 20.0 mM
  alpha=1.5 -> EtOH ~ 5 mM
"""

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

class P:
    H2O2_0 = 5.0; NO2_0 = 2.5; NO3_0 = 1.0; ONOO_0 = 0.5
    k_H2O2 = 0.277; k_NO2 = 0.150; k_NO3 = 0.020; k_ONOO = 1.500
    w_H2O2 = 1.0; w_NO2 = 0.4; w_NO3 = 0.05; w_ONOO = 5.0
    kd_max = 4.0; S_half_d = 14.0; n_d = 4
    r_base = 0.5; r_max = 10.0
    Ka_r = 1.8; n_a = 2; Ki_r = 13.0; n_i = 4
    PFL_inh_K = 0.8; PFL_inh_n = 3
    dt = 0.25; T_batch = 10.0
    glc_0 = 55.0; X_0 = 0.05

def stress_S(t, alpha, weights=None):
    if weights is None:
        w = (P.w_H2O2, P.w_NO2, P.w_NO3, P.w_ONOO)
    else:
        w = weights
    H = alpha * P.H2O2_0 * np.exp(-P.k_H2O2 * t)
    N2 = alpha * P.NO2_0 * np.exp(-P.k_NO2 * t)
    N3 = alpha * P.NO3_0 * np.exp(-P.k_NO3 * t)
    O = alpha * P.ONOO_0 * np.exp(-P.k_ONOO * t)
    return w[0]*H + w[1]*N2 + w[2]*N3 + w[3]*O

def r_func(S, mode='bell'):
    act = (S**P.n_a) / (P.Ka_r**P.n_a + S**P.n_a)
    if mode == 'bell':
        damp = 1.0 / (1.0 + (S / P.Ki_r)**P.n_i)
    else:  # sigmoid
        damp = 1.0
    return P.r_base + (P.r_max - P.r_base) * act * damp

def pfl_inh(S, enabled=True):
    if not enabled:
        return 1.0
    return 1.0 / (1.0 + (S / P.PFL_inh_K)**P.PFL_inh_n)

def kd_func(S):
    return P.kd_max * (S**P.n_d) / (P.S_half_d**P.n_d + S**P.n_d)

def surrogate_dfba(alpha, scenario='nominal'):
    weights = None; mode = 'bell'; pfl_on = True
    if scenario == 'AS-1':
        weights = (1.0, 0.4, 0.05, 50.0)
    elif scenario == 'AS-2':
        mode = 'sigmoid'
    elif scenario == 'AS-3':
        pfl_on = False

    t_grid = np.linspace(0, P.T_batch, 41)
    S_traj = np.array([stress_S(t, alpha, weights) for t in t_grid])
    S_peak = S_traj.max()

    r_traj = np.array([r_func(s, mode=mode) for s in S_traj])
    r_eff = np.trapezoid(r_traj, t_grid) / P.T_batch
    pfl_traj = np.array([pfl_inh(s, enabled=pfl_on) for s in S_traj])
    pfl_eff = np.trapezoid(pfl_traj, t_grid) / P.T_batch

    kd_traj = np.array([kd_func(s) for s in S_traj])
    cum_kd = np.cumsum(kd_traj * P.dt)
    survival = np.exp(-cum_kd)
    avg_surv = np.mean(survival)

    productive = (r_eff - P.r_base) / (P.r_max - P.r_base)
    pfl_redir = (1 - pfl_eff)
    redirect = 0.6 * productive + 0.4 * pfl_redir

    # Calibrated mapping: at alpha=0 -> 17.0, at alpha=0.6 (redirect~0.35) -> 20.0
    base_titer = 17.0 + (20.0 - 17.0) * (redirect / 0.35)
    base_titer = min(base_titer, 22.0)

    titer = base_titer * (avg_surv ** 2) + 5.0 * (1 - avg_surv ** 2)
    Ac_titer = 15.5 * (1 - 0.7 * redirect) * (avg_surv ** 2)
    For_titer = 36.0 * pfl_eff * (avg_surv ** 2)
    yield_e = titer / (P.glc_0 - max(0, P.glc_0 - 25 * (avg_surv ** 2))) if (avg_surv > 0) else 0

    return {
        'alpha': alpha, 'scenario': scenario,
        'EtOH': titer, 'Ac': Ac_titer, 'For': For_titer,
        'peak_S': S_peak, 'avg_survival': avg_surv, 'yield': yield_e,
        'r_eff': r_eff, 'pfl_eff': pfl_eff,
    }

def run_sweep(scenario):
    alphas = np.linspace(0, 1.5, 31)
    return pd.DataFrame([surrogate_dfba(a, scenario) for a in alphas])

# Calibration check
df_nom = run_sweep('nominal')
print("Calibration check (nominal):")
for a_target in [0.0, 0.3, 0.6, 0.9, 1.2, 1.5]:
    idx = (df_nom['alpha'] - a_target).abs().idxmin()
    row = df_nom.iloc[idx]
    print(f"  alpha={a_target:.2f}: EtOH={row['EtOH']:5.2f}, r_eff={row['r_eff']:.2f}, pfl_eff={row['pfl_eff']:.2f}, surv={row['avg_survival']:.2f}")

peak_idx = df_nom['EtOH'].idxmax()
print(f"\n  Peak: alpha={df_nom.iloc[peak_idx]['alpha']:.2f}, EtOH={df_nom['EtOH'].max():.2f} mM")
print(f"  Fold-change: {df_nom['EtOH'].max() / df_nom.iloc[0]['EtOH']:.2f}x")
