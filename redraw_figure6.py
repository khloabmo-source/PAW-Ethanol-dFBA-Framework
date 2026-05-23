"""
Re-draw Figure 6 (adversarial sensitivity) with NO text overlap.

Fix: the 'Nominal operating' annotation in panel (b) overlapped the legend.
Solution:
  - Move the nominal-optimum reference label to a clear position (below the dashed line, left side)
  - Reposition the legend to upper-left where there is empty space
  - Increase right-axis headroom so bars and labels do not collide
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'savefig.facecolor': 'white', 'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
})

df = pd.read_csv('/home/claude/adversarial_sweeps.csv')

COL_SC = {
    'nominal': '#1f77b4',
    'AS-1': '#d62728',
    'AS-2': '#9467bd',
    'AS-3': '#ff7f0e',
}
LABEL = {
    'nominal': 'Nominal (Figure 1 baseline)',
    'AS-1': r'AS-1: $w_{ONOO}$ × 10',
    'AS-2': 'AS-2: Sigmoid controller (no damping)',
    'AS-3': 'AS-3: PFL inhibition removed',
}

fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor='white')

# ============================================================
# Panel (a): dose-response — unchanged, it was fine
# ============================================================
ax = axes[0]
for sc in ['nominal', 'AS-1', 'AS-2', 'AS-3']:
    sub = df[df['scenario']==sc].sort_values('alpha')
    ls = '-' if sc == 'nominal' else '--'
    lw = 2.5 if sc == 'nominal' else 1.8
    ax.plot(sub['alpha'], sub['EtOH'], ls, color=COL_SC[sc], lw=lw,
            label=LABEL[sc], marker='o' if sc=='nominal' else 's',
            markersize=4, markeredgecolor='white', markeredgewidth=0.4)

ax.axvspan(0.4, 0.7, alpha=0.10, color='#ffd700', zorder=0)
# Place 'Nominal operating window' text in clear upper area, left of the curves' descent
ax.text(0.55, 23.3, 'Nominal\noperating\nwindow', ha='center', va='top', fontsize=8,
        color='#806000', fontweight='bold')

ax.set_xlabel(r'Dimensionless PAW dose ($\alpha_{PAW}$)', fontsize=11)
ax.set_ylabel('Final ethanol titer (mM)', fontsize=11)
ax.set_title('(a) Adversarial dose-response (WT/iJO1366)',
             fontsize=11, fontweight='bold', loc='left', pad=8)
ax.legend(fontsize=8.5, loc='lower left', framealpha=0.95, edgecolor='#666')
ax.grid(True, alpha=0.3, linestyle='--', lw=0.5)
ax.set_xlim(-0.05, 1.55)
ax.set_ylim(0, 25)

# ============================================================
# Panel (b): REDESIGNED to eliminate overlap
# ============================================================
ax2 = axes[1]
scenarios = ['nominal', 'AS-1', 'AS-2', 'AS-3']
x_labels = ['Nominal', 'AS-1', 'AS-2', 'AS-3']
x_pos = np.arange(len(scenarios))

folds = []
peak_alphas = []
for sc in scenarios:
    sub = df[df['scenario']==sc]
    fold = sub['EtOH'].max() / sub[sub['alpha']==0].iloc[0]['EtOH']
    folds.append(fold)
    peak_idx_local = sub['EtOH'].values.argmax()
    peak_alphas.append(sub['alpha'].values[peak_idx_local])

# LEFT axis: fold-change bars (filled)
bars1 = ax2.bar(x_pos - 0.2, folds, 0.38,
                color=[COL_SC[s] for s in scenarios],
                edgecolor='#333', linewidth=0.8, zorder=3)
ax2.set_ylabel('Hormetic fold-change', fontsize=11, color='#1a1a1a')
ax2.set_ylim(0.90, 1.30)   # extra headroom so labels don't hit the top
ax2.axhline(1.0, color='#666', ls=':', lw=1, zorder=1)

# RIGHT axis: peak position bars (hatched, hollow)
ax2b = ax2.twinx()
bars2 = ax2b.bar(x_pos + 0.2, peak_alphas, 0.38,
                 color='none', edgecolor=[COL_SC[s] for s in scenarios],
                 linewidth=1.8, hatch='///', zorder=3)
ax2b.set_ylabel(r'Peak position $\alpha_{PAW}^*$', fontsize=11, color='#1a1a1a')
ax2b.set_ylim(0, 0.85)     # extra headroom

# Nominal optimum reference line on the RIGHT axis (peak alpha = 0.60)
ax2b.axhline(0.60, color='#d4a017', ls='--', lw=1.2, alpha=0.8, zorder=2)
# Label it on the FAR LEFT, BELOW the line — clear of everything
ax2b.text(-0.45, 0.605, r'nominal $\alpha^*$=0.60', fontsize=7.5,
          color='#806000', ha='left', va='bottom', style='italic')

# Value labels on bars (fold-change) — above each filled bar
for i, fold in enumerate(folds):
    ax2.text(i - 0.2, fold + 0.008, f'{fold:.2f}', ha='center', fontsize=8.5,
             fontweight='bold', color=COL_SC[scenarios[i]], zorder=5)
# Value labels (peak alpha) — above each hatched bar
for i, pa in enumerate(peak_alphas):
    ax2b.text(i + 0.2, pa + 0.02, f'{pa:.2f}', ha='center', fontsize=8.5,
              fontweight='bold', color=COL_SC[scenarios[i]], zorder=5)

ax2.set_xticks(x_pos)
ax2.set_xticklabels(x_labels)
ax2.set_xlim(-0.6, 3.6)
ax2.set_title('(b) Impact on hormetic fold-change and peak position',
              fontsize=11, fontweight='bold', loc='left', pad=8)
ax2.grid(True, alpha=0.3, linestyle='--', lw=0.5, axis='y', zorder=0)

# LEGEND: place at top-center, horizontal, clear of all bars and the nominal label
legend_handles = [
    Patch(facecolor='#888', edgecolor='#333', label='Fold-change (left axis)'),
    Patch(facecolor='none', edgecolor='#333', hatch='///', label=r'Peak $\alpha^*$ (right axis)'),
]
ax2.legend(handles=legend_handles, fontsize=8.5, loc='upper center',
           ncol=2, framealpha=0.95, edgecolor='#666',
           bbox_to_anchor=(0.5, 1.0), columnspacing=1.2)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/Figure_6_adversarial_sensitivity.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved: Figure_6_adversarial_sensitivity.png (overlap fixed)")
