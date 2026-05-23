"""
Figure 7: Controller benchmark — single-enzyme vs redox-partitioning.

Three panels demonstrating the methodological advantages of Controller B:
  (a) Dose-response: enzyme-choice fragility of Controller A
  (b) Cross-GSMM concordance comparison
  (c) By-product suppression realism
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

df = pd.read_csv('/home/claude/controller_benchmark_sweeps.csv')

COL_B = '#1f77b4'       # redox-partitioning (ours)
COL_A_ALCD = '#d62728'  # single-enzyme ALCD2x
COL_A_ADHE = '#ff7f0e'  # single-enzyme ADHEr (fragility)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), facecolor='white')

# ============================================================
# Panel (a): Dose-response — enzyme-choice fragility
# ============================================================
ax = axes[0]

B = df[df['config']=='B_iJO'].sort_values('alpha')
A_alcd = df[df['config']=='A_ALCD_iJO'].sort_values('alpha')
A_adhe = df[df['config']=='A_ADHE_iJO'].sort_values('alpha')

ax.plot(B['alpha'], B['EtOH'], '-o', color=COL_B, lw=2.5, ms=4,
        markeredgecolor='white', markeredgewidth=0.4,
        label='B: redox-partitioning (this work)')
ax.plot(A_alcd['alpha'], A_alcd['EtOH'], '--s', color=COL_A_ALCD, lw=2, ms=4,
        markeredgecolor='white', markeredgewidth=0.4,
        label='A: single-enzyme (ALCD2x target)')
ax.plot(A_adhe['alpha'], A_adhe['EtOH'], ':^', color=COL_A_ADHE, lw=2, ms=4,
        markeredgecolor='white', markeredgewidth=0.4,
        label='A: single-enzyme (ADHEr target)')

# Highlight the fragility gap at the peak
ax.annotate('', xy=(0.6, A_alcd['EtOH'].max()), xytext=(0.6, A_adhe['EtOH'].max()),
            arrowprops=dict(arrowstyle='<->', color='#666', lw=1.2))
ax.text(0.68, (A_alcd['EtOH'].max() + A_adhe['EtOH'].max())/2,
        'enzyme-choice\nfragility\n(Δfold = 0.19)', fontsize=8,
        color='#444', va='center', fontstyle='italic')

ax.set_xlabel(r'PAW dose ($\alpha_{PAW}$)', fontsize=10.5)
ax.set_ylabel('Final ethanol titer (mM)', fontsize=10.5)
ax.set_title('(a) Enzyme-choice fragility', fontsize=11, fontweight='bold', loc='left', pad=8)
ax.legend(fontsize=8, loc='lower left', framealpha=0.95, edgecolor='#666')
ax.grid(True, alpha=0.3, linestyle='--', lw=0.5)
ax.set_xlim(-0.05, 1.55)
ax.set_ylim(0, 25)

# ============================================================
# Panel (b): Cross-GSMM concordance
# ============================================================
ax2 = axes[1]

# Fold-changes for each controller on each GSMM
def fold(config):
    sub = df[df['config']==config]
    return sub['EtOH'].max() / sub[sub['alpha']==0].iloc[0]['EtOH']

folds_data = {
    'Controller A\n(single-enzyme)': [fold('A_ALCD_iJO'), fold('A_ALCD_iML')],
    'Controller B\n(redox-partition)': [fold('B_iJO'), fold('B_iML')],
}

x_pos = np.arange(2)
width = 0.35

iJO_vals = [folds_data['Controller A\n(single-enzyme)'][0], folds_data['Controller B\n(redox-partition)'][0]]
iML_vals = [folds_data['Controller A\n(single-enzyme)'][1], folds_data['Controller B\n(redox-partition)'][1]]

bars_J = ax2.bar(x_pos - width/2, iJO_vals, width, color='#1f77b4',
                 edgecolor='#333', linewidth=0.8, label='iJO1366')
bars_M = ax2.bar(x_pos + width/2, iML_vals, width, color='#d62728',
                 edgecolor='#333', linewidth=0.8, label='iML1515')

# Annotate spread
for i, (vj, vm) in enumerate(zip(iJO_vals, iML_vals)):
    spread = abs(vj - vm)
    ax2.text(i - width/2, vj + 0.008, f'{vj:.2f}', ha='center', fontsize=8.5, color='#1f77b4', fontweight='bold')
    ax2.text(i + width/2, vm + 0.008, f'{vm:.2f}', ha='center', fontsize=8.5, color='#d62728', fontweight='bold')
    ax2.text(i, 0.93, f'spread\n{spread:.3f}', ha='center', fontsize=8.5,
             color='#1a5e1a' if spread < 0.05 else '#8b0000', fontweight='bold')

ax2.axhline(1.0, color='#666', ls=':', lw=1)
ax2.set_xticks(x_pos)
ax2.set_xticklabels(['Controller A\n(single-enzyme)', 'Controller B\n(redox-partition)'], fontsize=9)
ax2.set_ylabel('Hormetic fold-change', fontsize=10.5)
ax2.set_ylim(0.90, 1.40)
ax2.set_title('(b) Cross-GSMM concordance', fontsize=11, fontweight='bold', loc='left', pad=8)
ax2.legend(fontsize=9, loc='upper right', framealpha=0.95, edgecolor='#666')
ax2.grid(True, alpha=0.3, linestyle='--', lw=0.5, axis='y')

# ============================================================
# Panel (c): By-product suppression realism
# ============================================================
ax3 = axes[2]

def byproduct_at(config, species, a=0.6):
    sub = df[df['config']==config]
    idx = (sub['alpha'] - a).abs().idxmin()
    return sub.loc[idx, species]

ac_baseline = 15.5
for_baseline = 36.0

ac_supp_A = (ac_baseline - byproduct_at('A_ALCD_iJO', 'Ac')) / ac_baseline * 100
for_supp_A = (for_baseline - byproduct_at('A_ALCD_iJO', 'For')) / for_baseline * 100
ac_supp_B = (ac_baseline - byproduct_at('B_iJO', 'Ac')) / ac_baseline * 100
for_supp_B = (for_baseline - byproduct_at('B_iJO', 'For')) / for_baseline * 100

categories = ['Acetate', 'Formate']
x_pos = np.arange(2)
width = 0.35

A_supp = [ac_supp_A, for_supp_A]
B_supp = [ac_supp_B, for_supp_B]

bars_A = ax3.bar(x_pos - width/2, A_supp, width, color=COL_A_ALCD,
                 edgecolor='#333', linewidth=0.8, label='A: single-enzyme', alpha=0.85)
bars_B = ax3.bar(x_pos + width/2, B_supp, width, color=COL_B,
                 edgecolor='#333', linewidth=0.8, label='B: redox-partition', alpha=0.85)

for i, (va, vb) in enumerate(zip(A_supp, B_supp)):
    ax3.text(i - width/2, va + 1.5, f'{va:.0f}%', ha='center', fontsize=9, color=COL_A_ALCD, fontweight='bold')
    ax3.text(i + width/2, vb + 1.5, f'{vb:.0f}%', ha='center', fontsize=9, color=COL_B, fontweight='bold')

ax3.set_xticks(x_pos)
ax3.set_xticklabels(categories, fontsize=10)
ax3.set_ylabel('By-product suppression at α=0.6 (%)', fontsize=10.5)
ax3.set_ylim(0, 80)
ax3.set_title('(c) By-product suppression realism', fontsize=11, fontweight='bold', loc='left', pad=8)
ax3.legend(fontsize=9, loc='upper left', framealpha=0.95, edgecolor='#666')
ax3.grid(True, alpha=0.3, linestyle='--', lw=0.5, axis='y')

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/Figure_7_controller_benchmark.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved: Figure_7_controller_benchmark.png")
