"""
Graphical Abstract for Biochemical Engineering Journal.
BEJ format: 1328 x 531 px (w x h), ratio ~2.5:1, readable at 13 x 5 cm.
Horizontal 4-panel flow: chemistry -> controller -> dose-response -> message.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'savefig.facecolor': 'white', 'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
})

COL_WT = '#1f77b4'; COL_KO = '#2ca02c'; COL_LY = '#d62728'
COL_BOX = '#34495e'; COL_ACCENT = '#e67e22'

# 1328 x 531 px at 300 dpi -> 4.43 x 1.77 inches
fig = plt.figure(figsize=(13.28/2.54*1.0, 5.31/2.54*1.0), dpi=300)
# Use cm: width 13.28cm, height 5.31cm
fig.set_size_inches(13.28/2.54, 5.31/2.54)

gs = gridspec.GridSpec(1, 4, figure=fig, wspace=0.45,
                       left=0.04, right=0.985, top=0.82, bottom=0.16)

# ---------- Panel A: PAW chemistry -> stress ----------
axA = fig.add_subplot(gs[0, 0])
t = np.linspace(0, 10, 100)
species = [('H\u2082O\u2082', 5.0, 0.277, '#1f77b4'),
           ('NO\u2082\u207b', 2.5, 0.150, '#ff7f0e'),
           ('ONOO\u207b', 0.5, 1.5, '#d62728')]
for name, c0, k, col in species:
    axA.plot(t, c0*np.exp(-k*t), color=col, lw=1.3, label=name)
axA.set_title('PAW chemistry', fontsize=7, fontweight='bold', pad=3)
axA.set_xlabel('time (h)', fontsize=6)
axA.set_ylabel('mM', fontsize=6)
axA.legend(fontsize=4.5, loc='upper right', framealpha=0.9, handlelength=1)
axA.tick_params(labelsize=5)
axA.grid(alpha=0.25, lw=0.4)

# ---------- Panel B: bell-shaped controller ----------
axB = fig.add_subplot(gs[0, 1])
S = np.linspace(0, 20, 200)
r = 0.5 + (10-0.5)*(S**2/(1.8**2+S**2))*(1/(1+(S/13)**4))
axB.plot(S, r, color=COL_ACCENT, lw=1.8)
axB.fill_between(S, r, alpha=0.15, color=COL_ACCENT)
peak_i = np.argmax(r)
axB.plot(S[peak_i], r[peak_i], '*', color='#c0392b', ms=9)
axB.set_title('Redox controller', fontsize=7, fontweight='bold', pad=3)
axB.set_xlabel('stress S', fontsize=6)
axB.set_ylabel('EtOH/formate', fontsize=6)
axB.tick_params(labelsize=5)
axB.grid(alpha=0.25, lw=0.4)
axB.text(0.5, 0.92, 'bell-shaped', transform=axB.transAxes, fontsize=5,
         ha='left', style='italic', color=COL_ACCENT)

# ---------- Panel C: hormetic dose-response ----------
axC = fig.add_subplot(gs[0, 2])
a = np.linspace(0, 1.5, 50)
def dose(base, peak_gain):
    surv = np.exp(-np.cumsum(np.where(a>0.6, (a-0.6)*3, 0))*0.02)
    return base*(1 + peak_gain*np.exp(-((a-0.6)/0.35)**2))*surv
axC.plot(a, dose(17, 0.18), color=COL_WT, lw=1.5, label='WT')
axC.plot(a, dose(18, 0.18), color=COL_KO, lw=1.3, ls='--', label='KO11*')
axC.plot(a, 25.2*np.exp(-np.cumsum(np.where(a>0.6,(a-0.6)*3,0))*0.02),
         color=COL_LY, lw=1.3, ls=':', label='LY01*')
axC.axvspan(0.4, 0.7, alpha=0.12, color='gold')
axC.set_title('Hormetic window', fontsize=7, fontweight='bold', pad=3)
axC.set_xlabel(r'$\alpha_{PAW}$', fontsize=6)
axC.set_ylabel('EtOH (mM)', fontsize=6)
axC.legend(fontsize=4.5, loc='lower left', framealpha=0.9, handlelength=1.3)
axC.tick_params(labelsize=5)
axC.grid(alpha=0.25, lw=0.4)
axC.annotate('+18\u201328%', xy=(0.6, 20), xytext=(0.85, 22),
             fontsize=5.5, color='#c0392b', fontweight='bold')

# ---------- Panel D: central message ----------
axD = fig.add_subplot(gs[0, 3])
axD.axis('off')
box = FancyBboxPatch((0.05, 0.30), 0.9, 0.45, boxstyle="round,pad=0.03",
                     linewidth=1.3, edgecolor=COL_BOX, facecolor='#ecf0f1',
                     transform=axD.transAxes)
axD.add_patch(box)
axD.text(0.5, 0.66, 'PAW', ha='center', fontsize=10, fontweight='bold',
         color=COL_ACCENT, transform=axD.transAxes)
axD.text(0.5, 0.52, 'may phenocopy', ha='center', fontsize=5.5,
         style='italic', transform=axD.transAxes)
axD.text(0.5, 0.40, r'$\Delta pta$–$ackA$', ha='center', fontsize=9,
         fontweight='bold', color=COL_BOX, transform=axD.transAxes)
axD.text(0.5, 0.16, 'non-genetic\nengineering stimulus', ha='center',
         fontsize=5, color='#555', transform=axD.transAxes)
axD.set_title('Testable hypothesis', fontsize=7, fontweight='bold', pad=3)

# Flow arrows between panels (in figure coordinates)
for x in [0.265, 0.505, 0.745]:
    fig.patches.append(FancyArrowPatch((x, 0.49), (x+0.018, 0.49),
        transform=fig.transFigure, arrowstyle='-|>', mutation_scale=8,
        color='#888', lw=1.0))

# Overall title strip
fig.text(0.5, 0.95, 'Cross-model hybrid dFBA framework for PAW-stimulated '
         '$\\it{E.\\ coli}$ ethanol fermentation',
         ha='center', fontsize=7.5, fontweight='bold')

fig.savefig('/mnt/user-data/outputs/Graphical_Abstract_BEJ.png', dpi=300,
            facecolor='white', bbox_inches='tight', pad_inches=0.05)
plt.close()

# Verify dimensions
from PIL import Image
img = Image.open('/mnt/user-data/outputs/Graphical_Abstract_BEJ.png')
print(f"Saved Graphical_Abstract_BEJ.png: {img.size[0]} x {img.size[1]} px")
print(f"Ratio: {img.size[0]/img.size[1]:.2f} (target ~2.5)")
