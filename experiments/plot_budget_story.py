"""Create a source-mask -> admissible boundary -> count interval figure.

Optional figure dependencies: NumPy and Matplotlib. Numerical claims are read
from executed exact records; floats are used only to draw the figure.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[1]
charts = json.loads((ROOT / 'data/charts.json').read_text())
ch = next(x for x in charts if x['id'] == '4037abe9fea38e')
with gzip.open(ROOT / 'results/budgets/native_frontiers.json.gz', 'rt') as f:
    artifact = json.load(f)
model = next(x for x in artifact['models'] if x['id'] == ch['id'] and x['setting'] == 'baseline')
gene = next(x for x in model['genes'] if x['gene'] == 'PECAM1')
cell = next(x for x in gene['cells'] if x['cell'] == 153124)
record = next(x for x in cell['frontier'] if x['count'] == 1)
state = gene['enumeration']['states'][record['state_index']]
threshold = float(F(record['minimum_deformation']))
px = float(F(ch['pixel_size_um']))
mask = np.load(ROOT / 'data/native/roi_01/masks.npz')['cell_mask']
crop = mask[ch['row']:ch['row'] + ch['height'], ch['left']:ch['right']]
colored = np.zeros_like(crop, dtype=int)
for i, label in enumerate(ch['cell_labels']): colored[crop == label] = i + 1
source_t, source_x = zip(*[(float(F(t)), float(F(x))) for t, x in ch['baselines'][0]])
segments = [[float(F(x)) for x in row] for row in state['boundaries'][0]]
tt = [segments[0][0]] + [row[1] for row in segments]
xx = [segments[0][2]] + [row[3] for row in segments]
plt.rcParams.update({'font.size': 10.2, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.8), layout='constrained')
for ax in axes[0]:
    ax.imshow(colored, origin='lower', interpolation='nearest', aspect='equal',
              extent=[0, crop.shape[1] * px, -px / 2, (crop.shape[0] - .5) * px],
              cmap=ListedColormap(['#eeeeee', '#d8e8f4', '#f6e1c3']), vmin=0, vmax=2)
    for pt in ch['points']:
        special = pt['gene'] == 'PECAM1'
        ax.scatter(float(F(pt['x'])), float(F(pt['t'])), s=65 if special else 8,
                   marker='*' if special else '.', color='#a4262c' if special else '#555555', zorder=5)
    ax.plot(source_x, source_t, color='#333333', lw=1.3)
    ax.set(xlabel='Transverse coordinate (µm)', ylabel='Longitudinal (µm)', ylim=(0, source_t[-1]))
axes[0, 0].set_title('A  Source masks and RNA', loc='left', fontweight='bold')
axes[0, 1].set_title('B  A verified alternative', loc='left', fontweight='bold')
axes[0, 1].plot(xx, tt, color='#1665a8', lw=1.8)
t = np.unique(np.r_[tt, source_t])
axes[0, 1].fill_betweenx(t, np.interp(t, source_t, source_x), np.interp(t, tt, xx), color='#1665a8', alpha=.2)
axes[0, 0].text(.03, .96, 'Cell 153124', transform=axes[0, 0].transAxes, va='top', fontsize=8.5)
axes[0, 0].text(.97, .96, '153133', transform=axes[0, 0].transAxes, va='top', ha='right', fontsize=8.5)
ax = axes[1, 0]
ax.set_title('C  Exact count interval', loc='left', fontweight='bold')
ax.step([0, threshold, 4], [0, 1, 1], where='post', color='#1665a8', lw=2, label='Upper endpoint')
ax.plot([0, 4], [0, 0], color='#333333', lw=1.3, label='Lower endpoint')
ax.fill_between([threshold, 4], [0, 0], [1, 1], color='#1665a8', alpha=.16)
ax.axvline(threshold, ls=':', lw=1, color='#777777')
ax.set(xlim=(0, 4), ylim=(-.1, 1.25), yticks=[0, 1], xlabel='Area budget δ (µm²)', ylabel='PECAM1 count in 153124')
ax.legend(frameon=False, fontsize=8.5, loc='upper left')
ax = axes[1, 1]; ax.axis('off')
ax.set_title('D  Conditional interpretation', loc='left', fontweight='bold')
ax.text(0, .84, 'δ = 1 µm²:  [0, 0]\nObserved zero is preserved.\n\nδ = 4 µm²:  [0, 1]\nDetection can change.\n\nApplies within this local family.\nBudget remains uncalibrated.',
        transform=ax.transAxes, va='top', linespacing=1.5)
dest = ROOT / 'docs/budgets'; dest.mkdir(parents=True, exist_ok=True)
fig.savefig(dest / 'budget_story.png', dpi=220)
plt.close(fig)
print(dest / 'budget_story.png')
