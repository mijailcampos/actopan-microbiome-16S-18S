import os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('DATA', os.path.join(HERE, 'data'))  # inputs, see README
OUT = os.environ.get('OUT', os.path.join(HERE, 'out'))
os.makedirs(OUT, exist_ok=True)
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, HERE); import arial  # noqa: E402,F401
from matplotlib.colors import Normalize

# Figure 5 (LAJAR): bacterial genera of public health relevance, mean relative abundance (%)
# per sample type, read directly from the n=17 Kruskal-Wallis table (feces F1-F5 only).
CSV = os.path.join(DATA, 'Kruskal_Wallis_all_taxa_n17.csv')
OUT = os.path.join(OUT, 'Figure5')

rows = ['Vibrio', 'Aeromonas', 'Escherichia-Shigella', 'Arcobacter', 'Enterococcus']
col_keys = ['Water', 'Oyster', 'Feces']
cols = ['Water', 'Oyster', 'Human\nfeces']

means = {}
with open(CSV) as f:
    for r in csv.DictReader(f):
        if r['marker'] == '16S' and r['level'] == 'genus':
            g = r['taxon'].split('g__')[-1]
            if g in rows:
                means[g] = [float(r[f'mean_pct_{k}']) for k in col_keys]
assert set(means) == set(rows), means.keys()

pct = np.array([means[g] for g in rows])
matrix = np.log10(pct + 0.001)

fig, ax = plt.subplots(figsize=(7.2, 4.4))
norm = Normalize(vmin=-3, vmax=matrix.max())
im = ax.imshow(matrix, cmap='Reds', norm=norm, aspect='auto')

for i in range(len(rows)):
    for j in range(len(cols)):
        text_color = 'white' if matrix[i, j] > (norm.vmin + norm.vmax) / 2 else 'black'
        ax.text(j, i, f'{pct[i, j]:.2f}%', ha='center', va='center', color=text_color, fontsize=9)

ax.set_xticks(range(len(cols)))
ax.set_xticklabels(cols, fontweight='bold')
ax.set_yticks(range(len(rows)))
ax.set_yticklabels(rows, fontstyle='italic')
ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=2)
ax.tick_params(which='minor', length=0)
for spine in ax.spines.values():
    spine.set_visible(False)

cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('log$_{10}$(rel. abundance % + 0.001)')

for ext in ('png',):
    kw = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tiff' else {}
    plt.savefig(f'{OUT}.{ext}', dpi=600, bbox_inches='tight', **kw)
    print(f'Saved: {OUT}.{ext}')
for g, v in zip(rows, pct):
    print(g, ' '.join(f'{x:.2f}' for x in v))
