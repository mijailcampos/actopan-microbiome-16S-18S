import os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('DATA', os.path.join(HERE, 'data'))  # inputs, see README
OUT = os.environ.get('OUT', os.path.join(HERE, 'out'))
os.makedirs(OUT, exist_ok=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, HERE); import arial  # noqa: E402,F401
import numpy as np

with open(os.path.join(DATA, 'table.tsv')) as f:
    lines = f.read().splitlines()
header = lines[1].split('\t')
samples = header[1:]
data = {s: {} for s in samples}       # genus-level (best available label)
data_phy = {s: {} for s in samples}   # phylum-level
def best_label(taxstr):
    ranks = [p.strip() for p in taxstr.split(';')]
    for p in reversed(ranks):
        name = p.split('__', 1)[-1].strip('_') if '__' in p else p.strip('_')
        if name:
            prefix = p.split('__')[0] if '__' in p else ''
            rank_map = {'g':'', 'f':' (fam.)', 'o':' (ord.)', 'c':' (class)', 'p':' (phy.)', 'sk':'', 'd':''}
            return name + rank_map.get(prefix, '')
    return 'Unresolved'

def phylum_label(taxstr):
    ranks = [p.strip() for p in taxstr.split(';')]
    for p in ranks:
        if p.startswith('p__'):
            name = p[3:].strip('_')
            if name:
                return name
    return 'Unresolved (phylum)'

for line in lines[2:]:
    parts = line.split('\t')
    taxon = best_label(parts[0])
    phy = phylum_label(parts[0])
    vals = [float(x) for x in parts[1:]]
    for s, v in zip(samples, vals):
        if v > 0:
            data[s][taxon] = data[s].get(taxon, 0) + v
            data_phy[s][phy] = data_phy[s].get(phy, 0) + v

WATER = ['AguaP1_S252_L001','AguaP2_S264_L001','AguaP3_S276_L001','AguaP4_S288_L001','AguaP5_S203_L001']
FECES = ['HecesP6_S215_L001','HecesP7_S227_L001','HecesP8_S239_L001','HecesP9_S251_L001','HecesP10_S263_L001']
FECES_LABEL = {'HecesP6_S215_L001':'F1','HecesP7_S227_L001':'F2','HecesP8_S239_L001':'F3',
               'HecesP9_S251_L001':'F4','HecesP10_S263_L001':'F5'}
WATER_LABEL = {s: f'W{i+1}' for i, s in enumerate(WATER)}

def top_n_plus_other(sample_list, src, n=6):
    all_taxa = {}
    for s in sample_list:
        total = sum(src[s].values())
        for t, v in src[s].items():
            all_taxa[t] = all_taxa.get(t, 0) + v / total
    top = sorted(all_taxa, key=lambda t: -all_taxa[t])[:n]
    return top

def stacked_matrix(sample_list, src, top_taxa):
    mat = np.zeros((len(top_taxa) + 1, len(sample_list)))
    for j, s in enumerate(sample_list):
        total = sum(src[s].values())
        other = 100.0
        for i, t in enumerate(top_taxa):
            pct = 100 * src[s].get(t, 0) / total if total > 0 else 0
            mat[i, j] = pct
            other -= pct
        mat[-1, j] = max(other, 0)
    return mat

def build_color_map(top_a, top_b, pal_a, pal_b, shared_pal):
    """Shared taxa (present in both top-N lists) get a common color from shared_pal;
    taxa unique to each group keep their own palette."""
    shared = [t for t in top_a if t in top_b]
    color_map = {}
    for i, t in enumerate(shared):
        color_map[t] = shared_pal[i % len(shared_pal)]
    ia = 0
    for t in top_a:
        if t in color_map:
            continue
        color_map[t] = pal_a[ia % len(pal_a)]
        ia += 1
    ib = 0
    for t in top_b:
        if t in color_map:
            continue
        color_map[t] = pal_b[ib % len(pal_b)]
        ib += 1
    return color_map

def draw_panel(ax, sample_list, src, labels_map, title, color_map, others_color='#BAB0AC', n=6, depth_labels=False):
    top = top_n_plus_other(sample_list, src, n=n)
    mat = stacked_matrix(sample_list, src, top)
    labels = top + ['Others']
    x = np.arange(len(sample_list))
    bottom = np.zeros(len(sample_list))
    for i, lab in enumerate(labels):
        c = others_color if lab == 'Others' else color_map[lab]
        ax.bar(x, mat[i], bottom=bottom, color=c, label=lab, width=0.65)
        bottom += mat[i]
    ax.set_xticks(x)
    ax.set_xticklabels([labels_map[s] for s in sample_list])
    ax.set_ylabel('Relative abundance (%)')
    ax.set_title(title, fontsize=10.5, fontweight='bold')
    if depth_labels:
        depths = [sum(src[s].values()) for s in sample_list]
        ax.set_ylim(0, 118)
        for xi, d in zip(x, depths):
            ax.text(xi, 103, f'{int(d)}\nreads', ha='center', va='bottom', fontsize=7, color='#555')
    else:
        ax.set_ylim(0, 100)
    ax.legend(fontsize=6.8, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)

# Water and feces each get their own palette; taxa shared between the two (at a given
# taxonomic level) get a common color from a third, distinct palette instead.
palette_water = ['#264653', '#2A9D8F', '#E9C46A', '#F4A261', '#E76F51', '#8AB17D']
palette_feces = ['#3A0CA3', '#7209B7', '#F72585', '#4CC9F0', '#4361EE', '#560BAD']
palette_shared = ['#FFB703', '#FB8500', '#9E2A2B', '#540B0E']
OTHERS_COLOR = '#BAB0AC'

top_water_phy = top_n_plus_other(WATER, data_phy)
top_feces_phy = top_n_plus_other(FECES, data_phy)
top_water_gen = top_n_plus_other(WATER, data)
top_feces_gen = top_n_plus_other(FECES, data)

cmap_phy = build_color_map(top_water_phy, top_feces_phy, palette_water, palette_feces, palette_shared)
cmap_gen = build_color_map(top_water_gen, top_feces_gen, palette_water, palette_feces, palette_shared)

fig, axes = plt.subplots(2, 2, figsize=(12, 10.5))

draw_panel(axes[0,0], WATER, data_phy, WATER_LABEL,
           '(a) Water (n = 5), phylum level\n(top 6 + Others)', color_map=cmap_phy)
draw_panel(axes[0,1], WATER, data, WATER_LABEL,
           '(b) Water (n = 5), genus level\n(top 6 + Others)', color_map=cmap_gen)
draw_panel(axes[1,0], FECES, data_phy, FECES_LABEL,
           '(c) Feces (n = 5 independent family pools),\nphylum level (top 6 + Others)', color_map=cmap_phy, depth_labels=True)
draw_panel(axes[1,1], FECES, data, FECES_LABEL,
           '(d) Feces (n = 5 independent family pools),\ngenus level (top 6 + Others)', color_map=cmap_gen, depth_labels=True)

plt.tight_layout()
out = os.path.join(OUT, 'Figure4.png')
import os
os.makedirs(os.path.dirname(out), exist_ok=True)
plt.savefig(out, dpi=600, bbox_inches='tight')
print('Saved:', out)
print('Feces non-host depths:', {FECES_LABEL[s]: sum(data[s].values()) for s in FECES})
