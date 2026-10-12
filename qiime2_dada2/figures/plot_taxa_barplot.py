import sys, zipfile, tempfile
import biom
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import arial  # noqa: F401

def load_biom(qza_path):
    z = zipfile.ZipFile(qza_path)
    name = [n for n in z.namelist() if n.endswith('feature-table.biom')][0]
    with tempfile.NamedTemporaryFile(suffix='.biom', delete=False) as f:
        f.write(z.read(name))
        path = f.name
    return biom.load_table(path)

def read_metadata(path):
    d = {}
    with open(path) as f:
        lines = f.read().splitlines()
    for line in lines[2:]:
        sid, styp = line.split('\t')
        d[sid] = styp
    return d

def short_name(taxon, rank_prefix):
    # take the last non-empty rank token, strip prefix
    parts = [p.strip() for p in taxon.split(';')]
    for p in reversed(parts):
        p2 = p.strip('_')
        if p2 and not p2.endswith('__'):
            # strip rank prefix like g__ or p__
            if '__' in p2:
                p2 = p2.split('__', 1)[1]
            if p2:
                return p2
    return taxon

def group_means(biom_table, meta, groups_order=('Water', 'Oyster', 'Feces')):
    ids = list(biom_table.ids(axis='observation'))
    samp_ids = biom_table.ids(axis='sample')
    groups = {g: [] for g in groups_order}
    for sid in samp_ids:
        g = meta.get(sid)
        if g in groups:
            groups[g].append(sid)
    data = {}  # taxon -> {group: mean%}
    for taxon in ids:
        row = {}
        for g in groups_order:
            sids = groups[g]
            if not sids:
                row[g] = 0.0
                continue
            vals = [biom_table.get_value_by_ids(taxon, s) * 100 for s in sids]
            row[g] = float(np.mean(vals))
        data[taxon] = row
    return data, groups

def make_figure(marker, phylum_qza, genus_qza, meta_path, out_path, rank_prefix_phylum, rank_prefix_genus, title_suffix='', sanitary_genera=frozenset()):
    meta = read_metadata(meta_path)
    groups_order = ('Water', 'Oyster', 'Feces')
    colors = {'Water': '#2ca02c', 'Oyster': '#1f77b4', 'Feces': '#d62728'}

    # Panel A: phylum, top 8 + Others
    pt = load_biom(phylum_qza)
    pdata, groups = group_means(pt, meta, groups_order)
    overall = {t: np.mean(list(v.values())) for t, v in pdata.items()}
    phyla_sorted_all = sorted(pdata.keys(), key=lambda t: -overall[t])
    top8 = phyla_sorted_all[:8]
    rest = phyla_sorted_all[8:]
    phylum_rows = {}
    for t in top8:
        phylum_rows[short_name(t, rank_prefix_phylum)] = pdata[t]
    rest_row = {g: sum(pdata[t][g] for t in rest) for g in groups_order}
    phylum_rows['Others'] = rest_row
    phyla_labels = list(phylum_rows.keys())

    # Panel B: genus, top 10 + Others
    gt = load_biom(genus_qza)
    gdata, _ = group_means(gt, meta, groups_order)
    overall_g = {t: np.mean(list(v.values())) for t, v in gdata.items()}
    genus_sorted_all = sorted(gdata.keys(), key=lambda t: -overall_g[t])
    top10 = genus_sorted_all[:10]
    others = genus_sorted_all[10:]
    genus_rows = {}
    for t in top10:
        name = short_name(t, rank_prefix_genus)
        label = f'{name} *' if name in sanitary_genera else name
        genus_rows[label] = gdata[t]
    others_row = {g: sum(gdata[t][g] for t in others) for g in groups_order}
    genus_rows['Others'] = others_row
    genus_labels = list(genus_rows.keys())

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A stacked bars (bars = environments, segments = top 8 phyla + Others)
    ax = axes[0]
    bottoms = np.zeros(len(groups_order))
    cmap = plt.get_cmap('tab20')
    for i, label in enumerate(phyla_labels):
        color = '#888888' if label == 'Others' else cmap(i % 20)
        vals = np.array([phylum_rows[label][g] for g in groups_order])
        ax.bar(groups_order, vals, bottom=bottoms, label=label, color=color)
        bottoms += vals
    ax.set_ylabel('Relative abundance (%)')
    ax.set_title(f'(a) Phylum-level composition (top 8 + Others){title_suffix}')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, ncol=1)
    ax.set_ylim(0, 100)

    # Panel B stacked bars (top10 genera + Others)
    ax2 = axes[1]
    bottoms2 = np.zeros(len(groups_order))
    cmap2 = plt.get_cmap('tab20b')
    for i, label in enumerate(genus_labels):
        vals = np.array([genus_rows[label][g] for g in groups_order])
        color = '#888888' if label == 'Others' else cmap2(i % 20)
        ax2.bar(groups_order, vals, bottom=bottoms2, label=label, color=color)
        bottoms2 += vals
    ax2.set_ylabel('Relative abundance (%)')
    ax2.set_title(f'(b) Genus-level composition (top 10 + Others){title_suffix}')
    leg2 = ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, ncol=1)
    for text, label in zip(leg2.get_texts(), genus_labels):
        if label != 'Others':
            text.set_fontstyle('italic')
    ax2.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    print(f'Saved: {out_path}')
