import zipfile, tempfile, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import arial  # noqa: F401

def read_metadata(path):
    d = {}
    with open(path) as f:
        lines = f.read().splitlines()
    for line in lines[2:]:
        sid, styp = line.split('\t')
        d[sid] = styp
    return d

def read_pcoa(qza_path):
    z = zipfile.ZipFile(qza_path)
    name = [n for n in z.namelist() if n.endswith('data/ordination.txt')][0]
    text = z.read(name).decode()
    lines = text.splitlines()

    # parse "Proportion explained" section
    prop_idx = [i for i, l in enumerate(lines) if l.split('\t')[0] == 'Proportion explained'][0]
    prop_vals = [float(x) for x in lines[prop_idx + 1].split('\t') if x.strip()]

    # parse "Site" section (sample coordinates) -- header line itself is "Site\tN\tM"
    site_idx = [i for i, l in enumerate(lines) if l.split('\t')[0] == 'Site'][0]
    n_samples, n_axes = [int(x) for x in lines[site_idx].split('\t')[1:3]]
    coords = {}
    for i in range(n_samples):
        parts = lines[site_idx + 1 + i].split('\t')
        sid = parts[0]
        vals = [float(x) for x in parts[1:]]
        coords[sid] = vals
    return coords, prop_vals

def plot_pcoa(qza_path, meta_path, title, out_path, dpi=200):
    coords, prop = read_pcoa(qza_path)
    meta = read_metadata(meta_path)
    colors = {'Water': '#2ca02c', 'Oyster': '#1f77b4', 'Feces': '#d62728'}

    fig, ax = plt.subplots(figsize=(7, 6))
    plotted_groups = set()
    for sid, xy in coords.items():
        styp = meta.get(sid, 'Unknown')
        color = colors.get(styp, '#888888')
        label = styp if styp not in plotted_groups else None
        plotted_groups.add(styp)
        ax.scatter(xy[0], xy[1], color=color, s=80, edgecolor='black', linewidth=0.5, label=label)

    pc1_pct = prop[0] * 100
    pc2_pct = prop[1] * 100
    ax.set_xlabel(f'PC1 ({pc1_pct:.1f}%)')
    ax.set_ylabel(f'PC2 ({pc2_pct:.1f}%)')
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.legend(title='Sample type')
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight')
    print(f'Saved: {out_path} | grupos graficados: {sorted(plotted_groups)} | n total puntos: {len(coords)}')
