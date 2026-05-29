"""
Curvas de rarefacción para 16S y 18S.
Cada muestra = una curva coloreada por grupo (Oyster/Human feces/Water).
Banda sombreada = intervalo de confianza 95%.
Leyenda con una entrada por grupo.

Salidas:
  doctorado_16S/05_figures/fig_rarefaction_16S.{png,tiff,pdf}
  doctorado_18S/05_figures/fig_rarefaction_18S.{png,tiff,pdf}
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

# ── Tipografía global ─────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.size':        10,
    'axes.labelsize':   10,
    'xtick.labelsize':  9,
    'ytick.labelsize':  9,
    'legend.fontsize':  9,
    'legend.title_fontsize': 9,
})

# ── Constantes ────────────────────────────────────────────────────────────────
FIGSIZE   = (7.0, 5.8)
DPI_PNG   = 300
DPI_HIRES = 600

GROUP_COLORS = {
    'Oyster':      '#3A86FF',
    'Human feces': '#FF6B6B',
    'Water':       '#06D6A0',
}
LINE_ALPHA = 0.75
CI_ALPHA   = 0.13
LW         = 1.2


def sample_to_group(name):
    if 'Crassostrea' in name: return 'Oyster'
    if 'Heces'       in name: return 'Human feces'
    return 'Water'


def make_rarefaction_figure(rarefaction_path, out_prefix, title=None):
    df = pd.read_csv(rarefaction_path, sep='\t', index_col=0)

    # Columnas de media (0.03-NombreMuestra)
    mean_cols = [c for c in df.columns if c.startswith('0.03-')]
    samples   = [c.replace('0.03-', '') for c in mean_cols]

    fig, ax = plt.subplots(figsize=FIGSIZE, constrained_layout=True)

    seen_groups = set()
    for mean_col, sample in zip(mean_cols, samples):
        lci_col = f'lci-{sample}'
        hci_col = f'hci-{sample}'
        group   = sample_to_group(sample)
        color   = GROUP_COLORS[group]

        x   = df.index.values
        y   = df[mean_col].values
        lci = df[lci_col].values if lci_col in df.columns else y
        hci = df[hci_col].values if hci_col in df.columns else y

        # Recortar NaN al final
        valid = ~np.isnan(y)
        x, y, lci, hci = x[valid], y[valid], lci[valid], hci[valid]

        ax.plot(x, y, color=color, lw=LW, alpha=LINE_ALPHA, zorder=3)
        ax.fill_between(x, lci, hci, color=color, alpha=CI_ALPHA, zorder=2)
        seen_groups.add(group)

    # ── Ejes ──────────────────────────────────────────────────────────────
    ax.set_xlabel('Number of sequences sampled')
    ax.set_ylabel('Number of OTUs')
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f'{int(v):,}'
    ))
    ax.yaxis.set_major_locator(mticker.AutoLocator())
    ax.tick_params(axis='both', direction='out')
    ax.grid(axis='both', linestyle='--', alpha=0.4, linewidth=0.7, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.spines['left'].set_linewidth(0.8)

    # ── Leyenda (una entrada por grupo) ───────────────────────────────────
    group_order = ['Oyster', 'Human feces', 'Water']
    handles = [
        Line2D([0], [0], color=GROUP_COLORS[g], lw=2.0, label=g)
        for g in group_order if g in seen_groups
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        loc='lower right',
        handlelength=1.8,
        labelspacing=0.4,
    )

    # Colorear etiquetas de la leyenda con el color del grupo
    for text in ax.get_legend().get_texts():
        grp = text.get_text()
        if grp in GROUP_COLORS:
            text.set_color(GROUP_COLORS[grp])
            text.set_fontweight('bold')

    for ext, dpi in [('png', DPI_PNG), ('tiff', DPI_HIRES), ('pdf', DPI_HIRES)]:
        fig.savefig(f'{out_prefix}.{ext}', dpi=dpi)
    plt.close(fig)
    print(f'Guardado: {out_prefix}.{{png,tiff,pdf}}')


# ── 16S ───────────────────────────────────────────────────────────────────────
make_rarefaction_figure(
    rarefaction_path='example_data/16S/final.opti_mcc.groups.rarefaction',
    out_prefix='figures/16S/fig_rarefaction_16S',
)

# ── 18S ───────────────────────────────────────────────────────────────────────
make_rarefaction_figure(
    rarefaction_path='example_data/18S/final.opti_mcc.groups.rarefaction',
    out_prefix='figures/18S/fig_rarefaction_18S',
)
