"""
Genera las figuras de abundancia relativa a nivel de filo (combinadas, 3 ambientes)
para 16S y 18S con estética idéntica.

Salidas:
  doctorado_16S/05_figures/fig_filos_combinados_16S.{png,tiff,pdf}
  doctorado_18S/05_figures/fig_filos_combinados_18S.{png,tiff,pdf}
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms
from matplotlib.cm import get_cmap

# ── Estética unificada ────────────────────────────────────────────────────────
FIGSIZE      = (10.0, 5.8)
DPI_PNG      = 300
DPI_HIRES    = 600          # TIFF y PDF
CMAP_NAME    = 'tab20'
COLOR_OTROS  = 'lightgrey'
N_TOP        = 15
FONTSIZE_AX  = 10
FONTSIZE_TK  = 9
FONTSIZE_LEG = 9

AGUA_LABELS = {
    'AguaP1': 'OCW', 'AguaP2': 'HZW', 'AguaP3': 'GZW',
    'AguaP4': 'RMW', 'AguaP5': 'CZW',
}

GROUP_COLORS = {
    'Oyster':      '#3A86FF',
    'Human feces': '#FF6B6B',
    'Water':       '#06D6A0',
}


def get_sample_cols(df):
    """Extrae solo columnas de frecuencia relativa de muestras individuales."""
    pat = re.compile(r'^(CrassostreavirginicaG|HecesP|AguaP)\d+: rel\. freq\. \(%\)$')
    return [c for c in df.columns if pat.match(c)]


def rename_and_order(cols):
    """Devuelve dict de renombre y lista ordenada (Ostión → Heces → Agua)."""
    g_cols = sorted([c for c in cols if c.startswith('Crassostrea')],
                    key=lambda c: int(re.search(r'G(\d+)', c).group(1)))
    h_cols = sorted([c for c in cols if c.startswith('Heces')],
                    key=lambda c: int(re.search(r'P(\d+)', c).group(1)))
    a_cols = sorted([c for c in cols if c.startswith('Agua')],
                    key=lambda c: int(re.search(r'P(\d+)', c).group(1)))

    rename, order = {}, []
    for c in g_cols:
        n = re.search(r'G(\d+)', c).group(1)
        lbl = f'O{n}'; rename[c] = lbl; order.append(lbl)
    for i, c in enumerate(h_cols):
        lbl = f'F{i+1}'; rename[c] = lbl; order.append(lbl)
    for c in a_cols:
        key = re.match(r'^(AguaP\d+)', c).group(1)
        lbl = AGUA_LABELS.get(key, key); rename[c] = lbl; order.append(lbl)
    return rename, order


def add_group_bar(ax, sample_order):
    n_oyster = sum(1 for s in sample_order if re.match(r'^O\d+$', s))
    n_heces  = sum(1 for s in sample_order if re.match(r'^F\d+$', s))
    n_agua   = len(sample_order) - n_oyster - n_heces

    groups = [
        (-0.5,                                   n_oyster - 0.5,                    'Oyster'),
        (n_oyster - 0.5,                         n_oyster + n_heces - 0.5,          'Human feces'),
        (n_oyster + n_heces - 0.5,               n_oyster + n_heces + n_agua - 0.5, 'Water'),
    ]

    trans  = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    bar_y  = -0.22
    bar_h  = 0.07

    for x0, x1, label in groups:
        color = GROUP_COLORS[label]
        rect = mpatches.Rectangle(
            (x0, bar_y), x1 - x0, bar_h,
            transform=trans, facecolor=color,
            edgecolor='white', linewidth=0.8,
            clip_on=False, zorder=5,
        )
        ax.add_patch(rect)
        ax.text(
            (x0 + x1) / 2, bar_y + bar_h / 2, label,
            transform=trans, ha='center', va='center',
            fontsize=8, color='white', fontweight='bold',
            clip_on=False, zorder=6,
        )


def plot_filos(tsv_path, filo_col, filos_excluir, merge_labels, out_prefix):
    df = pd.read_csv(tsv_path, sep='\t')

    # Filtrar filos no deseados y unificar etiquetas
    df[filo_col] = df[filo_col].replace(merge_labels)
    df = df[~df[filo_col].isin(filos_excluir)].copy()

    # Extraer y renombrar columnas de muestras individuales
    rel_cols = get_sample_cols(df)
    rename_dict, sample_order = rename_and_order(rel_cols)

    # Agregar por filo, transponer (muestras = filas)
    df_agg = df.groupby(filo_col)[rel_cols].sum()
    df_agg.columns = [rename_dict[c] for c in df_agg.columns]
    df_agg = df_agg[sample_order]

    # Normalizar cada muestra al 100%
    total = df_agg.sum(axis=0).replace(0, np.nan)
    mat = (df_agg.div(total, axis=1).fillna(0) * 100).T   # muestras × filos

    # Top N filos por abundancia global
    top = mat.sum(axis=0).sort_values(ascending=False).head(N_TOP).index.tolist()
    mat_plot = mat[top].copy()
    mat_plot['Otros'] = (100 - mat_plot.sum(axis=1)).clip(lower=0)

    # Paleta unificada
    cmap = get_cmap(CMAP_NAME)
    colors = [cmap(i / N_TOP) for i in range(len(top))] + [COLOR_OTROS]

    # Figura
    fig, ax = plt.subplots(figsize=FIGSIZE)
    mat_plot[top + ['Otros']].plot(
        kind='bar', stacked=True, ax=ax,
        color=colors, width=0.8, linewidth=0
    )

    ax.set_ylabel('Relative abundance (%)', fontsize=FONTSIZE_AX)
    ax.set_xlabel('')
    ax.set_ylim(0, 100)
    ax.tick_params(axis='x', labelsize=FONTSIZE_TK, rotation=45)
    ax.tick_params(axis='y', labelsize=FONTSIZE_TK)
    ax.set_xticklabels(ax.get_xticklabels(), ha='right')
    ax.grid(axis='y', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles, labels,
        title='Phylum', title_fontsize=FONTSIZE_LEG,
        fontsize=FONTSIZE_LEG,
        bbox_to_anchor=(1.02, 1), loc='upper left',
        frameon=False,
    )

    add_group_bar(ax, sample_order)

    plt.tight_layout()
    fig.savefig(f'{out_prefix}.png',  dpi=DPI_PNG,   bbox_inches='tight')
    fig.savefig(f'{out_prefix}.tiff', dpi=DPI_HIRES, bbox_inches='tight')
    fig.savefig(f'{out_prefix}.pdf',  dpi=DPI_HIRES, bbox_inches='tight')
    plt.close(fig)
    print(f'Guardado: {out_prefix}.{{png,tiff,pdf}}')


# ── 16S ───────────────────────────────────────────────────────────────────────
plot_filos(
    tsv_path    = 'statistics/16S/CV_vs_HE_vs_WA_filo_16S_clear.tsv',
    filo_col    = 'Level_2',
    filos_excluir = ['unknown_unclassified', 'Cyanobacteria_Chloroplast'],
    merge_labels  = {},
    out_prefix  = 'figures/16S/fig_filos_combinados_16S',
)

# ── 18S ───────────────────────────────────────────────────────────────────────
plot_filos(
    tsv_path    = 'statistics/18S/WA_vs_HE_CV_filo_18S_clear.tsv',
    filo_col    = 'Level_2',
    filos_excluir = ['Eukaryota_unclassified', 'Archaeplastida'],
    merge_labels  = {'Cryptista:nucl': 'Cryptista'},
    out_prefix  = 'figures/18S/fig_filos_combinados_18S',
)
