"""
Figuras de dos paneles (horizontal): abundancia relativa promedio por ambiente
a nivel de Filo (A) y Género (B).

Panel B — Género: top 10 más abundantes globalmente + géneros de relevancia
sanitaria (marcados con * en la leyenda y colores cálidos).

Salidas:
  doctorado_16S/05_figures/fig_filos_generos_16S.{png,tiff,pdf}
  doctorado_18S/05_figures/fig_filos_generos_18S.{png,tiff,pdf}
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
    'font.family':           'sans-serif',
    'font.size':             9,
    'axes.titlesize':        9,
    'axes.labelsize':        9,
    'xtick.labelsize':       9,
    'ytick.labelsize':       9,
    'legend.fontsize':       8,
    'legend.title_fontsize': 8,
})

# ── Constantes ────────────────────────────────────────────────────────────────
FIGSIZE      = (14.0, 5.8)
DPI_PNG      = 300
DPI_HIRES    = 600
N_TOP_PHYLUM = 12
N_TOP_GENUS  = 10
COLOR_OTROS  = '#CFCFCF'

GROUP_ORDER  = ['Oyster', 'Human feces', 'Water']
GROUP_COLORS = {
    'Oyster':      '#3A86FF',
    'Human feces': '#FF6B6B',
    'Water':       '#06D6A0',
}

# Filo — Tableau-10 + 2 extensiones
PALETTE_PHYLUM = [
    '#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F',
    '#EDC948', '#B07AA1', '#FF9DA7', '#9C755F', '#BAB0AC',
    '#499894', '#86BCB6',
]

# Género — tonos fríos/neutros para los más abundantes
PALETTE_GENUS_ABUNDANT = [
    '#0072B2', '#56B4E9', '#009E73', '#44AA99', '#332288',
    '#88CCEE', '#117733', '#999933', '#DDCC77', '#4477AA',
]

# Género — tonos cálidos para géneros de relevancia sanitaria
PALETTE_GENUS_PATHOGEN = [
    '#CC3311', '#EE7733', '#EE3377', '#AA3377', '#BB5566',
]

# Géneros de relevancia sanitaria para 16S
PATHOGENS_16S = [
    'Vibrio', 'Escherichia_Shigella', 'Aeromonas', 'Arcobacter', 'Enterococcus',
]

# Géneros de relevancia sanitaria para 18S
PATHOGENS_18S = [
    'Cryptosporidium', 'Colpodella', 'Entamoeba',
]


# ── Utilidades ────────────────────────────────────────────────────────────────
def get_sample_cols(df):
    pat = re.compile(r'^(CrassostreavirginicaG|HecesP|AguaP)\d+: rel\. freq\. \(%\)$')
    return [c for c in df.columns if pat.match(c)]


def sample_to_group(col):
    if col.startswith('Crassostrea'): return 'Oyster'
    elif col.startswith('Heces'):     return 'Human feces'
    return 'Water'


def grouped_mean(df, taxa_col, rel_cols):
    """Devuelve DataFrame (ambiente × taxón) con promedios de abundancia relativa."""
    df_agg = df.groupby(taxa_col)[rel_cols].sum()
    total  = df_agg.sum(axis=0).replace(0, np.nan)
    mat    = (df_agg.div(total, axis=1).fillna(0) * 100).T
    mat['_group'] = [sample_to_group(c) for c in mat.index]
    return mat.groupby('_group').mean().loc[GROUP_ORDER]


# ── Panel individual ──────────────────────────────────────────────────────────
def plot_taxa_panel(ax, tsv_path, taxa_col, taxa_excluir, merge_labels,
                   legend_title, n_top, palette,
                   pathogens=None,
                   palette_pathogen=None,
                   italic_legend=False,
                   legend_fontsize=8):

    df = pd.read_csv(tsv_path, sep='\t')
    df[taxa_col] = df[taxa_col].replace(merge_labels)
    df = df[~df[taxa_col].isin(taxa_excluir)].copy()

    rel_cols  = get_sample_cols(df)
    mg        = grouped_mean(df, taxa_col, rel_cols)

    # ── Selección de taxa ──────────────────────────────────────────────────
    top_n = (mg.sum(axis=0)
               .sort_values(ascending=False)
               .head(n_top).index.tolist())

    if pathogens:
        valid_path  = [p for p in pathogens if p in mg.columns]
        extra_path  = [p for p in valid_path if p not in top_n]
        # Si un patógeno ya está en top_n, lo tratamos como patógeno (color cálido)
        all_taxa    = top_n + extra_path
        path_set    = set(valid_path)
    else:
        all_taxa = top_n
        path_set = set()

    # ── Colores ────────────────────────────────────────────────────────────
    a_idx = p_idx = 0
    colors = []
    for t in all_taxa:
        if t in path_set:
            colors.append(palette_pathogen[p_idx % len(palette_pathogen)])
            p_idx += 1
        else:
            colors.append(palette[a_idx % len(palette)])
            a_idx += 1
    colors.append(COLOR_OTROS)

    # ── Construcción de la matriz a plotear ────────────────────────────────
    mat_plot = mg[all_taxa].copy()
    mat_plot['Others'] = (100 - mat_plot.sum(axis=1)).clip(lower=0)

    # Renombrar patógenos con asterisco para la leyenda
    rename = {t: f'{t} *' if t in path_set else t for t in all_taxa}
    rename['Others'] = 'Others'
    mat_plot = mat_plot.rename(columns=rename)
    col_order = [rename[t] for t in all_taxa] + ['Others']

    mat_plot[col_order].plot(
        kind='bar', stacked=True, ax=ax,
        color=colors, width=0.55, linewidth=0, legend=False,
    )

    # ── Estilo de ejes ─────────────────────────────────────────────────────
    ax.set_ylabel('Relative abundance (%)')
    ax.set_xlabel('')
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(20))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(10))
    ax.tick_params(axis='x', rotation=0, length=0, pad=6)
    ax.tick_params(axis='y', which='both', direction='out')
    ax.grid(axis='y', linestyle='--', alpha=0.45, linewidth=0.7, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.spines['left'].set_linewidth(0.8)

    for tick, group in zip(ax.get_xticklabels(), GROUP_ORDER):
        tick.set_color(GROUP_COLORS[group])
        tick.set_fontweight('bold')

    # ── Leyenda ────────────────────────────────────────────────────────────
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors]
    labels  = col_order

    # Separador visual y nota de asterisco si hay patógenos
    if path_set:
        sep    = Line2D([0], [0], color='none')
        note   = Line2D([0], [0], color='none')
        handles += [sep, note]
        labels  += ['', '* sanitary relevance']

    ax.legend(
        handles, labels,
        title=legend_title,
        bbox_to_anchor=(1.02, 1), loc='upper left',
        frameon=False, borderaxespad=0,
        fontsize=legend_fontsize,
        title_fontsize=legend_fontsize,
        handlelength=1.0, handleheight=0.85,
        labelspacing=0.25, borderpad=0,
    )

    # Cursiva en nombres de taxón (solo si italic_legend=True)
    if italic_legend:
        _skip = {'Others', '', '* sanitary relevance'}
        for text in ax.get_legend().get_texts():
            if text.get_text() not in _skip:
                text.set_style('italic')


# ── Figura de dos paneles ─────────────────────────────────────────────────────
def make_two_panel_figure(panel_a, panel_b, out_prefix):
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=FIGSIZE, constrained_layout=True,
    )
    fig.get_layout_engine().set(wspace=0.04)

    plot_taxa_panel(
        ax=ax_a, n_top=N_TOP_PHYLUM,
        palette=PALETTE_PHYLUM, palette_pathogen=[],
        italic_legend=False, legend_fontsize=10,
        tsv_path=panel_a['tsv'], taxa_col=panel_a['col'],
        taxa_excluir=panel_a.get('excluir', []),
        merge_labels=panel_a.get('merge', {}),
        legend_title=panel_a['legend_title'],
    )
    plot_taxa_panel(
        ax=ax_b, n_top=N_TOP_GENUS,
        palette=PALETTE_GENUS_ABUNDANT,
        palette_pathogen=PALETTE_GENUS_PATHOGEN,
        pathogens=panel_b.get('pathogens'),
        italic_legend=True, legend_fontsize=9,
        tsv_path=panel_b['tsv'], taxa_col=panel_b['col'],
        taxa_excluir=panel_b.get('excluir', []),
        merge_labels=panel_b.get('merge', {}),
        legend_title=panel_b['legend_title'],
    )

    for ax, label in [(ax_a, 'A)'), (ax_b, 'B)')]:
        ax.text(-0.10, 1.04, label, transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom', ha='left')

    for ext, dpi in [('png', DPI_PNG), ('tiff', DPI_HIRES), ('pdf', DPI_HIRES)]:
        fig.savefig(f'{out_prefix}.{ext}', dpi=dpi)
    plt.close(fig)
    print(f'Guardado: {out_prefix}.{{png,tiff,pdf}}')


# ── 16S ───────────────────────────────────────────────────────────────────────
make_two_panel_figure(
    panel_a={
        'tsv':          'statistics/16S/CV_vs_HE_vs_WA_filo_16S_clear.tsv',
        'col':          'Level_2',
        'excluir':      ['unknown_unclassified', 'Cyanobacteria_Chloroplast'],
        'merge':        {},
        'legend_title': 'Phylum',
    },
    panel_b={
        'tsv':          'statistics/16S/CV_vs_HE_vs_WA_genero_16S_clear.tsv',
        'col':          'Level_6',
        'excluir':      [],
        'merge':        {},
        'legend_title': 'Genus',
        'pathogens':    PATHOGENS_16S,
    },
    out_prefix='figures/16S/fig_filos_generos_16S',
)

# ── 18S ───────────────────────────────────────────────────────────────────────
make_two_panel_figure(
    panel_a={
        'tsv':          'statistics/18S/WA_vs_HE_CV_filo_18S_clear.tsv',
        'col':          'Level_2',
        'excluir':      ['Eukaryota_unclassified', 'Archaeplastida'],
        'merge':        {'Cryptista:nucl': 'Cryptista'},
        'legend_title': 'Phylum',
    },
    panel_b={
        'tsv':          'statistics/18S/WA_vs_HE_CV_genero_18S_clear.tsv',
        'col':          'Level_8',
        'excluir':      [],
        'merge':        {'Chroomonas:nucl': 'Chroomonas'},
        'legend_title': 'Genus',
        'pathogens':    PATHOGENS_18S,
    },
    out_prefix='figures/18S/fig_filos_generos_18S',
)
