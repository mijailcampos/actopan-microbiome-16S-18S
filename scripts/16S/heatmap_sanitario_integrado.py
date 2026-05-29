"""
Heatmap integrado de taxones de relevancia sanitaria (16S + 18S).
Filas: Vibrio, Aeromonas, Escherichia/Shigella, Cryptosporidium, Entamoeba
Columnas: Water, Oyster, Human feces
Valores: log10(abundancia relativa promedio % + 0.001)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings("ignore")

# ── Rutas ─────────────────────────────────────────────────────────────────────
PATH_16S = "example_data/16S/final.opti_mcc.0.03.cons.tax.summary"
PATH_18S = "example_data/18S/final.opti_mcc.0.03.cons.tax.summary"
OUT_PDF  = "figures/16S/heatmap_sanitario_integrado.pdf"
OUT_TIFF = "figures/16S/heatmap_sanitario_integrado.tiff"
OUT_PNG  = "figures/16S/heatmap_sanitario_integrado.png"

# ── Agrupación de muestras ────────────────────────────────────────────────────
GROUPS = {
    "Water":        ["AguaP1", "AguaP2", "AguaP3", "AguaP4", "AguaP5"],
    "Oyster":       ["CrassostreavirginicaG1", "CrassostreavirginicaG10",
                     "CrassostreavirginicaG2", "CrassostreavirginicaG3",
                     "CrassostreavirginicaG4", "CrassostreavirginicaG5",
                     "CrassostreavirginicaG6", "CrassostreavirginicaG7",
                     "CrassostreavirginicaG8", "CrassostreavirginicaG9"],
    "Human\nfeces": ["HecesP6", "HecesP7", "HecesP8", "HecesP9",
                     "HecesP10", "HecesP11", "HecesP12", "HecesP13", "HecesP14"],
}

# (nombre_display, patrón, dataset)
TAXA_CONFIG = [
    ("Vibrio",               "^Vibrio$",              "16S"),
    ("Aeromonas",            "^Aeromonas$",           "16S"),
    ("Escherichia/Shigella", "^Escherichia_Shigella$", "16S"),
    ("Cryptosporidium",      "^Cryptosporidium$",     "18S"),
    ("Entamoeba",            "^Entamoeba$",           "18S"),
]

PSEUDO = 0.001   # constante para log10(x + PSEUDO)

# ── Helpers ───────────────────────────────────────────────────────────────────
def load_relative_abundance(path):
    df = pd.read_csv(path, sep="\t")
    root = df[df["taxon"] == "Root"].iloc[0]
    sample_cols = [c for c in df.columns
                   if c not in ("taxlevel", "rankID", "taxon", "daughterlevels", "total")]
    totals = root[sample_cols].astype(float)
    rel = df.copy()
    rel[sample_cols] = df[sample_cols].astype(float).div(totals) * 100
    return rel, sample_cols


def get_taxon_row(df, pattern, sample_cols):
    rows = df[df["taxon"].str.match(pattern)]
    if rows.empty:
        return pd.Series(np.zeros(len(sample_cols)), index=sample_cols)
    exact = rows[rows["taxon"].str.fullmatch(pattern.strip("^$"))]
    return (exact if not exact.empty else rows).iloc[0][sample_cols]


# ── Cargar y construir matriz de abundancias crudas (%) ───────────────────────
rel16, cols16 = load_relative_abundance(PATH_16S)
rel18, cols18 = load_relative_abundance(PATH_18S)

env_labels  = list(GROUPS.keys())
taxa_labels = [t[0] for t in TAXA_CONFIG]
n_taxa, n_env = len(TAXA_CONFIG), len(env_labels)

raw = np.zeros((n_taxa, n_env))
for i, (_, pattern, origin) in enumerate(TAXA_CONFIG):
    rel, scols = (rel16, cols16) if origin == "16S" else (rel18, cols18)
    vals = get_taxon_row(rel, pattern, scols)
    for j, members in enumerate(GROUPS.values()):
        valid = [m for m in members if m in vals.index]
        raw[i, j] = vals[valid].mean() if valid else 0.0

# Transformación log10
log_data = np.log10(raw + PSEUDO)

print("Abundancias relativas promedio (%):")
print(pd.DataFrame(raw, index=taxa_labels, columns=env_labels).round(3).to_string())
print(f"\nRango log10: {log_data.min():.2f} → {log_data.max():.2f}")

# ── Paleta blanco → rojo oscuro ───────────────────────────────────────────────
cmap = LinearSegmentedColormap.from_list(
    "wr", ["#FFFFFF", "#FEE5D9", "#FC9272", "#EF3B2C", "#99000D"], N=256
)

# ── Figura ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7.0, 5.8))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

vmin, vmax = log_data.min(), log_data.max()
im = ax.imshow(log_data, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)

# Etiquetas de ejes
ax.set_xticks(range(n_env))
ax.set_xticklabels(env_labels, fontsize=10, fontweight="bold", color="#111111")
ax.set_yticks(range(n_taxa))
ax.set_yticklabels(taxa_labels, fontsize=10, color="#111111", style="italic")
ax.tick_params(which="both", length=0)
ax.xaxis.set_tick_params(pad=6)

# Anotación de valores crudos en celdas
threshold = vmin + (vmax - vmin) * 0.62
for i in range(n_taxa):
    for j in range(n_env):
        v = raw[i, j]
        lv = log_data[i, j]
        txt = f"{v:.2f}%" if v >= 0.01 else "0%"
        fc  = "white" if lv > threshold else "#222222"
        ax.text(j, i, txt, ha="center", va="center",
                fontsize=8, color=fc)

# Separador discreto entre bacterias y protozoarios
ax.axhline(2.5, color="#888888", linewidth=0.8, linestyle=(0, (4, 3)))

# Cuadrícula de celdas muy sutil
ax.set_xticks(np.arange(-0.5, n_env, 1), minor=True)
ax.set_yticks(np.arange(-0.5, n_taxa, 1), minor=True)
ax.grid(which="minor", color="#DDDDDD", linewidth=0.5)
ax.tick_params(which="minor", length=0)

# Sin bordes
for sp in ax.spines.values():
    sp.set_visible(False)

# Anotaciones laterales Bacteria / Protozoa
bkw = dict(xycoords="axes fraction", fontsize=8, fontstyle="italic",
           color="#444444", ha="center", annotation_clip=False)
ax.annotate("Bacteria", xy=(-0.36, 0.72), rotation=90, va="center", **bkw)
ax.annotate("Protozoa", xy=(-0.36, 0.18), rotation=90, va="center", **bkw)

# Corchetes
ckw = dict(xycoords="axes fraction", textcoords="axes fraction",
           arrowprops=dict(arrowstyle="-", color="#777777", lw=0.9),
           annotation_clip=False)
ax.annotate("", xy=(-0.27, 1.01), xytext=(-0.27, 0.43), **ckw)
ax.annotate("", xy=(-0.27, 0.35), xytext=(-0.27, -0.01), **ckw)

# Título compacto dentro de la figura
ax.set_title("Health-relevant taxa · 16S + 18S", fontsize=9,
             color="#444444", pad=8, loc="center")

# Colorbar etiquetada en escala log10 con ticks significativos
cbar = fig.colorbar(im, ax=ax, fraction=0.044, pad=0.04)
cbar.set_label("log₁₀(rel. abundance % + 0.001)", fontsize=8, color="#333333")
cbar.ax.tick_params(labelsize=7.5)
cbar.outline.set_visible(False)

# Ticks de la barra en valores de log10 redondos dentro del rango
log_ticks = [t for t in [-3, -2, -1, 0, 1, 2] if vmin <= t <= vmax]
cbar.set_ticks(log_ticks)
cbar.set_ticklabels([f"{t}" for t in log_ticks])

plt.tight_layout(pad=0.8)

# ── Exportar ──────────────────────────────────────────────────────────────────
fig.savefig(OUT_PDF,  dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_TIFF, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PNG,  dpi=300, bbox_inches="tight", facecolor="white")

print(f"\nFiguras guardadas en:")
print(f"  {OUT_PDF}")
print(f"  {OUT_TIFF}")
print(f"  {OUT_PNG}")
