"""
PCoA (Bray-Curtis) – 18S microbiota
Genera figura lista para publicación en revista científica.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from scipy.spatial.distance import braycurtis

# ── Rutas ────────────────────────────────────────────────────────────────────
SHARED   = "example_data/18S/final.opti_mcc.shared"
META     = "metadata/metadatos_18S.txt"
OUT_PNG  = "figures/18S/fig_PCoA_18S_pub.png"
OUT_TIFF = "figures/18S/fig_PCoA_18S_pub.tiff"

# ── Colores y etiquetas ───────────────────────────────────────────────────────
PALETTE = {"WA": "#4472C4", "CV": "#E8856A", "HE": "#4CAF50"}
LABELS  = {"WA": "WA", "CV": "CV", "HE": "HE"}


# ── 1. Cargar OTU table ───────────────────────────────────────────────────────
shared = pd.read_csv(SHARED, sep="\t", index_col="Group")
shared = shared.drop(columns=["label", "numOtus"], errors="ignore")

# Normalizar a abundancia relativa
otu = shared.div(shared.sum(axis=1), axis=0)

# ── 2. Metadatos ──────────────────────────────────────────────────────────────
meta = pd.read_csv(META, sep="\t", index_col="SampleID", usecols=[0, 1],
                   names=["SampleID", "Type"], header=0)
meta.index = meta.index.str.strip()

# Alinear
samples = otu.index.intersection(meta.index)
otu  = otu.loc[samples]
meta = meta.loc[samples]

# ── 3. Distancia Bray-Curtis ──────────────────────────────────────────────────
n = len(samples)
D = np.zeros((n, n))
arr = otu.values
for i in range(n):
    for j in range(i + 1, n):
        d = braycurtis(arr[i], arr[j])
        D[i, j] = d
        D[j, i] = d

# ── 4. PCoA (classical MDS / eigendecomposición) ─────────────────────────────
def pcoa(dist_matrix):
    A = -0.5 * (dist_matrix ** 2)
    row_means  = A.mean(axis=1, keepdims=True)
    col_means  = A.mean(axis=0, keepdims=True)
    grand_mean = A.mean()
    H = A - row_means - col_means + grand_mean

    eigenvalues, eigenvectors = np.linalg.eigh(H)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues  = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    pos = eigenvalues > 1e-10
    ev_pos  = eigenvalues[pos]
    vec_pos = eigenvectors[:, pos]

    scores    = vec_pos * np.sqrt(ev_pos)
    explained = ev_pos / ev_pos.sum() * 100
    return scores, explained

scores, explained = pcoa(D)
pc1_pct = explained[0]
pc2_pct = explained[1]

pc_df = pd.DataFrame(
    {"PC1": scores[:, 0], "PC2": scores[:, 1]},
    index=samples
)
pc_df["Type"] = meta.loc[samples, "Type"].values

# Orientar ejes para consistencia con figura original (CV a la izquierda de PC1)
cv_mean_pc1 = pc_df.loc[pc_df["Type"] == "CV", "PC1"].mean()
wa_mean_pc1 = pc_df.loc[pc_df["Type"] == "WA", "PC1"].mean()
if cv_mean_pc1 > wa_mean_pc1:          # si CV queda a la derecha, invertir PC1
    pc_df["PC1"] = -pc_df["PC1"]
cv_mean_pc2 = pc_df.loc[pc_df["Type"] == "CV", "PC2"].mean()
wa_mean_pc2 = pc_df.loc[pc_df["Type"] == "WA", "PC2"].mean()
if cv_mean_pc2 > wa_mean_pc2:          # si CV queda arriba de WA, invertir PC2
    pc_df["PC2"] = -pc_df["PC2"]

# ── 5. Elipse de confianza (95 %) ─────────────────────────────────────────────
def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """Dibuja elipse de confianza basada en covarianza."""
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    order = eigenvalues.argsort()[::-1]
    eigenvalues  = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(eigenvalues)
    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width, height=height, angle=angle,
        **kwargs
    )
    ax.add_patch(ellipse)


# ── 6. Figura ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8.0, 5.5))

# Fondo y cuadrícula ultra-tenue
ax.set_facecolor("white")
fig.patch.set_facecolor("white")
ax.grid(True, color="#EFEFEF", linewidth=0.4, zorder=0)
ax.set_axisbelow(True)

# Solo ejes X e Y (sin spine superior ni derecho)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Orden de grupos para la leyenda
group_order = ["WA", "CV", "HE"]

# Dibujar elipses primero (capa inferior)
for grp in group_order:
    sub = pc_df[pc_df["Type"] == grp]
    if len(sub) < 3:
        continue
    confidence_ellipse(
        sub["PC1"].values, sub["PC2"].values, ax,
        n_std=2.0,
        facecolor=PALETTE[grp], alpha=0.10,
        edgecolor=PALETTE[grp], linewidth=1.8,
        linestyle="-", zorder=1, clip_on=True
    )

# Dibujar puntos encima (capa superior)
for grp in group_order:
    sub = pc_df[pc_df["Type"] == grp]
    ax.scatter(
        sub["PC1"], sub["PC2"],
        color=PALETTE[grp],
        s=70, zorder=3,
        edgecolors="none",
        label=LABELS[grp]
    )

# Ejes — 9 pt en figura final (~6.5 cm ancho de columna doble)
ax.set_xlabel(f"PC1 ({pc1_pct:.2f}%)", fontsize=9)
ax.set_ylabel(f"PC2 ({pc2_pct:.2f}%)", fontsize=9)
ax.tick_params(axis="both", labelsize=8)

# Límites automáticos con margen sobre los puntos (no sobre la elipse)
pad = 0.08
x_vals = pc_df["PC1"].values
y_vals = pc_df["PC2"].values
x_range = x_vals.max() - x_vals.min()
y_range = y_vals.max() - y_vals.min()
ax.set_xlim(x_vals.min() - pad * x_range, x_vals.max() + pad * x_range)
ax.set_ylim(y_vals.min() - pad * y_range, y_vals.max() + pad * y_range)

# Líneas de cuadrante en 0 (divisiones visibles)
ax.axhline(0, color="#999999", linewidth=0.8, zorder=1)
ax.axvline(0, color="#999999", linewidth=0.8, zorder=1)

# Spines delgados
for spine in ax.spines.values():
    spine.set_linewidth(0.9)

# Sin título (va como leyenda de figura en el manuscrito)
ax.set_title("")

# Leyenda externa
legend = ax.legend(
    title="Tipo de muestra",
    title_fontsize=9,
    fontsize=8.5,
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    borderaxespad=0,
    frameon=False,
    handletextpad=0.4,
    labelspacing=0.6
)
legend.get_title().set_fontweight("bold")

plt.tight_layout(rect=[0, 0, 0.82, 1])

# ── 7. Exportar ───────────────────────────────────────────────────────────────
fig.savefig(OUT_PNG,  dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_TIFF, dpi=600, bbox_inches="tight", facecolor="white")
print(f"PC1 = {pc1_pct:.2f}%  |  PC2 = {pc2_pct:.2f}%")
print(f"Guardado: {OUT_PNG}")
print(f"Guardado: {OUT_TIFF}")
