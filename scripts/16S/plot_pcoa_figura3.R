library(tidyverse)
library(vegan)    # vegdist()

# Run from the repo root: Rscript scripts/16S/plot_pcoa_figura3.R
shared <- read.delim(
  "example_data/16S/final.opti_mcc.shared",
  header   = TRUE,
  row.names = 2,
  sep      = "\t"
)
shared <- shared[, -c(1, 2)]   # quitar columnas 'label' y 'numOtus'

bc_dist <- vegdist(shared, method = "bray")

# ── 2. PCoA con cmdscale ──────────────────────────────────────────────────────
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)

eig     <- pcoa_res$eig
eig_pos <- sum(eig[eig > 0])
pct     <- eig / eig_pos * 100

pc1_pct <- round(pct[1], 1)
pc2_pct <- round(pct[2], 1)
cat(sprintf("PC1: %.1f%%   PC2: %.1f%%\n", pc1_pct, pc2_pct))

# ── 3. Coordenadas + metadatos ────────────────────────────────────────────────
coords <- as.data.frame(pcoa_res$points)
colnames(coords) <- c("PC1", "PC2")
coords$muestra   <- rownames(coords)

coords <- coords %>%
  mutate(
    grupo = case_when(
      grepl("^Agua",        muestra) ~ "WA",
      grepl("^Crassostrea", muestra) ~ "CV",
      grepl("^Heces",       muestra) ~ "HE"
    ),
    grupo = factor(grupo, levels = c("WA", "CV", "HE"))
  )

# ── 4. Colores y etiquetas ────────────────────────────────────────────────────
colores <- c("WA" = "#4393C3", "CV" = "#F4A582", "HE" = "#4DAC26")

# ── 5. Gráfica ────────────────────────────────────────────────────────────────
p <- ggplot(coords, aes(x = PC1, y = PC2, color = grupo, fill = grupo)) +

  # Elipses: relleno muy suave, borde delgado para no tapar puntos
  stat_ellipse(
    geom        = "polygon",
    level       = 0.95,
    alpha       = 0.07,
    linewidth   = 0.55,
    show.legend = FALSE
  ) +

  # Puntos
  geom_point(size = 3.2, alpha = 0.92, stroke = 0.25) +

  labs(
    x     = paste0("PC1 (", pc1_pct, "%)"),
    y     = paste0("PC2 (", pc2_pct, "%)"),
    color = "Tipo de muestra",
    fill  = "Tipo de muestra"
  ) +

  scale_color_manual(values = colores) +
  scale_fill_manual( values = colores) +

  theme_classic(base_size = 12) +
  theme(
    axis.text        = element_text(size = 10, color = "black"),
    axis.title       = element_text(size = 11),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 10),
    legend.position  = "right",
    legend.key       = element_rect(fill = "white", color = NA),
    # Cuadrícula muy suave — solo referencia visual
    panel.grid.major = element_line(color = "grey96", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    # Fondo blanco puro en toda la figura
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black", linewidth = 0.5)
  )

# ── 6. Guardar (PNG 600 dpi + TIFF 600 dpi + PDF vectorial) ──────────────────
ggsave("figures/16S/figura3_pcoa.png",  p, width = 6.5, height = 5, dpi = 600)
ggsave("figures/16S/figura3_pcoa.tiff", p, width = 6.5, height = 5, dpi = 600,
       compression = "lzw")
ggsave("figures/16S/figura3_pcoa.pdf",  p, width = 6.5, height = 5)
cat("Figura guardada: figura3_pcoa.png / figura3_pcoa.tiff / figura3_pcoa.pdf\n")
