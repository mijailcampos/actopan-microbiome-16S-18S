library(tidyverse)
library(vegan)
library(ggplot2)
library(scales)

# Run from the repo root: Rscript scripts/18S/fig_filos_por_ambiente.R
tax_file <- "example_data/18S/final.opti_mcc.0.03.cons.tax.summary"
out_dir  <- "figures/18S"

# ── 1. Leer tabla de taxonomía Mothur ─────────────────────────────────────────
tax_raw <- read.table(tax_file, header = TRUE, sep = "\t",
                      check.names = FALSE, stringsAsFactors = FALSE)

# ── 2. Filtrar nivel phylum (taxlevel == 2) ────────────────────────────────────
phylum_raw <- tax_raw %>% filter(taxlevel == 2)

# Columnas de muestras (todo después de "total")
sample_cols <- c(
  "AguaP1", "AguaP2", "AguaP3", "AguaP4", "AguaP5",
  "CrassostreavirginicaG2", "CrassostreavirginicaG3", "CrassostreavirginicaG4",
  "CrassostreavirginicaG5", "CrassostreavirginicaG6", "CrassostreavirginicaG7",
  "CrassostreavirginicaG8",
  "HecesP10", "HecesP11", "HecesP12", "HecesP13", "HecesP14",
  "HecesP6",  "HecesP7",  "HecesP8",  "HecesP9"
)

water_cols <- c("AguaP1", "AguaP2", "AguaP3", "AguaP4", "AguaP5")
cv_cols    <- c("CrassostreavirginicaG2", "CrassostreavirginicaG3",
                "CrassostreavirginicaG4", "CrassostreavirginicaG5",
                "CrassostreavirginicaG6", "CrassostreavirginicaG7",
                "CrassostreavirginicaG8")
he_cols    <- c("HecesP10", "HecesP11", "HecesP12", "HecesP13", "HecesP14",
                "HecesP6",  "HecesP7",  "HecesP8",  "HecesP9")

# ── 3. Abundancia relativa por muestra con vegan::decostand ───────────────────
counts_mat <- as.matrix(phylum_raw[, sample_cols])
rownames(counts_mat) <- phylum_raw$taxon

# decostand "total" normaliza cada columna (muestra) a proporciones
rel_mat <- decostand(t(counts_mat), method = "total") * 100  # muestras × filos
rel_mat <- t(rel_mat)                                        # filos × muestras

# ── 4. Limpiar y fusionar nombres de taxones ──────────────────────────────────
# Renombrar en el data.frame (no en el matrix) para evitar duplicados en rownames
rel_df <- as.data.frame(rel_mat) %>%
  rownames_to_column("taxon") %>%
  mutate(taxon = case_when(
    taxon == "Eukaryota_unclassified" ~ "Sin clasificar",
    taxon == "Cryptista:nucl"         ~ "Cryptista",
    TRUE                              ~ taxon
  )) %>%
  group_by(taxon) %>%
  summarise(across(everything(), sum), .groups = "drop")

# ── 5. Promedio por ambiente ───────────────────────────────────────────────────
rel_df <- rel_df %>%
  mutate(
    Agua              = rowMeans(across(all_of(water_cols))),
    `C. virginica`    = rowMeans(across(all_of(cv_cols))),
    `Heces humanas`   = rowMeans(across(all_of(he_cols))),
    total_mean        = (Agua + `C. virginica` + `Heces humanas`) / 3
  )

# ── 6. Seleccionar top 6 filos + "Otros" ──────────────────────────────────────
n_top <- 6
top_taxa <- rel_df %>%
  arrange(desc(total_mean)) %>%
  slice_head(n = n_top) %>%
  pull(taxon)

plot_df <- rel_df %>%
  mutate(taxon_plot = if_else(taxon %in% top_taxa, taxon, "Otros")) %>%
  group_by(taxon_plot) %>%
  summarise(
    Agua            = sum(Agua),
    `C. virginica`  = sum(`C. virginica`),
    `Heces humanas` = sum(`Heces humanas`),
    .groups = "drop"
  ) %>%
  mutate(total_mean = (Agua + `C. virginica` + `Heces humanas`) / 3)

# Orden: top phyla desc, "Otros" siempre al final (base de la barra)
ordered_taxa <- plot_df %>%
  filter(taxon_plot != "Otros") %>%
  arrange(desc(total_mean)) %>%
  pull(taxon_plot)
ordered_taxa <- c("Otros", ordered_taxa)   # Otros en base, más abundante arriba

# ── 7. Formato largo para ggplot ───────────────────────────────────────────────
long_df <- plot_df %>%
  select(taxon_plot, Agua, `C. virginica`, `Heces humanas`) %>%
  pivot_longer(
    cols      = c("Agua", "C. virginica", "Heces humanas"),
    names_to  = "Ambiente",
    values_to = "Abundancia"
  ) %>%
  mutate(
    taxon_plot = factor(taxon_plot, levels = ordered_taxa),
    Ambiente   = factor(Ambiente,
                        levels = c("Agua", "C. virginica", "Heces humanas"))
  )

# ── 8. Paleta de colores contrastantes ────────────────────────────────────────
pal_all <- c(
  "TSAR"           = "#2C4A8C",   # azul marino
  "Archaeplastida" = "#56B4E9",   # azul claro
  "Cryptista"      = "#3CBFB5",   # verde azulado
  "Obazoa"         = "#1D6B3D",   # verde oscuro
  "Amoebozoa"      = "#D4820A",   # ámbar/naranja
  "Haptista"       = "#F5C518",   # amarillo
  "Sin clasificar" = "#9E9E9E",   # gris medio
  "Otros"          = "#D3D3D3"    # gris claro
)

colors_used <- pal_all[levels(long_df$taxon_plot)]

# Etiqueta legenda ordenada de más a menos abundante (Otros al final)
legend_labels <- rev(ordered_taxa)   # el geom_bar apila de abajo arriba

# ── 9. Etiquetas eje X (cursiva para Crassostrea) ─────────────────────────────
x_labels <- c(
  "Agua",
  expression(italic("Crassostrea virginica")),
  "Heces humanas"
)

# ── 10. Figura ────────────────────────────────────────────────────────────────
p <- ggplot(long_df,
            aes(x = Ambiente, y = Abundancia, fill = taxon_plot)) +
  geom_bar(stat = "identity", width = 0.65,
           color = "white", linewidth = 0.25) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 101),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(
    values = colors_used,
    breaks = c(ordered_taxa[-1], "Otros"),  # más abundante primero, Otros al final
    name   = "Filo"
  ) +
  labs(
    x = "Ambiente",
    y = "Abundancia relativa (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # Fondo completamente blanco
    panel.background   = element_rect(fill = "white"),
    plot.background    = element_rect(fill = "white", color = NA),
    # Cuadrícula mínima solo horizontal
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.7, fill = NA),
    # Ejes — tamaño legible a formato revista
    axis.text.x        = element_text(size = 12, color = "black", vjust = 1),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.title.x       = element_text(size = 14, face = "bold", margin = margin(t = 9)),
    axis.title.y       = element_text(size = 14, face = "bold", margin = margin(r = 9)),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    # Leyenda
    legend.title       = element_text(size = 13, face = "bold"),
    legend.text        = element_text(size = 12),
    legend.key.size    = unit(0.6, "cm"),
    legend.spacing.y   = unit(0.15, "cm"),
    legend.position    = "right",
    legend.background  = element_rect(fill = "white", color = NA)
  )

print(p)

# ── 11. Exportar alta resolución ──────────────────────────────────────────────
ggsave(file.path(out_dir, "fig18S_filos_ambiente_resumen.pdf"),  p, width = 8, height = 5.5, dpi = 600)
ggsave(file.path(out_dir, "fig18S_filos_ambiente_resumen.tiff"), p, width = 8, height = 5.5, dpi = 600, compression = "lzw")
ggsave(file.path(out_dir, "fig18S_filos_ambiente_resumen.png"),  p, width = 8, height = 5.5, dpi = 600)

message("✓ Figura exportada en figures/18S/")

# ── 12. Tabla de medias para verificación ─────────────────────────────────────
cat("\nAbundancia media (%) por ambiente y filo:\n")
print(
  plot_df %>%
    filter(taxon_plot != "Otros") %>%
    arrange(desc(total_mean)) %>%
    mutate(across(where(is.numeric), \(x) round(x, 2))) %>%
    select(taxon_plot, Agua, `C. virginica`, `Heces humanas`, total_mean),
  n = Inf
)
