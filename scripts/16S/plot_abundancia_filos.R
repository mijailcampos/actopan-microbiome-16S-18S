library(tidyverse)

# Run from the repo root: Rscript scripts/16S/plot_abundancia_filos.R
tax <- read.delim(
  "example_data/16S/final.opti_mcc.0.03.cons.tax.summary",
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE
)

# R usa los números de fila como row.names automáticamente; no hay que eliminar columna

# Limpiar comillas en nombres de taxones
tax$taxon <- gsub('"', '', tax$taxon)

# ── 2. Definir columnas por grupo ─────────────────────────────────────────────
agua_cols <- c("AguaP1", "AguaP2", "AguaP3", "AguaP4", "AguaP5")

ostion_cols <- c(
  "CrassostreavirginicaG1",  "CrassostreavirginicaG2",
  "CrassostreavirginicaG3",  "CrassostreavirginicaG4",
  "CrassostreavirginicaG5",  "CrassostreavirginicaG6",
  "CrassostreavirginicaG7",  "CrassostreavirginicaG8",
  "CrassostreavirginicaG9",  "CrassostreavirginicaG10"
)

heces_cols <- c(
  "HecesP6",  "HecesP7",  "HecesP8",  "HecesP9",
  "HecesP10", "HecesP11", "HecesP12", "HecesP13", "HecesP14"
)

# ── 3. Filtrar nivel de filo (taxlevel == 2) ──────────────────────────────────
phyla <- tax %>% filter(taxlevel == 2)

# ── 4. Abundancia relativa por muestra → promedio por grupo ──────────────────
# Calcula % por columna (muestra) antes de promediar; así cada muestra pesa igual
rel_por_muestra <- function(df_filas, cols) {
  mat        <- as.matrix(df_filas[, cols])
  col_totals <- colSums(mat)
  col_totals[col_totals == 0] <- 1        # evitar división por cero
  sweep(mat, 2, col_totals, "/") * 100
}

agua_rel   <- rel_por_muestra(phyla, agua_cols)
ostion_rel <- rel_por_muestra(phyla, ostion_cols)
heces_rel  <- rel_por_muestra(phyla, heces_cols)

df_grupos <- tibble(
  taxon           = phyla$taxon,
  Agua            = rowMeans(agua_rel),
  Ostión          = rowMeans(ostion_rel),
  `Heces humanas` = rowMeans(heces_rel)
)

# ── 5. Top 7 filos por promedio global; resto → "Otros" ───────────────────────
n_top <- 7

df_grupos <- df_grupos %>%
  mutate(media_global = (Agua + Ostión + `Heces humanas`) / 3) %>%
  arrange(desc(media_global))

# Excluir filos "unclassified" del top (no son filos reales)
df_clasificados <- df_grupos %>%
  filter(!grepl("unclassified|_unclassified", taxon, ignore.case = TRUE))

top_phyla <- df_clasificados$taxon[seq_len(min(n_top, nrow(df_clasificados)))]

df_plot <- df_grupos %>%
  mutate(filo = if_else(taxon %in% top_phyla, taxon, "Otros")) %>%
  group_by(filo) %>%
  summarise(
    Agua            = sum(Agua),
    Ostión          = sum(Ostión),
    `Heces humanas` = sum(`Heces humanas`),
    .groups         = "drop"
  )

# ── 6. Formato largo ──────────────────────────────────────────────────────────
df_long <- df_plot %>%
  pivot_longer(
    cols      = c(Agua, Ostión, `Heces humanas`),
    names_to  = "Grupo",
    values_to = "Abundancia"
  )

# Orden del eje X (historia ambiental → hospedero → humano)
df_long$Grupo <- factor(df_long$Grupo,
                         levels = c("Agua", "Ostión", "Heces humanas"))

# Orden de apilado: más abundante abajo, Otros arriba
# top_phyla ya está ordenado de mayor a menor por media_global
filo_order        <- c(top_phyla, "Otros")
df_long$filo      <- factor(df_long$filo, levels = filo_order)

# ── 7. Paleta de colores ──────────────────────────────────────────────────────
# Colores específicos para los filos más comunes en 16S
colores_filos <- c(
  "Proteobacteria"    = "#2166AC",
  "Pseudomonadota"    = "#2166AC",
  "Firmicutes"        = "#4DAC26",
  "Bacillota"         = "#4DAC26",
  "Bacteroidota"      = "#F46D43",
  "Bacteroidetes"     = "#F46D43",
  "Actinobacteriota"  = "#762A83",
  "Actinobacteria"    = "#762A83",
  "Fusobacteriota"    = "#D6604D",
  "Fusobacteria"      = "#D6604D",
  "Cyanobacteria"     = "#74ADD1",
  "Verrucomicrobiota" = "#A6D96A",
  "Cyanobacteria"              = "#74ADD1",
  "Cyanobacteria_Chloroplast"  = "#74ADD1",
  "Verrucomicrobia"            = "#A6D96A",
  "Verrucomicrobiota"          = "#A6D96A",
  "Campilobacterota"           = "#E6AB02",
  "Spirochaetota"              = "#E78AC3",
  "Planctomycetota"            = "#66C2A5",
  "Chloroflexi"                = "#FC8D62",
  "Otros"                      = "#BDBDBD"
)

# Colores extra para filos no reconocidos
extra_pool <- c("#E66101", "#998EC3", "#FDB863", "#B2ABD2",
                "#5E4FA2", "#3288BD", "#ABDDA4", "#FDAE61")

palette_final <- character(0)
extra_i       <- 1L
for (f in filo_order) {
  if (f %in% names(colores_filos)) {
    palette_final[f] <- colores_filos[f]
  } else {
    palette_final[f] <- extra_pool[extra_i]
    extra_i          <- extra_i + 1L
  }
}

# ── 8. Gráfica ────────────────────────────────────────────────────────────────
p <- ggplot(df_long, aes(x = Grupo, y = Abundancia, fill = filo)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = palette_final,
    name   = "Filo",
    # Leyenda: orden de arriba a abajo igual que la barra de abajo a arriba
    guide  = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(
    limits = c(0, 100.01),   # .01 para que 100 no quede cortado
    expand = c(0, 0),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    x = NULL,
    y = "Abundancia relativa (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x        = element_text(size = 13, color = "black"),
    axis.text.y        = element_text(size = 11, color = "black"),
    axis.title.y       = element_text(size = 12),
    axis.line          = element_line(color = "black"),
    legend.position    = "right",
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.5, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA)
  )

# ── 9. Guardar ────────────────────────────────────────────────────────────────
ggsave("figures/16S/abundancia_relativa_filos.pdf", p, width = 7, height = 6, dpi = 300)
ggsave("figures/16S/abundancia_relativa_filos.png", p, width = 7, height = 6, dpi = 300)

cat("\nTop", n_top, "filos identificados:\n")
cat(paste0("  ", seq_along(top_phyla), ". ", top_phyla, "\n"), sep = "")
cat("\nFigura guardada como:\n")
cat("  abundancia_relativa_filos.pdf\n")
cat("  abundancia_relativa_filos.png\n")
