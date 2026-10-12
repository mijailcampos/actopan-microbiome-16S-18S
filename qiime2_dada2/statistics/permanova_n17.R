suppressMessages(library(vegan))

# Usage: Rscript permanova_n17.R [path/to/braycurtis-n17-export/distance-matrix.tsv]
args <- commandArgs(trailingOnly = TRUE)
dm_path <- if (length(args) > 0) args[1] else "16S/reanalysis_fecal_n5/braycurtis-n17-export/distance-matrix.tsv"
dm <- read.table(dm_path,
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
dm <- as.dist(as.matrix(dm))

ids <- labels(dm)
sample_type <- ifelse(grepl("^Agua", ids), "Water",
               ifelse(grepl("^Crassostrea", ids), "Oyster",
               ifelse(grepl("^Heces", ids), "Feces", NA)))
meta <- data.frame(sample_type = factor(sample_type), row.names = ids)

cat("n por grupo:\n"); print(table(meta$sample_type))

set.seed(123)
ad <- adonis2(dm ~ sample_type, data = meta, permutations = 999)
cat("\n=== PERMANOVA (adonis2, n=17, corregido) ===\n")
print(ad)

bd <- betadisper(dm, meta$sample_type)
set.seed(123)
pt <- permutest(bd, permutations = 999)
cat("\n=== PERMDISP (betadisper, n=17, corregido) ===\n")
print(pt)
cat("\nDistancias promedio al centroide por grupo:\n")
print(tapply(bd$distances, meta$sample_type, mean))
