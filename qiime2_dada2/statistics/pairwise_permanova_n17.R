suppressMessages(library(vegan))

# Usage: Rscript pairwise_permanova_n17.R [path/to/braycurtis-n17-export/distance-matrix.tsv]
args <- commandArgs(trailingOnly = TRUE)
dm_path <- if (length(args) > 0) args[1] else "16S/reanalysis_fecal_n5/braycurtis-n17-export/distance-matrix.tsv"
dm <- read.table(dm_path,
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
dm <- as.matrix(dm)
ids <- rownames(dm)
sample_type <- ifelse(grepl("^Agua", ids), "Water",
               ifelse(grepl("^Crassostrea", ids), "Oyster",
               ifelse(grepl("^Heces", ids), "Feces", NA)))
names(sample_type) <- ids

pairwise <- function(g1, g2) {
  keep <- sample_type %in% c(g1, g2)
  sub_dm <- as.dist(dm[keep, keep])
  meta <- data.frame(sample_type = factor(sample_type[keep]))
  set.seed(123)
  ad <- adonis2(sub_dm ~ sample_type, data = meta, permutations = 999)
  n <- sum(keep)
  list(n = n, F = ad$F[1], R2 = ad$R2[1], p = ad$`Pr(>F)`[1], df = paste0(ad$Df[1], ",", ad$Df[2]))
}

pairs <- list(c("Feces","Oyster"), c("Feces","Water"), c("Oyster","Water"))
results <- lapply(pairs, function(p) pairwise(p[1], p[2]))
names(results) <- sapply(pairs, function(p) paste(p[1], "vs", p[2]))

pvals <- sapply(results, function(r) r$p)
qvals <- p.adjust(pvals, method = "BH")

cat("n por grupo global:\n"); print(table(sample_type))
cat("\n")
for (i in seq_along(results)) {
  r <- results[[i]]
  cat(sprintf("%-18s n=%d  df=%s  F=%.3f  R2=%.3f  p=%.3f  q=%.4f\n",
              names(results)[i], r$n, r$df, r$F, r$R2, r$p, qvals[i]))
}
