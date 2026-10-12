# QIIME 2 / DADA2 analysis (manuscript version)

This folder contains the code used for the analyses reported in the manuscript
*Relationship between the microbiota composition of the American oyster
(Crassostrea virginica) intestine, river water, and human feces in Chachalacas,
Veracruz, Mexico*. It supersedes the earlier mothur/OptiClust workflow kept in the
rest of this repository.

The QIIME 2 commands were reconstructed from the provenance recorded inside the
artifacts used for the manuscript, so every parameter below is the one that was
actually run. No sequence data, feature tables or sample coordinates are
distributed here; raw reads are available at NCBI BioProject
[PRJNA1471323](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1471323).

## Software

| Tool | Version |
|------|---------|
| QIIME 2 | 2024.10.1 |
| q2-dada2 / DADA2 | 2024.10.0 / 1.30.0 |
| cutadapt | 4.9 |
| R / vegan | 4.3 / 2.7-1 (identical results with 2.6-8) |
| Python | 3.10+ (numpy, pandas, scipy, statsmodels, biom-format, matplotlib) |

## Layout

```
qiime2_dada2/
├── pipeline/
│   ├── 01_references.sh                    # SILVA 138 (16S) and PR2 (18S) classifiers, primer-trimmed
│   ├── 02_amplicon_pipeline.sh             # import → cutadapt → DADA2 → taxonomy → diversity → composition
│   └── 03_independent_fecal_pools_16S.sh   # 16S beta diversity with independent fecal pools (F1–F5)
├── statistics/
│   ├── permanova_n17.R                     # PERMANOVA (adonis2) + PERMDISP, 999 permutations
│   ├── pairwise_permanova_n17.R            # pairwise PERMANOVA, Benjamini–Hochberg
│   ├── kruskal_wallis_n17.py               # Kruskal–Wallis per taxon (phylum, genus), BH-FDR
│   ├── shared_genera.py                    # genera present in ≥50 % of samples of every group
│   └── classify_ambiguous.py               # genus labels not resolved to a named genus
└── figures/
    ├── make_figure1_map.py                 # Figure 1, sampling map
    ├── plot_fig2_n5.py                     # Figure 2, 16S composition
    ├── plot_fig3_pcoa_600dpi.py            # Figure 3, 16S PCoA (Bray–Curtis, genus level)
    ├── plot_fig4_18S_v2.py                 # Figure 4, 18S composition (host-filtered)
    ├── plot_fig5_heatmap_600dpi.py         # Figure 5, taxa of sanitary relevance
    └── arial.py, plot_pcoa.py, plot_taxa_barplot.py   # shared helpers
```

## Key parameters

| Step | 16S rRNA (V3–V4) | 18S rRNA (V4) |
|------|------------------|---------------|
| Primers | 341F CCTACGGGNGGCWGCAG / 805R GACTACHVGGGTATCTAATCC | CYGCGGTAATTCCAGCTC / AYGGTATCTRATCRTCTTYG |
| cutadapt | error rate 0.1, overlap 3, discard untrimmed | same |
| DADA2 truncation (F/R) | 260 / 200 | 270 / 230 |
| DADA2 other | maxEE 2/2, truncQ 2, min overlap 12, consensus chimera removal | same |
| Reference | SILVA 138 SSURef NR99 | PR2 (8 ranks) |
| Classifier | naive Bayes, primer-extracted reads, confidence 0.7 | same |
| Host filtering | — | Mollusca removed (`taxa filter-table`) |
| Rarefaction depth | 5,421 | 986 (host-filtered) |
| Collapse level (phylum / genus) | 2 / 6 | 3 / 7 |

Alpha diversity: Chao1, Shannon and inverse Simpson (1 / `dominance_vector`).
Beta diversity: Bray–Curtis on genus-level tables rarefied to the depths above.

Fecal pools F6–F9 partition the same 20 participants by demographic group, so
they are not independent of F1–F5. Community-level comparisons use F1–F5 only
(Water n = 5, Oyster n = 7, Feces n = 5).

## Running

```bash
# 1. classifiers (needs silva-138-99-seqs.qza / silva-138-99-tax.qza and PR2 FASTA + taxonomy)
bash pipeline/01_references.sh

# 2. per marker; each marker directory needs manifest.tsv and sample-metadata.tsv
bash pipeline/02_amplicon_pipeline.sh 16S
bash pipeline/02_amplicon_pipeline.sh 18S
bash pipeline/03_independent_fecal_pools_16S.sh

# 3. statistics
Rscript statistics/permanova_n17.R 16S/reanalysis_fecal_n5/braycurtis-n17-export/distance-matrix.tsv
Rscript statistics/pairwise_permanova_n17.R 16S/reanalysis_fecal_n5/braycurtis-n17-export/distance-matrix.tsv
BASE=. python statistics/kruskal_wallis_n17.py

# 4. figures: inputs are read from $DATA (default figures/data/), output to $OUT
#    Arial TTF files are expected in $ARIAL_DIR (default figures/font/ttf/)
DATA=path/to/inputs OUT=figures/out python figures/plot_fig2_n5.py
```

Figure inputs expected in `$DATA`:

| Figure | Files |
|--------|-------|
| 2 | `phylum-table-relfreq.qza`, `genus-table-relfreq.qza`, `sample-metadata-n20.tsv` (16S) |
| 3 | `braycurtis-pcoa-n17.qza`, `sample-metadata.tsv` (16S) |
| 4 | `table.tsv` (`genus-table-nohost.qza` exported to TSV with taxonomy, 18S) |
| 5 | `Kruskal_Wallis_all_taxa_n17.csv` (output of `kruskal_wallis_n17.py`) |

Figure 1 reads `sites.tsv` and basemap files from its own folder; these are not
distributed. Running the scripts on the original inputs reproduces Figures 1–5
of the manuscript pixel for pixel.
