# Microbiome Analysis — Actopan River Oyster, Water and Human Feces (16S & 18S rRNA)

Characterization of prokaryotic (16S rRNA) and eukaryotic (18S rRNA) microbial communities in *Crassostrea virginica* (Eastern oyster), river water, and human fecal samples collected at Chachalacas, Veracruz, Mexico.

> **Manuscript analyses:** the results reported in the manuscript were produced with the QIIME 2 / DADA2 (ASV) workflow in [`qiime2_dada2/`](qiime2_dada2/). The mothur / OptiClust workflow described below is the earlier version of the analysis and is kept for reference.

## Study Overview

| Feature | 16S rRNA | 18S rRNA |
|---------|----------|----------|
| Target community | Prokaryotes (bacteria) | Eukaryotes (protists, fungi) |
| Reference database | SILVA v138 | PR2 v5.1.0 |
| Clustering | OptiClust (97% similarity) | OptiClust (97% similarity) |
| Sample types | Water (WA), Oyster (CV), Feces (HE) | Water (WA), Oyster (CV), Feces (HE) |

## Raw Sequencing Data

Raw FASTQ files are deposited at NCBI under:

**BioProject: [PRJNA1471323](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1471323)**

## Repository Structure

```
actopan-microbiome-16S-18S/
├── notebooks/
│   ├── 16S_analisis_microbiota.ipynb   # Prokaryote community analysis & figures
│   └── 18S_analisis_microbiota.ipynb   # Eukaryote community analysis & figures
├── scripts/
│   ├── 16S/                            # Python and R scripts for 16S analysis
│   │   ├── generar_tabla_diversidad_alfa.py   # Alpha diversity indices
│   │   ├── fig_rarefaction.py                 # Rarefaction curves
│   │   ├── plot_pcoa_figura3.R                # PCoA Bray-Curtis (vegan)
│   │   ├── plot_abundancia_filos.R            # Phylum barplot (tidyverse)
│   │   ├── fig_filos_combinados.py            # Phylum composition figure
│   │   ├── fig_filos_generos_combinados.py    # Genus composition figure
│   │   └── heatmap_sanitario_integrado.py     # Sanitary taxa heatmap
│   └── 18S/                            # Python and R scripts for 18S analysis
│       ├── fig_pcoa_18S.py                    # PCoA Bray-Curtis (scipy)
│       └── fig_filos_por_ambiente.R           # Phylum by environment (vegan)
├── statistics/
│   ├── 16S/                            # Differential abundance and community comparison results
│   └── 18S/                            # Differential abundance and community comparison results
├── figures/
│   ├── 16S/                            # Publication figures (PNG, 300–600 dpi)
│   └── 18S/                            # Publication figures (PNG, 300–600 dpi)
├── metadata/
│   ├── metadatos_16S.txt               # Sample IDs and group assignments
│   └── metadatos_18S.txt
├── example_data/
│   ├── 16S/                            # Mothur processed outputs + alpha diversity tables
│   └── 18S/                            # Mothur processed outputs
├── requirements.txt
└── .gitignore
```

## How to Reproduce the Analyses

### 1. Clone the repository

```bash
git clone https://github.com/mijailcampos/actopan-microbiome-16S-18S.git
cd actopan-microbiome-16S-18S
```

### 2. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 3. Install R dependencies

```r
install.packages(c("tidyverse", "vegan", "ggplot2", "scales"))
```

### 4. Run the notebooks

```bash
jupyter notebook notebooks/16S_analisis_microbiota.ipynb
jupyter notebook notebooks/18S_analisis_microbiota.ipynb
```

The notebooks read input data from `statistics/` and `example_data/` using relative paths. Run from the repo root directory.

### 5. Run individual scripts

All scripts must be run from the **repo root**:

```bash
# Python
python scripts/16S/generar_tabla_diversidad_alfa.py
python scripts/16S/fig_rarefaction.py
python scripts/18S/fig_pcoa_18S.py

# R
Rscript scripts/16S/plot_pcoa_figura3.R
Rscript scripts/16S/plot_abundancia_filos.R
Rscript scripts/18S/fig_filos_por_ambiente.R
```

## Statistical Analyses

Pairwise community comparisons between sample groups (CV, WA, HE) were performed at phylum and genus levels using Welch's t-test and White's non-parametric t-test with Benjamini–Hochberg FDR correction. Results are provided as tab-separated tables in `statistics/16S/` and `statistics/18S/`.

Beta diversity was assessed using Bray–Curtis dissimilarity with PCoA ordination. PERMANOVA results for 18S are in `statistics/18S/PERMANOVA_18S.tsv`.

Alpha diversity indices (Chao1, Shannon, Inverse Simpson) were calculated in Mothur after rarefaction.

## Software Versions

| Software | Version |
|----------|---------|
| Mothur | 1.48.0 |
| Python | 3.10 |
| pandas | ≥ 1.5 |
| matplotlib | ≥ 3.6 |
| numpy | ≥ 1.23 |
| scipy | ≥ 1.9 |
| R | 4.x |
| vegan | 2.6 |
| tidyverse | 2.x |

## Sample Information

- **WA** — River water samples, Actopan River, Chachalacas, Veracruz, Mexico (n = 5)
- **CV** — *Crassostrea virginica* tissue samples (n = 10)
- **HE** — Human fecal samples (n = 9)

Sample IDs in metadata files are anonymous laboratory codes. No personal or identifying information is included.

## Citation

> [Citation will be added upon publication]

## License

See [LICENSE](LICENSE) for details.
