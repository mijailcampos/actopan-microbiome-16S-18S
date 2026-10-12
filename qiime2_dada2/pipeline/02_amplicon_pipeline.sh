#!/usr/bin/env bash
# QIIME 2 / DADA2 processing of one marker (16S or 18S).
# Commands and parameters reconstructed from the provenance of the artifacts
# used in the manuscript (QIIME 2 2024.10.1, q2-dada2 2024.10.0, DADA2 1.30.0,
# cutadapt 4.9).
#
# Usage:  bash 02_amplicon_pipeline.sh 16S
#         bash 02_amplicon_pipeline.sh 18S
#
# Expects, inside the marker directory:
#   manifest.tsv         PairedEndFastqManifestPhred33V2 (FASTQ from PRJNA1471323)
#   sample-metadata.tsv  see ../metadata/
set -euo pipefail

MARKER=${1:?usage: $0 16S|18S}
REFS=${REFS:-$(pwd)/refs}
THREADS=${THREADS:-32}
cd "$MARKER"

case "$MARKER" in
  16S)
    FWD=CCTACGGGNGGCWGCAG;  REV=GACTACHVGGGTATCTAATCC   # 341F / 805R (V3-V4)
    TRUNC_F=260; TRUNC_R=200
    CLASSIFIER=$REFS/16S/classifier-16S-v3v4.qza
    DEPTH=5421
    PHYLUM_LEVEL=2; GENUS_LEVEL=6                       # SILVA ranks
    ;;
  18S)
    FWD=CYGCGGTAATTCCAGCTC; REV=AYGGTATCTRATCRTCTTYG    # V4
    TRUNC_F=270; TRUNC_R=230
    CLASSIFIER=$REFS/18S/classifier-18S-v4.qza
    DEPTH=986                                           # host-filtered table
    PHYLUM_LEVEL=3; GENUS_LEVEL=7                       # PR2 ranks
    ;;
  *) echo "unknown marker: $MARKER" >&2; exit 1 ;;
esac

# 1. Import
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path manifest.tsv \
  --output-path demux.qza

# 2. Primer removal (reads without primer are discarded)
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f "$FWD" --p-front-r "$REV" \
  --p-error-rate 0.1 --p-overlap 3 --p-discard-untrimmed \
  --o-trimmed-sequences demux-trimmed.qza

# 3. ASV inference
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-trimmed.qza \
  --p-trunc-len-f "$TRUNC_F" --p-trunc-len-r "$TRUNC_R" \
  --p-max-ee-f 2 --p-max-ee-r 2 --p-trunc-q 2 --p-min-overlap 12 \
  --p-pooling-method independent --p-chimera-method consensus \
  --p-n-threads "$THREADS" \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

# 4. Taxonomy
qiime feature-classifier classify-sklearn \
  --i-reads rep-seqs.qza \
  --i-classifier "$CLASSIFIER" \
  --p-confidence 0.7 --p-n-jobs 16 \
  --o-classification taxonomy.qza

# 5. 18S only: remove host (oyster) reads before any downstream analysis
TABLE=table.qza; SUFFIX=""
if [ "$MARKER" = 18S ]; then
  qiime taxa filter-table \
    --i-table table.qza --i-taxonomy taxonomy.qza \
    --p-exclude Mollusca --p-mode contains \
    --o-filtered-table table-nohost.qza
  TABLE=table-nohost.qza; SUFFIX="-nohost"
fi

# 6. Alpha and beta diversity at a fixed depth (non-phylogenetic core metrics)
qiime diversity core-metrics \
  --i-table "$TABLE" \
  --p-sampling-depth "$DEPTH" \
  --m-metadata-file sample-metadata.tsv \
  --output-dir "core-metrics-results$SUFFIX"
# Inverse Simpson reported in the manuscript = 1 / dominance_vector.

# 7. Rarefaction curves (10 steps up to 20,000 reads)
qiime diversity alpha-rarefaction \
  --i-table "$TABLE" \
  --p-max-depth 20000 --p-steps 10 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization "alpha-rarefaction$SUFFIX.qzv"

# 8. Taxonomic composition (phylum and genus)
qiime taxa collapse --i-table "$TABLE" --i-taxonomy taxonomy.qza \
  --p-level "$PHYLUM_LEVEL" --o-collapsed-table "phylum-table$SUFFIX.qza"
qiime feature-table relative-frequency --i-table "phylum-table$SUFFIX.qza" \
  --o-relative-frequency-table "phylum-table$SUFFIX-relfreq.qza"

qiime taxa collapse --i-table "$TABLE" --i-taxonomy taxonomy.qza \
  --p-level "$GENUS_LEVEL" --o-collapsed-table "genus-table$SUFFIX.qza"
qiime feature-table relative-frequency --i-table "genus-table$SUFFIX.qza" \
  --o-relative-frequency-table "genus-table$SUFFIX-relfreq.qza"

# 9. Genus-level Bray-Curtis on the rarefied genus table
qiime feature-table rarefy --i-table "genus-table$SUFFIX.qza" \
  --p-sampling-depth "$DEPTH" \
  --o-rarefied-table "genus-table$SUFFIX-rarefied.qza"
qiime diversity beta --i-table "genus-table$SUFFIX-rarefied.qza" \
  --p-metric braycurtis \
  --o-distance-matrix "genus-braycurtis-distance$SUFFIX.qza"
qiime diversity pcoa --i-distance-matrix "genus-braycurtis-distance$SUFFIX.qza" \
  --o-pcoa "genus-braycurtis-pcoa$SUFFIX.qza"
