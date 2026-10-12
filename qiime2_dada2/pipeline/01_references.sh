#!/usr/bin/env bash
# Build primer-specific naive Bayes classifiers for 16S (SILVA 138) and 18S (PR2).
# Commands and parameters reconstructed from the QIIME 2 provenance of the
# taxonomy.qza artifacts used in the manuscript (QIIME 2 2024.10.1).
set -euo pipefail

REFS=${REFS:-refs}
THREADS=${THREADS:-16}
mkdir -p "$REFS/16S" "$REFS/18S"

# ---------------------------------------------------------------------------
# 16S: SILVA 138 SSURef NR99 (pre-formatted RESCRIPt release distributed by
# QIIME 2 Data Resources: silva-138-99-seqs.qza / silva-138-99-tax.qza).
# Provenance of those files: get-silva-data (v138, SSURef_NR99, species labels)
# -> parse-silva-taxonomy -> cull-seqs (5 degenerates, homopolymer 8)
# -> filter-seqs-length-by-taxon (Archaea 900, Bacteria 1200, Eukaryota 1400)
# -> dereplicate (mode uniq).
# ---------------------------------------------------------------------------
qiime feature-classifier extract-reads \
  --i-sequences "$REFS/16S/silva-138-99-seqs.qza" \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-identity 0.8 --p-min-length 50 --p-max-length 0 \
  --p-read-orientation both \
  --p-n-jobs "$THREADS" \
  --o-reads "$REFS/16S/ref-seqs-v3v4.qza"

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "$REFS/16S/ref-seqs-v3v4.qza" \
  --i-reference-taxonomy "$REFS/16S/silva-138-99-tax.qza" \
  --o-classifier "$REFS/16S/classifier-16S-v3v4.qza"

# ---------------------------------------------------------------------------
# 18S: PR2 (8 ranks), imported from FASTA + QIIME-formatted taxonomy TSV.
# ---------------------------------------------------------------------------
qiime tools import --type 'FeatureData[Sequence]' \
  --input-path "$REFS/18S/pr2_seqs.fasta" \
  --output-path "$REFS/18S/pr2-seqs.qza"

qiime tools import --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat \
  --input-path "$REFS/18S/pr2_tax_q2.tsv" \
  --output-path "$REFS/18S/pr2-tax.qza"

qiime feature-classifier extract-reads \
  --i-sequences "$REFS/18S/pr2-seqs.qza" \
  --p-f-primer CYGCGGTAATTCCAGCTC \
  --p-r-primer AYGGTATCTRATCRTCTTYG \
  --p-identity 0.8 --p-min-length 50 --p-max-length 0 \
  --p-read-orientation both \
  --p-n-jobs "$THREADS" \
  --o-reads "$REFS/18S/ref-seqs-v4.qza"

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "$REFS/18S/ref-seqs-v4.qza" \
  --i-reference-taxonomy "$REFS/18S/pr2-tax.qza" \
  --o-classifier "$REFS/18S/classifier-18S-v4.qza"
