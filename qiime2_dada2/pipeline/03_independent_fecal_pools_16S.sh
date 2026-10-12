#!/usr/bin/env bash
# 16S community-level comparisons with independent fecal pools only (F1-F5).
# Fecal pools F6-F9 partition the same 20 participants by demographic group,
# so they are not independent of F1-F5 and are excluded to avoid
# pseudoreplication (Water n = 5, Oyster n = 7, Feces n = 5).
# Reconstructed from the provenance of braycurtis-pcoa-n17.qza.
set -euo pipefail
cd 16S
mkdir -p reanalysis_fecal_n5

qiime feature-table filter-samples \
  --i-table genus-table-rarefied.qza \
  --m-metadata-file ../metadata/exclude-heces-demografico.tsv \
  --p-exclude-ids \
  --o-filtered-table reanalysis_fecal_n5/genus-table-rarefied-n17.qza

qiime diversity beta \
  --i-table reanalysis_fecal_n5/genus-table-rarefied-n17.qza \
  --p-metric braycurtis \
  --o-distance-matrix reanalysis_fecal_n5/braycurtis-n17.qza

qiime diversity pcoa \
  --i-distance-matrix reanalysis_fecal_n5/braycurtis-n17.qza \
  --o-pcoa reanalysis_fecal_n5/braycurtis-pcoa-n17.qza

for a in genus-table-rarefied-n17 braycurtis-n17; do
  qiime tools export --input-path "reanalysis_fecal_n5/$a.qza" \
    --output-path "reanalysis_fecal_n5/$a-export"
done
