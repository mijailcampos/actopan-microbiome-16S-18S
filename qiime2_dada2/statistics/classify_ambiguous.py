import re

import sys
# Usage: python classify_ambiguous.py [path/to/genus-table-rarefied-n17-export/feature-table.tsv]
path = sys.argv[1] if len(sys.argv) > 1 else "16S/reanalysis_fecal_n5/genus-table-rarefied-n17-export/feature-table.tsv"
with open(path) as f:
    lines = f.read().splitlines()

genera = []
for line in lines[2:]:
    taxon = line.split('\t')[0]
    if 'g__' in taxon:
        g = taxon.split('g__')[-1].strip()
        if g:
            genera.append(g)
genera = sorted(set(genera), key=lambda x: x.lower())

EXACT_EXCLUDE = {
    'uncultured', 'uncultured_bacterium', 'unidentified', 'unclassified', 'unknown',
    'ambiguous_taxa', 'incertae_sedis', 'metagenome', 'environmental_sample', 'sericytochromatia'
}
ORGANELLE = {'chloroplast', 'mitochondria'}
PHYLUM_PREFIXES = ('bacteroidetes_', 'proteobacteria_', 'firmicutes_', 'actinobacteria_',
                    'verrucomicrobia_', 'chloroflexi_', 'cyanobacteria_', 'planctomycetes_')

def is_ambiguous(g):
    low = g.lower()
    if low in EXACT_EXCLUDE or low in ORGANELLE:
        return 'placeholder/uncultured explicito'
    if low.endswith('aceae'):
        return 'rango familia (-aceae)'
    if low.endswith('ales'):
        return 'rango orden (-ales)'
    if '_group' in low or '_clade' in low or '_lineage' in low or 'subgroup' in low:
        return 'clado sin nombre de genero ("_group/_clade/_lineage/subgroup")'
    if 'ucg' in low or 'ucg-' in low:
        return 'codigo "UCG" (uncultured genus, SILVA)'
    if low.startswith(PHYLUM_PREFIXES):
        return 'prefijo de filo + codigo de clado'
    first_token = re.split(r'[_(]', g)[0]
    stripped = first_token.replace('-', '')
    if not stripped.isalpha() or len(stripped) < 4:
        return 'codigo alfanumerico de clon ambiental (no nombre latinizado)'
    return None

excluded = {}
kept = []
for g in genera:
    reason = is_ambiguous(g)
    if reason:
        excluded[g] = reason
    else:
        kept.append(g)

print(f"Total etiquetas: {len(genera)}")
print(f"Excluidas (no identificadas a nivel genero): {len(excluded)}")
print(f"Mantenidas como genero identificado: {len(kept)}")
print("\n== EXCLUIDAS (revisar) ==")
for g in sorted(excluded, key=lambda x: x.lower()):
    print(f"  {g:45s} -> {excluded[g]}")
