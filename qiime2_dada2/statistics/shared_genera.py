import csv, collections

import sys
# Usage: python shared_genera.py [path/to/genus-table-rarefied-n17-export/feature-table.tsv]
path = sys.argv[1] if len(sys.argv) > 1 else "16S/reanalysis_fecal_n5/genus-table-rarefied-n17-export/feature-table.tsv"

with open(path) as f:
    lines = f.read().splitlines()

header = lines[1].split('\t')  # skip "# Constructed from biom file"
samples = header[1:]

def group_of(s):
    if s.startswith('Agua'): return 'Water'
    if s.startswith('Crassostreavirginica'): return 'Oyster'
    if s.startswith('Heces'): return 'Feces'
    return 'Unknown'

groups = [group_of(s) for s in samples]
group_counts = collections.Counter(groups)
print("Muestras por grupo:", dict(group_counts))

thresholds = {'Water': 3, 'Oyster': 4, 'Feces': 3}  # >=50% redondeado arriba, de 5,7,5

genus_presence = collections.defaultdict(lambda: collections.Counter())

for line in lines[2:]:
    parts = line.split('\t')
    taxon = parts[0]
    counts = [float(x) for x in parts[1:]]
    if 'g__' in taxon:
        genus = taxon.split('g__')[-1].strip()
        if genus == '' or genus.startswith('uncultured') is False:
            pass
    else:
        genus = None
    if not genus:
        continue
    for grp, c in zip(groups, counts):
        if c > 0:
            genus_presence[genus][grp] += 1

shared = []
not_shared_sanitary = []
sanitary = ['Escherichia-Shigella', 'Aeromonas', 'Vibrio', 'Arcobacter', 'Enterococcus']

for genus, pres in genus_presence.items():
    ok = all(pres.get(g, 0) >= thresholds[g] for g in thresholds)
    if ok:
        shared.append(genus)

print(f"\nTotal generos con al menos 1 lectura en algun lado: {len(genus_presence)}")
print(f"Total generos SHARED bajo el criterio >=50% en los 3 grupos: {len(shared)}")

print("\n== Taxones sanitarios: cumplen o no el criterio ==")
for g in sanitary:
    pres = genus_presence.get(g, collections.Counter())
    w, o, f = pres.get('Water',0), pres.get('Oyster',0), pres.get('Feces',0)
    ok = w>=thresholds['Water'] and o>=thresholds['Oyster'] and f>=thresholds['Feces']
    print(f"  {g:22s} Water {w}/5  Oyster {o}/7  Feces {f}/5   -> {'SHARED' if ok else 'no shared'}")
