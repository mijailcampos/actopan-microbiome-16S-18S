import os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('DATA', os.path.join(HERE, 'data'))  # inputs, see README
OUT = os.environ.get('OUT', os.path.join(HERE, 'out'))
os.makedirs(OUT, exist_ok=True)
import sys
sys.path.insert(0, HERE)
from plot_pcoa import plot_pcoa

# Figure 3 (LAJAR): same data/script as Figura3_PCoA_16S_corregida.png (n=17, feces F1-F5), 600 dpi
base = os.path.join(DATA, '')
for ext in ('png',):
    plot_pcoa(
        base + 'braycurtis-pcoa-n17.qza',
        os.path.join(DATA, 'sample-metadata.tsv'),
        'PCoA (Bray-Curtis) - 16S rRNA microbiota, genus level - n=17 (feces F1-F5)',
        os.path.join(OUT, 'Figure3.png'),
        dpi=600,
    )
