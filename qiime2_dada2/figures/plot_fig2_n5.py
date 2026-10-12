import os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('DATA', os.path.join(HERE, 'data'))  # inputs, see README
OUT = os.environ.get('OUT', os.path.join(HERE, 'out'))
os.makedirs(OUT, exist_ok=True)
# Figura 2 corregida: composicion 16S con heces solo F1-F5 (P6-P10), ostion 10 pools (relfreq, no rarificado)
import sys, re
sys.path.insert(0, HERE)
import plot_taxa_barplot as ptb

EXACT = {"uncultured","uncultured_bacterium","unidentified","unclassified","unknown","ambiguous_taxa",
         "incertae_sedis","metagenome","environmental_sample"}
RANK = {"f":"fam.","o":"ord.","c":"class","p":"phylum"}

def ok(name):
    low = name.lower()
    return bool(name) and low not in EXACT and "uncultured" not in low

def label(taxon, rank_prefix):
    parts = [p.strip() for p in taxon.split(";")]
    d = {}
    for p in parts:
        if "__" in p:
            k, v = p.split("__", 1)
            d[k] = v.strip()
    if rank_prefix == "g__":
        g = d.get("g", "")
        if ok(g) and not g.lower().endswith(("aceae", "ales")):
            return g
        for k in ("f", "o", "c", "p"):
            v = d.get(k, "")
            if ok(v):
                return f"{v} ({RANK[k]})"
        return "Unassigned"
    return ptb.short_name(taxon, rank_prefix)

ptb.short_name_orig = ptb.short_name
_orig = ptb.short_name
def short_name(taxon, rank_prefix):
    return label(taxon, rank_prefix) if rank_prefix == "g__" else _orig(taxon, rank_prefix)
ptb.short_name = short_name

D = DATA
ptb.make_figure("16S rRNA (prokaryotic)",
    f"{D}/phylum-table-relfreq.qza", f"{D}/genus-table-relfreq.qza",
    f"{D}/sample-metadata-n20.tsv",
    os.path.join(OUT, 'Figure2.png'),
    rank_prefix_phylum="p__", rank_prefix_genus="g__",
    sanitary_genera={"Vibrio","Escherichia-Shigella","Aeromonas","Arcobacter","Enterococcus"})

# Verificacion: valores de heces y conteo de etiquetas de genero
meta = ptb.read_metadata(f"{D}/sample-metadata-n20.tsv")
from collections import Counter
print("n por grupo:", Counter(meta.values()))
for q, pref, n in ((f"{D}/phylum-table-relfreq.qza","p__",8),(f"{D}/genus-table-relfreq.qza","g__",12)):
    t = ptb.load_biom(q)
    t = t.filter(list(meta.keys()), axis="sample", inplace=False)
    data, _ = ptb.group_means(t, meta)
    top = sorted(data, key=lambda k: -sum(data[k].values())/3)[:n]
    for k in top:
        print(pref, f"{short_name(k,pref):35s}", {g: round(v,2) for g,v in data[k].items()})
    if pref == "g__":
        nz = [k for k in data if sum(data[k].values())>0]
        print("etiquetas de genero con lecturas (n=20):", len(nz))
