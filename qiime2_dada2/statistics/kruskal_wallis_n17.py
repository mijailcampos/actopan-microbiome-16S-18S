#!/usr/bin/env python3
import subprocess, os, tempfile
import pandas as pd
import numpy as np
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
import biom

# BASE = directory containing the 16S/ outputs of pipeline/02 and 03
BASE = os.environ.get("BASE", ".")
OUT = f"{BASE}/16S/reanalysis_fecal_n5"

EXCLUDE_IDS = {"HecesP11_S169_L001", "HecesP12_S181_L001", "HecesP13_S98_L001", "HecesP14_S110_L001"}

JOBS = [
    ("16S", "phylum", f"{BASE}/16S/phylum-table-relfreq.qza", f"{BASE}/16S/sample-metadata.tsv", 1.0),
    ("16S", "genus",  f"{BASE}/16S/genus-table-relfreq.qza",  f"{BASE}/16S/sample-metadata.tsv", 2.0),
]

def export_biom_to_df(qza_path, tmpdir):
    outdir = os.path.join(tmpdir, "export")
    subprocess.run([os.environ.get("QIIME", "qiime"), "tools", "export", "--input-path", qza_path, "--output-path", outdir],
                    check=True, capture_output=True, text=True)
    biom_path = os.path.join(outdir, "feature-table.biom")
    table = biom.load_table(biom_path)
    df = table.to_dataframe(dense=True)
    return df.T

def load_metadata(meta_path):
    md = pd.read_csv(meta_path, sep="\t")
    md = md[md.iloc[:, 0] != "#q2:types"]
    md = md.set_index(md.columns[0])
    return md

all_results = []
for marker, level, qza, meta_path, min_pct in JOBS:
    with tempfile.TemporaryDirectory() as tmpdir:
        df = export_biom_to_df(qza, tmpdir)
    df = df.drop(index=[i for i in EXCLUDE_IDS if i in df.index])
    md = load_metadata(meta_path)
    df = df.join(md["sample-type"], how="inner")
    groups = df["sample-type"].unique().tolist()
    print(f"[{marker}-{level}] n por grupo tras excluir F6-F9:", df["sample-type"].value_counts().to_dict())
    feature_cols = [c for c in df.columns if c != "sample-type"]

    rows = []
    for feat in feature_cols:
        group_vals = [df.loc[df["sample-type"] == g, feat].values for g in groups]
        group_vals = [v for v in group_vals if len(v) > 0]
        if len(group_vals) < 2:
            continue
        all_vals = np.concatenate(group_vals)
        if np.all(all_vals == all_vals[0]):
            stat, p = np.nan, 1.0
        else:
            try:
                stat, p = kruskal(*group_vals)
            except ValueError:
                stat, p = np.nan, 1.0
        means = {g: df.loc[df["sample-type"] == g, feat].mean() * 100 for g in groups}
        sds = {g: df.loc[df["sample-type"] == g, feat].std() * 100 for g in groups}
        max_mean = max(means.values())
        rows.append({"marker": marker, "level": level, "taxon": feat, "H": stat, "p": p,
                     **{f"mean_pct_{g}": means[g] for g in groups},
                     **{f"sd_pct_{g}": sds[g] for g in groups},
                     "max_mean_pct": max_mean})

    res = pd.DataFrame(rows)
    valid = res.dropna(subset=["p"]).copy()
    reject, p_adj, _, _ = multipletests(valid["p"], method="fdr_bh")
    valid["p_adj"] = p_adj
    res = res.merge(valid[["taxon", "p_adj"]], on="taxon", how="left")
    res["min_pct_threshold"] = min_pct
    res["pass_filter"] = (res["p_adj"] < 0.05) & (res["max_mean_pct"] >= min_pct)
    all_results.append(res)

final = pd.concat(all_results, ignore_index=True)
final = final.sort_values(["marker", "level", "p_adj"])
out_all = f"{OUT}/Kruskal_Wallis_all_taxa_n17.csv"
out_sig = f"{OUT}/Kruskal_Wallis_significant_taxa_n17.csv"
final.to_csv(out_all, index=False)
final[final["pass_filter"]].to_csv(out_sig, index=False)
print("\nWrote:", out_all)
print("Wrote:", out_sig)
print("\n== Significativos (p_adj<0.05 y sobre umbral de abundancia) ==")
print(final[final["pass_filter"]][["marker","level","taxon","p_adj","max_mean_pct"]].to_string(index=False))
