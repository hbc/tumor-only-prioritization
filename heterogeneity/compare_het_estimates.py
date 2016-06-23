#!/usr/bin/env python
"""Compare measures of heterogeneity from multiple estimation methods.
"""
import csv
import glob
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def main(het_dir, purecn_file, tcga_file):
    bubbletree_df = bubbletree_parse(het_dir)
    purecn_df = purecn_parse(purecn_file)
    want_ids = set(bubbletree_df["sample"].unique()) | set(purecn_df["sample"].unique())
    tcga_df = tcga_parse(tcga_file, want_ids)

    df = pd.merge(bubbletree_df, purecn_df, on="sample", how="outer")
    df = pd.merge(df, tcga_df, on="sample", how="outer")
    print df.head()

    sns.despine()
    sns.set(style="white")
    g = sns.PairGrid(df)
    g.map_lower(sns.regplot)
    g.map_upper(plt.scatter)
    g.map_diag(sns.kdeplot, legend=False)
    plt.subplots_adjust(top=0.95)
    g.fig.suptitle("Heterogeneity method comparison: BubbleTree, PureCN, ABSOLUTE, ESTIMATE", size=16)
    g.fig.savefig("plots/het-comparison.pdf")
    g.fig.savefig("plots/detailed/het-comparison.png")
    df.to_csv("plots/het-comparison.csv", index_label="sample")

def bubbletree_parse(het_dir):
    out = {"sample": [], "bubbletree": []}
    for fname in glob.glob(os.path.join(het_dir, "*", "bubbletree", "*_prevalence.txt")):
        with open(fname) as in_handle:
            reader = csv.reader(in_handle)
            reader.next()  # header
            sample, purity, prev, ploidy = reader.next()[:4]
            if not ((purity == "1" and ploidy == "2") or (purity == "1.0" and ploidy == "2.0")):
                name = parse_tcga_name(sample)
                if name not in out["sample"]:
                    out["sample"].append(name)
                    out["bubbletree"].append(float(purity))
    return pd.DataFrame(out)

def purecn_parse(in_file):
    out = {"sample": [], "purecn": []}
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        reader.next()  # header
        for sample, purity in (xs[:2] for xs in reader):
            name = parse_tcga_name(sample)
            if name not in out["sample"]:
                out["sample"].append(name)
                out["purecn"].append(float(purity))
    return pd.DataFrame(out)

def tcga_parse(in_file, want_ids):
    to_use = ["estimate", "absolute", "cpe"]  # ["ihc", "lump"]
    out = {"sample": []}
    for k in to_use:
        out[k] = []
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        while 1:
            header = reader.next()
            if header[0].startswith("Sample"):
                break
        header = [x.lower() for x in header]
        for line in reader:
            val = dict(zip(header, line))
            if val["sample id"] in want_ids:
                out["sample"].append(val["sample id"])
                for k in to_use:
                    if val[k] == "NaN":
                        out[k].append(None)
                    else:
                        out[k].append(float(val[k]))
    return pd.DataFrame(out)

def parse_tcga_name(name):
    if name.startswith("TCGA"):
        parts = name.split("-")
        name = "-".join(parts[:4])
    return name

if __name__ == "__main__":
    main(*sys.argv[1:])
