#!/usr/bin/env python
"""Plot cohort sample grouping and comparisons from pre-calculated heterogeneity.
"""
import collections
import os
import sys

from bcbio import utils

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main(freq_files):
    af_key = "adjust_af"
    #af_key = "af"
    plot_dir = utils.safe_makedir("plots")
    df = pd.concat([pd.read_csv(f) for f in freq_files])
    df = df[df[af_key] < 1.0]
    fishplot_by_groups(df, plot_dir, af_key)
    plot_by_groups(df, plot_dir, af_key)
    plot_by_genes(df, plot_dir, af_key)

def _round_freq(freq):
    #return int(round(freq * 100.0, -1))
    return int(round(freq * 100.0 / 20.0) * 20)

def fishplot_by_groups(df, plot_dir, af_key):
    """Group samples into clones by frequency changes for plotting.
    """
    for (cohort, group), cur_df in df.groupby(["cohort", "group"]):
        clones = collections.defaultdict(int)
        for pos, pos_df in cur_df.groupby(["pos"]):
            pre = None
            post = None
            for _, row in pos_df.iterrows():
                known = row["known"]
                if "pre" in row["group_class"].lower() or "primary" in row["group_class"].lower():
                    pre = _round_freq(row[af_key])
                elif "post" in row["group_class"].lower() or "metastatic" in row["group_class"].lower():
                    post = _round_freq(row[af_key])
                else:
                    raise ValueError("Unexpected group: %s" % row["group_class"])
            clones[(pre, post, known)] += 1
        print cohort, group
        import pprint
        pprint.pprint(dict(clones))

def plot_by_groups(df, plot_dir, af_key):
    """Plot allele frequencies of grouped/paired samples.
    """
    out_file = os.path.join(plot_dir, "cohort-group-af-comparison.pdf")
    df["sample_label"] = df.apply(lambda row: "%s\n%s" % (row["group_class"], row["sample"]), axis=1)
    sns.despine()
    sns.set(style="white")
    with PdfPages(out_file) as pdf_out:
        for (cohort, group), cur_df in df.groupby(["cohort", "group"]):
            labels = sorted(list(cur_df["sample_label"].unique()))
            labels.reverse()
            cur_df["sample_label"].categories = labels
            g = sns.violinplot(x=af_key, y="sample_label", data=cur_df, inner=None, bw=.1)
            #sns.swarmplot(x=af_key, y="sample_label", data=cur_df, color="w", alpha=.5)
            try:
                group = int(group)
            except ValueError:
                pass
            g.set_title("%s: %s" % (cohort, group))
            g = _af_violinplot_shared(g)
            pdf_out.savefig(g.figure)
            plt.clf()
    return out_file

def plot_by_genes(df, plot_dir, af_key):
    """Plot allele frequencies of known cancer genes in primary, relapse status
    """
    out_file = os.path.join(plot_dir, "driver-af-comparison.pdf")
    df = df[pd.notnull(df["known"])]
    with PdfPages(out_file) as pdf_out:
        for cohort, cohort_df in df.groupby(["cohort"]):
            labels = sorted(list(cohort_df["status"].unique()))
            labels.reverse()
            cohort_df["status"].categories = labels
            g = sns.violinplot(x=af_key, y="status", data=cohort_df, inner=None)
            g.set_title("%s -- %s cancer genes" % (cohort, len(cohort_df["known"].unique())))
            g = _af_violinplot_shared(g)
            pdf_out.savefig(g.figure)
            plt.clf()
        for cohort, cohort_df in df.groupby(["cohort"]):
            for gene, gene_df in cohort_df.groupby(["known"]):
                if len(gene_df["status"].unique()) > 1 and len(gene_df) > 10:
                    gene_df["sample_label"] = gene_df.apply(
                        lambda row: "%s\n(%s variants)" %
                        (row["status"],
                         len(gene_df[gene_df["status"] == row["status"]])),
                        axis=1)
                    labels = list(gene_df["sample_label"].unique())
                    labels.reverse()
                    gene_df["sample_label"].categories = labels
                    g = sns.violinplot(x=af_key, y="sample_label", data=gene_df, inner=None, bw=.1)
                    sns.swarmplot(x=af_key, y="sample_label", data=gene_df, color="w", alpha=.5)
                    g.set_title("%s -- %s" % (cohort, gene))
                    g = _af_violinplot_shared(g)
                    pdf_out.savefig(g.figure)
                    plt.clf()
    return out_file

def _af_violinplot_shared(g):
    g.set_xlim(0, 1.0)
    g.set_xlabel("Adjusted allele frequency (copy number and purity)")
    g.tick_params(axis="y", which="major", labelsize=8)
    g.set_ylabel("")
    g.figure.tight_layout()
    return g

if __name__ == "__main__":
    main(sys.argv[1:])
