#!/usr/bin/env python
"""Plot cohort sample grouping and comparisons from pre-calculated heterogeneity.
"""
import os
import sys

from bcbio import utils

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main(freq_file):
    plot_dir = utils.safe_makedir("plots")
    df = pd.read_csv(freq_file)
    plot_by_groups(df, plot_dir)
    plot_by_genes(df, plot_dir)

def plot_by_groups(df, plot_dir):
    """Plot allele frequencies of grouped/paired samples.
    """
    out_file = os.path.join(plot_dir, "cohort-groups.pdf")
    df["sample_label"] = df.apply(lambda row: "%s\n%s" % (row["sample"], row["group_class"]), axis=1)
    sns.despine()
    sns.set(style="white")
    with PdfPages(out_file) as pdf_out:
        for (cohort, group), cur_df in df.groupby(["cohort", "group"]):
            g = sns.violinplot(x="af", y="sample_label", data=cur_df, inner=None)
            try:
                group = int(group)
            except ValueError:
                pass
            g.set_title("%s: %s" % (cohort, group))
            g = _af_violinplot_shared(g)
            pdf_out.savefig(g.figure)
            plt.clf()
    return out_file

def plot_by_genes(df, plot_dir):
    """Plot allele frequencies of known cancer genes in primary, relapse status
    """
    out_file = os.path.join(plot_dir, "driver-af-dist.pdf")
    df = df[pd.notnull(df["known"])]
    with PdfPages(out_file) as pdf_out:
        for cohort, cohort_df in df.groupby(["cohort"]):
            g = sns.violinplot(x="af", y="status", data=cohort_df, inner=None)
            g.set_title("%s -- %s cancer genes, %s samples" % (cohort,
                                                               len(cohort_df["known"].unique()),
                                                               len(cohort_df["sample"].unique())))
            g = _af_violinplot_shared(g)
            pdf_out.savefig(g.figure)
            plt.clf()
        for cohort, cohort_df in df.groupby(["cohort"]):
            for gene, gene_df in cohort_df.groupby(["known"]):
                if len(gene_df["status"].unique()) > 1:
                    g = sns.violinplot(x="af", y="status", data=gene_df, inner=None)
                    g.set_title("%s -- %s, %s samples" % (cohort, gene, len(gene_df["sample"].unique())))
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
    main(*sys.argv[1:])
