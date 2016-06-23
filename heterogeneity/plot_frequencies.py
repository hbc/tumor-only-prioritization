#!/usr/bin/env python
"""Plot cohort sample grouping and comparisons from pre-calculated heterogeneity.
"""
import collections
import imp
import os
import sys

from bcbio import utils

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main(input_files):
    if len(input_files) > 0 and input_files[0].endswith(".py"):
        config_file = input_files[0]
        freq_files = input_files[1:]
    else:
        config_file = None
        freq_files = input_files
    config = imp.load_source("config", config_file) if config_file else None
    plot_dir = utils.safe_makedir("plots")
    df = pd.concat([pd.read_csv(f) for f in freq_files])
    organize_samples(df)
    plot_af_adjustments(df, plot_dir, config)

    af_key = "adjust_af"
    # af_key = "af"
    df = df[df[af_key] < 1.0]
    fishplot_by_groups(df, plot_dir, af_key, config)
    plot_by_groups(df, plot_dir, af_key, config)
    plot_by_genes(df, plot_dir, af_key, config)

def organize_samples(df):
    out_file = "db/samples-by-cohort.csv"
    df = df[["sample", "cohort", "status", "group", "group_class"]].drop_duplicates()
    df.to_csv(out_file, index_label="sample")

def _round_freq(freq):
    #return int(round(freq * 100.0, -1))
    return int(round(freq * 100.0 / 20.0) * 20)

def fishplot_by_groups(df, plot_dir, af_key, config):
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
        try:
            group = int(group)
        except ValueError:
            pass
        if not config or(cohort, group) in config.group_detailed:
            print cohort, group
            import pprint
            pprint.pprint(dict(clones))

def plot_by_groups(df, plot_dir, af_key, config):
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
            if config and (cohort, group) in config.group_detailed:
                out_dir = utils.safe_makedir(os.path.join(plot_dir, "detailed"))
                out_file = os.path.join(out_dir, "group-%s-%s.png" % (cohort, group))
                g.figure.savefig(out_file)
            plt.clf()
    return out_file

def plot_af_adjustments(df, plot_dir, config):
    """Plot results of adjusting allele frequencies in groups
    """
    out_file = os.path.join(plot_dir, "summarize-af-adjustment.pdf")
    df["sample_label"] = df.apply(lambda row: "%s\n%s" % (row["group_class"], row["sample"]), axis=1)
    sns.despine()
    sns.set(style="white")
    with PdfPages(out_file) as pdf_out:
        for (cohort, group), cur_df in df.groupby(["cohort", "group"]):
            labels = sorted(list(cur_df["sample_label"].unique()))
            labels.reverse()
            cur_df["sample_label"].categories = labels
            cur_df = cur_df[["sample_label", "adjust_af", "af"]]
            cur_df = pd.melt(cur_df, id_vars=["sample_label"], value_vars=["adjust_af", "af"],
                             var_name="Allele frequency")
            cur_df["Allele frequency"] = cur_df["Allele frequency"].map({"adjust_af": "Purity and copy adjusted",
                                                                         "af": "Raw"})
            g = sns.violinplot(x="value", y="sample_label", data=cur_df, inner=None, bw=.1,
                               hue="Allele frequency", split=True)
            try:
                group = int(group)
            except ValueError:
                pass
            g.set_title("%s: %s" % (cohort, group))
            g.set_xlabel("Allele frequency")
            g.set_xlim(0, 1.0)
            g.tick_params(axis="y", which="major", labelsize=8)
            g.set_ylabel("")
            g.figure.tight_layout()
            pdf_out.savefig(g.figure)
            if config and (cohort, group) in config.group_detailed:
                out_dir = utils.safe_makedir(os.path.join(plot_dir, "detailed"))
                out_file = os.path.join(out_dir, "af-adjustment-%s-%s.png" % (cohort, group))
                g.figure.savefig(out_file)
            plt.clf()
    return out_file

def plot_by_genes(df, plot_dir, af_key, config):
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
                    if config and (cohort, gene) in config.driver_detailed:
                        out_dir = utils.safe_makedir(os.path.join(plot_dir, "detailed"))
                        out_file = os.path.join(out_dir, "driver-%s-%s.png" % (cohort, gene))
                        g.figure.savefig(out_file)
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
