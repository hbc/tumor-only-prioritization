#!/usr/bin/env python
"""Run a tumor-only BubbleTree analysis using VarDict and Seq2c inputs.
"""
import collections
import csv
import glob
import gzip
import math
import os
import subprocess
import sys

import pandas as pd
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from bcbio import utils
from bcbio.heterogeneity import bubbletree
from bcbio.variation import vcfutils

def main(vardict_flat, seq2c_input, annotate_conf, priority_file, name):
    var_samples = read_samples(vardict_flat)
    cnv_samples = read_samples(seq2c_input)
    freq_out_file = "%s-dlbcl-af-dist.pdf" % name
    sample_freqs = collections.defaultdict(dict)
    known_freqs = {}
    with PdfPages(freq_out_file) as pdf_out:
        for i, sample in enumerate(sorted(list(var_samples & cnv_samples))):
            sample_vcf = annotate_vcf(sample_variants_as_vcf(sample, vardict_flat), annotate_conf)
            seq2c_cns = sample_seq2c_as_cns(sample, seq2c_input)
            known_genes, known_positions = prioritize_variants(sample_vcf, seq2c_cns, priority_file)
            print sample, known_genes
            vrn_info = {"vrn_file": sample_vcf, "variantcaller": "vardict"}
            calls_by_name = {"seq2c": {"cns": seq2c_cns,
                                       "variantcaller": "seq2c"}}
            data = {"dirs": {"work": os.getcwd()},
                    "rgnames": {"sample": sample}}
            somatic_info = vcfutils.PairedData(None, sample, None, None, None, None, data, None)
            bubbletree_out = glob.glob(os.path.join("heterogeneity", sample, "bubbletree",
                                                    "*_prevalence.txt"))
            if len(bubbletree_out) == 0:
                bubbletree.run(vrn_info, calls_by_name, somatic_info, do_plots=False)
                bubbletree_out = glob.glob(os.path.join("heterogeneity", sample, "bubbletree",
                                                        "*_prevalence.txt"))
            freqs, sample_known_freqs = plot_frequencies(sample, sample_vcf, seq2c_cns, known_positions,
                                                         bubbletree_out[0], pdf_out)
            known_freqs[sample] = sample_known_freqs
            for driver in known_genes:
                sample_freqs[driver][sample] = freqs

    for driver, cur_sample_freqs in sample_freqs.items():
        if len(cur_sample_freqs) > 1:
            print driver, len(cur_sample_freqs)
            comparison_out_file = "%s-af-comparison-%s.pdf" % (name, driver.split(":")[-1])
            comparison_plot(driver, cur_sample_freqs, known_freqs, comparison_out_file)

def comparison_plot(driver, sample_freqs, known_freqs, out_file):
    metrics = {"sample": [], "af": []}
    known_metrics = {"sample": [], "af": []}
    for sample, freqs in sorted(sample_freqs.items()):
        sample_label = sample.split("_")[0].replace("AZ-", "")
        metrics["sample"].extend([sample_label] * len(freqs))
        metrics["af"].extend(freqs)
        for kfreq in known_freqs[sample]:
            known_metrics["sample"].append(sample_label)
            known_metrics["af"].append(kfreq)

    df = pd.DataFrame(metrics)
    df_known = pd.DataFrame(known_metrics)
    sns.despine()
    sns.set(style="white")
    g = sns.violinplot(x="af", y="sample", data=df, inner=None)
    if len(df_known) > 0:
        sns.swarmplot(x="af", y="sample", data=df_known, color="w", alpha=.5)
    g.set_title(driver)
    g.set_xlim(0, 1.0)
    g.figure.savefig(out_file)
    plt.clf()

def plot_frequencies(sample, sample_vcf, seq2c_cns, known_positions,
                     bubbletree_out, pdf_out):
    """Plot non-germline frequencies, adjusted by purity and copy number.
    """
    freqs = []
    known_freqs = []
    copy_adjust = cns_to_relative_copy(seq2c_cns)
    purity = parse_bubbletree(bubbletree_out)
    for rec in pysam.VariantFile(sample_vcf):
        if not bubbletree.is_population_germline(rec):
            cur_copy_adjust = 1.0
            for chrom, start, end, cadj in copy_adjust:
                if chrom == "chr%s" % rec.chrom and rec.start >= start and rec.start < end:
                    cur_copy_adjust = cadj
                    break
            baf = float(rec.info["AF"][0] if isinstance(rec.info["AF"], (tuple, list)) else rec.info["AF"])
            af = (baf * cur_copy_adjust) / purity
            # larger than one variants are likely germline since contamination
            # adjustment gives a non-sense frequency
            if af < 1:
                freqs.append(af)
                if (rec.chrom, rec.start) in known_positions:
                    known_freqs.append(af)

    sns.despine()
    sns.set(style="white")
    g = sns.distplot(freqs, kde=False, rug=True, bins=20)
    g.set_title("%s: purity %0.1f%%" % (sample, purity * 100.0))
    g.set_xlim(0, 1.0)
    g.set_xlabel("Adjusted allele frequency (copy number and purity)")
    pdf_out.savefig(g.figure)
    plt.clf()
    return freqs, known_freqs

def prioritize_variants(vcf_file, cns_file, priority_file):
    """Find known cancer associated genes in small variants and CNVs.
    """
    known, known_positions = _prioritize_vcf(vcf_file, priority_file)
    known += _prioritize_cns(cns_file, priority_file)
    if known:
        if len(known) > 7:
            known = ["multiple"]
        else:
            known = list(set(known))
    else:
        known = ["unclassified"]
    return known, known_positions

def _prioritize_vcf(in_file, priority_file):
    out_file = "%s-known.vcf.gz" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        cmd = ["bcbio-prioritize", "known", "-i", in_file, "-k", priority_file, "-o", out_file]
        subprocess.check_call(cmd)
    known = []
    known_positions = set([])
    if vcfutils.vcf_has_variants(out_file):
        for rec in pysam.VariantFile(out_file):
            if not bubbletree.is_population_germline(rec):
                known.extend(rec.info["KNOWN"])
                known_positions.add((rec.chrom, rec.start))
    return known, known_positions

def _chr_sort(region):
    chrom, start, end = region[:3]
    chrom = chrom.replace("chr", "")
    try:
        chrom = int(chrom)
    except ValueError:
        pass
    return (chrom, int(start), int(end))

def _prioritize_cns(in_file, priority_file):
    out_file = "%s-known.bed.gz" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        bed_file = "%s.bed" % utils.splitext_plus(in_file)[0]
        if not utils.file_uptodate(bed_file, in_file):
            cn_changes = []
            with open(in_file) as in_handle:
                in_handle.next()  # header
                for chrom, start, end, gene, log2, _ in (x.split("\t") for x in in_handle):
                    cn = math.pow(2, float(log2)) * 2.0
                    if cn < 1 or cn > 4:
                        cn_changes.append((chrom.replace("chr", ""), start, end, "%s_%.1f" % (gene, cn)))
            with open(bed_file, "w") as out_handle:
                for chrom, start, end, gene in sorted(cn_changes, key=_chr_sort):
                    out_handle.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, gene))
        cmd = ["bcbio-prioritize", "known", "-i", bed_file, "-k", priority_file, "-o", out_file]
        subprocess.check_call(cmd)
    known = []
    with gzip.open(out_file) as in_handle:
        for line in in_handle:
            match = line.split("\t")[-1]
            if match.strip():
                known.append(match.strip())
    return known

def cns_to_relative_copy(in_cns):
    out = []
    with open(in_cns) as in_handle:
        in_handle.next()  # header
        for chrom, start, end, gene, log2, probes in (l.split("\t") for l in in_handle):
            out.append((chrom, int(start), int(end), math.pow(2, float(log2))))
    return out

def parse_bubbletree(in_file):
    with open(in_file) as in_handle:
        in_handle.next()  # header
        return float(in_handle.next().split(",")[1])

def sample_seq2c_as_cns(sample, seq2c_input):
    """Extract seq2c sample output as a CNS formatted file.
    """
    out_file = "%s.cns" % sample
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle, dialect="excel-tab")
        header = ["chromosome", "start", "end", "gene", "log2", "probes"]
        writer.writerow(header)
        with open(seq2c_input) as in_handle:
            header = in_handle.readline().strip().split("\t")
            for line in in_handle:
                parts = line.split("\t")
                if parts[0] == sample:
                    cur = dict(zip(header, parts))
                    writer.writerow([cur["Chr"], cur["Start"], cur["End"], cur["Gene"], cur["Log2ratio"],
                                     cur["Total_Seg"]])
    return out_file

def sample_variants_as_vcf(sample, vardict_flat):
    """Retrieve variants for the input sample
    """
    out_file = "%s.vcf" % sample
    with open(out_file, "w") as out_handle:
        out_handle.write("##fileformat=VCFv4.2\n")
        out_handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        out_handle.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        out_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
        with open(vardict_flat) as in_handle:
            header = in_handle.next().strip().split("\t")
            for line in in_handle:
                parts = line.strip().split("\t")
                if parts[0] == sample:
                    cur = dict(zip(header, parts))
                    if len(cur["Ref"]) == 1 and len(cur["Alt"]) == 1:
                        chrom = cur["Chr"].replace("chr", "")
                        try:
                            chrom = int(chrom)
                        except ValueError:
                            continue
                        curout = [str(chrom), cur["Start"], ".", cur["Ref"], cur["Alt"], ".", "PASS",
                                  "AF=%s;DP=%s" % (cur["AlleleFreq"], cur["Depth"]), "GT", "0/1"]
                        out_handle.write("\t".join(curout) + "\n")
    return out_file

def annotate_vcf(in_file, conf_file):
    out_file = "%s-annotate.vcf.gz" % os.path.splitext(in_file)[0]
    if not utils.file_exists(out_file):
        cmd = "vcfanno {conf_file} {in_file} | bgzip -c > {out_file}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    vcfutils.bgzip_and_index(out_file)
    return out_file

def read_samples(in_file):
    """Read sample names from VarDict or Seq2c flat output.
    """
    samples = set([])
    with open(in_file) as in_handle:
        in_handle.next()  # header
        for line in in_handle:
            samples.add(line.split()[0])
    return samples

if __name__ == "__main__":
    main(*sys.argv[1:])
