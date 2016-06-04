#!/usr/bin/env python
"""Run a tumor-only BubbleTree analysis using VarDict and Seq2c inputs.
"""
import csv
import glob
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

def main(vardict_flat, seq2c_input, annotate_conf):
    var_samples = read_samples(vardict_flat)
    cnv_samples = read_samples(seq2c_input)
    out_file = "dlbcl-breast-lung-af-dist.pdf"
    sample_freqs = {}
    with PdfPages(out_file) as pdf_out:
        for sample in sorted(list(var_samples & cnv_samples)):
            print sample
            sample_vcf = annotate_vcf(sample_variants_as_vcf(sample, vardict_flat), annotate_conf)
            seq2c_cns = sample_seq2c_as_cns(sample, seq2c_input)
            vrn_info = {"vrn_file": sample_vcf, "variantcaller": "vardict"}
            calls_by_name = {"seq2c": {"cns": seq2c_cns,
                                       "variantcaller": "seq2c"}}
            data = {"dirs": {"work": os.getcwd()},
                    "rgnames": {"sample": sample}}
            somatic_info = vcfutils.PairedData(None, sample, None, None, None, None, data, None)
            bubbletree.run(vrn_info, calls_by_name, somatic_info, do_plots=False)
            freqs = plot_frequencies(sample, sample_vcf, seq2c_cns,
                                     glob.glob(os.path.join("heterogeneity", sample, "bubbletree",
                                                            "*_prevalence.txt"))[0],
                                     pdf_out)
            sample_freqs[sample] = freqs
    comparison_plot(sample_freqs)

def comparison_plot(sample_freqs):
    out_file = "breast-af-comparison.pdf"
    metrics = {"sample": [], "af": []}
    for sample, freqs in sample_freqs.items():
        sample = sample.split("_")[0].replace("AZ-", "")
        metrics["sample"].extend([sample] * len(freqs))
        metrics["af"].extend(freqs)

    df = pd.DataFrame(metrics)
    sns.despine()
    sns.set(style="white")
    g = sns.swarmplot(x="af", y="sample", data=df)
    g.set_xlim(0, 1.0)
    g.figure.savefig(out_file)

def plot_frequencies(sample, sample_vcf, seq2c_cns, bubbletree_out, pdf_out):
    """Plot non-germline frequencies, adjusted by purity and copy number.
    """
    freqs = []
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

    sns.despine()
    sns.set(style="white")
    g = sns.distplot(freqs, kde=False, rug=True, bins=20)
    g.set_title("%s: purity %0.1f%%" % (sample, purity * 100.0))
    g.set_xlim(0, 1.0)
    g.set_xlabel("Adjusted allele frequency (copy number and purity)")
    pdf_out.savefig(g.figure)
    plt.clf()
    return freqs

def cns_to_relative_copy(in_cns):
    out = []
    with open(in_cns) as in_handle:
        in_handle.next()  # header
        for chrom, start, end, gene, log2, probes in (l.split("\t") for l in in_handle):
            out.append((chrom, int(start), int(end), math.pow(2,float(log2))))
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
        out_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
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
                                  "AF=%s;DP=%s" % (cur["AlleleFreq"], cur["Depth"])]
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
