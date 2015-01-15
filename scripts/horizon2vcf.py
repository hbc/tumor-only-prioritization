#!/usr/bin/env python
"""Convert Horizon known data into an output VCF file for validation.

Usage:
  horizon2vcf.py <excel file> <reference file>
"""
import os
import sys

import pandas as pd
from pyfaidx import Fasta

def main(in_file, ref_file):
    out_file = "%s.vcf" % os.path.splitext(in_file)[:1]
    df = pd.read_excel(in_file, header=None)
    df.columns = ["coords", "blank", "gene", "event", "freq", "copies1", "copies2",
                  "cytoband", "in_pulldown"]
    for i, row in df.iterrows():
        chrom, start = row["coords"].split()[:2]
        df.loc[i, "chrom"] = int(chrom.replace("chr", ""))
        df.loc[i, "start"] = int(start)
    df = df.sort(["chrom", "start"])
    ref_fasta = Fasta(ref_file)
    with open(out_file, "w") as out_handle:
        h = ["##fileformat=VCFv4.1",
             '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '#' + "\t".join("CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO FORMAT Horizon".split())]
        out_handle.write("\n".join(h) + "\n")
        for i, row in df.iterrows():
            if row["in_pulldown"].lower() == "yes":
                chrom = str(int(row["chrom"]))
                start = int(row["start"])
                ref = ref_fasta[chrom][start - 1:start].seq
                parts = [chrom, str(start), ".", ref, "N", ".", ".", "VAF=%.3f" % row["freq"],
                         "GT", "0/1"]
                out_handle.write("\t".join(parts) + "\n")

if __name__ == "__main__":
    main(*sys.argv[1:])
