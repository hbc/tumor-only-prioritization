"""Convert validated COLO-829 calls from Excel spreadsheets into VCF.

http://www.nature.com/nature/journal/v463/n7278/full/nature08658.html

Usage:
  colo2vcf.py <nature08658-s2.xls> <nature08658-s3.xls> <fasta ref file>
"""
import sys

import pandas as pd
from pyfaidx import Fasta

def main(snp_file, indel_file, ref_file):
    ref_fasta = Fasta(ref_file)
    out_file = "colo289-validated.vcf"
    snpdf = pd.read_excel(snp_file)
    snpdf = snpdf[snpdf["Validation_status"].isin(["SOMATIC", "COSMIC", "REAL, STATUS UNKNOWN"])]
    snpdf = calculate_call(snpdf)
    snpdf = snpdf[["Chromosome", "Position", "Reference", "Mutant", "call", "Validation_status"]]
    snpdf.columns = ["chr", "pos", "ref", "alt", "call", "valstatus"]

    idf = pd.read_excel(indel_file)
    idf = idf[idf["Validation_status"].isin(["SOMATIC", "REAL SOMATIC UNKNOWN"])]
    idf = add_ref_alt(idf, ref_fasta)
    idf = idf[["Chromosome", "Start", "ref", "alt", "call", "Validation_status"]]
    idf.columns = ["chr", "pos", "ref", "alt", "call", "valstatus"]
    df = pd.concat([snpdf, idf])
    df = df.sort(["chr", "pos"])

    with open(out_file, "w") as out_handle:
        h = ["##fileformat=VCFv4.1",
             '##INFO=<ID=VS,Number=1,Type=String,Description="Validation status">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '#' + "\t".join("CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO FORMAT COLO-829".split())]
        out_handle.write("\n".join(h) + "\n")
        for i, row in df.iterrows():
            parts = [str(row["chr"]), str(int(row["pos"])), ".", row["ref"], row["alt"], ".", ".",
                     "VS=%s" % row["valstatus"].split()[0].replace(",", ""), "GT", row["call"]]
            out_handle.write("\t".join(parts) + "\n")

def add_ref_alt(df, ref_fasta):
    for i, row in df.iterrows():
        chrom = row.Chromosome
        if row.Chromosome == 23:
            chrom = "X"
        start = row.Start - 1
        ref = ref_fasta[chrom][start:start + 1].seq
        if not ref:
            ref = "N"
        if row.Type.lower() == "insertion":
            alt = ref + row.Sequence
        else:
            alt = ref
            ref = ref + row.Sequence
        df.loc[i, "Chromosome"] = chrom
        df.loc[i, "ref"] = ref
        df.loc[i, "alt"] = alt
        df.loc[i, "call"] = "1/1"
    return df

def calculate_call(df):
    def get_call(row):
        z = row.Validation_zygosity if row.Validation_zygosity else row.Zygosity
        return "0/1" if not pd.isnull(z) and z.lower() == "het" else "1/1"
    df["call"] = df.apply(get_call, axis=1)
    return df

if __name__ == "__main__":
    main(*sys.argv[1:])
