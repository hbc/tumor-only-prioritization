#!/usr/bin/env python
"""Evaluate concordance of tumor-only prioritized variants versus tumor/normal calls.

Usage:
  evaluate_tumoronly.py config-inputs.yaml

To improve standard calling on:

2 232458818 in HB2 -- rejected and priority filtered by exac, 1000g
~10% in East Asian populations
http://exac.broadinstitute.org/variant/2-232458818-C-T
"""
import collections
import csv
import itertools
import os
import sys

from tabulate import tabulate
import toolz as tz
import vcf
import yaml

from bcbio import utils
from bcbio.provenance import do

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    confirmed = read_confirmed(tz.get_in(["inputs", "confirmed"], config),
                               tz.get_in(["inputs", "skip"], config))
    evals = [evaluate_tumoronly(b, c, config, confirmed.get(b, []))
             for b, c in itertools.product(tz.get_in(["inputs", "batches"], config),
                                           tz.get_in(["inputs", "callers"], config))]
    header = ["sample", "caller", "confirmed", "enrichment", "additional", "filtered"]
    print tabulate(evals, header, tablefmt="orgtbl")

def evaluate_tumoronly(batch, caller, config, confirmed):
    name = tz.get_in(["inputs", "name_base"], config) % batch
    n_file = os.path.join(tz.get_in(["dirs", "final"], config), name, "%s-%s.vcf.gz" % (name, caller))
    t_file = os.path.join(tz.get_in(["dirs", "final"], config), name + "-t", "%s-t-%s.vcf.gz" % (name, caller))

    n_file = os.path.abspath(n_file)
    t_file = os.path.abspath(t_file)
    region_file = os.path.abspath(tz.get_in(["resources", "region_file"], config))

    cmp_dir = os.path.abspath(utils.safe_makedir(os.path.join(tz.get_in(["dirs", "out"], config),
                                                              "%s-%s" % (batch, caller))))
    work_dir = os.path.join(cmp_dir, "%s-evalwork" % utils.splitext_plus(os.path.basename(t_file))[0],
                            "work")
    with utils.chdir(cmp_dir):
        out_file = os.path.join(work_dir, "eval-config-grading.yaml")
        if not utils.file_exists(out_file):
            cmd = ["java", "-Xms500m", "-Xmx2000m", "-XX:+UseSerialGC",
                   "-jar", tz.get_in(["resources", "bcbio_variation"], config),
                   "variant-utils", "comparetwo", t_file, n_file, tz.get_in(["resources", "ref_file"], config),
                   region_file, "-s", name]
            do.run(cmd, "bcbio.variation comparison")
    with open(out_file) as in_handle:
        stats = yaml.load(in_handle)[0]
    col = check_overlap(confirmed, os.path.join(work_dir, "%s-eval-ref-concordance.vcf" % name))
    dol = check_overlap(confirmed, os.path.join(work_dir, "%s-eval-ref-discordance-annotate.vcf" % name))
    print "--", batch, caller
    print "Confirmed: total %s, concordant %s, discordant %s" % (len(confirmed), col, dol)
    con = stats["concordant"]["concordant"]["total"]
    missing = sum(stats["discordant"]["snp"]["missing"].values() +
                  tz.get_in(["discordant", "indel", "missing"], stats, {}).values())
    extra = sum(stats["discordant"]["snp"]["extra"].values() +
                tz.get_in(["discordant", "indel", "extra"], stats, {}).values())
    print "Total: concordant %s, missing %s, extra %s" % (con, missing, extra)
    enrichment, filtered = _find_filtered(t_file, extra)
    return [name, caller, "%s / %s" % (col, len(confirmed)), enrichment, extra, filtered]

def _find_filtered(fname, extra):
    """Identify the filtered inputs in the original VCF file.
    """
    filtered = 0
    with utils.open_gzipsafe(fname) as in_handle:
        for rec in vcf.Reader(in_handle, fname):
            if "LowPriority" in rec.FILTER:
                filtered += 1
    enrichment = "%sx" % (int((extra + filtered) / float(extra)))
    return enrichment, filtered

def check_overlap(confirmed, call_file):
    confirmed = set(confirmed)
    reader = vcf.Reader(filename=call_file)
    ol = 0
    for rec in reader:
        key = (rec.CHROM, int(rec.POS) - 1, rec.REF, ",".join(str(x) for x in rec.ALT))
        if key in confirmed:
            confirmed.remove(key)
            ol += 1
    #print "Remaining", list(confirmed)
    return ol

def read_confirmed(in_file, to_skip):
    to_skip = set([tuple(xs) for xs in to_skip])
    out = collections.defaultdict(list)
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        header = reader.next()
        for chrom, start, _, ref, alt, _, sample in (xs[:7] for xs in reader):
            sample = int(sample.replace("HB", ""))
            key = (chrom.replace("chr", ""), int(start) - 1, ref, alt)
            if key not in to_skip:
                out[sample].append(key)
    return dict(out)

if __name__ == "__main__":
    main(*sys.argv[1:])
