#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: script.py <GFF file> <FASTA file>")
    sys.exit(1)

gff_file = sys.argv[1]
fasta_file = sys.argv[2]

# === Load reference genome ===
ref = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# === Collect CDS regions per transcript ===
cds_regions = defaultdict(list)

with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        chrom, _, feature, start, end, _, strand, _, attr = fields
        if feature != "CDS":
            continue
        start, end = int(start), int(end)
        transcript_id = None
        for x in attr.split(";"):
            if "Parent=" in x:
                transcript_id = x.split("Parent=")[1].split(",")[0]
        if transcript_id:
            cds_regions[(chrom, strand, transcript_id)].append((start, end))

# === Codon table ===
codon_table = CodonTable.unambiguous_dna_by_id[1]

# === Output containers ===
four_fold_sites = set()
other_coding_sites = set()

# === Process each transcript ===
for (chrom, strand, tid), regions in cds_regions.items():
    regions = sorted(regions, key=lambda x: x[0])
    seq = "".join(str(ref[chrom].seq[start - 1:end]) for start, end in regions)
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())

    for codon_start in range(0, len(seq) - 2, 3):
        codon = seq[codon_start:codon_start+3]
        if "N" in codon or len(codon) != 3:
            continue
        try:
            aa = codon_table.forward_table[codon]
        except KeyError:
            continue

        for i in range(3):
            original_base = codon[i]
            degenerate = True
            for alt in "ACGT":
                if alt == original_base:
                    continue
                test_codon = codon[:i] + alt + codon[i+1:]
                try:
                    alt_aa = codon_table.forward_table[test_codon]
                except KeyError:
                    degenerate = False
                    break
                if alt_aa != aa:
                    degenerate = False
                    break

            # === Map codon position to genomic position ===
            cds_pos = codon_start + i
            remaining = cds_pos
            if strand == "+":
                for start, end in regions:
                    region_len = end - start + 1
                    if remaining < region_len:
                        pos = start + remaining
                        break
                    remaining -= region_len
            else:  # reverse strand
                for start, end in reversed(regions):
                    region_len = end - start + 1
                    if remaining < region_len:
                        pos = end - remaining
                        break
                    remaining -= region_len

            site_str = f"{chrom}:{pos}"
            if degenerate:
                four_fold_sites.add(site_str)
            else:
                other_coding_sites.add(site_str)

# === Write output ===
with open("sites_4D.txt", "w") as f:
    for s in sorted(four_fold_sites):
        f.write(s + "\n")

with open("sites_otherCDS.txt", "w") as f:
    for s in sorted(other_coding_sites - four_fold_sites):
        f.write(s + "\n")
