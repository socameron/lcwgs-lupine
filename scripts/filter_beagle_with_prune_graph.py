#!/usr/bin/env python3

import sys
import gzip

# --- Input arguments ---
beagle_input = sys.argv[1]
pruned_lists = sys.argv[2:-1]
output_beagle = sys.argv[-1]

# --- Function to normalize SNP names ---
# NOTE: the .beagle.gz file uses "_" as a separator. e.g: Scaffold_1__1_contigs__length_29266999_10000086
# Whereas, the prune .txt file uses ":" as a separator. e.g: Scaffold_1__1_contigs__length_29266999:10000086

def normalize(id_str):
    # split on colon and re-join with underscore
    chrom, pos = id_str.split(":", 1)
    return f"{chrom}_{pos}"


# --- Step 1: Load all pruned SNPs into a set ---
pruned_snps = set()
for fn in pruned_lists:
    with open(fn) as f:
        for line in f:
            pruned_snps.add(normalize(line.strip()))

# --- Step 2: Filter Beagle file ---
kept = 0
with gzip.open(beagle_input, "rt") as r, gzip.open(output_beagle, "wt") as w:
    for i, line in enumerate(r):
        if i == 0:
            if not line.lower().startswith("marker"):
                sys.exit("Error: Unexpected header")
            w.write(line)
        else:
            snp = line.split()[0]
            if snp in pruned_snps:
                w.write(line)
                kept += 1

if kept == 0:
    sys.exit("Error: no SNPs retained (ID format mismatch?)")

# --- Step 3: Final report ---
print(f"Kept {kept} SNPs out of {len(pruned_snps)} listed", file=sys.stderr)
print(f"Wrote pruned Beagle file to {output_beagle}", file=sys.stderr)
