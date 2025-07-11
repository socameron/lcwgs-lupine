#!/usr/bin/env python3

import sys
import gzip

# --- Input arguments ---
pruned_lists = sys.argv[1:-1]  # All input .txt files
output_beagle = sys.argv[-1]  # Final output .beagle.gz file

# --- Step 1: Load all pruned SNPs into a set ---
pruned_snps = set()
for file in pruned_lists:
    with open(file, 'r') as f:
        for line in f:
            pruned_snps.add(line.strip())

# --- Step 2: Read stdin (or the original beagle file) ---
# NOTE: We assume you are running this on the *combined* unfiltered Beagle file
# If you're not, you may need to modify this to open a hardcoded or parameterized file
# For this example, we will just read from stdin for clarity

# TO ADAPT: Replace this with your original Beagle file path if needed
beagle_input = "results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.beagle.gz"

with gzip.open(beagle_input, "rt") as infile, gzip.open(output_beagle, "wt") as outfile:
    for line_num, line in enumerate(infile):
        if line_num == 0:
            outfile.write(line)  # Keep the header
        else:
            snp_id = line.split()[0]  # Marker ID (e.g., Scaffold_1__..._12345)
            if snp_id in pruned_snps:
                outfile.write(line)
