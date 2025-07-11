#!/usr/bin/env python3
import sys
import gzip
import math
import argparse

def convert_line(line):
    # Split the line on whitespace.
    parts = line.strip().split()
    if len(parts) < 3:
        return None  # not enough columns
    marker = parts[0]
    
    # Map allele numeric codes to nucleotide letters
    allele_map = {"0": "A", "1": "C", "2": "G", "3": "T"}
    allele1 = allele_map.get(parts[1], parts[1])
    allele2 = allele_map.get(parts[2], parts[2])
    
    # Split the marker by underscore.
    # For example:
    # "Scaffold_11__2_contigs__length_31776681_10200"
    tokens = marker.split("_")
    if len(tokens) < 2:
        chrom = marker
        pos = "NA"
    else:
        chrom = tokens[0] + "_" + tokens[1]
        pos = tokens[-1]
    
    # Convert each likelihood from column 4 onward to a phred score.
    phred_scores = []
    for p_str in parts[3:]:
        try:
            p = float(p_str)
        except ValueError:
            p = 0.0
        # Convert probability to phred. Use a high value for p==0.
        if p == 0.0:
            phred = 100
        else:
            phred = int(-10 * math.log10(p) + 0.5)
        phred_scores.append(str(phred))
    
    # Build the final output line (tab-delimited)
    output_columns = [chrom, marker, pos, allele1, allele2] + phred_scores
    return "\t".join(output_columns)

def main():
    parser = argparse.ArgumentParser(
        description="Convert ANGSD .beagle.gz file to RZooROH-compatible GL format"
    )
    parser.add_argument("input", help="Input .beagle.gz file")
    parser.add_argument("output", help="Output text file")
    args = parser.parse_args()

    with gzip.open(args.input, "rt") as infile, open(args.output, "w") as outfile:
        # Skip the header line
        header = infile.readline()
        for line in infile:
            conv_line = convert_line(line)
            if conv_line is not None:
                outfile.write(conv_line + "\n")

if __name__ == "__main__":
    main()
