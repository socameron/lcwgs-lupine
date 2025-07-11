#!/usr/bin/env python3

import gzip
import sys

def extract_key(line):
    parts = line.strip().split()
    return (parts[0], int(parts[1]))  # chrom (scaffold), position

def main(input_files, output_file):
    data = []

    print(f"Reading and collecting lines from {len(input_files)} files...")

    # Extract header from the first file
    with gzip.open(input_files[0], "rt") as f:
        header = f.readline()

    for file in input_files:
        with gzip.open(file, "rt") as f:
            first_line = f.readline()  # skip header
            if first_line.strip() != header.strip():
                print(f"Warning: header mismatch in {file}", file=sys.stderr)
            for line in f:
                if line.strip():
                    key = extract_key(line)
                    data.append((key, line))

    print(f"Sorting {len(data)} lines...")
    data.sort()

    print(f"Writing to {output_file}...")
    with gzip.open(output_file, "wt") as out:
        out.write(header)
        for _, line in data:
            out.write(line)

    print("Done.")

if __name__ == "__main__":
    input_list = sys.argv[1:-1]
    output_path = sys.argv[-1]
    main(input_list, output_path)
