import gzip
import sys

# Function to categorize SNP as transition or transversion
def categorize_snp(major, minor):
    transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
    transversions = [('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'),
                     ('G', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'G')]

    if (major, minor) in transitions:
        return 'transition'
    elif (major, minor) in transversions:
        return 'transversion'
    else:
        return 'unknown'

def main(mafs_file, output_file):
    transitions_count = 0
    transversions_count = 0

    # Open and read the .mafs.gz file
    with gzip.open(mafs_file, 'rt') as f:
        # Skip the header
        next(f)
        for line in f:
            cols = line.strip().split()
            major = cols[2]  # Major allele
            minor = cols[3]  # Minor allele

            category = categorize_snp(major, minor)
            if category == 'transition':
                transitions_count += 1
            elif category == 'transversion':
                transversions_count += 1

    # Calculate the Ts/Tv ratio
    if transversions_count > 0:
        ts_tv_ratio = transitions_count / transversions_count
    else:
        ts_tv_ratio = None

    # Output results
    with open(output_file, 'w') as out:
        out.write(f"Transitions: {transitions_count}\n")
        out.write(f"Transversions: {transversions_count}\n")
        if ts_tv_ratio is not None:
            out.write(f"Ts/Tv Ratio: {ts_tv_ratio:.2f}\n")
        else:
            out.write("No transversions found, Ts/Tv ratio is undefined.\n")

if __name__ == "__main__":
    # Get command-line arguments: mafs_file and output_file
    mafs_file = sys.argv[1]
    output_file = sys.argv[2]
    main(mafs_file, output_file)
