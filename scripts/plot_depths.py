import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
import sys

# Function to read and aggregate depth from a single gzipped file in chunks
def aggregate_depth(file):
    all_depths = []
    with gzip.open(file, 'rt') as f:  # Open in text mode for reading
        for chunk in pd.read_csv(f, sep='\t', chunksize=1000000, dtype={'depth': 'int32'}, usecols=['depth']):
            all_depths.extend(chunk['depth'].tolist())
    return all_depths

# Function to calculate the bottom 10% and top 2% cutoffs
def calculate_cutoffs(depths):
    bottom_10_cutoff = np.percentile(depths, 10)
    top_2_cutoff = np.percentile(depths, 98)
    return bottom_10_cutoff, top_2_cutoff

# Main script
if __name__ == "__main__":
    depth_file = sys.argv[1]  # The first argument is the input file
    output_file = sys.argv[2]  # The second argument is the output file

    # Aggregate depths
    aggregated_depths = aggregate_depth(depth_file)

    # Calculate cutoffs
    bottom_10_cutoff, top_2_cutoff = calculate_cutoffs(aggregated_depths)

    # Plotting the histogram with increased bin size and limited x-axis
    plt.figure(figsize=(10, 6))
    plt.hist(aggregated_depths, bins=100, color='skyblue', edgecolor='black', range=(0, 1000))  # Adjust bin size and range
    plt.axvline(bottom_10_cutoff, color='red', linestyle='--', label=f'Bottom 10% cutoff: {bottom_10_cutoff:.2f}')
    plt.axvline(top_2_cutoff, color='red', linestyle='--', label=f'Top 2% cutoff: {top_2_cutoff:.2f}')
    plt.title('Histogram of Depth for Scaffold')
    plt.xlabel('Total Site Depth')
    plt.ylabel('Number of Sites')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend()
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_file)
