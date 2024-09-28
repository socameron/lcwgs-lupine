import pandas as pd
import matplotlib.pyplot as plt
import sys
import gzip

# Function to read and aggregate depth from multiple gzipped files
def aggregate_depth(files):
    all_depths = []
    for file in files:
        with gzip.open(file, 'rt') as f:  # Open in text mode for reading
            df = pd.read_csv(f, sep='\t')
            all_depths.extend(df['total_depth'].tolist())
    return all_depths

# Main script
if __name__ == "__main__":
    depth_files = sys.argv[1:-1]  # All arguments except the last one are input files
    output_file = sys.argv[-1]  # The last argument is the output file

    # Aggregate depths
    aggregated_depths = aggregate_depth(depth_files)

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(aggregated_depths, bins=50, color='skyblue', edgecolor='black')
    plt.title('Histogram of Aggregate Depths Across All Scaffolds')
    plt.xlabel('Total Site Depth')
    plt.ylabel('Number of Sites')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_file)