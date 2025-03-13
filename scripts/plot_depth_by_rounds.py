import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import gzip
from glob import glob

def load_depth_summary(file_list, sequencing_round):
    """Loads depth data from multiple scaffolds and computes median depth per scaffold."""
    scaffold_depths = []

    for filepath in file_list:
        scaffold_name = filepath.split("/")[-1].replace("_depth_est.txt.gz", "")  # Extract scaffold name
        depths = []

        with gzip.open(filepath, "rt") as f:
            for line in f:
                fields = line.strip().split()
                depth = int(fields[2])
                depths.append(depth)

        if depths:
            median_depth = pd.Series(depths).median()  # Compute median depth for this scaffold
            scaffold_depths.append({"Scaffold": scaffold_name, "MedianDepth": median_depth, "Sequencing Round": sequencing_round})

    return pd.DataFrame(scaffold_depths)

def main(first_round_files, second_round_files, output_plot):
    """Loads depth data, plots boxplot of depth per sequencing round, and saves the figure."""
    # Convert space-separated input file paths into lists
    first_round_files = first_round_files.split()
    second_round_files = second_round_files.split()

    # Load depth summaries for all scaffolds
    first_round_depth = load_depth_summary(first_round_files, "First Round")
    second_round_depth = load_depth_summary(second_round_files, "Second Round")

    # Combine both rounds
    depth_data = pd.concat([first_round_depth, second_round_depth])

    # Plot depth distributions
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="Sequencing Round", y="MedianDepth", data=depth_data)
    sns.stripplot(x="Sequencing Round", y="MedianDepth", data=depth_data, color=".3", alpha=0.5, jitter=True)  # Add individual points
    plt.xlabel("Sequencing Round")
    plt.ylabel("Median Depth per Scaffold")
    plt.title("Depth Distribution by Sequencing Round")
    plt.grid(True)

    # Save plot
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    print(f"Plot saved to {output_plot}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python plot_depths_by_round.py <first_round_depth_files> <second_round_depth_files> <output_plot>")
        sys.exit(1)

    first_round_files = sys.argv[1]
    second_round_files = sys.argv[2]
    output_plot = sys.argv[3]

    main(first_round_files, second_round_files, output_plot)
