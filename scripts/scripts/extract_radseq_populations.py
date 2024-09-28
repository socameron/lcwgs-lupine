# extract_radseq_populations.py

import re
import sys

def extract_populations(sample_file):
    populations = {}
    with open(sample_file, 'r') as file:
        for line in file:
            sample_path = line.strip()
            sample_name = sample_path.split('/')[-1]
            match = re.match(r'sorted\.trim\.([A-Z0-9-]+)-.*\.p', sample_name)
            if match:
                population = match.group(1)
                if population not in populations:
                    populations[population] = []
                populations[population].append(sample_name)
    return populations

if __name__ == "__main__":
    sample_file = sys.argv[1]
    populations = extract_populations(sample_file)
    for population, samples in populations.items():
        with open(f"{population}.txt", "w") as f:
            for sample in samples:
                f.write(sample + "\n")
    # Print populations for Snakemake
    for population in populations:
        print(population)
