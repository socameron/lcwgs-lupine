import pandas as pd

# Define the mapping of sample prefixes (old levels) to population codes (new levels)
prefix_to_population = {
    "ARCA": "E.MI.3",
    "RLPLV": "E.ON.8",
    "PPP": "E.ON.9",
    "HPW": "E.ON.5",
    "LW-KBS": "E.ON.12",
    "IDNP-MW": "C.IN.2",
    "LCTGP": "C.OH.1",
    "MFNP": "C.MI.2",
    "NCBW-FOT": "C.ON.1",
    "SWCP": "C.ON.3"
}

# Read the list of BAM file paths from a text file (assuming you saved it as "bam_list.txt")
with open("data/lists/hap2/all_populations_realign_hap2.txt", "r") as file:
    bam_samples = [line.strip() for line in file]

# Extract sample information
sample_data = []
for bam in bam_samples:
    # Extract the sample name from the path (e.g., "SWCP-25" from "results/bam_realign/hap2/SWCP-25_hap2_realign.bam")
    sample_name = bam.split("/")[-1].split("_hap2_realign.bam")[0]  # Isolate the sample name
    
    # Split the sample name into parts using "-" and reconstruct the prefix if it contains more than one part
    parts = sample_name.split("-")
    sample_prefix = "-".join(parts[:2]) if parts[0] + "-" + parts[1] in prefix_to_population else parts[0]
    
    # Map the sample prefix to the population code
    population_code = prefix_to_population.get(sample_prefix, "Unknown")
    
    # Store the sample name and its corresponding population code
    sample_data.append([population_code, sample_name])

# Create a DataFrame from the sample data
df = pd.DataFrame(sample_data, columns=["PopulationCode", "Sample"])

# Save the matched samples to a tab-delimited .info file
df.to_csv("data/lists/Batch1_popln_to_samples.info", sep="\t", index=False, header=False)

# Print the result
print(df)
