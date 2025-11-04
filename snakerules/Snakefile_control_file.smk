###########################
##    WILDCARDS SET-UP   ##
###########################

import os
import shlex
import glob
import subprocess
import pandas as pd
import yaml
import re
import itertools
from pathlib import Path






#### Functions to set up samples list in config from different sequencing batches

# In summary:
  #{sample} refers to the sample name received from the sequencer. It includes replicates and other identifying info.
  #{sample_prefix} refers to a trimmed sample name and merges replicates.

# Rule to update the configuration file
# The python script creates a list of sample names and their corresponding .fastq files for both reads R1 and R2
rule update_config:
  output:
    "snakeprofile/config.yaml"
  shell:
    "python scripts/extract_lcWGS_samples.py"

# Function to get samples from the config file
def get_samples(batch):
    with open("snakeprofile/config.yaml", 'r') as file:
        config = yaml.safe_load(file)
    return config['batches'][batch]

# Function to get batches from the config file
def get_batches():
    with open("snakeprofile/config.yaml", 'r') as file:
        config = yaml.safe_load(file)
    return list(config['batches'].keys())

# Function to get .fastq filenames from all 3 batches from the config file
def get_fastq_paths(batch, sample, read):
    with open("snakeprofile/config.yaml", 'r') as file:
        config = yaml.safe_load(file)
    if batch in config['batches'] and sample in config['batches'][batch]:
        return config['batches'][batch][sample].get(read, f"data/{batch}/{sample}_L001_{read}_001.fastq.gz")
    # Fallback for batch_1 and batch_2:
    return f"data/{batch}/{sample}_L001_{read}_001.fastq.gz"

# Ensure that the config is updated before any other rule runs
# Get all batches from the config
batches = get_batches()

# Get samples for each batch
samples_batch_1 = get_samples('batch_1')
samples_batch_2 = get_samples('batch_2')
samples_batch_3 = get_samples('batch_3')

# Try using this for expand in rule all: 
  #  [f"results/bam_raw/hap2/{batch}/{sample}_hap2.bam"
  #   for batch in get_batches()
  #   for sample in list(get_samples(batch).keys())],

# Get basename of the fastq files
def get_basename(batch, sample, read):
    fastq_path = get_fastq_paths(batch, sample, read)
    if fastq_path:
        return os.path.basename(fastq_path).replace(".fastq.gz", "")
    return f"{sample}_{read}"  # Fallback in case of missing file


# Create list for fastqc outputs
expected_fastqc_outputs = []

for batch in get_batches():
    for sample in get_samples(batch).keys():
        # Compute the basenames from the full FASTQ paths.
        # This should always return a non-None string.
        r1_base = get_basename(batch, sample, "R1")
        r2_base = get_basename(batch, sample, "R2")

        expected_fastqc_outputs.extend([
            f"results/fastqc_raw/{batch}/{r1_base}_fastqc.html",
            f"results/fastqc_raw/{batch}/{r1_base}_fastqc.zip",
            f"results/fastqc_raw/{batch}/{r2_base}_fastqc.html",
            f"results/fastqc_raw/{batch}/{r2_base}_fastqc.zip"
        ])


#### Functions to set up lists of sample prefixes for merging .BAM files from different sequencing batches

# List of haplotypes
haps = ["1", "2"]

# List of populations
POPULATIONS = ["HPW", "IDNP-MW", "LCTGP", "MFNP", "PPP", "RLPLV", "SWCP", "APB", "BSL", "BSNA", "CPB", "FMB", "GRAY", "NBWA", "NGP", "PBBT", "RSB", "UWA"]

def extract_sample_prefix(sample_name):
    for pop in POPULATIONS:
        if sample_name.startswith(pop + "-"):
            # Get the portion before the first underscore.
            prefix = sample_name.split("_")[0]
            # If the prefix includes a replicate suffix, remove it.
            if "-rep" in prefix:
                prefix = prefix.split("-rep")[0]
            return prefix
    return sample_name

# Find replicates across BAM files as inputs into merge_replicates rule
def find_replicates(sample_prefix, hap):
    """
    Finds all BAM files corresponding to a given sample (using the unique sample_prefix)
    across batches (batch_1, batch_2, batch_3). The function expects the BAM file names to
    start exactly with sample_prefix and optionally a "-rep<number>" (for replicates), 
    immediately followed by an underscore and then "hap{hap}".
    """
    # Broad glob pattern to capture candidates
    pattern = f"results/bam_raw/hap{hap}/batch*/{sample_prefix}*_hap{hap}.bam"
    files = glob.glob(pattern, recursive=True)
    
    # Compile a regex to match:
    #   ^                 --> start of the filename
    #   sample_prefix     --> exact sample prefix (escaped)
    #   (?!\d)           --> not immediately followed by a digit (so "PPP-1" won't match "PPP-10")
    #   (?:-rep\d+)?     --> optionally match a replicate indicator like "-rep2"
    #   _hap             --> then an underscore and "hap..."
    regex = re.compile(r'^' + re.escape(sample_prefix) + r'(?:-rep\d+)?_')

    
    # Filter files based on the regex match on the base name
    filtered_files = [f for f in files if regex.search(os.path.basename(f))]
    
    # Write a debug file for inspection
    debug_path = f"debug/debug_find_replicates_{sample_prefix}_hap{hap}.txt"
    os.makedirs(os.path.dirname(debug_path), exist_ok=True)
    with open(debug_path, "w") as f:
        f.write(f"Pattern: {pattern}\n")
        if filtered_files:
            f.write("Found files:\n")
            for file in filtered_files:
                f.write(f"{file}\n")
        else:
            f.write("No files found.\n")
    
    return filtered_files

# Get all unique sample prefixes by extracting from samples listed in config file
def list_sample_prefixes():
    samples_batch_1 = get_samples('batch_1')  # returns a dict of {sample_name: {...}}
    samples_batch_2 = get_samples('batch_2')
    samples_batch_3 = get_samples('batch_3')
    # Combine sample names from all batches
    all_samples = set(
        list(samples_batch_1.keys()) +
        list(samples_batch_2.keys()) +
        list(samples_batch_3.keys())
    )
    # Now extract the prefix from each sample name using your extract_sample_prefix function.
    sample_prefixes = {extract_sample_prefix(sample) for sample in all_samples}
    return list(sample_prefixes)

# List of sample prefixes
sample_prefixes = list_sample_prefixes()
   
# Use sample prefixes and designate to a particular population
def get_population_sample_prefixes(population):
    # Get all sample prefixes from the existing list_sample_prefixes function
    all_sample_prefixes = list_sample_prefixes()
    # Filter the sample prefixes based on the population
    population_sample_prefixes = [prefix for prefix in all_sample_prefixes if prefix.startswith(population)]
    return population_sample_prefixes

#### Functions to create scaffold names list per scaffold
# Get first 24 scaffold names to later estimate LD per scaffold
# Not using other scaffolds as they don't quite align to the expected 24 chromosomes
# Use the actual python script below to create it.

rule get_24_scaffold_names:
  input:
    "data/reference/lupinehap1.fasta",
    "data/reference/lupinehap2.fasta"
  output:
    "results/scaffolds/hap1_24_scaffolds.txt",
    "results/scaffolds/hap2_24_scaffolds.txt"
  shell:
    """
    python scripts/extract_24_scaffold_names_by_hap.py
    """

rule get_all_scaffold_names:
  input:
    "data/reference/lupinehap1.fasta",
    "data/reference/lupinehap2.fasta"
  output:
    "results/scaffolds/hap1_all_scaffolds.txt",
    "results/scaffolds/hap2_all_scaffolds.txt"
  shell:
    """
    python scripts/extract_all_scaffold_names_by_hap.py
    """

with open("results/scaffolds/hap1_24_scaffolds.txt", "r") as file:
    HAP1SCAFFOLDS = [line.strip() for line in file.readlines()]

with open("results/scaffolds/hap2_24_scaffolds.txt", "r") as file:
    HAP2SCAFFOLDS = [line.strip() for line in file.readlines()]

# Split scaffold names by comma and create a list
HAP1SCAFFOLDS = [name.strip() for name in HAP1SCAFFOLDS]
HAP2SCAFFOLDS = [name.strip() for name in HAP2SCAFFOLDS]
HAP1SCAFFOLD_PREFIXES = [s.split("__")[0] for s in HAP1SCAFFOLDS]
HAP2SCAFFOLD_PREFIXES = [s.split("__")[0] for s in HAP2SCAFFOLDS]

# Function give full scaffold name given prefix of scaffold
def map_prefix_to_full_scaffold(prefix, hap_type):
    scaffold_list = HAP1SCAFFOLDS if hap_type == 1 else HAP2SCAFFOLDS
    for scaffold in scaffold_list:
        if scaffold.startswith(prefix):
            return scaffold
    return None  
# Return None or an appropriate default if not found



#### Bunch of parameters for using ngsParalog, Fst, ANGSD, etc

# Function to check if bams exist
def bams_to_fix():
    pops = ["HPW","LCTGP","IDNP-MW","RLPLV","PPP","MFNP","SWCP","APB","BSL","BSNA","CPB","FMB","GRAY","NBWA","NGP","PBBT","RSB","UWA"]
    bams = []
    for pop in pops:
        bl = f"data/lists/hap2/{pop}_clipped_hap2.txt"
        if not os.path.exists(bl): continue
        with open(bl) as f:
            for line in f:
                p = line.strip()
                if p: bams.append(p)
    return bams



###########################
## BEGINNING OF WORKFLOW ##
###########################


# Expand the final files using the updated configuration
rule all:
  input:
    expand("results/plots/hap2/depths/all_rounds/{hap2scaffold_prefix}_depth_histogram.png", hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES),
    "data/lists/hap2/all_populations_realign_hap2.txt"

    # inbreeding coefficents
    #expand("results/ngsF/hap2/by_popln/{population}_ngsF_inbreeding_final.lrt", population = POPULATIONS),
    #"results/plots/hap2/ngsF/ngsF-HMM_inbreeding_coeff.tiff"






    # currentNe estK
    #expand("results/currentNe2/hap2/{population}/{population}_estK_currentNe2_OUTPUT.txt", population=POPULATIONS),
    #expand("results/currentNe2/hap2/{population}/{population}_estK_currentNe2_mix_OUTPUT.txt", population=POPULATIONS),
    # fixed-k (now includes k=0)
    #expand("results/currentNe2/hap2/{population}/{population}_k{k}_currentNe2_OUTPUT.txt", population=POPULATIONS, k=K_GRID),
    #expand("results/currentNe2/hap2/{population}/{population}_k{k}_currentNe2_mix_OUTPUT.txt", population=POPULATIONS, k=K_GRID),
    # -x
    #expand("results/currentNe2/hap2/{population}/{population}_x_currentNe2_OUTPUT.txt", population=POPULATIONS)


include: "snakerules/00_lcwgs_dataprep.smk"
include: "snakerules/01_ngsRelate_relatedness.smk"
include: "snakerules/02_ngsF-HMM_individual_inbreeding.smk"
include: "snakerules/03_RZooROH_runs_of_homozygosity.smk"
include: "snakerules/04_PCAngsd_population_structure.smk"
include: "snakerules/05_ANGSD_diversity_stats.smk"
include: "snakerules/06_ANGSD_Fst_and_IBD.smk"
include: "snakerules/07_ngsLD_linkage_disequilibrium.smk"
include: "snakerules/08_currentNe2_contemporary_effective_popln_size.smk"
include: "snakerules/09_fastDFE_distribution_fitness_effects.smk"
