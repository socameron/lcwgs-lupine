import os
import glob
import yaml
import re

def extract_sample_name(filename):
    base_filename = os.path.basename(filename)

    # Ignore checksum (.md5) files
    if base_filename.endswith(".md5"):
        return None, None

    # Batch_1 & Batch_2: Extract sample name based on "_L001"
    parts = base_filename.split('_')
    for i, part in enumerate(parts):
        if 'L001' in part:
            clean_name = '_'.join(parts[:i])  # Example: HPW-10_2-2695941_S13
            return clean_name, clean_name  # Keep the same for compatibility

    # Batch_3: Extract both number and letter-based sample names
    match = re.search(r'NEBNext_dual_i5_(\w+\d*)\.(\w+-\d+(?:-rep\d+)?)', base_filename)
    if match:
        clean_name = match.group(2)  # Extracts the sample name, e.g., `NBWA-1`, `PBBT-3`
        return clean_name, clean_name  # ✅ Store only clean sample name

    return None, None


batch_directories = ["data/batch_1", "data/batch_2", "data/batch_3"]
batch_samples = {}

for directory in batch_directories:
    batch_name = os.path.basename(directory)
    fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz"))
    samples = {}

    for filename in fastq_files:
        clean_name, _ = extract_sample_name(filename)
        if clean_name:
            if clean_name not in samples:
                samples[clean_name] = {"R1": None, "R2": None}

            # Assign correct R1/R2 files
            if "_R1" in filename:
                samples[clean_name]["R1"] = filename
            elif "_R2" in filename:
                samples[clean_name]["R2"] = filename

    batch_samples[batch_name] = samples

config_path = "snakeprofile/config.yaml"

if os.path.exists(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
else:
    config = {}

if 'batches' not in config:
    config['batches'] = {}

updated = False
for batch_name, samples in batch_samples.items():
    if batch_name not in config['batches'] or config['batches'][batch_name] != samples:
        config['batches'][batch_name] = samples
        updated = True

if updated:
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    print(f"Updated {config_path} with batch sample names.")
else:
    print(f"No changes made to {config_path}.")
