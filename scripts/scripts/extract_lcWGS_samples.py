import os
import glob
import yaml

def extract_sample_name(filename):
    base_filename = os.path.basename(filename)
    parts = base_filename.split('_')
    for i, part in enumerate(parts):
        if 'L001' in part:
            return '_'.join(parts[:i])
    return None

batch_directories = ["data/batch_1", "data/batch_2"]
batch_samples = {}

for directory in batch_directories:
    batch_name = os.path.basename(directory)
    fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz"))
    samples = set()
    for filename in fastq_files:
        sample_name = extract_sample_name(filename)
        if sample_name:
            samples.add(sample_name)
    batch_samples[batch_name] = sorted(list(samples))

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
