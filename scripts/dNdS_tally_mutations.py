import os
import pandas as pd
import argparse

# Define the directory containing .subst files
parser = argparse.ArgumentParser(description="Process input directory and output file")
parser.add_argument("--input_dir", required=True, help="Path to the input directory with .subst files")
parser.add_argument("--output_file", required=True, help="Path to the output Excel file")
args = parser.parse_args()

input_directory = args.input_dir
output_file = args.output_file

# Create a list to hold the data
data = []
synon_sum = 0
nonsyn_sum = 0
missen_sum = 0  # Variable to count MISSEN mutations

# Loop through all files in the directory
for filename in os.listdir(input_directory):
    if filename.endswith(".subst"):
        file_path = os.path.join(input_directory, filename)
        with open(file_path, 'r') as file:
            content = file.read().strip()
            # Split the content into individual columns based on the colon separator
            columns = []
            for part in content.split(', '):
                split_part = part.split(': ')
                if len(split_part) == 2:
                    columns.append(split_part[1])
                else:
                    columns.append(split_part[0])  # Add the part itself if it doesn’t follow "key: value" format
            # Add the filename as the first column
            data.append([filename] + columns)

            # Add to the appropriate sum based on mutation type
            frequency = float(columns[3])  # Frequency is the 4th column (index 3)
            if columns[-1] == 'SYNON':
                synon_sum += frequency
            elif columns[-1] == 'NONSYN':
                nonsyn_sum += frequency
            elif columns[-1] == 'MISSEN':  # New condition for MISSEN mutations
                missen_sum += frequency

# Define column names
column_names = ['Filename', 'Position', 'Ref', 'Alts', 'Frequency', 'Original Codon', 'Alternative Codon', 'AA Change', 'Change Type']

# Convert the data to a DataFrame
df = pd.DataFrame(data, columns=column_names)

# Ensure the Frequency column is numeric
df['Frequency'] = pd.to_numeric(df['Frequency'])

# Save the DataFrame to an CSV file
df.to_csv(output_file, index=False)


# Print the sums of frequencies
print(f"Sum of frequencies for SYNON: {synon_sum}")
print(f"Sum of frequencies for NONSYN: {nonsyn_sum}")
print(f"Sum of frequencies for MISSEN: {missen_sum}")

print(f"Excel file created at: {output_file}")
