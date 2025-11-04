def extract_scaffolds(fasta_file, output_file):
    """Extract all scaffold names from a FASTA file."""
    scaffolds = []

    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                scaffold_name = line.split()[0][1:]  # Get header, remove '>'
                scaffolds.append(scaffold_name)

    # Write all scaffold names to the output file
    with open(output_file, "w") as out_file:
        for scaffold in scaffolds:
            out_file.write(f"{scaffold}\n")

# Example usage
fasta_file_hap1 = "data/reference/hap1/lupinehap1.fasta"
fasta_file_hap2 = "data/reference/hap2/lupinehap2.fasta"
output_file_hap1 = "results/scaffolds/hap1_all_scaffolds.txt"
output_file_hap2 = "results/scaffolds/hap2_all_scaffolds.txt"

extract_scaffolds(fasta_file_hap1, output_file_hap1)
extract_scaffolds(fasta_file_hap2, output_file_hap2)
