def extract_scaffolds(fasta_file, output_file, max_scaffolds=24):
    scaffolds = []
    
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                scaffold_name = line.split()[0][1:]  # Extract the scaffold name
                scaffolds.append(scaffold_name)
                
                # Stop if we have enough scaffolds
                if len(scaffolds) >= max_scaffolds:
                    break
    
    with open(output_file, "w") as out_file:
        for scaffold in scaffolds:
            out_file.write(f"{scaffold}\n")

# Example usage
fasta_file_hap1 = "data/reference/hap1/lupinehap1.fasta"
fasta_file_hap2 = "data/reference/hap2/lupinehap2.fasta"
output_file_hap1 = "results/scaffolds/hap1_scaffolds.txt"
output_file_hap2 = "results/scaffolds/hap2_scaffolds.txt"

extract_scaffolds(fasta_file_hap1, output_file_hap1)
extract_scaffolds(fasta_file_hap2, output_file_hap2)