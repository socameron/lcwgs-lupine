from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Check transcripts for dN/dS analysis")
parser.add_argument("--input_file", required=True, help="Path to the input FASTA file")
parser.add_argument("--output_file_with_stops", required=True, help="Path to output file for transcripts with premature stop codons")
parser.add_argument("--output_file_without_stops", required=True, help="Path to output file for transcripts without premature stop codons")
parser.add_argument("--output_file_protein", required=True, help="Path to output file for translated proteins")
parser.add_argument("--output_file_saved_transcripts", required=True, help="Path to output file for saved transcripts")
args = parser.parse_args()

# Use the provided arguments in place of hard-coded file paths
input_file = args.input_file
output_file_with_stops = args.output_file_with_stops
output_file_without_stops = args.output_file_without_stops
output_file_protein = args.output_file_protein
output_file_saved_transcripts = args.output_file_saved_transcripts

# Define stop codons for the standard genetic code
stop_codons = ["TAA", "TAG", "TGA"]

# Initialize lists to store transcripts with and without premature stop codons
premature_stop_transcripts = []
no_premature_stop_transcripts = []
translated_protein_sequences = []
saved_transcripts = []
premature_stop_count = 0

# Process each sequence in the FASTA file
for record in SeqIO.parse(input_file, "fasta"):
    # Get the sequence of the transcript
    seq = record.seq
    found_stop = False
    
    # Iterate through codons, stopping before the final codon
    for i in range(0, len(seq) - 3, 3):
        codon = str(seq[i:i+3])
        if codon in stop_codons:
            # Found a premature stop codon
            premature_stop_transcripts.append(record)
            premature_stop_count += 1
            found_stop = True
            break
    
    # If no premature stop was found, add to the no_premature_stop_transcripts list
    if not found_stop:
        no_premature_stop_transcripts.append(record)
        
        # Translate the full transcript sequence into protein
        protein_seq = seq.translate(to_stop=False)  # Translate the entire sequence
        
        # Check for conditions: no internal stops and starts with 'M'
        if '*' not in protein_seq[:-1] and protein_seq[0] == 'M':
            # Save the protein sequence and corresponding transcript
            protein_record = record[:]
            protein_record.seq = protein_seq
            translated_protein_sequences.append(protein_record)
            saved_transcripts.append(record)

# Write the output FASTA files with line width 0 (single line for sequences)
with open(output_file_with_stops, "w") as output_handle_with_stops:
    writer = FastaWriter(output_handle_with_stops, wrap=None)
    writer.write_file(premature_stop_transcripts)

with open(output_file_without_stops, "w") as output_handle_without_stops:
    writer = FastaWriter(output_handle_without_stops, wrap=None)
    writer.write_file(no_premature_stop_transcripts)

with open(output_file_protein, "w") as output_handle_protein:
    writer = FastaWriter(output_handle_protein, wrap=None)
    writer.write_file(translated_protein_sequences)

with open(output_file_saved_transcripts, "w") as output_handle_saved_transcripts:
    writer = FastaWriter(output_handle_saved_transcripts, wrap=None)
    writer.write_file(saved_transcripts)

# Print the results
print(f"Transcripts with premature stop codons have been saved to '{output_file_with_stops}'.")
print(f"Number of transcripts with premature stop codons: {premature_stop_count}")
print(f"Transcripts without premature stop codons have been saved to '{output_file_without_stops}'.")
print(f"Translated protein sequences (starting with 'M' and with no internal stops) have been saved to '{output_file_protein}'.")
print(f"Corresponding transcripts have been saved to '{output_file_saved_transcripts}'.")
print(f"Number of valid translated proteins and corresponding transcripts: {len(translated_protein_sequences)}")
