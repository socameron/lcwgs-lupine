import concurrent.futures
import os
import logging
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Process output directory and input files")
parser.add_argument("--output_dir", required=True, help="Path to the output directory")
parser.add_argument("--transcripts", required=True, help="Path to the transcripts FASTA file")
parser.add_argument("--transcript_proteins", required=True, help="Path to the transcript proteins FASTA file")
parser.add_argument("--vcf", required=True, help="Path to the VCF file")
args = parser.parse_args()

# Define paths from command-line arguments
OUTPUT_DIR = args.output_dir
transcripts_file = args.transcripts
proteins_file = args.transcript_proteins
vcf_file = args.vcf

# Create the output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Using output directory: {OUTPUT_DIR}")
print(f"Transcripts file: {transcripts_file}")
print(f"Proteins file: {proteins_file}")
print(f"VCF file: {vcf_file}")

# Configure logging
logging.basicConfig(filename=os.path.join(OUTPUT_DIR, 'script_debug.log'), level=logging.DEBUG, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Global dictionary to keep track of .fa file counts for each transcript_id
fa_file_counts = {}

# Global dictionary for tracking protein variant frequencies
protein_variant_frequencies = {}
# Global dictionary for tracking the total counts of alleles processed for each transcript variant
total_alleles_count = {}

def read_fasta(filepath):
    sequences = {}
    with open(filepath, 'r') as file:
        sequence_id = None
        sequence_lines = []
        for line in file:
            if line.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] = ''.join(sequence_lines).replace('\n', '')
                sequence_id = line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(line.strip())
        if sequence_id:
            sequences[sequence_id] = ''.join(sequence_lines).replace('\n', '')
    return sequences

def parse_vcf_line(line):
    """Parse a single VCF line and return relevant info."""
    parts = line.strip().split('\t')
    chrom, pos, _, ref, alts, _, _, _, format, *samples = parts
    return chrom, int(pos), ref, alts.split(','), format, samples

def translate_seq(nucleotide_seq):
    """Translate nucleotide sequence to amino acid sequence using the standard genetic code."""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein_seq = ''
    for i in range(0, len(nucleotide_seq) - 2, 3):
        codon = nucleotide_seq[i:i + 3]
        protein_seq += codon_table.get(codon, '?')
    return protein_seq

def process_variant(chrom, pos, ref, alt_alleles, format, samples, transcript_id, transcript_seq, protein_seq):
    logging.debug(f"Processing variant: chrom={chrom}, pos={pos}, ref={ref}, alt_alleles={alt_alleles}")
    
    # Determine if the position is polymorphic by scanning all genotypes
    alleles = set()
    for sample in samples:
        genotypes = sample.split(':')[0].replace('|', '/').split('/')
        alleles.update(genotypes)
    
    if len(alleles) <= 1:  # Not polymorphic
        return
    
    alt_allele_present = False
    alt_allele_count = 0
    total_genotypes = -2

    for sample in samples:
        genotypes = sample.split(':')[0].replace('|', '/').split('/')
        alt_allele_count += genotypes.count('1')
        total_genotypes += len([g for g in genotypes if g != '.'])

    variant_frequency = alt_allele_count / ((total_genotypes) * 2) if total_genotypes > 0 else 0

    if alt_allele_count == 0:  # No alternative alleles present
        return

    # Continue with the original processing
    fa_file_path = os.path.join(OUTPUT_DIR, f"{transcript_id}.fa")
    with open(fa_file_path, 'w') as fa_file:
        fa_file.write(f">{transcript_id}\n{transcript_seq}\n")

    for alt_index, alt in enumerate(alt_alleles):
        if alt == '*':
            continue
        alt_nucleotide_seq = transcript_seq[:pos - 1] + alt + transcript_seq[pos:]
        alt_protein_seq = translate_seq(alt_nucleotide_seq)

        # Calculate the protein variant title
        fa_file_title = f"{transcript_id}_{pos}_{ref}>{alt}"

        if fa_file_title not in protein_variant_frequencies:
            protein_variant_frequencies[fa_file_title] = alt_allele_count
            total_alleles_count[fa_file_title] = total_genotypes * 2
        else:
            protein_variant_frequencies[fa_file_title] += alt_allele_count
            total_alleles_count[fa_file_title] += total_genotypes * 2
        
        # Calculate the original and alternative amino acids and their positions
        ref_codon_start = (pos - 1) - (pos - 1) % 3
        ref_codon = transcript_seq[ref_codon_start:ref_codon_start + 3]
        alt_codon = alt_nucleotide_seq[ref_codon_start:ref_codon_start + 3]
        ref_aa = translate_seq(ref_codon)
        alt_aa = translate_seq(alt_codon)
        protein_pos = (ref_codon_start // 3) + 1
        aa_change = f"{alt_aa}{protein_pos}{ref_aa}"

        # Determine if the change is synonymous, non-synonymous, or a stop codon (MISSEN)
        if ref_aa == alt_aa:
            change_type = "SYNON"
        elif alt_aa == '*':
            change_type = "MISSEN"
        else:
            change_type = "NONSYN"

        # Write .subst file
        subst_file_path = os.path.join(OUTPUT_DIR, f"{transcript_id}_{pos}.subst")
        
        logging.debug(f"Writing .fa file: {fa_file_path}")
        logging.debug(f"Writing .subst file: {subst_file_path}")
        
        with open(subst_file_path, 'w') as subst_file:
            subst_file.write(f"Position: {pos}, Ref: {ref}, Alts: {', '.join(alt_alleles)}, Frequency: {2*variant_frequency:.4f}, Original Codon: {ref_codon}, Alternative Codon: {alt_codon}, AA Change: {aa_change}, {change_type}\n")


def process_variants(vcf_file, transcripts, proteins):
    with open(vcf_file, 'r') as file, concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        futures = []
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines
            chrom, pos, ref, alt_alleles, format, samples = parse_vcf_line(line)
            
            transcript_id = chrom
            transcript_seq = transcripts.get(transcript_id, '')
            protein_seq = proteins.get(transcript_id, '')
            
            if transcript_seq and protein_seq:
                futures.append(executor.submit(process_variant, chrom, pos, ref, alt_alleles, format, samples, transcript_id, transcript_seq, protein_seq))
        
        concurrent.futures.wait(futures)

# Load your data
transcripts = read_fasta(transcripts_file)
proteins = read_fasta(proteins_file)

# Process variants
process_variants(vcf_file, transcripts, proteins)

logging.info("Script completed successfully.")
