import os
import glob

# Define a function to extract the sample name from a fastq.gz filename
def extract_sample_name(filename):
    # Get the base filename without the directory
    base_filename = os.path.basename(filename)
    
    # Split the filename by underscores
    parts = base_filename.split('_')
    
    # Find the part that contains 'L001' and get the sample name before it
    for i, part in enumerate(parts):
        if 'L001' in part:
            return '_'.join(parts[:i])

# Defining wildcard for haplotype and sequence batch #
wildcard_constraints:
    hap = r'1|2',
    batch = r'1|2'

# Directory containing the fastq.gz files 
# Change 1 to {batch}
directory = "data/batch_1"

# Get a list of all fastq.gz files in the directory
fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz"))

# Initialize a set to store unique sample names
SAMPLES = set()

# Iterate over the fastq.gz files and extract sample names
for filename in fastq_files:
    sample_name = extract_sample_name(filename)
    SAMPLES.add(sample_name)

# Convert the set to a sorted list for consistent order
SAMPLES = sorted(list(SAMPLES))
POPULATIONS = set(sample.split('-')[0] for sample in SAMPLES)








###########################
## BEGINNING OF WORKFLOW ##
###########################



# final files desired (can change)
rule all:
  input:
    expand("results/angsd/hap1/ngs_prep_{population}", population=POPULATIONS)



# Trim adapter ends off each sequence file using Trimmomatic 
# When batch 2 is available, use wildcard {batch}
rule trim_reads:
  input:
    r1="data/batch_{batch}/{sample}_L001_R1_001.fastq.gz",
    r2="data/batch_{batch}/{sample}_L001_R2_001.fastq.gz",
  output:
    r1="results/trimmed/{sample}_R1.fastq.gz",
    r1_unp="results/trimmed/{sample}_R1_unpaired.fastq.gz",
    r2="results/trimmed/{sample}_R2.fastq.gz",
    r2_unp="results/trimmed/{sample}_R2_unpaired.fastq.gz"
  log:
    "results/logs/trim_reads/{sample}.log"
  envmodules:
    "trimmomatic/0.39"
  params:
    adapters="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa"
  shell:
    "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE {input.r1} {input.r2} "
    "{output.r1} {output.r1_unp} {output.r2} {output.r2_unp} "
    "ILLUMINACLIP:{params.adapters}:2:30:10 "
    "LEADING:3 "
    "TRAILING:3 "
    "SLIDINGWINDOW:4:15 "
    "MINLEN:36 2> {log}"



# Unzip trimmed files as fastqc 0.12.0 cannot properly read compressed files
rule unzip_files:
  input:
    zipped_r1="results/batch_{batch}/trimmed/{sample}_R1.fastq.gz",
    zipped_r2="results/batch_{batch}trimmed/{sample}_R2.fastq.gz"
  output: 
    unzipped_r1="results/batch_{batch}/trimmed_unzip/{sample}_R1.fastq",
    unzipped_r2="results/batch_{batch}/trimmed_unzip/{sample}_R2.fastq"
  shell:
    "gunzip -c {input.zipped_r1} > {output.unzipped_r1} && gunzip -c {input.zipped_r2} > {output.unzipped_r2}"



# Run FastQC per each trimmed sequence file
rule fastqc:
  input:
    fastq_r1="results/trimmed_unzip/batch_{batch}/{sample}_R1.fastq",
    fastq_r2="results/trimmed_unzip/batch_{batch}/{sample}_R2.fastq"
  output:
    html_report_r1="results/fastqc/batch_{batch}/{sample}_R1_fastqc.html",
    zip_report_r1="results/fastqc/batch_{batch}/{sample}_R1_fastqc.zip",
    html_report_r2="results/fastqc/batch_{batch}/{sample}_R2_fastqc.html",
    zip_report_r2="results/fastqc/batch_{batch}/{sample}_R2_fastqc.zip"
  envmodules:
    "fastqc/0.12.0"
  shell:
    "fastqc {input.fastq_r1} {input.fastq_r2} --outdir results/fastqc"



# Create an aggregated FastQC report using MultiQC.
# Note that we create separate MultiQC reports for both batch 1 and 2
rule multiqc_raw:
  input:
    fastqc_dir="results/batch_{batch}/fastqc"
  output:
    html_report="results/multiqc/multiqc_raw_batch_{batch}/multiqc_raw.html"
  params:
    output_dir="results/multiqc/multiqc_raw_batch_{batch}"
  shell:
    "multiqc -o {params.output_dir} {input.fastqc_dir} "



# Creating faidx index for reference genome
rule faidx_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta",
  output:
    "data/reference/hap{hap}/lupinehap{hap}.fasta.fai",
  log:
    "results/logs/refgen/lupinehap{hap}_faidx.log",
  envmodules:
    "samtools/1.17"
  shell:
    "samtools faidx {input} 2> {log} "



# Rules for indexing reference genomes (haplotypes 1 and 2)
rule index_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
  log:
    "results/logs/refgen/lupinehap{hap}_bwa_index.log"
  envmodules:
    "bwa/0.7.17"
  shell:
    "bwa index {input} 2> {log}"




# Rules for creating dictionary files
rule create_dict:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.dict"
  log:
    "results/logs/refgen/hap{hap}_dict.log"
  envmodules:
    "samtools/1.17"
  shell:
    "samtools dict {input} > {output} 2> {log}"



# Mapping/Aligning reads to haplotypes
rule map_reads:
  input:
    r1="results/trimmed/{sample}_R1.fastq.gz",
    r2="results/trimmed/{sample}_R2.fastq.gz",
    genome="data/reference/hap{hap}/lupinehap{hap}.fasta",
    idx=multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  output:
    "results/bam/hap{hap}/{sample}_hap{hap}.bam"
  log:
    "results/logs/map_reads/hap{hap}/{sample}_hap{hap}.log"
  benchmark:
    "results/benchmarks/map_reads/{sample}_hap{hap}.bmk"
  envmodules:
    "bwa/0.7.17",
    "samtools/1.17"
  threads: 8 
  params:
    RG="-R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA' "
    # Excluded information on flowcell, lane, and barcode index b/c not thought necessary
  shell:
    " (bwa mem {params.RG} -t {threads} {input.genome} {input.r1} {input.r2} | "
    " samtools view -u | "
    " samtools sort - > {output}) 2> {log}"



# Mark duplicates 
rule mark_duplicates:
  input:
    "results/bam/hap{hap}/{sample}_hap{hap}.bam"
  output:
    bam="results/bam_mkdup/hap{hap}/{sample}_hap{hap}_mkdup.bam",
    bai="results/bam_mkdup/hap{hap}/{sample}_hap{hap}_mkdup.bai",
    metrics="results/qc/mkdup_metrics/{sample}_hap{hap}.metrics"
  log:
    "results/logs/mark_duplicates/hap{hap}/{sample}_hap{hap}.log"
  envmodules:
    "gatk/4.2.5.0"
  shell:
    " gatk MarkDuplicates  "
    "  --CREATE_INDEX "
    "  -I {input} "
    "  -O {output.bam} "  
    "  -M {output.metrics} "
    "  2> {log} "



# Remove duplicates using PICARD (in GATK)
# Separate rule because we want to keep original bam files
rule remove_duplicates:
  input:
    bam="results/bam_mkdup/hap{hap}/{sample}_hap{hap}_mkdup.bam"
  output:
    filtered_bam="results/bam_filtered/hap{hap}/{sample}_hap{hap}_filtered.bam",
    filtered_metrics="results/qc/filtered_metrics/{sample}_hap{hap}.metrics"
  log:
    "results/logs/remove_duplicates/hap{hap}/{sample}_hap{hap}_remove_duplicates.log"
  params:
    java_opts="-Xmx8g"
  envmodules:
    "gatk/4.2.5.0"
  threads: 4
  shell:
    """
    gatk MarkDuplicates \
    --java-options "{params.java_opts}" \
    --REMOVE_DUPLICATES \
    --INPUT {input.bam} \
    --OUTPUT {output.filtered_bam} \
    --METRICS_FILE {output.filtered_metrics} \
    2> {log}
    """



# Create MultiQC report on duplicate metrics
rule multiqc_filtered:
  input:
    filtered_metrics=expand("results/qc/filtered_metrics/{sample}_hap{hap}.metrics", sample=SAMPLES, hap=[1, 2])
  output:
    html_report="results/multiqc/multiqc_filtered_batch_1/multiqc_filtered.html",
  params:
    output_dir="results/multiqc/multiqc_filtered_batch_1"
  shell:
    "multiqc -f -o {params.output_dir} {input.filtered_metrics}"



# Call SNPs using angsd
# Necessary so only running ngsParalog on variant sites; reduces computational time\
# Note: we do NOT filter marked duplicates (-remove_bads set to 0); set minimum base quality to 0, minMaf to 0, minMapQ to 0, and minimum individuals 1
rule ngs_prep_hap1:
  wildcard_constraints:
    population="|".join(POPULATIONS)
  input:
    bam=lambda wildcards: expand("results/bam_mkdup/hap1/{sample}_hap1_mkdup.bam", sample=[s for s in SAMPLES if s.startswith(wildcards.population)])
  output:
    angsd_output="results/angsd/hap1/ngs_prep_{population}"
  params:
    ref="data/reference/hap1/lupinehap1.fasta",
    threads=10
  log:
    "results/logs/angsd/angsd_ngs_prep_{population}.log"
  envmodules:
    "angsd/0.939"
  threads: 10
  shell:
    """
    angsd -b {input.bam} \
    -ref {params.ref} \
    -out {output.angsd_output} \
    -remove_bads 0 \
    -GL 2 \ # Genotype likelihoods based on Bayesian GATK models. Could swap for '1' for SAMtool models.
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMapQ 0 \
    -minQ 0 \
    -minInd 1 \ 
    -minMaf 0 \ 
    -threads {params.threads} \
    2> {log}
    """

# Create likelihood ratio of mismapping reads using ngsParalog
# Requires mpileup files - so will need to convert using samtools
rule ngs_paralog_hap1:
  input: 
    bam_mkdup=expand("results/bam_mkdup/hap1/{sample}_hap1_mkdup.bam", sample=[s for s in SAMPLES if s.startswith('HPW')])
  output:
    paralog_output="results/ngs_paralog/hap1_paralog_HPW.lr"
  log:
    "results/logs/ngs_paralog/hap1_paralog_HPW.log"
  envmodules:
    "samtools/1.17"
  shell: #pipe ngsParalog into samtools
    """ 
    samtools mpileup {input.bam_mkdup} -q 0 -Q 0 --ff DUP 2>> {log} | \
    /home/socamero/ngsParalog/ngsParalog calcLR -infile - -outfile {output.paralog_output} \
    -minQ 20 -minind 27 -mincov 1 \
    2>> {log}
    """


# If using filtered data, insert: results/filtered_bam/hap{hap}/{sample}_hap{hap}_filtered.bam
# Otherwise, can remove marked duplicates using -remove_bads
# Convert BAM files to .beagle file using angsd
rule angsd_hap1:
  input:
    bam=expand("results/bam_mkdup/hap1/{sample}_hap1_mkdup.bam", sample=SAMPLES)
  output:
    angsd_output="results/angsd/hap1/angsd_hap1"
  params:
    ref="data/reference/hap1/lupinehap1.fasta",
    threads=10
  log:
    "results/logs/angsd/angsd_hap1.log"
  envmodules:
    "angsd/0.939"
  threads: 10
  shell:
    """
    angsd -b {input.bam} \
    -ref {params.ref} \
    -out {output.angsd_output} \
    -remove_bads 1 \
    -GL 2 \ # Genotype likelihoods based on Bayesian GATK models. Could swap for '1' for SAMtool models.
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMapQ 30 \
    -minQ 20 \
    -minInd 1 \
    -minMaf 0.05 \
    -doGlf 2 \
    -doPost 1 \
    -doGeno 32 \  # Beagle output format
    -threads {params.threads} \
    2> {log}
    """



rule angsd_hap2:
  input:
    bam=expand("results/bam_mkdup/hap2/{sample}_hap2_mkdup.bam", sample=SAMPLES)
  output:
    angsd_output="results/angsd/hap2/angsd_hap2"
  params:
    ref="data/reference/hap2/lupinehap2.fasta",
    threads=10
  log:
    "results/logs/angsd/angsd_hap2.log"
  envmodules:
    "angsd/0.939"
  threads: 10
  shell:
    """
    angsd -b {input.bam} \
    -ref {params.ref} \
    -out {output.angsd_output} \
    -remove_bads 1 \
    -GL 2 \ # Genotype likelihoods based on Bayesian GATK models. Could swap for '1' for SAMtool models.
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMapQ 30 \
    -minQ 20 \
    -minInd 1 \
    -minMaf 0.05 \
    -doGlf 2 \
    -doPost 1 \
    -doGeno 32 \ # Beagle output format
    -threads {params.threads} \
    2> {log}
    """






  


## UNUSED RULES ##


# Calling genotype likelihoods (GLs) for variant calling using bam files w/ marked duplicates
rule make_gvcfs:
  input:
    bam="results/mkdup/hap{hap}/{sample}_hap{hap}_filtered.bam",
    bai="results/mkdup/hap{hap}/{sample}_hap{hap}_filtered.bai",
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    idx="data/reference/hap{hap}/lupinehap{hap}.dict",
    fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    gvcf="results/gvcf/hap{hap}/{sample}_hap{hap}.g.vcf.gz",
    idx="results/gvcf/hap{hap}/{sample}_hap{hap}.g.vcf.gz.tbi",
  params:
    java_opts="-Xmx4g",
    tempdir="/tmp/bam"
  log:
    "results/logs/make_gvcfs/{sample}_hap{hap}.log"
  envmodules:
    "gatk/4.2.5.0",
    "StdEnv/2020"
  threads: 4
  shell:
    " gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
    " --native-pair-hmm-threads {threads} "
    " -R {input.ref} "
    " -I {input.bam} "
    " -O {output.gvcf} "
    " --tmp-dir {params.tempdir} "
    " -ERC GVCF > {log} 2> {log} "



# create GenomicsDB workspace directory to store data for efficient joint genotyping
# need separate rules for each haplotype
rule hap1_import_genomics_db:
  input:
    ref="data/reference/hap1/lupinehap1.fasta"
  output:
    gdb="results/genomics_db/hap1"
  params:
    java_opts="-Xmx4g",
    tempdir="/tmp/genomics_DB/hap1",
    sample_map="data/genomicsDB/lupine_sample_map_hap1.txt"
  log:
    "results/logs/import_genomics_db/hap1.log"
  envmodules:
    "gatk/4.2.5.0",
    "StdEnv/2020"
  threads: 4
  shell:
    """
    gatk --java-options "{params.java_opts}" GenomicsDBImport 
    --genomicsdb-shared-posixfs-optimizations 
    --genomicsdb-workspace-path {output.gdb} 
    --sample-name-map {params.sample_map}
    --L genome 
    -R {input.ref} 
    --tmp-dir {params.tempdir} 
    &> {log}
    """
    

rule hap2_import_genomics_db:
  input:
    ref="data/reference/hap2/lupinehap2.fasta"
  output:
    gdb="results/genomics_db/hap2"
  params:
    java_opts="-Xmx4g",
    tempdir="/tmp/genomics_DB/hap2",
    sample_map="data/genomicsDB/lupine_sample_map_hap2.txt"
  log:
    "results/logs/import_genomics_db/hap2.log"
  envmodules:
    "gatk/4.2.5.0",
    "StdEnv/2020"
  threads: 4
  shell:
    """
    gatk --java-options "{params.java_opts}" GenomicsDBImport 
    --genomicsdb-shared-posixfs-optimizations 
    --genomicsdb-workspace-path {output.gdb} 
    --sample-name-map {params.sample_map}
    --L genome 
    -R {input.ref} 
    --tmp-dir {params.tempdir} 
    &> {log}
    """