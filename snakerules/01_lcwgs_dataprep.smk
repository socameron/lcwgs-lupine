##########################
#### DATA PREPARATION ####
##########################

# Trim adapter ends off each sequence file using Trimmomatic
rule trim_reads:
  input:
    r1=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R1"),
    r2=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R2")
  output:
    r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    r1_unp="results/fastq_trimmed/{batch}/{sample}_R1_unpaired.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz",
    r2_unp="results/fastq_trimmed/{batch}/{sample}_R2_unpaired.fastq.gz"
  log:
    "results/logs/trim_reads/{batch}/{sample}.log"
  envmodules:
    "trimmomatic/0.39"
  params:
    adapters="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa"
  shell:
    "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE {input.r1} {input.r2} "
    "{output.r1} {output.r1_unp} {output.r2} {output.r2_unp} "
    "ILLUMINACLIP:{params.adapters}:2:30:10:1 " #previously had ':True' but this failed even though there's documentation of this online
    "LEADING:3 "
    "TRAILING:3 "
    "SLIDINGWINDOW:4:15 " # removes low quality bases
    "MINLEN:36 2> {log}"


# Unzip trimmed files as fastqc 0.12.0 cannot properly read compressed files. 
# Tried with 0.12.1 and it works!
rule unzip_files:
  input:
    zipped_r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    zipped_r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output: 
    unzipped_r1="results/fastq_trimmed_unzip/{batch}/{sample}_R1.fastq",
    unzipped_r2="results/fastq_trimmed_unzip/{batch}/{sample}_R2.fastq"
  shell:
    "gunzip -c {input.zipped_r1} > {output.unzipped_r1} && gunzip -c {input.zipped_r2} > {output.unzipped_r2}"

# Run FastQC per each trimmed sequence file
# Attempting zipped files since using fastqc/0.12.1
rule fastqc:
  input:
    fastq_r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    fastq_r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output:
    html_report_r1="results/fastqc/{batch}/{sample}_R1_fastqc.html",
    zip_report_r1="results/fastqc/{batch}/{sample}_R1_fastqc.zip",
    html_report_r2="results/fastqc/{batch}/{sample}_R2_fastqc.html",
    zip_report_r2="results/fastqc/{batch}/{sample}_R2_fastqc.zip"
  log:
    path="results/logs/fastQC/{batch}/{sample}.log"
  envmodules:
    "fastqc/0.12.1"
  shell:
    "fastqc {input.fastq_r1} {input.fastq_r2} --outdir results/fastqc/{wildcards.batch} 2> {log.path}"

# A little different from above because using the raw fastq paths as well as deriving the filenames to get expected fastqc outputs
# Call in rule all with: expected_fastqc_outputs (this is a predetermined list from above)
rule fastqc_raw:
  input:
    fastq_r1=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R1"),
    fastq_r2=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R2")
  output:
    html_report_r1="results/fastqc_raw/{batch}/{sample}_R1_fastqc.html",
    zip_report_r1="results/fastqc_raw/{batch}/{sample}_R1_fastqc.zip",
    html_report_r2="results/fastqc_raw/{batch}/{sample}_R2_fastqc.html",
    zip_report_r2="results/fastqc_raw/{batch}/{sample}_R2_fastqc.zip"
  params:
    # Compute the default names produced by fastqc from the full input filenames.
    r1_default=lambda wildcards: os.path.basename(get_fastq_paths(wc.batch, wc.sample, "R1")).replace(".fastq.gz", ""),
    r2_default=lambda wildcards: os.path.basename(get_fastq_paths(wc.batch, wc.sample, "R2")).replace(".fastq.gz", "")
  log:
    path="results/logs/fastQC_raw/{batch}/{sample}.log"
  envmodules:
    "fastqc/0.12.1"
  shell:
    """
    # Run FastQC: it will produce files named like:
    #   {params.r1_default}_fastqc.html and {params.r1_default}_fastqc.zip
    #   for read 1, and similarly for read 2.
    fastqc {input.fastq_r1} {input.fastq_r2} --outdir results/fastqc_raw/{wildcards.batch} 2> {log.path}
    
    # Rename the outputs to use your config sample name (the shortened version)
    mv results/fastqc_raw/{wildcards.batch}/{params.r1_default}_fastqc.html {output.html_report_r1}
    mv results/fastqc_raw/{wildcards.batch}/{params.r1_default}_fastqc.zip  {output.zip_report_r1}
    mv results/fastqc_raw/{wildcards.batch}/{params.r2_default}_fastqc.html {output.html_report_r2}
    mv results/fastqc_raw/{wildcards.batch}/{params.r2_default}_fastqc.zip  {output.zip_report_r2}
    """

# Create an aggregated FastQC report using MultiQC.
# Note that we create separate MultiQC reports for batch 1 to 3; and .fastqc files must exist prior to calling it on snakemake
rule multiqc_trimmed:
  input:
    fastqc_dir="results/fastqc/{batch}"
  output:
    html_report="results/multiqc/multiqc_report_trimmed_{batch}.html"
  log:
    path="results/logs/multiqc/{batch}.log"
  params:
    fastqc_dir="results/fastqc/{batch}"
  shell:
    "multiqc -n {output.html_report} {input.fastqc_dir} 2> {log.path}"

rule multiqc_raw:
  input:
    fastqc_dir="results/fastqc_raw/{batch}"
  output:
    html_report="results/multiqc_raw/multiqc_report_raw_{batch}.html"
  log:
    path="results/logs/multiqc_raw/{batch}.log"
  params:
    fastqc_dir="results/fastqc_raw/{batch}"
  shell:
    "multiqc -n {output.html_report} {input.fastqc_dir} 2> {log.path}"

# Creating faidx index for reference genome
rule faidx_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  log:
    "results/logs/refgen/lupinehap{hap}_faidx.log"
  envmodules:
    "samtools/1.20"
  shell:
    "samtools faidx {input} 2> {log}"

# Rules for indexing reference genomes (haplotypes 1 and 2)
rule index_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  log:
    "results/logs/refgen/lupinehap{hap}_bwa_index.log"
  envmodules:
    "bwa/0.7.18"
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
    "samtools/1.20"
  shell:
    "samtools dict {input} > {output} 2> {log}"

# Removing reads smaller than 70bp so 'bwa mem' works better
# Note for a read to be kept, it must be greater than 70bp in BOTH read 1 and 2.
# There are some tools available in envmodules 'seqkit', but I found none that were compatible
rule filter_short_reads_with_seqkit:
  input:
    r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output:
    r1_filtered="results/fastq_filtered/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_filtered="results/fastq_filtered/{batch}/{sample}_R2_filtered.fastq.gz"
  envmodules:
    "StdEnv/2020",
    "seqkit/2.3.1"
  log:
    "results/logs/filter_fastq_70bp/{batch}/{sample}_filter_fastq.log"
  shell:
    """
    seqkit seq -m 70 {input.r1} | gzip > {output.r1_filtered}
    seqkit seq -m 70 {input.r2} | gzip > {output.r2_filtered}
    """

# Pair filtered reads back together
# Note: naming convention stays the same! Just different folder
rule pair_filtered_reads:
  input:
    r1_filtered="results/fastq_filtered/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_filtered="results/fastq_filtered/{batch}/{sample}_R2_filtered.fastq.gz"
  output:
    r1_paired="results/fastq_paired/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_paired="results/fastq_paired/{batch}/{sample}_R2_filtered.fastq.gz"
  envmodules:
    "StdEnv/2020",
    "seqkit/2.3.1"
  log:
    "results/logs/pair_fastq_70bp/{batch}/{sample}_pair_fastq.log"
  shell:
    """
    seqkit pair \
    -1 {input.r1_filtered} \
    -2 {input.r2_filtered} \
    -O results/fastq_paired/{wildcards.batch}
    """

# NOTE: Prior to mapping, some people like to merge fastqs from the same individual/library and remove PCR duplicates prior to mapping using SuperDeduper from HTStream (last used 2020)
# This might be helpful in reducing heterozygote excess
# I have opted not to do this yet as it may be outdated.

# Mapping/Aligning reads to reference haplotypes
rule map_reads:
  input:
    r1="results/fastq_paired/{batch}/{sample}_R1_filtered.fastq.gz",
    r2="results/fastq_paired/{batch}/{sample}_R2_filtered.fastq.gz",
    genome="data/reference/hap{hap}/lupinehap{hap}.fasta",
    idx=multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  output:
    "results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.bam"
  log:
    "results/logs/map_reads/hap{hap}/{batch}/{sample}_hap{hap}.log"
  envmodules:
    "bwa/0.7.18",
    "samtools/1.20"
  threads: 12 
  params:
    # Biological sample (constant across libraries)
    SM=lambda wc: extract_sample_prefix(wc.sample),
    # Per-library identifiers (unique per {sample} string)
    ID=lambda wc: wc.sample,
    LB=lambda wc: extract_sample_prefix(wc.sample) + "_LB",
    PU=lambda wc: f"{wc.batch}.{wc.sample}"
  shell:
    r"""
    RG='@RG\tID:{params.ID}\tSM:{params.SM}\tLB:{params.LB}\tPL:ILLUMINA\tPU:{params.PU}'

    (bwa mem -R "$RG" -t {threads} {input.genome} {input.r1} {input.r2} \
      | samtools view -u \
      | samtools sort -o {output}) 2> {log}
    """

# ID = unique ID | LB = library | PL = platform | PU = platform unit | SM = sample


# To call 'raw' .bam files, we use the get_fastq_paths function in previous rules because the naming conventions are different between batch_1/batch_2 and batch_3. 
# In rule all, you can call:
    #[f"results/bam_raw/hap2/{batch}/{sample}_hap2.bam"
     #for batch in get_batches()
     #for sample in list(get_samples(batch).keys())],

# Sometimes the metadata is randomly not changed (e.g ComputeCanada suddenly replaces all your new files with backups)
rule fix_rg:
  input:
    bam = "results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.bam"
  output:
    bam = "results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.rg.bam",
    bai = "results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.rg.bam.bai"
  log:
    "results/logs/fix_rg/hap{hap}/{batch}/{sample}.log"
  envmodules:
    "samtools/1.20"
  threads: 2
  params:
    SM = lambda wc: extract_sample_prefix(wc.sample),           # e.g. "LCTGP-49"
    ID = lambda wc: wc.sample,                                   # keep lane-unique ID
    LB = lambda wc: extract_sample_prefix(wc.sample) + "_LB",
    PU = lambda wc: f"{wc.batch}.{wc.sample}"
  shell:
    r"""
    samtools addreplacerg \
      -r ID:{params.ID} -r SM:{params.SM} -r LB:{params.LB} -r PL:ILLUMINA -r PU:{params.PU} \
      -o {output.bam} {input.bam} 2> {log}
    samtools index -f {output.bam} 2>> {log}
    """



##############################
#### DATA QUALITY CONTROL ####
##############################

#### DATA QUALITY CONTROL ####

# Steps in QC :
  # 1. Mark PCR and optical duplicates
  # 2. Produce HWE and SFS plots to evaluate mapping problems
  # 3. Identify paralogous regions causing mapping problems
      # - Requires ngsParalog and 1st run of ANGSD
      # - ANGSD: SNP call and input SNPs (aka polymorphic sites) into ngsParalog
      # - ngsParalog: probablistic call of paralogous regions
      # - After ngsParalog, calculate p-values based on chi-sq df = 1 with Benjamini-Hochberg correction
  # 4. Verify heterozygote excess is reduced
      # - Re-run ANGSD excluding paralogous regions and PCR duplicates
      # - Visualize SFS

## STEP 1: MERGE REPLICATES & MARK PCR & OPTICAL DUPLICATES 


# merge replicates if found, otherwise rename and move to hap{hap}/merged/
# NEW: added rm results/bam_raw/hap{hap}/merged/*rep2*.bam because the file remained in the folder. 
rule merge_replicates:
  input:
    lambda wildcards: find_replicates(wildcards.sample_prefix, wildcards.hap)
  output:
    "results/bam_raw/hap{hap}/merged/{sample_prefix}_hap{hap}.bam"
  log:
    "results/logs/merge_replicates/hap{hap}/{sample_prefix}_hap{hap}.log"
  envmodules:
    "samtools/1.17"
  shell:
    """
    module load StdEnv/2020
    module load samtools/1.17
    echo "Input files: {input}" >> {log}
    echo "Output file: {output}" >> {log}
    
    # If no input files are found, exit with an error.
    if [ -z "$(echo {input} | tr -d '[:space:]')" ]; then
      echo "No files found for {wildcards.sample_prefix} in hap{wildcards.hap}" >> {log}
      exit 1
    # If there's only one file, just copy it.
    elif [ $(echo {input} | wc -w) -eq 1 ]; then
      echo "Single file found for {wildcards.sample_prefix} in hap{wildcards.hap}, copying to merged folder." >> {log}
      cp {input} {output}
    else
      echo "Multiple files found for {wildcards.sample_prefix} in hap{wildcards.hap}. Merging..." >> {log}
      samtools merge -f -o {output} {input} 2>> {log}
    fi
    sync

    rm -f results/bam_raw/hap{wildcards.hap}/merged/*rep2*.bam
    """




# Marking and removing PCR duplicates + index
# Note that: /bam_mkdup/ is where marked duplicates are marked AND removed.
# marked duplicates here are now removed (updated March 26, 2025)
rule mark_remove_duplicates:
  input:
    raw_bam="results/bam_raw/hap{hap}/merged/{sample_prefix}_hap{hap}.bam"
  output:
    bam="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bam",
    bai="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bai",
    metrics="results/qc/mkdup_metrics/{sample_prefix}_hap{hap}.metrics"
  log:
    "results/logs/mark_remove_duplicates/hap{hap}/{sample_prefix}_hap{hap}.log"
  envmodules:
    "gatk/4.4.0.0"
  shell:
    """
    gatk MarkDuplicates \
    --CREATE_INDEX \
    -I {input.raw_bam} \
    -O {output.bam} \
    -M {output.metrics} \
    --REMOVE_DUPLICATES true \
    2> {log}
    """


# Clip overlapping reads
# Previous we loaded packages nixpkgs/16.09 intel/2018.3 to run bamutil/1.0.14 but these might be depreciated with new changes to the Supercomputer
# We use samtools/1.18 because it is compatible with the output of bamutil/1.0.14
rule clip_overlapping_reads:
  input:
    bam="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bam"
  output:
    clipped_bam="results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bam",
    clipped_index="results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bai"
  log:
    "results/logs/clip_overlap/hap{hap}/{sample_prefix}_clip_overlapping_reads.log"
  shell:
    """
    module --force purge
    module load StdEnv/2020
    module load bamutil/1.0.14
    bam clipOverlap --in {input.bam} --out {output.clipped_bam} 2> {log} || true
    module --force purge
    module load StdEnv/2023 samtools/1.18
    samtools index -b {output.clipped_bam} -o {output.clipped_index} --threads 4
    module --force purge
    """

# After downstream analyses (admixture), it appears that LCTGP-19 and SWCP-19 are mislabeled and need to swap names.
# The first and second raw reads (fastq files) for both SWCP-19 and LCTGP-19 still match. It was just mislabelled at the sequencing step!
# NOTE: The metadata in the RG group is incorrect still and hasn't been updated!
rule rename_specific_bam_files:
  input:
    lctgp_bam="results/bam_clipped/hap{hap}/LCTGP-19_hap{hap}_clipped.bam",
    lctgp_bai="results/bam_clipped/hap{hap}/LCTGP-19_hap{hap}_clipped.bai",
    swcp_bam="results/bam_clipped/hap{hap}/SWCP-19_hap{hap}_clipped.bam",
    swcp_bai="results/bam_clipped/hap{hap}/SWCP-19_hap{hap}_clipped.bai"
  output:
    checkpoint="results/checkpoints/hap{hap}/rename_specific_files_checkpoint.txt"
  log:
    "results/logs/rename_specific_files/hap{hap}/rename_specific_files.log"
  shell:
    """
    # Temporary files to avoid overwriting
    tmp_lctgp_bam="results/bam_clipped/hap{wildcards.hap}/LCTGP-19_temp.bam"
    tmp_lctgp_bai="results/bam_clipped/hap{wildcards.hap}/LCTGP-19_temp.bai"
    tmp_swcp_bam="results/bam_clipped/hap{wildcards.hap}/SWCP-19_temp.bam"
    tmp_swcp_bai="results/bam_clipped/hap{wildcards.hap}/SWCP-19_temp.bai"

    # Copy the files to temporary locations
    cp {input.lctgp_bam} $tmp_lctgp_bam
    cp {input.lctgp_bai} $tmp_lctgp_bai
    cp {input.swcp_bam} $tmp_swcp_bam
    cp {input.swcp_bai} $tmp_swcp_bai

    # Rename LCTGP-19 to SWCP-19 and vice versa using the temporary files
    mv $tmp_lctgp_bam results/bam_clipped/hap{wildcards.hap}/SWCP-19_hap{wildcards.hap}_clipped.bam
    mv $tmp_lctgp_bai results/bam_clipped/hap{wildcards.hap}/SWCP-19_hap{wildcards.hap}_clipped.bai
    mv $tmp_swcp_bam results/bam_clipped/hap{wildcards.hap}/LCTGP-19_hap{wildcards.hap}_clipped.bam
    mv $tmp_swcp_bai results/bam_clipped/hap{wildcards.hap}/LCTGP-19_hap{wildcards.hap}_clipped.bai

    # Remove the temporary files
    rm -f $tmp_lctgp_bam $tmp_lctgp_bai $tmp_swcp_bam $tmp_swcp_bai

    # Create a checkpoint file to indicate renaming is complete
    echo "LCTGP-19 and SWCP-19 file names successfully swapped and renamed!" > {output.checkpoint}
    """


# Create bam list per population for entry into realign indels 
# Set for haplotype 2
rule generate_clipped_bam_list_per_population:
  input:
    expand("results/bam_clipped/hap2/{sample_prefix}_hap2_clipped.bam", sample_prefix=sample_prefixes),
    checkpoint="results/checkpoints/hap2/rename_specific_files_checkpoint.txt"
  output:
    "data/lists/hap2/{population}_clipped_hap2.list"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")


# Create interval list of indels
# We use apptainer, since I've had issues running gatk/3.8 with simply module load on the DRAC clusters
# MUST PRE-INSTALL gatk/3.8 container in login node
# apptainer pull gatk3.sif docker://broadinstitute/gatk3:3.8-1
rule indel_list:
  input:
    bam_list = "data/lists/hap{hap}/{population}_clipped_hap{hap}.list",
    reference = "data/reference/hap{hap}/lupinehap{hap}.fasta",
  output:
    intervals = "data/lists/hap{hap}/{population}_hap{hap}_indels.intervals",
  log:
    "results/logs/indel_list/hap{hap}/{population}_hap{hap}_indel_list.log",
  threads: 4
  envmodules:
    "apptainer/1.3.5",
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output.intervals})" "$(dirname {log})"

    # Resolve your host scratch root (e.g., /home/socamero/links/scratch)
    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    if [[ ! -d "$HOST_SCRATCH" ]]; then
        echo "Host scratch path not found: $HOST_SCRATCH" >&2
        exit 1
    fi

    # Bind current project dir (for relative paths) AND map host scratch -> /links/scratch inside container
    env -u JAVA_TOOL_OPTIONS \
    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      /home/socamero/gatk3.sif \
      java -Xms2g -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R {input.reference} \
        -I {input.bam_list} \
        -o {output.intervals} \
        -drf BadMate \
        -nt {threads} \
      2> {log}
    """

rule realign_indels:
  input:
    bam = "results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bam",
    ref = "data/reference/hap{hap}/lupinehap{hap}.fasta",
    intervals = lambda wildcards: "data/lists/hap" + wildcards.hap + "/" + next(pop for pop in POPULATIONS if pop in wildcards.sample_prefix) + f"_hap{wildcards.hap}_indels.intervals",
  output:
    realigned_bam = "results/bam_realign/hap{hap}/{sample_prefix}_hap{hap}_realign.bam",
  log:
    "results/logs/realign_indels/hap{hap}/{sample_prefix}_hap{hap}_realign_indels.log",
  envmodules:
    "apptainer/1.3.5",
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {input.intervals})" "$(dirname {log})"

    # Resolve your host scratch root (e.g., /home/socamero/links/scratch)
    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    if [[ ! -d "$HOST_SCRATCH" ]]; then
        echo "Host scratch path not found: $HOST_SCRATCH" >&2
        exit 1
    fi

    # Bind your current project dir (for relative paths) + /links/scratch (for any absolute paths in lists)
    env -u JAVA_TOOL_OPTIONS \
    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      /home/socamero/gatk3.sif \
      java -Xms2g -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R {input.ref} \
        -I {input.bam} \
        -targetIntervals {input.intervals} \
        -o {output.realigned_bam} \
        -drf BadMate \
      &> {log}
    """

# Rule to create .txt file of BAM files for each haplotype  includes all populations
# Only applies to hap2
rule generate_bam_list_realign:
  input:
    expand("results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam", sample_prefix=sample_prefixes)
  output:
    "data/lists/hap2/all_realign_hap2.txt"
  run:
    bam_files = input
    output_file = output[0]

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")



rule generate_bam_list_clipped:
  input:
    expand("results/bam_clipped/hap2/{sample_prefix}_hap2_clipped.bam", sample_prefix=sample_prefixes)
  output:
    "data/lists/hap2/all_clipped_hap2.txt"
  run:
    bam_files = input
    output_file = output[0]

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")



# Merge bam lists into batches 1+2 and 3 because of different sequencing platforms
# Only applies to hap2
rule generate_bam_list_per_round:
  input:
    expand("results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam", sample_prefix=sample_prefixes)
  output:
    first_round="data/lists/hap2/first_round_hap2.txt",
    second_round="data/lists/hap2/second_round_hap2.txt"
  run:
    first_round_bams = []
    second_round_bams = []

    # Define sequencing rounds
    first_round_populations = {"HPW", "IDNP-MW", "LCTGP", "MFNP", "PPP", "RLPLV", "SWCP"}
    
    # Extract population name dynamically
    def extract_population(sample_prefix):
        for pop in POPULATIONS:
            if sample_prefix.startswith(pop):
                return pop  # Return the correct population name
        return None  # If no match is found, which should not happen

    for bam_file in input:
        sample_prefix = bam_file.split("/")[-1].split("_hap")[0]  # Extract sample prefix
        population = extract_population(sample_prefix)  # Get correct population

        if population in first_round_populations:
            first_round_bams.append(bam_file)
        else:
            second_round_bams.append(bam_file)

    # Write BAM lists for each sequencing round
    with open(output.first_round, "w") as f:
        for bam in first_round_bams:
            f.write(f"{bam}\n")

    with open(output.second_round, "w") as f:
        for bam in second_round_bams:
            f.write(f"{bam}\n")


# Merge all bams to estimate per sample depth
rule merge_bams:
  input:
    bam_list="data/lists/hap{hap}/all_realign_hap{hap}.txt"
  output:
    merged_bam="results/bam_merge/hap{hap}/merged_hap{hap}.bam"
  log:
    "results/logs/merge_bams/merge_bam_{hap}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools merge -b {input.bam_list} {output.merged_bam} -@ 12 2> {log}
    """

rule index_bams:
  input:
    merged_bam="results/bam_merge/hap{hap}/merged_hap{hap}.bam"
  output:
    merged_bai="results/bam_merge/hap{hap}/merged_hap{hap}.bam.bai"
  log:
    "results/logs/merged_bams/index_bai_{hap}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools index {input.merged_bam} {output.merged_bai} 2> {log}
    """

ruleorder: 
    estimate_depth_round_1_per_scaffold > estimate_depth_hap2
  
ruleorder:
    estimate_depth_round_2_per_scaffold > estimate_depth_hap2

# Estimate depth per bp per sample
# -d flag : set max-depth to infinity
# -A flag : don't skip anamalous read pairs marked w/ FLAG but missing properly paired flag set
# -q flag : minimum mapping quality
# -Q flag : minimum base quality
# -r flag : specify regions for pileup; needs indexed BAM files
# -a flag : all positions including unused reference sequences
# -f flag : faidx-indexed reference fasta file
# -ff flag : exclude SECONDARY, QCFAIL, DUP
# Note: previously used mpileup with all 189 .bam files, but now merged them as one file because the positions are different across samples with flag -s
# Using ngsQC tools 'bamstats', created by Dr. Tyler Linderoth

# Estimating depth per scaffold!
rule estimate_depth_hap1:
  input:
    bam="results/bam_merge/hap1/merged_hap1.bam",
    bai="results/bam_merge/hap1/merged_hap1.bam.bai",
    ref="data/reference/hap1/lupinehap1.fasta"
  output:
    depth="results/depths/hap1/all_rounds/{hap1scaffold_prefix}_depth_est.txt.gz"
  params:
    hap1scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap1scaffold_prefix, 1)
  envmodules:
    "samtools/1.20"
  shell:
    """
    bamstats {input.bam} -r {params.hap1scaffold}\
    -A -d 77000000 -q 0 -Q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP\
    -s -aa -f {input.ref} | gzip > {output.depth}
    """

rule estimate_depth_hap2:
  input:
    bam="results/bam_merge/hap2/merged_hap2.bam",
    bai="results/bam_merge/hap2/merged_hap2.bam.bai",
    ref="data/reference/hap2/lupinehap2.fasta"
  output:
    depth="results/depths/hap2/all_rounds/{hap2scaffold_prefix}_depth_est.txt.gz"
  params:
    hap2scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap2scaffold_prefix, 2)
  envmodules:
    "samtools/1.20"
  shell:
    """
    bamstats {input.bam} -r {params.hap2scaffold} \
    -A -d 77000000 -q 0 -Q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP \
    -s -aa -f {input.ref} | gzip > {output.depth}
    """


# Merge bams per sequencing round (and platform)
rule merge_bams_round_1:
  input:
    first_round_bams="data/lists/hap{hap}/first_round_hap{hap}.txt"
  output:
    first_round_merged="results/bam_merge/hap{hap}/first_round_hap{hap}.bam",
    first_round_bai="results/bam_merge/hap{hap}/first_round_hap{hap}.bam.bai"
  log:
    first_round_log="results/logs/merge_bams/first_round_hap{hap}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools merge -b {input.first_round_bams} {output.first_round_merged} -@ 12 2> {log.first_round_log}
    samtools index {output.first_round_merged} {output.first_round_bai} 2>> {log.first_round_log}
    """

rule merge_bams_round_2:
  input:
    second_round_bams="data/lists/hap{hap}/second_round_hap{hap}.txt"
  output:
    second_round_merged="results/bam_merge/hap{hap}/second_round_hap{hap}.bam",
    second_round_bai="results/bam_merge/hap{hap}/second_round_hap{hap}.bam.bai"
  log:
    second_round_log="results/logs/merge_bams/second_round_hap{hap}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools merge -b {input.second_round_bams} {output.second_round_merged} -@ 12 2> {log.second_round_log}
    samtools index {output.second_round_merged} {output.second_round_bai} 2>> {log.second_round_log}
    """


# Only applies to hap2 
rule estimate_depth_round_1_per_scaffold:
  input:
    first_round_bam="results/bam_merge/hap2/first_round_hap2.bam",
    first_round_bai="results/bam_merge/hap2/first_round_hap2.bam.bai",
    ref="data/reference/hap2/lupinehap2.fasta"
  output:
    first_round_depth="results/depths/hap2/round_1/{hap2scaffold_prefix}_depth_est.txt.gz"
  params:
    hap2scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap2scaffold_prefix, 2)
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Estimate depth for first sequencing round per scaffold
    bamstats {input.first_round_bam} -r {params.hap2scaffold} \
    -A -d 77000000 -q 0 -Q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP \
    -s -aa -f {input.ref} | gzip > {output.first_round_depth}
    """

# Only applies to hap2 
rule estimate_depth_round_2_per_scaffold:
  input:
    second_round_bam="results/bam_merge/hap2/second_round_hap2.bam",
    second_round_bai="results/bam_merge/hap2/second_round_hap2.bam.bai",
    ref="data/reference/hap2/lupinehap2.fasta"
  output:
    second_round_depth="results/depths/hap2/round_2/{hap2scaffold_prefix}_depth_est.txt.gz"
  params:
    hap2scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap2scaffold_prefix, 2)
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Estimate depth for second sequencing round per scaffold
    bamstats {input.second_round_bam} -r {params.hap2scaffold} \
    -A -d 77000000 -q 0 -Q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP \
    -s -aa -f {input.ref} | gzip > {output.second_round_depth}
    """



# Plot aggregated depths per sequencing rounds/platform
rule plot_depth_per_round:
  input:
    first_round_depth=expand("results/depths/hap2/round_1/{hap2scaffold_prefix}_depth_est.txt.gz", hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES),
    second_round_depth=expand("results/depths/hap2/round_2/{hap2scaffold_prefix}_depth_est.txt.gz", hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    plot="results/plots/hap2/depths/depth_by_round_boxplot.png"
  log:
    "results/logs/plot_depth_per_round.log"
  shell:
    """
    python scripts/plot_depths_by_round.py {input.first_round_depth} {input.second_round_depth} {output.plot} > {log} 2>&1
    """



# Plot aggregate depths per population
rule plot_depth_per_scaffold_hap1:
  input:
    depth_file="results/depths/hap1/all_rounds/{hap1scaffold_prefix}_depth_est.txt.gz"
  output:
    plot="results/plots/hap1/depths/all_rounds/{hap1scaffold_prefix}_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths.py {input.depth_file} {output.plot}
    """

rule plot_depth_per_scaffold_hap2:
  input:
    depth_file="results/depths/hap2/all_rounds/{hap2scaffold_prefix}_depth_est.txt.gz"
  output:
    plot="results/plots/hap2/depths/all_rounds/{hap2scaffold_prefix}_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths.py {input.depth_file} {output.plot}
    """

rule plot_aggregate_depths_hap1:
  input:
    aggregated_depth=expand("results/depths/hap1/all_rounds/{hap1scaffold_prefix}_depth_est.txt.gz", hap1scaffold_prefix=HAP1SCAFFOLD_PREFIXES)
  output:
    plot="results/plots/hap1/depths/hap1_aggregate_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths_WG.py {input.aggregated_depth} {output.plot}
    """

rule plot_aggregate_depths_hap2:
  input:
    aggregated_depth=expand("results/depths/hap2/all_rounds/{hap2scaffold_prefix}_depth_est.txt.gz", hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    plot="results/plots/hap2/depths/hap2_aggregate_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths_WG.py {input.aggregated_depth} {output.plot}
    """



## STEP 2: Quality control check of HWE and SFS

# Rule to create .txt file of BAM files per population
rule generate_bam_hap2_list_per_population:
  input:
    expand("results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam", sample_prefix=sample_prefixes),
  output:
    "data/lists/hap2/{population}_realign_hap2.txt"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")



rule generate_bam_hap1_list_per_population:
  input:
    expand("results/bam_realign/hap1/{sample_prefix}_hap1_realign.bam", sample_prefix=sample_prefixes),
  output:
    "data/lists/hap1/{population}_realign_hap1.txt"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap1_" in bam_file:
                output.write(f"{bam_file}\n")


# ANGSD by population:  To calculate SFS without ngsParalog filtering
# We try different depth levels to get a better distribution of HWE. 

rule angsd_SFS_control_by_population_on_all_sites:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt"
  output:
    arg_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.arg",
    mafs_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.mafs.gz",
    depth_sample="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.depthSample",
    depth_global="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.depthGlobal",
    saf_1="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.idx",
    saf_2="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    min_MAF=lambda wildcards: wildcards.min_MAF,
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
    #min_depth=lambda wildcards: wildcards.min_depth,
    #max_depth=lambda wildcards: wildcards.max_depth
  log:
    "results/logs/angsd/hap{hap}/raw/by_popln/angsd_raw_sites_control_hap{hap}_{population}_minMAF{min_MAF}.log"
  envmodules:
    "angsd/0.940"
  threads: 12
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 1\
    -C 50\
    -minMapQ 30\
    -minQ 20\
    -minInd {params.minInd}\
    -minMAF {params.min_MAF}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doSaf 1\
    -anc {params.ref}\
    &> {log}
    """

# -minMaf 0.01 makes everything with a steep SFS

# Rule for generating global SFS with different depth settings
rule global_SFS_control_by_population:
  input:
    saf_idx="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.idx"
  output:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.sfs"
  log:
    sfs1="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.log"
  threads: 40
  shell:
    """
    winsfs {input.saf_idx} -t {threads} --seed 1 -v > {output.sfs} 2> {log.sfs1}
    """


rule fold_global_SFS_control_by_population:
  input:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.sfs"
  output:
    sfs_folded="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_folded_minMAF{min_MAF}.sfs"
  log:
    sfs2="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_folded_minMAF{min_MAF}.log"
  threads: 40
  shell:
    """
    winsfs view --fold {input.sfs} -v > {output.sfs_folded} 2> {log.sfs2}
    """


# Rule for plotting the folded SFS with different depth settings
rule global_SFS_control_by_population_plots:
  input:
    unfolded_sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.sfs",
    folded_sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_folded_minMAF{min_MAF}.sfs"
  output:
    unfolded_plot="results/plots/hap{hap}/SFS/{population}_globalSFS_control_unfolded_minMAF{min_MAF}.png",
    folded_plot="results/plots/hap{hap}/SFS/{population}_globalSFS_control_folded_minMAF{min_MAF}.png"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    """
    Rscript scripts/SFS_1D_graph.R {input.folded_sfs} {output.folded_plot}
    Rscript scripts/SFS_1D_graph.R {input.unfolded_sfs} {output.unfolded_plot}
    """


# Running ANGSD HWE analysis with different depth settings
rule angsd_HWE_by_population_on_control_SNPs:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.arg",
    mafs_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.mafs.gz",
    hwe_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.hwe.gz",
    depth_sample="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.depthSample",
    depth_global="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.depthGlobal"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    min_depth=lambda wildcards: wildcards.min_depth,
    max_depth=lambda wildcards: wildcards.max_depth
  log:
    "results/logs/angsd/hap{hap}/raw/by_popln/angsd_raw_SNPs_control_hap{hap}_{population}_min{min_depth}_max{max_depth}.log"
  envmodules:
    "angsd/0.940"
  threads: 12
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 1\
    -C 50\
    -minMapQ 30\
    -minQ 20 \
    -minMaf 0.01\
    -minHWEpval 0.01\
    -setMinDepth {params.min_depth}\
    -setMaxDepth {params.max_depth}\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    &> {log}
    """


# Unzipping the HWE output with different depth settings
rule unzip_hwe_control:
  input:
    zipped_hwe="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.hwe.gz"
  output:
    unzip_hwe="results/angsd/hap{hap}/raw/by_popln/{population}_hwe_raw_SNPs_control_min{min_depth}_max{max_depth}.lr"
  shell:
    """
    zcat {input.zipped_hwe} > {output.unzip_hwe}
    """

# Generating HWE histograms with different depth settings
rule hwe_histogram_control:
  input:
    lr_file="results/angsd/hap{hap}/raw/by_popln/{population}_hwe_raw_SNPs_control_min{min_depth}_max{max_depth}.lr"
  output:
    plot="results/plots/hap{hap}/HWE/{population}_hwe_raw_SNPs_control_min{min_depth}_max{max_depth}.png"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/hetero_excess_header.R {input.lr_file} {output.plot}"



## STEP 3: IDENTIFY PARALOGOUS REGIONS 

rule generate_clipped_bam_hap2_list_per_population:
  input:
    expand("results/bam_clipped/hap2/{sample_prefix}_hap2_clipped.bam", sample_prefix=sample_prefixes),
    checkpoint="results/checkpoints/hap2/rename_specific_files_checkpoint.txt"
  output:
    "data/lists/hap2/{population}_clipped_hap2.txt"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")

rule generate_clipped_bam_hap2_list_per_population:
  input:
    expand("results/bam_clipped/hap2/{sample_prefix}_hap2_clipped.bam", sample_prefix=sample_prefixes),
    checkpoint="results/checkpoints/hap2/rename_specific_files_checkpoint.txt"
  output:
    "data/lists/hap2/{population}_clipped_hap2.txt"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")



# Call SNPs with liberal rules for input to ngsParalog
# NOTE: -SNP_pvalue 0.05 so very liberal SNP calls
# NOTE: NO filters applied to gather as much as data as possible as input for ngsParalog, except we put missing data to 20% of individuals
rule angsd_raw_SNP:
  input:
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.txt"
  output:
    arg_file="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_clipped.arg",
    mafs_file="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_clipped.mafs.gz",
    hwe_file="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_clipped.hwe.gz",
    depth_sample="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_clipped.depthSample",
    depth_global="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_clipped.depthGlobal"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_clipped",
    hap="{hap}",
    population="{population}",
    minInd=lambda wildcards, input: max(1, int(0.2 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap{hap}/raw/ngsParalog_input/angsd_SNP_raw_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 12
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -C 50\
    -GL 1\
    -SNP_pval 0.1\
    -minMapQ 30\
    -minQ 20\
    -minInd {params.minInd} \
    -baq 2\
    -only_proper_pairs 1\
    -P {threads}\
    -doHWE 1\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    &> {log}
    """



# Create BED files so including only SNPs into ngsParalog
# NOTE: BED files indexed at 0bp because SAMtools to create pileup requires 0bp index
rule convert_mafs_to_bed:
  input:
    mafs_gz="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_clipped.mafs.gz"
  output:
    bed_file="results/bed/hap{hap}/raw_SNPs/{population}_raw_SNPs_clipped.BED"
  shell:
   """
   gunzip -c {input.mafs_gz} | awk 'NR>1 {{print $1, $2 - 1, $2}}' > {output.bed_file}
   dos2unix {output.bed_file}  # Add this line to convert line endings
   """


rule ngsParalog_hap2:
  input:
    bam_ngsPara=lambda wildcards: expand("results/bam_clipped/hap2/{sample}_hap2_clipped.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    ref="data/reference/hap2/lupinehap2.fasta",
    bed_file="results/bed/hap2/raw_SNPs/{population}_raw_SNPs_clipped.BED"
  output:
    paralog_output="results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr"
  log:
    "results/logs/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.log"
  params:
    hap2scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap2scaffold_prefix, 2),
    minInd=lambda wildcards, input: max(1, int(0.8 * len(input.bam_ngsPara)))
  envmodules:
    "samtools/1.20"
  shell:
    """
    rm -f {output.paralog_output} #remove existing output file if it exists
    touch {output.paralog_output}
    samtools mpileup {input.bam_ngsPara} -A -d 77000000 -q 0 -Q 0 --ff UNMAP,QCFAIL,DUP \
    -l {input.bed_file} -r {params.hap2scaffold} -f {input.ref} 2>> {log} | \
    /home/socamero/ngsParalog/ngsParalog calcLR -infile - -outfile {output.paralog_output} -allow_overwrite 1 \
    -minQ 20 -minind {params.minInd} -mincov 1 \
    -runinfo 1 \
    2>> {log} || true
    """


rule concatenate_paralog_hap2:
  input:
    scaffold_files=lambda wildcards: expand("results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{{population}}-{hap2scaffold_prefix}.lr", population=wildcards.population, hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    paralog_final="results/ngs_paralog/hap2/concat/{population}_paralog_clipped_WGS.lr"
  log:
    "results/logs/ngs_paralog/hap2/concat/{population}_paralog_clipped_WGS.log"
  shell:
    """
    cat {input.scaffold_files} >> {output.paralog_final}
    """


# Print summary of calcLR quantile ranges
rule quantile_summary:
  input:
    "results/ngs_paralog/hap{hap}/concat/{population}_paralog_clipped_WGS.lr"
  output:
    "results/ngs_paralog/hap{hap}/quantile_summary/{population}_calcLR_quantile_summary.txt"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/calcLR_quantile_summary.R {input} {output}"


# Identify false positives from ngsParalog (grab only true Paralogs) and filter out
# NOTE: R script indexes positions back to 1bp start
rule ngsParalog_false_pos:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_clipped_WGS.lr"
  output:
    deviant_snps="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected.BED",
    deviant_snps_bp1="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_clipped_bp1_BH_corrected.lr"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/ngs_paralog_false_pos.R {input.lr_file} {output.deviant_snps} {output.deviant_snps_bp1}
    """


# Identify false positives for varying levels of Benjamini Hochberg critical values
rule ngsParalog_false_pos_BH:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_clipped_WGS.lr"
  output:
    deviant_snps_bp1_BH="results/bed/hap{hap}/deviant_SNPs/BH_correction/{population}_deviant_SNPs_bp1_clipped_BH{BH_VAR}.lr"
  params:
    BH_VAR="{BH_VAR}"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/ngs_paralog_false_pos_BH.R {input.lr_file} {output.deviant_snps_bp1_BH} {params.BH_VAR}
    """


# Create Manhattan Plots [NOTE: Specific Scaffold # otherwise it'll take forever for whole genome!] 
rule ngsParalog_manhattan:
  input:
    lr_file="results/ngs_paralog/hap{hap}/by_popln/{population}_scaffolds/{population}-Scaffold_1.lr"
  output:
    plot="results/plots/hap{hap}/ngsParalog/by_popln/{population}_manhattan_plot_clipped_scaffold_1.png"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/ngs_paralog_graphs.R {input.lr_file} {output.plot}"


## dupHMM - Discover paralogous regions rather than individual sites ##
# Requires calcLR .lr files and depth of coverage .tsv files 

# Convert ngsParalog calcLR outputs (.lr files) to BED format
# This is NOT BH corrected. This is just here for future use
rule convert_calcLR_output_to_BED:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_clipped_WGS.lr",
  output:
    bed="results/bed/hap{hap}/deviant_SNPs_clipped/{population}_deviant_SNPs_calcLR.BED"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/convert_calcLR_to_bed.R {input.lr_file} {output.bed}"


# Generate average depth of coverage per site for dupHMM
# Only estimate depth at sites deemed 'paralogous' from raw ngsParalog (pre BH filter)
# Use pre BH filter because take all data as input!
# Only applies to hap2 (need rule for hap1)
rule estimate_average_coverage_per_population:
  input:
    bam_files=lambda wildcards: expand("results/bam_clipped/hap2/{sample}_hap2_clipped.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    raw_lr="results/ngs_paralog/hap2/concat/{population}_paralog_clipped_WGS.lr"
  output:
    bed="results/bed/hap2/deviant_SNPs/{population}_deviant_SNPs_calcLR.BED",
    avg_cov_file="results/coverage/hap2/{population}_average_deviant_SNP_coverage.tsv"
  threads: 2
  log:
    "results/logs/coverage/hap2/{population}_average_coverage.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Convert .lr file to BED format (0-based start)
    awk '{{print $1 "\\t" ($2-1) "\\t" $2}}' {input.raw_lr} > {output.bed}

    # Calculate the average coverage per position for all BAM files together using the calcLR BED file
    samtools depth -q 0 -Q 0 -J -a -b {output.bed} {input.bam_files} | \
    awk '{{cov[$1"\\t"$2]+=$3; count[$1"\\t"$2]++}} END {{for (pos in cov) print pos, cov[pos] / count[pos]}}' | \
    sort -k1,1V -k2,2n > {output.avg_cov_file}
    """

# Print summary statistics
rule estimate_coverage_statistics:
  input:
    avg_cov_file="results/coverage/hap{hap}/{population}_average_deviant_SNP_coverage.tsv"
  output:
    stats="results/coverage/hap{hap}/{population}_coverage_stats.txt"
  shell:
    """
    python scripts/estimate_coverage_stats.py {input.avg_cov_file} > {output.stats}
    """


# Reduce deviant SNP calls (aka calcLR LRs) to the number of known sites with depth data.
# ANGSD estimates depth globally and per individual, but NOT per site (see -doDepth 1)
rule filter_paralog_by_coverage:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_clipped_WGS.lr",
    cov_file="results/coverage/hap{hap}/{population}_average_deviant_SNP_coverage.tsv"
  output:
    filtered_lr="results/bed/hap{hap}/deviant_SNPs/{population}_coverage_only_paralog_clipped_WGS.lr"
  log:
    "results/logs/filter_lr/hap{hap}/{population}_filtered_deviant_SNPs_clipped.log"
  shell:
    """
    # Filter the .lr file based on the coverage file, ensuring no empty lines
    awk 'NR==FNR && $0!="" {{cov[$1" "$2]=1; next}} $0!="" && cov[$1" "$2] {{print}}' {input.cov_file} {input.lr_file} > {output.filtered_lr}
    """

# dupHMM_setup is run across the entire genome Whereas dupHMM_run is across each scaffold. 
rule dupHMM_setup:
  input:
    lr_file="results/bed/hap{hap}/deviant_SNPs/{population}_coverage_only_paralog_clipped_WGS.lr",
    cov_file="results/coverage/hap{hap}/{population}_average_deviant_SNP_coverage.tsv"
  output:
    param_output="results/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_clipped.par"
  params:
    r_script_path="/home/socamero/ngsParalog/dupHMM.R",
    name_output="results/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_clipped"
  log:
    "results/logs/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_clipped_setup.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript {params.r_script_path} \
    --lrfile {input.lr_file} \
    --outfile {params.name_output} \
    --covfile {input.cov_file} \
    --lrquantile 0.975 \
    --paramOnly 1 \
    2> {log}
    """


# Now generate coverage per scaffold and per population
# We need to input coverage into dupHMM!
rule estimate_coverage_per_population_and_scaffold_hap1:
  input:
    bam_files=lambda wildcards: expand("results/bam_clipped/hap1/{sample}_hap1_clipped.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    raw_lr="results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.lr"
  output:
    bed="results/bed/hap1/deviant_SNPs/{population}_scaffolds/{population}_{hap1scaffold_prefix}.BED",
    avg_cov_file="results/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.tsv"
  log:
    "results/logs/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.log"
  threads: 2
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Check if the .lr file is empty and create a BED file accordingly
    if [ -s {input.raw_lr} ]; then
      # Convert non-empty .lr file to BED format (0-based start)
      awk '{{print $1 "\\t" ($2-1) "\\t" $2}}' {input.raw_lr} > {output.bed}
    else
      # Create an empty BED file if .lr file is empty
      touch {output.bed}
    fi

    # If the BED file is not empty, calculate the average coverage
    if [ -s {output.bed} ]; then
      samtools depth -@ 2 -q 0 -Q 0 -J -a -b {output.bed} {input.bam_files} | \
      awk '{{cov[$1"\\t"$2]+=$3; count[$1"\\t"$2]++}} END {{for (pos in cov) print pos, cov[pos] / count[pos]}}' | \
      sort -k1,1V -k2,2n > {output.avg_cov_file}
    else
      # Create an empty average coverage file if BED file is empty
      touch {output.avg_cov_file}
    fi
    """


rule estimate_coverage_per_population_and_scaffold_hap2:
  input:
    bam_files=lambda wildcards: expand("results/bam_clipped/hap2/{sample}_hap2_clipped.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    raw_lr="results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr"
  output:
    bed="results/bed/hap2/deviant_SNPs/{population}_scaffolds/{population}_{hap2scaffold_prefix}.BED",
    avg_cov_file="results/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.tsv"
  log:
    "results/logs/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.log"
  threads: 2
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Check if the .lr file is empty and create a BED file accordingly
    if [ -s {input.raw_lr} ]; then
      # Convert non-empty .lr file to BED format (0-based start)
      awk '{{print $1 "\\t" ($2-1) "\\t" $2}}' {input.raw_lr} > {output.bed}
    else
      # Create an empty BED file if .lr file is empty
      touch {output.bed}
    fi

    # If the BED file is not empty, calculate the average coverage
    if [ -s {output.bed} ]; then
      samtools depth -@ 2 -q 0 -Q 0 -J -a -b {output.bed} {input.bam_files} | \
      awk '{{cov[$1"\\t"$2]+=$3; count[$1"\\t"$2]++}} END {{for (pos in cov) print pos, cov[pos] / count[pos]}}' | \
      sort -k1,1V -k2,2n > {output.avg_cov_file}
    else
      # Create an empty average coverage file if BED file is empty
      touch {output.avg_cov_file}
    fi
    """

# Filter .lr files by which there is depth data
# This is similar as before, but we need it by scaffold since...
# dupHMM setup is run across the entire genome ! Where as dupHMM run is run across each scaffold instead. 
rule filter_scaffold_lr_by_coverage_hap1:
  input:
    lr_file="results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.lr",
    cov_file="results/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.tsv"
  output:
    filtered_lr="results/bed/hap1/deviant_SNPs/{population}_scaffolds/{population}_{hap1scaffold_prefix}_coverage_only_paralog_clipped.lr"
  log:
    "results/logs/filter_lr/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_filtered_deviant_SNPs_clipped.log"
  shell:
    """
    # Filter the .lr file based on the coverage file, ensuring no empty lines
    awk 'NR==FNR && $0!="" {{cov[$1" "$2]=1; next}} $0!="" && cov[$1" "$2] {{print}}' {input.cov_file} {input.lr_file} > {output.filtered_lr}

    # If the filtered file is empty, create an empty file
    if [ ! -s {output.filtered_lr} ]; then
      touch {output.filtered_lr}
    fi
    """


rule filter_scaffold_lr_by_coverage_hap2:
  input:
    lr_file="results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr",
    cov_file="results/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.tsv"
  output:
    filtered_lr="results/bed/hap2/deviant_SNPs/{population}_scaffolds/{population}_{hap2scaffold_prefix}_coverage_only_paralog_clipped.lr"
  log:
    "results/logs/filter_lr/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_filtered_deviant_SNPs_clipped.log"
  shell:
    """
    # Filter the .lr file based on the coverage file, ensuring no empty lines
    awk 'NR==FNR && $0!="" {{cov[$1" "$2]=1; next}} $0!="" && cov[$1" "$2] {{print}}' {input.cov_file} {input.lr_file} > {output.filtered_lr}

    # If the filtered file is empty, create an empty file
    if [ ! -s {output.filtered_lr} ]; then
      touch {output.filtered_lr}
    fi
    """

# Run dupHMM with beginning estimated parameters
# adjust --lrquantile based on manhattan plot! 
# adjust --maxcoverage based on .tsv files! Seems like some coverages are HIGH! ~50.. To be slightly conservative, choosing 25 (see coverage stats)
# NOTE: Some scaffolds only have 1 SNP with coverage and/or paralog data and can't be used in dupHMM so we skip these. 
rule dupHMM_run_hap1:
  input:
    lr_file = "results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.lr",
    cov_file = "results/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.tsv",
    parameters = "results/ngs_paralog/hap1/dupHMM/{population}_dupHMM_clipped.par"
  output:
    paralog_region = "results/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{population}-{hap1scaffold_prefix}_dupHMM_clipped_run.rf"
  params:
    r_script_path = "/home/socamero/ngsParalog/dupHMM.R",
    name_output = "results/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{population}-{hap1scaffold_prefix}_dupHMM_clipped_run"
  log:
    "results/logs/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{population}-{hap1scaffold_prefix}_dupHMM_clipped_run.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    if [ -s {input.lr_file} ] && [ $(cat {input.lr_file} | wc -l) -gt 1 ]; then
        Rscript {params.r_script_path} \
        --lrfile {input.lr_file} \
        --outfile {params.name_output} \
        --covfile {input.cov_file} \
        --lrquantile 0.97 \
        --maxcoverage 25 \
        --paramfile {input.parameters} \
        2> {log} || touch {output.paralog_region}
    else
        echo "Skipping .lr file due to insufficient data (empty or only one row): {input.lr_file}" >> {log}
        touch {output.paralog_region}
    fi
    """

# NOTE: coverage is lower using hap2
rule dupHMM_run_hap2:
  input:
    lr_file = "results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr",
    cov_file = "results/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.tsv",
    parameters = "results/ngs_paralog/hap2/dupHMM/{population}_dupHMM_clipped.par"
  output:
    paralog_region = "results/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{population}-{hap2scaffold_prefix}_dupHMM_clipped_run.rf"
  params:
    r_script_path = "/home/socamero/ngsParalog/dupHMM.R",
    name_output = "results/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{population}-{hap2scaffold_prefix}_dupHMM_clipped_run"
  log:
    "results/logs/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{population}-{hap2scaffold_prefix}_dupHMM_clipped_run.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    if [ -s {input.lr_file} ] && [ $(cat {input.lr_file} | wc -l) -gt 1 ]; then
        Rscript {params.r_script_path} \
        --lrfile {input.lr_file} \
        --outfile {params.name_output} \
        --covfile {input.cov_file} \
        --lrquantile 0.97 \
        --maxcoverage 20 \
        --paramfile {input.parameters} \
        2> {log} || touch {output.paralog_region}
    else
        echo "Skipping .lr file due to insufficient data (empty or only one row): {input.lr_file}" >> {log}
        touch {output.paralog_region}
    fi
    """


# Combine all dupHMM outputs together into one file per population
rule concatenate_dupHMM_hap1:
  input:
    dupHMM_files=lambda wildcards: expand("results/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{{population}}-{hap1scaffold_prefix}_dupHMM_clipped_run.rf", population=wildcards.population, hap1scaffold_prefix=HAP1SCAFFOLD_PREFIXES)
  output:
    dupHMM_final="results/ngs_paralog/hap1/dupHMM/{population}_dupHMM_clipped_WGS.lr"
  log:
    "results/logs/ngs_paralog/hap1/dupHMM/{population}_dupHMM_clipped_WGS.log"
  shell:
    """
    cat {input.dupHMM_files} >> {output.dupHMM_final}
    """


rule concatenate_dupHMM_hap2:
  input:
    dupHMM_files=lambda wildcards: expand("results/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{{population}}-{hap2scaffold_prefix}_dupHMM_clipped_run.rf", population=wildcards.population, hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    dupHMM_final="results/ngs_paralog/hap2/dupHMM/{population}_dupHMM_clipped_WGS.lr"
  log:
    "results/logs/ngs_paralog/hap2/dupHMM/{population}_dupHMM_clipped_WGS.log"
  shell:
    """
    cat {input.dupHMM_files} >> {output.dupHMM_final}
    """


# Convert dupHMM outputs (.lr files) to BED format
rule convert_dupHMM_output_to_BED:
  input:
    lr_file="results/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_clipped_WGS.lr"
  output:
    bed="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM.BED"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/convert_dupHMM_to_bed.R {input.lr_file} {output.bed}"

    
## STEP 4: VERIFY PARALOGS REDUCED -> look at SFS with RE-RUN ANGSD

# Print out all sites (monomorphic or not) to later filter out paralogous regions
rule angsd_raw_sites_by_popln:
  input:
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.txt"
  output:
    all_sites_gz="results/angsd/hap{hap}/raw/by_popln/{population}_all_sites.pos.gz",
    all_sites_arg="results/angsd/hap{hap}/raw/by_popln/{population}_all_sites.arg",
    all_sites_bed="results/bed/hap{hap}/raw_sites/{population}_all_sites.BED"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    angsd_out="results/angsd/hap{hap}/raw/by_popln/{population}_all_sites"
  log:
    angsd_log="results/logs/angsd/hap{hap}/raw/by_popln/{population}_all_sites.log"
  threads: 12
  envmodules:
    "angsd/0.940"
  shell:
    """
    # Extract all sites across the genome
    angsd -bam {input.bam_list} -ref {params.ref} -out {params.angsd_out} \
    -doCounts 1 -dumpCounts 1 -P {threads} &> {log.angsd_log}

    # Convert the ANGSD output to a BED format file
    gunzip -c {params.angsd_out}.pos.gz | awk 'NR > 1 {{print $1, $2-1, $2}}' > {output.all_sites_bed}
    """


# Format all BED files accordingly and make sure they're in UNIX format + tab delimited
rule process_all_sites_bed:
  input:
    "results/bed/hap{hap}/raw_sites/{population}_all_sites.BED"
  output:
    "results/bed/hap{hap}/raw_sites/by_popln/{population}_all_sites_processed.BED"
  shell:
    """
    dos2unix {input}
    awk 'NF >= 3 {{print $1"\\t"$2"\\t"$3}}' {input} > {output}

    echo "Lines before processing in {input}:"
    wc -l {input}
    echo "Lines after processing in {output}:"
    wc -l {output}
    """

rule process_calcLR_deviant_snps:
  input:
    "results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected.BED"
  output:
    "results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected_processed.BED"
  shell:
    """
    dos2unix {input}
    awk 'NF >= 3 {{print $1"\\t"$2"\\t"$3}}' {input} > {output}

    echo "Lines before processing in {input}:"
    wc -l {input}
    echo "Lines after processing in {output}:"
    wc -l {output}
    """

rule process_dupHMM_deviant_regions:
  input:
    "results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM.BED"
  output:
    "results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM_processed.BED"
  shell:
    """
    dos2unix {input}
    awk 'NF >= 3 {{print $1"\\t"$2"\\t"$3}}' {input} > {output}

    echo "Lines before processing in {input}:"
    wc -l {input}
    echo "Lines after processing in {output}:"
    wc -l {output}
    """


# Filter out deviant SNPs from all known sites
rule filter_all_sites_by_popln_calcLR:
  input:
    processed_all_sites_bed="results/bed/hap{hap}/raw_sites/by_popln/{population}_all_sites_processed.BED",
    processed_calcLR_snps="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected_processed.BED",
  output:
    filtered_calcLR_bed="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.BED",
    filtered_calcLR_txt="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.txt"
  log:
    calcLR_log="results/logs/bedtools/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_sites.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out deviant calcLR SNPs using bedtools and save to BED
    bedtools subtract -a {input.processed_all_sites_bed} -b {input.processed_calcLR_snps} > {output.filtered_calcLR_bed} 2> {log.calcLR_log}

    # Convert the filtered BED files to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3}}' {output.filtered_calcLR_bed} > {output.filtered_calcLR_txt}
    """


rule filter_all_sites_by_popln_dupHMM:
  input:
    processed_dupHMM_regs="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM_processed.BED",
    filtered_calcLR_bed="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.BED"
  output:
    filtered_dupHMM_bed="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.BED",
    filtered_dupHMM_txt="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt"
  log:
    dupHMM_log="results/logs/bedtools/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out dupHMM regions from the calcLR BED file
    bedtools subtract -a {input.filtered_calcLR_bed} -b {input.processed_dupHMM_regs} > {output.filtered_dupHMM_bed} 2> {log.dupHMM_log}

    # Convert the filtered BED files to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3}}' {output.filtered_dupHMM_bed} > {output.filtered_dupHMM_txt}
    """


# Index filtered all sites BED file
rule index_all_sites_by_popln_calcLR:
  input: 
    canonical_calcLR_sites="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.txt"
  output: 
    calcLR_index="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.txt.bin"
  envmodules:
    "angsd/0.940"
  shell: 
    "angsd sites index {input.canonical_calcLR_sites}"


rule index_all_sites_by_popln_dupHMM:
  input:
    canonical_dupHMM_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt"
  output:
    dupHMM_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  envmodules:
    "angsd/0.940"
  shell:
    "angsd sites index {input.canonical_dupHMM_sites}"


# ANGSD by population: To calculate SFS (check if heterozygote excess reduced) with filtered sites by ngsParalog
# Previously attempted just calcLR outputs but this did not remove paralogs. We now incorporate dupHMM and calcLR
rule angsd_SFS_by_population_on_all_sites:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  output:
    arg_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.arg",
    mafs_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.mafs.gz",
    depth_sample="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.depthSample",
    depth_global="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.depthGlobal",
    saf_1="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.idx",
    saf_2="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap{hap}/canonical/by_popln/angsd_canonical_sites_dupHMM_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 12
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 1\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd {params.minInd}\
    -minMaf 0.01\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doSaf 1\
    -anc {params.ref}\
    &> {log}
    """


# FILTERED SFS: Optimize and calculate SFS with folded spectra (--fold) since ancestral state unknown
# Can bootstrap to get confidence intervals
# We use 'winsfs' rather than 'realSFS' because realSFS is computational heavy (>100GB RAM required) plus less accurate
# Based on winsfs github, we allocate ~150GB RAM for ~1B sites (extracted from gzip and wc -l)
# We set seed at 1 for reproducibility
rule global_SFS_by_population:
  input:
    saf_idx="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.idx"
  output:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM.sfs"
  log:
    sfs1="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM.log"
  threads: 40
  shell:
    """
    winsfs {input.saf_idx} -t {threads} --seed 1 -v > {output.sfs} 2> {log.sfs1}
    """


rule fold_global_SFS_by_population:
  input:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM.sfs"
  output:
    sfs_folded="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM_folded.sfs"
  log:
    sfs2="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM_folded.log"
  threads: 40
  shell:
    """
    winsfs view --fold {input.sfs} -v > {output.sfs_folded} 2> {log.sfs2}
    """


# Create SFS Plots on filtered data
rule global_SFS_by_population_plots:
  input:
    sfs_file="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM_folded.sfs"
  output:
    plot="results/plots/hap{hap}/SFS/{population}_globalSFS_dupHMM_folded.png"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/SFS_1D_graph.R {input.sfs_file} {output.plot}"


# To check F distribution on filtered SNPs using dupHMM and calcLR
rule angsd_HWE_by_population_on_dupHMM_SNPs:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.arg",
    mafs_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.hwe.gz",
    depth_sample="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.depthSample",
    depth_global="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.depthGlobal",
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    # Use ~80% of smaples as minimum number of individuals
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap{hap}/canonical/by_popln/angsd_canonical_SNPs_dupHMM_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 12
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 1\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd {params.minInd}\
    -minHWEpval 0.01\
    -minMaf 0.01\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    &> {log}
    """


# Unzip hwe to extract desired variables (i.e F values)
rule unzip_hwe_dupHMM:
  input:
    zipped_hwe="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.hwe.gz"
  output:
    unzip_hwe="results/angsd/hap{hap}/canonical/by_popln/{population}_hwe_canonical_SNPs_dupHMM.lr"
  shell:
    """
    zcat {input.zipped_hwe} > {output.unzip_hwe}
    """


# Create histogram of F values for canonical SNPs (excludes header)
rule hwe_histogram_dupHMM:
  input:
    lr_file="results/angsd/hap{hap}/canonical/by_popln/{population}_hwe_canonical_SNPs_dupHMM.lr"
  output:
    plot="results/plots/hap{hap}/HWE/{population}_hwe_canonical_SNPs_dupHMM.png"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/hetero_excess_header.R {input.lr_file} {output.plot}"


## ESTIMATE DEPTH USING RepAdapt Scripts (credit Gabriele Nocchi)
# See: https://github.com/RepAdapt/snp_calling_simple_bcftools_slurm/blob/main/08_cnv_2.sh

# The gene annotation is based on haplotype 2, so we use haplotype2
rule prepare_depth_files:
  input:
    fai="data/reference/hap2/lupinehap2.fasta.fai",
    fasta="data/reference/hap2/lupinehap2.fasta",
    gff="data/annotation/MCG3698_Lupinus_perennis.annotation.gff"
  output:
    genome_bed="data/reference/hap2/lupinus_perennis_genome.bed",
    windows_bed="data/reference/hap2/lupinus_perennis_windows.bed",
    windows_list="data/reference/hap2/lupinus_perennis_windows.list",
    genes_bed="data/reference/hap2/lupinus_perennis_genes.bed",
    genes_list="data/reference/hap2/lupinus_perennis_genes.list"
  shell:
    """
    # Create a genome file for bedtools (chromosome and length)
    awk '{{print $1"\\t"$2}}' {input.fai} > {output.genome_bed}
    
    # Create a BED file of 5000 bp windows using the FAI file
    awk -v w=5000 '{{chr = $1; chr_len = $2;
      for (start = 0; start < chr_len; start += w) {{
          end = ((start + w) < chr_len ? (start + w) : chr_len);
          print chr "\\t" start "\\t" end;
      }}
    }}' {input.fai} > {output.windows_bed}
    
    # Create a sorted list of window locations
    awk -F "\\t" '{{print $1":"$2"-"$3}}' {output.windows_bed} | sort -k1,1 > {output.windows_list}
    
    # Create a BED file for each gene using the GFF file
    awk '$3 == "gene" {{print $1"\\t"$4"\\t"$5}}' {input.gff} | uniq > {output.genes_bed}
    
    # Sort the gene BED file based on the order of the reference (from the FAI)
    cut -f1 {input.fai} | while read chr; do 
        awk -v chr="$chr" '$1 == chr {{print $0}}' {output.genes_bed} | sort -k2,2n; 
    done > genes.sorted.bed
    mv genes.sorted.bed {output.genes_bed}
    
    # Create a sorted list of gene locations
    awk -F "\\t" '{{print $1":"$2"-"$3}}' {output.genes_bed} | sort -k1,1 > {output.genes_list}
    """

rule estimate_depth_RepAdapt:
  input:
    bam="results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam"
  output:
    temp_depth="results/depths/RepAdapt_temp/{sample_prefix}.depth"
  log:
    "results/logs/create_temp_depth/{sample_prefix}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools depth --reference {input.fasta} -aa {input.bam} > {output.temp_depth}
    """

rule estimate_depth_RepAdapt_stats:
  input:
    temp_depth="results/depths/RepAdapt_temp/{sample_prefix}.depth",
    genome_bed="data/reference/hap2/lupinus_perennis_genome.bed",
    windows_bed="data/reference/hap2/lupinus_perennis_windows.bed",
    windows_list="data/reference/hap2/lupinus_perennis_windows.list",
    genes_bed="data/reference/hap2/lupinus_perennis_genes.bed",
    genes_list="data/reference/hap2/lupinus_perennis_genes.list"
  output:
    wg="results/depths/RepAdapt_method/{sample_prefix}-wg.txt",
    genes_sorted="results/depths/RepAdapt_method/{sample_prefix}-genes.sorted.tsv",
    windows_sorted="results/depths/RepAdapt_method/{sample_prefix}-windows.sorted.tsv"
  params:
    temp_window="results/depths/RepAdapt_temp/{sample_prefix}-windows.tsv",
    temp_genes="results/depths/RepAdapt_temp/{sample_prefix}-genes.tsv"
  log:
    "results/logs/estimate_depth_RepAdapt/{sample_prefix}.log"
  envmodules:
    "samtools/1.20",
    "bedtools/2.31.0"
  shell:
    """
    set +o pipefail #we force the rule because we tested it line by line and it works, just somehow not in one rule

    # Gene depth analysis: compute the mean depth per gene
    awk '{{print $1"\\t"$2"\\t"$2"\\t"$3}}' {input.temp_depth} | bedtools map -a {input.genes_bed} -b stdin -c 4 -o mean -null 0 -g {input.genome_bed} | awk -F "\\t" '{{print $1":"$2"-"$3"\\t"$4}}' | sort -k1,1 > {params.temp_genes} || true

    join -a 1 -e 0 -o '1.1 2.2' -t $'\\t' {input.genes_list} {params.temp_genes} > {output.genes_sorted} || true

    # Window depth analysis: compute the mean depth per window
    awk '{{print $1"\\t"$2"\\t"$2"\\t"$3}}' {input.temp_depth} | bedtools map -a {input.windows_bed} -b stdin -c 4 -o mean -null 0 -g {input.genome_bed} | awk -F "\\t" '{{print $1":"$2"-"$3"\\t"$4}}' | sort -k1,1 > {params.temp_window} || true

    join -a 1 -e 0 -o '1.1 2.2' -t $'\\t' {input.windows_list} {params.temp_window} > {output.windows_sorted} || true

    # Overall genome depth (average depth across all positions)
    awk '{{sum += $3; count++}} END {{if (count > 0) print sum/count; else print "No data"}}' {input.temp_depth} > {output.wg}

    # Cleanup temporary files
    rm -f {input.temp_depth} {params.temp_genes} {params.temp_window}
    """


rule combine_depth_RepAdapt:
  input:
    wg=expand("results/depths/RepAdapt_method/{sample_prefix}-wg.txt", sample_prefix=sample_prefixes),
    genes=expand("results/depths/RepAdapt_method/{sample_prefix}-genes.sorted.tsv", sample_prefix=sample_prefixes),
    windows=expand("results/depths/RepAdapt_method/{sample_prefix}-windows.sorted.tsv", sample_prefix=sample_prefixes)
  output:
    combined_windows="results/depths/RepAdapt_method/combined_windows.tsv",
    combined_genes="results/depths/RepAdapt_method/combined_genes.tsv",
    combined_wg="results/depths/RepAdapt_method/combined_wg.tsv"
  params:
    depth_header="results/depths/RepAdapt_temp/depthheader.txt",
    samples_list="results/depths/RepAdapt_temp/samples.txt",
    genes_temp="results/depths/RepAdapt_temp/combined-genes.temp",
    windows_temp="results/depths/RepAdapt_temp/combined-windows.temp",
    # Create a newline-separated string of sample names from your Python variable.
    sample_names="\n".join(sample_prefixes)
  log:
    "results/logs/estimate_depth_RepAdapt/combine_depths.log"
  shell:
    """
    set -o pipefail

    # Ensure temporary directory exists.
    mkdir -p results/depths/RepAdapt_temp

    # Create a file with the sample names using the Python-provided list.
    echo -e "{params.sample_names}" > {params.samples_list}

    # Create a header using the sample names.
    # Replace newlines with tabs for the header.
    echo -e "location\\t$(cat {params.samples_list} | tr '\n' '\\t')" > {params.depth_header}

    # Combine window depth results:
    while read samp; do 
      cut -f2 results/depths/RepAdapt_method/${{samp}}-windows.sorted.tsv > results/depths/RepAdapt_method/${{samp}}-windows.depthcol; 
    done < {params.samples_list}
    first_sample=$(head -n1 {params.samples_list})
    paste results/depths/RepAdapt_method/${{first_sample}}-windows.sorted.tsv $(tail -n +2 {params.samples_list} | sed 's/.*/results\/depths\/RepAdapt_method\/&-windows.depthcol/') > {params.windows_temp}
    cat {params.depth_header} {params.windows_temp} > {output.combined_windows}

    # Combine gene depth results:
    while read samp; do 
      cut -f2 results/depths/RepAdapt_method/${{samp}}-genes.sorted.tsv > results/depths/RepAdapt_method/${{samp}}-genes.depthcol; 
    done < {params.samples_list}
    first_sample=$(head -n1 {params.samples_list})
    paste results/depths/RepAdapt_method/${{first_sample}}-genes.sorted.tsv $(tail -n +2 {params.samples_list} | sed 's/.*/results\/depths\/RepAdapt_method\/&-genes.depthcol/') > {params.genes_temp}
    cat {params.depth_header} {params.genes_temp} > {output.combined_genes}

    # Combine whole-genome depth results:
    while read samp; do 
      echo -e "${{samp}}\t$(cat results/depths/RepAdapt_method/${{samp}}-wg.txt)"; 
    done < {params.samples_list} > {output.combined_wg}
    """



#

rule generate_bam_list_all_populations:
  input:
    expand("results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam", sample_prefix=sample_prefixes)
  output:
    "data/lists/hap2/all_populations_clipped_hap2.txt"
  run:
    bam_files = input
    output_file = output[0]

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")


rule combine_population_calcLR_bed_files:
  input:
    lambda wildcards: expand("results/bed/hap2/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected.BED", population=POPULATIONS)
  output:
    "results/bed/hap2/deviant_sites/hap2_combined_deviant_SNPs_realign_BH_correction.BED"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Combine all BED files for a specific haplotype and remove duplicates
    bedtools merge -i <(sort -k1,1 -k2,2n $(echo {{' '.join(input)}})) > {output}
    """


rule combine_population_dupHMM_bed_files:
  input:
    lambda wildcards: expand("results/bed/hap2/deviant_SNPs/{population}_deviant_SNPs_dupHMM.BED", population=POPULATIONS)
  output:
    "results/bed/hap2/deviant_sites/hap2_combined_dupHMM_regions.BED"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Combine all BED files for a specific haplotype and remove duplicates
    bedtools merge -i <(sort -k1,1 -k2,2n $(echo {{' '.join(input)}})) > {output}
    """


# Extract all known sites from ALL populations. This is to create list sites and later filter from paralogs
# Parallelize across scaffolds as it would take too long for all samples
rule angsd_raw_sites_all_poplns:
  input:
    bam_list="data/lists/hap2/all_populations_clipped_hap2.txt"
  output:
    all_sites_gz="results/angsd/hap2/raw/all_poplns/{hap2scaffold}_all_sites.pos.gz",
    all_sites_arg="results/angsd/hap2/raw/all_poplns/{hap2scaffold}_all_sites.arg"
  params:
    ref="data/reference/hap2/lupinehap2.fasta",
    angsd_out="results/angsd/hap2/raw/all_poplns/{hap2scaffold}_all_sites"
  log:
    angsd_log="results/logs/angsd/hap2/raw/all_poplns/{hap2scaffold}_all_sites.log"
  threads: 12
  envmodules:
    "angsd/0.940"
  shell:
    """
    # Extract all sites across the genome for the given scaffold from all populations
    angsd -bam {input.bam_list} \
          -ref {params.ref} \
          -out {params.angsd_out} \
          -doCounts 1 -dumpCounts 1 \
          -P {threads} \
          -r {wildcards.hap2scaffold} \
          &> {log.angsd_log}
    """

# Convert each scaffold's sites to BED format
rule convert_raw_sites_scaffold:
  input:
    all_sites_gz="results/angsd/hap2/raw/all_poplns/{hap2scaffold}_all_sites.pos.gz"
  output:
    bed="results/bed/hap2/raw_sites/by_scaffold/{hap2scaffold}_all_sites.BED"
  threads: 4
  shell:
    """
    # Convert the ANGSD output to BED format for the scaffold
    gzip -cd {input.all_sites_gz} | awk 'NR>1 && NF>=3 {{print $1"\t"$2-1"\t"$2}}' > {output.bed}
    dos2unix {output.bed}
    """

# Combine scaffold BED files 
rule combine_raw_sites_scaffolds:
  input:
    beds=expand("results/bed/hap2/raw_sites/by_scaffold/{hap2scaffold}_all_sites.BED", hap2scaffold=HAP2SCAFFOLDS)
  output:
    all_sites_bed="results/bed/hap2/raw_sites/all_poplns/all_sites.BED"
  threads: 4
  shell:
    """    
    # Conmbine BED format file and check for completeness of data
    cat {input.beds} > {output.all_sites_bed}
    dos2unix {output.all_sites_bed}
    """


rule filter_all_sites_all_populations_calcLR:
  input:
    all_sites_bed="results/bed/hap2/raw_sites/all_poplns/all_sites.BED",
    deviant_snps="results/bed/hap2/deviant_sites/hap2_combined_deviant_SNPs_realign_BH_correction.BED"
  output:
    filtered_sites_bed="results/bed/hap2/canonical_sites/filtered_calcLR/calcLR_filtered_sites.BED",
    filtered_sites_txt="results/bed/hap2/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt"
  log:
    calcLR_log="results/logs/bedtools/hap2/canonical_sites/filtered_calcLR/calcLR_filtered_sites.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out deviant sites for all populations using bedtools
    bedtools subtract -a {input.all_sites_bed} -b {input.deviant_snps} > {output.filtered_sites_bed} \
    2> {log.calcLR_log}

    # Convert the filtered BED file to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3}}' {output.filtered_sites_bed} > {output.filtered_sites_txt}
    """


rule filter_all_sites_all_populations_dupHMM:
  input:
    filtered_calcLR_bed="results/bed/hap2/canonical_sites/filtered_calcLR/calcLR_filtered_sites.BED",
    dupHMM_sites="results/bed/hap2/deviant_sites/hap2_combined_dupHMM_regions.BED"
  output:
    filtered_dupHMM_bed="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.BED",
    filtered_dupHMM_txt="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt"
  log:
    dupHMM_log="results/logs/bedtools/hap2/canonical_sites/filtered_dupHMM/dupHMM_filtered_sites.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out deviant sites for all populations using bedtools
    bedtools subtract -a {input.filtered_calcLR_bed} -b {input.dupHMM_sites} > {output.filtered_dupHMM_bed} \
    2> {log.dupHMM_log}

    # Convert the filtered BED file to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3, $3 + 1}}' {output.filtered_dupHMM_bed} > {output.filtered_dupHMM_txt}
    """


rule index_all_sites_all_popln_calcLR:
  input: 
    calcLR_sites="results/bed/hap2/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt"
  output: 
    bin_index="results/bed/hap2/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt.bin",
    idx_index="results/bed/hap2/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt.idx"
  envmodules:
    "angsd/0.940"
  shell: 
    """
    angsd sites index {input.calcLR_sites}
    """


rule index_all_sites_all_popln_dupHMM:
  input: 
    dupHMM_sites="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt"
  output: 
    bin_index="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.bin",
    idx_index="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.idx"
  envmodules:
    "angsd/0.940"
  shell: 
    """
    angsd sites index {input.dupHMM_sites}
    """


# END OF SCRIPT
# For specific analyses, Snakerules in scripts