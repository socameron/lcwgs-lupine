## ANALYSIS 3: Runs of Homozygosity (ROH) to look at longer regions of inbreeding across the genome
# PROGRAM 1: RZooRoH | https://cran.r-project.org/web/packages/RZooRoH/vignettes/zooroh-vignette.pdf

# These functions are necessary for bcftools/roh because it calls for individual .vcf files within each popln. 
# Function to expand VCF file paths
def generate_vcf_files():
    vcf_files = []
    for population in POPULATIONS:
        sample_prefixes = get_population_sample_prefixes(population)
        vcf_files.extend(
            [f"results/angsd/hap2/canonical/bcfROH_input/{population}/{prefix}_canonical_SNPs.vcf.gz" for prefix in sample_prefixes]
        )
    return vcf_files

# Function to generate annotated BCF file paths
def generate_annotated_vcfs():
    annotated_vcfs = []
    for population in POPULATIONS:
        sample_prefixes = get_population_sample_prefixes(population)
        annotated_vcfs.extend(
            [f"results/bcftools/hap2/roh_analysis/{population}/{prefix}_annotated.vcf.gz" for prefix in sample_prefixes]
        )
    return annotated_vcfs

# Use *generate_vcf_files() or *generate_annotated_vcfs()

# We call for .beagle.gz and .bcf files
# .beagle is based on genotype likelihoods (GL) and .bcf -> .vcf *can* create genotype probabilities, both of which are acceptable for RZooROH
# For this analysis, we use GLs from .beagle.gz. It is also zipped and takes less space. 
# REVIEW IF .BEAGLE.GZ IS BASED ON GL IF -DoGeno 1 (samples Genotype from posterior) 
rule angsd_for_RZooROH:
  input: 
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.arg",
    mafs_file="results/angsd/hap{hap}/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.hwe.gz",
    depth_sample="results/angsd/hap{hap}/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.depthSample",
    depth_global="results/angsd/hap{hap}/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.depthGlobal",
    beagle="results/angsd/hap{hap}/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.beagle.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/RZooROH_input/{population}/{population}_canonical_SNPs",
    scaffolds="results/scaffolds/hap{hap}_24_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap{hap}/canonical/RZooROH_input/{population}_angsd_canonical_SNPs.log"
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
    -setMaxDepth 4000\
    -minMapQ 30\
    -minQ 20\
    -minInd {params.minInd}\
    -minMaf 0.01\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doPost 2\
    -doCounts 1\
    -doGlf 2\
    -doDepth 1\
    &> {log}
    """


rule RZooROH_get_sample_names:
  input:
    bam_list="data/lists/hap2/{population}_clipped_hap2.txt"
  output:
    sample_names="data/lists/hap2/{population}_sample_names.txt"
  shell:
    r"""
    awk -F'/' '{{ f=$NF; sub(/_hap2_clipped\.bam$/, "", f); print f }}' {input.bam_list} > {output.sample_names}
    """


rule RZooROH_extract_allele_freq:
  input:
    mafs="results/angsd/hap2/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.mafs.gz"
  output:
    allelefreq="results/angsd/hap2/canonical/RZooROH_input/{population}/{population}_allelefreq.txt"
  shell:
    """
    # Skip the header line and print only the sixth column (adjust the column number as required)
    zcat {input.mafs} | awk 'NR>1 {{print $6}}' > {output.allelefreq}
    """


# Convert .beagle.gz files in Oxford Gen format for input into RZooROH
rule RZooROH_prep:
  input:
    beagle_file="results/angsd/hap2/canonical/RZooROH_input/{population}/{population}_canonical_SNPs.beagle.gz",
    script="scripts/convert_beagle_to_rzooroh.py"
  output:
    zoo_txt="results/angsd/hap2/canonical/RZooROH_input/{population}/{population}_Zoo_format.txt"
  envmodules:
    "python/3.11.5"
  shell:
    "python {input.script} {input.beagle_file} {output.zoo_txt}"


# Rscript to run RZooRoH
rule RZooROH_analysis:
  input:
    zoo_txt      = "results/angsd/hap2/canonical/RZooROH_input/{population}/{population}_Zoo_format.txt",
    sample_names = "data/lists/hap2/{population}_sample_names.txt",
    allele_freq  = "results/angsd/hap2/canonical/RZooROH_input/{population}/{population}_allelefreq.txt"
  output:
    placeholder = "results/RZooROH/hap2/RZooROH_analysis/{population}/{population}_placeholder.txt"
  log:
    "results/logs/RZooROH/hap2/{population}_RZooROH_analysis.log"
  envmodules:
    "r/4.4.0"
  params:
    outbase = "results/RZooROH/hap2/RZooROH_analysis/{population}/{population}_summary"
  threads: 4
  shell:
    """
    Rscript scripts/estimate_RZooRoH.R \
      {input.zoo_txt} \
      {input.sample_names} \
      {input.allele_freq} \
      {params.outbase} \
      {wildcards.population} \
      {threads}
    touch {output.placeholder}
    """


# Plot RZooRoH results!
# We plot F (recent) vs F (ancient) with the K value cutoff of 512 for 'recent' inbreeding and anything bigger than K=512 (or smaller fragment) as ancient inbreeding
# We use model 6a from every RZooRoH population run since it was significant!
rule plot_RZooROH_results:
  input:
    realized = expand("results/RZooROH/hap2/RZooROH_analysis/{population}/{population}_realized_MixKR_6a.csv", population = POPULATIONS),
    meta     = "data/lists/hap2/all_popln_geo_coord.csv",
    script   = "scripts/plot_RZooROH.R"
  params:
    recent_thr = 64
  output:
    pdf        = "results/plots/hap2/RZooROH/compare_inbreeding_6a.pdf",
    png_recent = "results/plots/hap2/RZooROH/compare_inbreeding_6a_Frecent.png",
    png_ancient= "results/plots/hap2/RZooROH/compare_inbreeding_6a_Fancient.png",
    png_scatter= "results/plots/hap2/RZooROH/compare_inbreeding_6a_scatter.png"
  threads: 2
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript {input.script} \
      {input.realized} \
      {input.meta} \
      {params.recent_thr} \
      {output.pdf} \
      results/plots/hap2/RZooROH/compare_inbreeding_6a
    """