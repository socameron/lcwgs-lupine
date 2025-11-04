## ANALYSIS 1 : Relatedness Structure or basically identify potential clones
# PROGRAM: ngsRelate | https://github.com/ANGSD/NgsRelate

# NOTE: set -minMaf to 0.05
# Previously, did not notice much difference in SFS plots between 0.01 and 0.0001
# ANGSD with -doGlf 3 to prepare input for ngsRelate
# Estimate allele frequencies, genotype likelihoods, and call SNPs
rule angsd_for_ngsRelate:
  input:
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.arg",
    mafs_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.hwe.gz",
    glf_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.glf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM",
    scaffolds="results/scaffolds/hap{hap}_24_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list)))),
    minDepth=lambda wildcards: get_depth_bounds(wildcards.population)[0],
    maxDepth=lambda wildcards: get_depth_bounds(wildcards.population)[1]
  log:
    "results/logs/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/angsd_canonical_SNPs_dupHMM_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
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
    -setMinDepth {params.minDepth}\
    -setMaxDepth {params.maxDepth}\
    -minMapQ 30\
    -minQ 20\
    -minInd {params.minInd}\
    #-minHWEpval 0.01\
    -minMaf 0.05\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    #-doHWE 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doCounts 1\
    -doGlf 3\
    -doPost 2\
    &> {log}
    """
  

rule ngsRelate_prep:
  input:
    mafs_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.mafs.gz",
  output:
    freq_file="results/ngsRelate/hap{hap}/{population}_freq"
  envmodules:
    "gcc",
    "htslib"
  log:
    "results/logs/ngsRelate/hap{hap}/{population}_freq_extraction.log"
  shell:
    """
    mkdir -p $(dirname {output.freq_file})
    zcat {input.mafs_file} | cut -f6 | sed 1d > {output.freq_file} 2>{log}
    """


# Identify potential clones using ngsRelate
rule ngsRelate_analysis:
  input:
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    glf_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.glf.gz",
    freq_file="results/ngsRelate/hap{hap}/{population}_freq"
  output:
    ngsrelate_output="results/ngsRelate/hap{hap}/{population}_ngsrelate.out"
  params:
    pop_size=lambda wildcards: sum(1 for _ in open(f"data/lists/hap{wildcards.hap}/{wildcards.population}_clipped_hap{wildcards.hap}.txt"))
  log:
    "results/logs/ngsRelate/hap{hap}/{population}_ngsrelate_analysis.log"
  threads: 4
  shell:
    """
    ngsRelate -g {input.glf_file} \
              -n {params.pop_size} \
              -f {input.freq_file} \
              -O {output.ngsrelate_output} \
              &> {log}
    """


rule plot_ngsRelate_output:
  input:
    ngsrelate_out=expand("results/ngsRelate/hap2/{population}_ngsrelate.out", population=POPULATIONS)
  output:
    boxplot="results/plots/hap2/ngsRelate/rab_boxplot.png",
    boxplot_nozeroes="results/plots/hap2/ngsRelate/rab_boxplot_no_zeroes.png"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_ngsRelate.R
    """


rule detect_clones:
  input:
    ngsrelate_out="results/ngsRelate/hap{hap}/{population}_ngsrelate.out"
  output:
    clone_summary="results/ngsRelate/hap{hap}/{population}_clone_summary.txt"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/detect_clones.R {input.ngsrelate_out} {output.clone_summary}"


rule ngsRelate_analysis_inbreeding:
  input:
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    glf_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.glf.gz",
    freq_file="results/ngsRelate/hap{hap}/{population}_freq"
  output:
    ngsrelate_output="results/ngsRelate/hap{hap}/{population}_inbreeding.out"
  params:
    pop_size=lambda wildcards: sum(1 for _ in open(f"data/lists/hap{wildcards.hap}/{wildcards.population}_clipped_hap{wildcards.hap}.txt"))
  log:
    "results/logs/ngsRelate/hap{hap}/{population}_ngsrelate_inbreeding_analysis.log"
  threads: 4
  shell:
    """
    ngsRelate -g {input.glf_file} \
              -n {params.pop_size} \
              -F 1\
              -f {input.freq_file} \
              -O {output.ngsrelate_output} \
              &> {log}
    """

rule plot_ngsRelate_inbreeding:
  input:
    ngsrelate_out=expand("results/ngsRelate/hap2/{population}_inbreeding.out", population=POPULATIONS)
  output:
    boxplot="results/plots/hap2/ngsRelate/F_boxplot.png",
    boxplot_nozeroes="results/plots/hap2/ngsRelate/F_boxplot_no_zeroes.png"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_ngsRelate_inbreeding.R
    """