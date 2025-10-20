## ANALYSIS 2: Individual-level inbreeding coefficients
# Program: ngsF and ngsF-HMM | https://github.com/fgvieira/ngsF-HMM

# ANGSD with -doGlf 3 to prepare input for ngsF
# Estimate allele frequencies, genotype likelihoods, call SNPs, and get site frequency spectra
# NOTE: ngsF requires variable sites only for the genotype likelihood input
rule angsd_for_ngsF:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  output:
    arg_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.arg",
    mafs_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.hwe.gz",
    glf_gz="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.gz",
    glf_pos="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.pos.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap{hap}/canonical/ngsF_input/by_popln/angsd_canonical_SNPs_hap{hap}_{population}.log"
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
    -setMinDepth 25\
    -setMaxDepth 3500\
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
    -doCounts 1\
    -anc {params.ref}\
    -doGlf 3\
    &> {log}
    """


# Grab number of variable sites (SNPs) per population
rule count_variable_sites:
  input:
    mafs_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.mafs.gz"
  output:
    site_count="results/ngsF/hap{hap}/by_popln/{population}_SNPs_count.txt"
  shell:
    """
    zcat {input.mafs_file} | wc -l | awk '{{print $1-1}}' > {output.site_count}
    """


# Estimate inbreeding coefficients
# For low-coverage data, set min_epsilon for lower threshold b/w 1e-5 and 1e-9 so algorithm keeps exploring before stopping
rule ngsF_analysis:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    GL3="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.gz",
    SNP_count="results/ngsF/hap{hap}/by_popln/{population}_SNPs_count.txt"
  output:
    ngsF_est="results/ngsF/hap{hap}/by_popln/{population}_ngsF_inbreeding_final.lrt"
  params:
    ngsF_output_base="results/ngsF/hap{hap}/by_popln/{population}_ngsF_inbreeding_iter",
    pop_size=lambda wildcards: sum(1 for _ in open(f"data/lists/hap{wildcards.hap}/{wildcards.population}_realign_hap{wildcards.hap}.txt")),
    n_sites=lambda wildcards: int(open(f"results/ngsF/hap{wildcards.hap}/by_popln/{wildcards.population}_SNPs_count.txt").read().strip()),
    iterations=5  # Number of iterations
  log:
    "results/logs/ngsF/hap{hap}/by_popln/{population}_inbreeding_estimate.log"
  threads: 4
  shell:
    """
    # Run the first iteration without the --init_values parameter
    echo "Running ngsF iteration 1"
    zcat {input.GL3} |\
    ngsF --glf -\
         --n_threads {threads}\
         --calc_LRT 1\
         --out {params.ngsF_output_base}_1\
         --n_ind {params.pop_size}\
         --n_sites {params.n_sites}\
         --init_values r\
         --min_epsilon 1e-7\
         &>> {log}

    # Loop over the remaining iterations
    for iter in $(seq 2 {params.iterations}); do
        echo "Running ngsF iteration $iter"

        # Use the .pars file from the previous iteration as the initial values for the current iteration
        zcat {input.GL3} |\
        ngsF --glf -\
             --n_threads {threads}\
             --calc_LRT 1\
             --out {params.ngsF_output_base}_$iter\
             --n_ind {params.pop_size}\
             --n_sites {params.n_sites}\
             --init_values {params.ngsF_output_base}_$((iter - 1)).pars\
             --min_epsilon 1e-7\
             &>> {log}
    done

    # Copy the final iteration output to the expected result
    cp {params.ngsF_output_base}_{params.iterations}.lrt {output.ngsF_est}
    """


# Use ngsF-HMM which uses a 2-step hidden-markov model
rule ngsF_HMM_analysis:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    GL3="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.gz",
    GL3_pos="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.pos.gz",
    SNP_count="results/ngsF/hap{hap}/by_popln/{population}_SNPs_count.txt"
  output:
    ngsFHMM_est="results/ngsF/hap{hap}/by_popln/{population}_ngsF-HMM_inbreeding.indF",
    GL3_unzipped="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf",
    GL3_pos_unzipped="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.pos"
  params:
    ngsFHMM_output_base="results/ngsF/hap{hap}/by_popln/{population}_ngsF-HMM_inbreeding",
    pop_size=lambda wildcards: sum(1 for _ in open(f"data/lists/hap{wildcards.hap}/{wildcards.population}_realign_hap{wildcards.hap}.txt")),
    n_sites=lambda wildcards: int(open(f"results/ngsF/hap{wildcards.hap}/by_popln/{wildcards.population}_SNPs_count.txt").read().strip())
  log:
    "results/logs/ngsF/hap{hap}/by_popln/{population}_inbreeding_HMM_estimate.log"
  threads: 4
  shell:
    """
    zcat {input.GL3} > {output.GL3_unzipped}
    zcat {input.GL3_pos} > {output.GL3_pos_unzipped}
    ngsF-HMM --geno {output.GL3_unzipped}\
         --n_threads {threads}\
         --pos {output.GL3_pos_unzipped}\
         --out {params.ngsFHMM_output_base}\
         --n_ind {params.pop_size}\
         --n_sites {params.n_sites}\
         --freq r\
         --indF r\
         --loglkl\
         --min_epsilon 1e-7\
         --seed 12345\
         --log 1\
         &>> {log}
    """

rule plot_ngsF_HMM:
  input:
    ngsFHMM_est=expand("results/ngsF/hap2/by_popln/{population}_ngsF-HMM_inbreeding.indF", population=POPULATIONS)
  output:
    plot="results/plots/hap2/ngsF/ngsF-HMM_inbreeding_coeff.tiff"
  log:
    "results/logs/ngsF/hap2/plot_ngsF_HMM_metrics.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_ngsF_HMM.R 
    """
