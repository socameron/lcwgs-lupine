## ANALYSIS 2: Individual-level inbreeding coefficients
# Program: ngsF https://github.com/fgvieira/ngsF and ngsF-HMM | https://github.com/fgvieira/ngsF-HMM


# ANGSD with -doGlf 3 to prepare input for ngsF
# Estimate genotype likelihoods (GL 1 with SAMtools algorithm), estimate allele frequencies (-doMaf 2 with unknown minor allele), filter out monomorphic sites, and print GLs to Beagle
# choose which fractions you want to run
rule angsd_for_ngsF:
  input:
    bam_list = "data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    canonical_sites = "results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index = "results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  output:
    mafs    = "results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.mafs.gz",
    arg     = "results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.arg",
    glf_pos = "results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.pos.gz",
    glf     = "results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.gz"
  params:
    ref = "data/reference/hap{hap}/lupinehap{hap}.fasta",
    filename = lambda wc: f"results/angsd/hap{wc.hap}/canonical/ngsF_input/by_popln/{wc.population}_canonical_SNPs_minInd{wc.minInd_prop}",
    scaffolds = "results/scaffolds/hap{hap}_24_scaffolds.txt",
    minInd = lambda wc, input: max(
      1,
      int(float(wc.minInd_prop) * sum(1 for _ in open(input.bam_list)))
    ),
    minDepth = lambda wc: get_depth_bounds(wc.population)[0],
    maxDepth = lambda wc: get_depth_bounds(wc.population)[1]
  log:
    "results/logs/angsd/hap{hap}/canonical/ngsF_input/by_popln/angsd_SNP_discover_hap{hap}_{population}_minInd{minInd_prop}.log"
  threads: 12
  envmodules:
    "angsd/0.940"
  shell:
    """
    angsd -bam {input.bam_list} \
      -ref {params.ref} \
      -anc {params.ref} \
      -out {params.filename} \
      -rf {params.scaffolds} \
      -sites {input.canonical_sites} \
      -GL 1 -C 50 -baq 2 -only_proper_pairs 1 \
      -minMapQ 30 -minQ 20 \
      -setMinDepth {params.minDepth} -setMaxDepth {params.maxDepth} \
      -minInd {params.minInd} \
      -doMajorMinor 1 \
      -doMaf 2 \
      -SNP_pval 1e-6 \
      -doGlf 3 \
      -doCounts 1 \
      -nThreads {threads} &> {log}
    """



# Grab number of variable sites (SNPs) per population
rule count_variable_sites:
  input:
    mafs_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.mafs.gz"
  output:
    site_count="results/ngsF/hap{hap}/by_popln/{population}_SNPs_count_minInd{minInd_prop}.txt"
  shell:
    """
    zcat {input.mafs_file} | wc -l | awk '{{print $1-1}}' > {output.site_count}
    """


# Estimate inbreeding coefficients
# For low-coverage data, set min_epsilon for lower threshold b/w 1e-5 and 1e-9 so algorithm keeps exploring before stopping
rule ngsF_analysis:
  input:
    bam_list   = "data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    GL3        = "results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.gz",
    SNP_count  = "results/ngsF/hap{hap}/by_popln/{population}_SNPs_count_minInd{minInd_prop}.txt"
  output:
    ngsF_est   = "results/ngsF/hap{hap}/by_popln/{population}_minInd{minInd_prop}_ngsF_inbreeding_final.lrt"
  params:
    ngsF_output_base = "results/ngsF/hap{hap}/by_popln/{population}_minInd{minInd_prop}_ngsF_inbreeding_iter",
    pop_size         = lambda wc: sum(1 for _ in open(f"data/lists/hap{wc.hap}/{wc.population}_clipped_hap{wc.hap}.txt")),
    iterations       = 5
  log:
    "results/logs/ngsF/hap{hap}/by_popln/{population}_minInd{minInd_prop}_inbreeding_estimate.log"
  threads: 4
  shell:
    """
    N_SITES=$(cat {input.SNP_count} | tr -d '\r')
    if [ "${{N_SITES}}" -le 0 ]; then
      echo "No SNPs to analyze (N_SITES=${{N_SITES}}). Exiting." >> {log}
      : > {output.ngsF_est}
      exit 0
    fi

    echo "Running ngsF for {wildcards.population} hap{wildcards.hap} minInd={wildcards.minInd_prop} with N_SITES=${{N_SITES}} and N_IND={params.pop_size}" >> {log}

    # Iteration 1: random init
    zcat {input.GL3} | \
    ngsF --glf - \
         --n_threads {threads} \
         --calc_LRT 1 \
         --out {params.ngsF_output_base}_1 \
         --n_ind {params.pop_size} \
         --n_sites ${{N_SITES}} \
         --init_values r \
         --min_epsilon 1e-6 \
         &>> {log}

    # Subsequent iterations: warm start from previous .pars
    for iter in $(seq 2 {params.iterations}); do
      prev=$((iter - 1))
      echo "Running ngsF iteration $iter" >> {log}
      zcat {input.GL3} | \
      ngsF --glf - \
           --n_threads {threads} \
           --calc_LRT 1 \
           --out {params.ngsF_output_base}_$iter \
           --n_ind {params.pop_size} \
           --n_sites ${{N_SITES}} \
           --init_values {params.ngsF_output_base}_$prev.pars \
           --min_epsilon 1e-7 \
           &>> {log}
    done

    cp {params.ngsF_output_base}_{params.iterations}.lrt {output.ngsF_est}
    """




# ANGSD calls for ngsF-HMM, which requires GL but not only on variable sites -> remove SNP_pval
rule angsd_for_ngsF_HMM:
  input:
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  output:
    glf_gz="results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.gz",
    glf_pos="results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.pos.gz",
    counts="results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.mafs.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    filename=lambda wildcards: f"results/angsd/hap{wildcards.hap}/canonical/ngsF_HMM_input/by_popln/{wildcards.population}_canonical_SNPs_minInd{wildcards.minInd_prop}",
    scaffolds="results/scaffolds/hap{hap}_24_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(float(wildcards.minInd_prop) * sum(1 for _ in open(input.bam_list)))),
    minDepth=lambda wildcards: get_depth_bounds(wildcards.population)[0],
    maxDepth=lambda wildcards: get_depth_bounds(wildcards.population)[1]
  log:
    "results/logs/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/angsd_canonical_SNPs_hap{hap}_{population}_minInd{minInd_prop}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list} \
      -ref {params.ref} \
      -anc {params.ref} \
      -out {params.filename} \
      -rf {params.scaffolds} \
      -sites {input.canonical_sites} \
      -GL 1 -C 50 -baq 2 -only_proper_pairs 1 \
      -minMapQ 30 -minQ 20 \
      -setMinDepth {params.minDepth} -setMaxDepth {params.maxDepth} \
      -minInd {params.minInd} \
      -doMajorMinor 1 \
      -doMaf 2 \
      -doGlf 3 \
      -doCounts 1 \
      -nThreads {threads} &> {log}
    """


# Grab number of sites per population
rule count_ngsF_HMM_sites:
  input:
    glf_pos="results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.pos.gz"
  output:
    site_count="results/ngsF_HMM/hap{hap}/by_popln/{population}_minInd{minInd_prop}_sites_count.txt"
  shell:
    """
    zcat {input.glf_pos} | wc -l > {output.site_count}
    """


# Use ngsF-HMM which uses a 2-step hidden-markov model\
rule ngsF_HMM_analysis:
  input:
    bam_list   = "data/lists/hap{hap}/{population}_clipped_hap{hap}.txt",
    GL3        = "results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.gz",
    GL3_pos    = "results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.pos.gz",
    site_count = "results/ngsF_HMM/hap{hap}/by_popln/{population}_minInd{minInd_prop}_sites_count.txt"
  output:
    ngsFHMM_est      = "results/ngsF_HMM/hap{hap}/by_popln/{population}_minInd{minInd_prop}_ngsF-HMM_inbreeding.indF",
    GL3_unzipped     = "results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf",
    GL3_pos_unzipped = "results/angsd/hap{hap}/canonical/ngsF_HMM_input/by_popln/{population}_canonical_SNPs_minInd{minInd_prop}.glf.pos"
  params:
    ngsFHMM_output_base = "results/ngsF_HMM/hap{hap}/by_popln/{population}_minInd{minInd_prop}_ngsF-HMM_inbreeding",
    pop_size = lambda wc: sum(1 for _ in open(f"data/lists/hap{wc.hap}/{wc.population}_clipped_hap{wc.hap}.txt"))
  log:
    "results/logs/ngsF_HMM/hap{hap}/by_popln/{population}_minInd{minInd_prop}_inbreeding_HMM_estimate.log"
  threads: 8
  shell:
    """
    if [ "$(cat {input.site_count})" -le 0 ]; then
      echo "No sites for {wildcards.population} hap{wildcards.hap} minInd={wildcards.minInd_prop}" >> {log}
      : > {output.ngsFHMM_est}
      : > {output.GL3_unzipped}
      : > {output.GL3_pos_unzipped}
      exit 0
    fi

    zcat {input.GL3}     > {output.GL3_unzipped}
    zcat {input.GL3_pos} > {output.GL3_pos_unzipped}

    ngsF-HMM \
      --geno {output.GL3_unzipped} \
      --n_threads {threads} \
      --pos {output.GL3_pos_unzipped} \
      --out {params.ngsFHMM_output_base} \
      --n_ind {params.pop_size} \
      --n_sites "$(cat {input.site_count})" \
      --freq r \
      --indF r \
      --loglkl \
      --min_epsilon 1e-6 \
      --seed 12345 \
      --log 1 \
      &>> {log}
    """

rule plot_ngsF_HMM:
  input:
    ngsFHMM_est=expand("results/ngsF_HMM/hap2/by_popln/{population}_minInd{minInd_prop}_ngsF-HMM_inbreeding.indF", population=POPULATIONS, minInd_prop = MININD_PROPS)
  output:
    plot="results/plots/hap2/ngsF_HMM/ngsF-HMM_inbreeding_coeff.tiff"
  log:
    "results/logs/ngsF_HMM/hap2/plot_ngsF_HMM_metrics.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_ngsF_HMM.R 
    """
