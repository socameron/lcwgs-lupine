## ANALYSIS 7: Linkage Disequilibrium Decay
# NOTE: this step is necessary for LD pruning prior to the GEA analysis and PCA
# PROGRAM: ngsLD

# We use more slightly stringent rules to call for genotype likelihoods -> minMaf = 0.01 (remove more rare alleles)
# But also liberal for minInd at 50%, and no min/max depth, and no HWE filter for -minHWEpval 0.05
# We use raw genotype likelihoods (GLs) than genotype probabilities (GPs)
# May want to consider exploring minMaf = 0.001 (allow more rare alleles)

rule angsd_for_ngsLD:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_canonical_sites.arg",
    mafs_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_canonical_sites.mafs.gz",
    glf_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_canonical_sites.beagle.gz",
    hwe_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_canonical_sites.hwe.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_canonical_sites",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.5 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap{hap}/canonical/ngsLD/by_popln/angsd_canonical_sites_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list} \
    -ref {params.ref} \
    -out {params.file_name} \
    -remove_bads 1 \
    -rf {params.scaffolds} \
    -GL 1 \
    -C 50 \
    -sites {input.canonical_sites} \
    -minMapQ 30 \
    -minQ 20 \
    -minInd {params.minInd} \
    -minMaf 0.001 \
    -SNP_pval 1e-6 \
    -baq 2 \
    -only_proper_pairs 1 \
    -nThreads {threads} \
    -doMajorMinor 1 \
    -doHWE 1 \
    -doMaf 1 \
    -doGlf 2 \
    -doCounts 1 \
    -doDepth 1 \
    &> {log}
    """

# Linkage disequilibrium analysis
# Make sure to differ max_kb_dist to see differences!
rule ngsLD_analysis:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    mafs_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_canonical_sites.mafs.gz",
    beagle_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_canonical_sites.beagle.gz"
  output:
    GL_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_GL.beagle.gz",
    pos_file="results/angsd/hap{hap}/canonical/ngsLD/by_popln/{population}_positions.pos.gz",
    ld_output="results/ngsLD/hap{hap}/{population}_ld_output.ld"
  params:
    pop_size=lambda wildcards: sum(1 for _ in open(f"data/lists/hap{wildcards.hap}/{wildcards.population}_realign_hap{wildcards.hap}.txt")),
    n_sites=lambda wildcards, input: sum(1 for line in gzip.open(input.mafs_file, "rt")) - 1
  log:
    "results/logs/ngsLD/hap{hap}/{population}_ld_estimation.log"
  threads: 4
  shell:
    """
    # Remove the header and first three columns of the Beagle file
    zcat {input.beagle_file} | awk 'NR>1 {{for (i=4; i<=NF; i++) printf "%s%s", $i, (i==NF?"\\n":"\\t")}}' | gzip > {output.GL_file}

    # Extract chromosome and position from the MAF file (no header)
    zcat {input.mafs_file} | awk 'NR>1 {{print $1"\\t"$2}}' | gzip > {output.pos_file}


    ngsLD --geno {output.GL_file} \
    --probs 1 \
    --n_ind {params.pop_size} \
    --n_sites {params.n_sites} \
    --pos {output.pos_file} \
    --max_kb_dist 100 \
    --out {output.ld_output} \
    --n_threads {threads} \
    --extend_out \
    --verbose 1 \
    &> {log}
    """


# identify linkage for pruning (necessary for by_popln) analyses
# can add --weight_filter to filter edge weights, but need 4th column on weights
# use pruned_list = SNPs remaining after pruning; pruned_SNPs = SNPs that were removed 
rule prune_graph_by_popln:
  input:
    ld_input="results/ngsLD/hap{hap}/{population}_ld_output.ld"
  output:
    tsv_file="results/ngsLD/hap{hap}/by_popln/{population}_ld.tsv",
    pruned_list="results/ngsLD/hap{hap}/by_popln/{population}_pruned_list.txt",
    pruned_snps="results/ngsLD/hap{hap}/by_popln/{population}_pruned_snps.txt"
  params:
    prune_graph="/home/socamero/prune_graph/target/release/prune_graph"
  log:
    "results/logs/ngsLD/hap{hap}/by_popln/{population}_prune_graph.log"
  shell:
    """
    # Extract snp1, snp2, r2 from ngsLD output (assumes tab-delimited with header)
    awk 'NR>1 {{print $1"\\t"$2"\\t"$7}}' {input.ld_input} > {output.tsv_file}

    # Run prune_graph to get list of pruned SNPs
    {params.prune_graph} --in {output.tsv_file} \
    --out {output.pruned_list} \
    --out-excl {output.pruned_snps} \
     &> {log}
    """


# Sample a fraction of the LD because there's about 1.5TB of data in total across all populations!
rule sample_ld_values:
  input:
    ld = "results/ngsLD/hap2/{population}_ld_output.ld"
  output:
    sample = "results/ngsLD/hap2/by_popln/{population}_ld_random_sample.txt"
  params:
    frac= 0.001
  shell:
    """ 
    # Skip header, sample column 7 (r2) with prob=frac
    awk -v f={params.frac} 'NR>1 {{ if (rand() < f) print $7 }}' {input.ld} > {output.sample}
    """

rule sample_ld_pairs:
  input:
    ld="results/ngsLD/hap2/{population}_ld_output.ld"
  output:
    sample="results/ngsLD/hap2/by_popln/{population}_ld_pairs_sample.txt"
  params:
    frac=0.001
  shell:
    r"""
    # NR>1 skip header; $3=dist(bp), $7=r2    awk -F '\t' -v f={params.frac} 'NR>1 {{ if (rand()<f) printf "%d\t%f\n", $3, $7 }}' {input.ld} > {output.sample}
    """


# NOTE: At some point this graph needs to incorporate the Hill & Weir (1988) model for expected LD, which accounts for popln size
rule plot_ld_distribution:
  input:
    random_r2=expand("results/ngsLD/hap2/by_popln/{population}_ld_random_sample.txt", population = POPULATIONS)
  output:
    plot="results/plots/hap2/ngsLD/ld_distribution_by_popln.pdf"
  params:
    script="scripts/plot_LD_distribution.R"
  log:
    "results/logs/ngsLD/hap2/plot_ld_distribution.log"
  threads: 2
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript {params.script} {input.random_r2} {output.plot} 2> {log}
    """

rule plot_ld_decay:
  input:
    random_r2_dist=expand("results/ngsLD/hap2/by_popln/{population}_ld_pairs_sample.txt", population = POPULATIONS)
  output:
    plot="results/plots/hap2/ngsLD/ld_decay_with_distance.pdf"
  params:
    script="scripts/plot_LD_decay.R"
  log:
    "results/logs/ngsLD/hap2/plot_ld_decay.log"
  threads: 2
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript {params.script} {input.random_r2_dist} {output.plot} 2> {log}
    """

# We use the Hill and Weir (1988) model and account for population size
# NOTE - This script will *attempt* to plot all of the LD data, but it might be massive depending on the species
# For Lupine, we do NOT use this script. 
# --vanilla for a clean R environment, --quiet to remove verbosity
rule LD_decay_by_popln:
  input:
    bam_list="data/lists/hap2/{population}_clipped_hap2.txt",
    block_script="/home/socamero/ngsLD/scripts/fit_LDdecay.R",
    ld_input="results/ngsLD/hap2/{population}_ld_output.ld"
  output:
    plot="results/plots/hap2/ngsLD/{population}_plot.pdf"
  params:
    pop_size=lambda wildcards: sum(1 for _ in open(f"data/lists/hap2/{wildcards.population}_clipped_hap2.txt"))
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript --vanilla --quiet {input.block_script} \
    --ld_files {input.ld_input} \
    --out {output.plot} \
    --n_ind {params.pop_size}
    """

# Indicate 1. scaffold name 2. starting position 3. end position
rule LD_blocks_by_popln:
  input:
    block_script="/home/socamero/ngsLD/scripts/LD_blocks.sh",
    ld_input="results/ngsLD/hap2/{population}_ld_output.ld"
  output:
    plot="results/plots/hap2/ngsLD/{population}_500k-2000k.pdf"
  shell:
    """
    cat {input.ld_input} | bash {input.block_script} Scaffold_1__1_contigs__length_29266999 500000 2000000

    mv *.pdf {wildcards.population}_500k-2000k.pdf
    mv {wildcards.population}_500k-2000k.pdf /home/socamero/scratch/{output.plot}
    """