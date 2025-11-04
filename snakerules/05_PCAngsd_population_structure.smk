## ANALYSIS 4: PCA for population structure + Admixture
# PROGRAM: PCAngsd | https://github.com/Rosemeis/pcangsd
# PROGRAM 2: EMU | https://github.com/Rosemeis/emu
  # EMU takes into account of non-random missing data long the genome. Since I previously masked paralogs, this could be these regions. 

# Population structure uses all data from all populations
# Therefore, this requires consolidating all .BAM files as well the .BED files from ngsParalog for filtering

# This analysis excludes pruning for LD
# Estimate SAF, HWE, GL with SNPs on entire population
rule angsd_for_PCAngsd:
  input:
    bam_list="data/lists/hap2/all_populations_clipped_hap2.txt",
    canonical_sites="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt",
    bin_index="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.bin",
    idx_index="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.idx"
  output:
    arg_file="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.arg",
    mafs_file="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.mafs.gz",
    hwe_file="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.hwe.gz",
    depth_sample="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.depthSample",
    depth_global="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.depthGlobal",
    saf_1="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.saf.idx",
    saf_2="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.saf.pos.gz",
    saf_3="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.saf.gz",
    beagle="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.beagle.gz"
  params:
    ref="data/reference/hap2/lupinehap2.fasta",
    file_name="results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs",
    scaffolds="results/scaffolds/hap2_24_scaffolds.txt"
  log:
    "results/logs/angsd/hap2/canonical/pcangsd_input/angsd_{hap2scaffold}_canonical_SNPs_all_poplns.log"
  envmodules:
    "angsd/0.940"
  threads: 12
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -r {wildcards.hap2scaffold}\
    -GL 1\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 4000\
    -minMapQ 30\
    -minQ 20\
    -minInd 100\
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
    -doSaf 1\
    -anc {params.ref}\
    -doGlf 2\
    -doPost 2\
    &> {log}
    """

# Combine ANGSD beagle outputs files 
rule combine_GLs_all_popln_scaffolds:
  input:
    GL_files=expand("results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.beagle.gz", hap2scaffold=HAP2SCAFFOLDS)
  output:
    all_sites_GL="results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.beagle.gz"
  params:
    script="scripts/combine_and_sort_beagle.py"
  threads: 4
  envmodules:
    "python/3.13.2"
  shell:
    """
    python {params.script} {input.GL_files} {output.all_sites_GL}
    """


# Combine ANGSD mafs output files
rule combine_mafs_all_popln_scaffolds:
  input:
    POS_files=expand("results/angsd/hap2/canonical/pcangsd_input/all_poplns_{hap2scaffold}_canonical_SNPs.mafs.gz", hap2scaffold=HAP2SCAFFOLDS)
  output:
    all_sites_POS="results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.mafs.gz"
  params:
    script="scripts/combine_and_sort_mafs.py"
  threads: 4
  envmodules:
    "python/3.13.2"
  shell:
    """
    python {params.script} {input.POS_files} {output.all_sites_POS}
    """


# Estimate LD across all samples for LD pruning for each scaffold (to parallelize)
rule ngsLD_all_poplns:
  input:
    bam_list="data/lists/hap2/all_populations_clipped_hap2.txt",
    mafs_file="results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.mafs.gz",
    beagle_file="results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.beagle.gz"
  output:
    GL_file="results/angsd/hap2/canonical/ngsLD/all_poplns_GL_{hap2scaffold_prefix}.beagle.gz",
    pos_file="results/angsd/hap2/canonical/ngsLD/all_poplns_positions_{hap2scaffold_prefix}.pos.gz",
    ld_output="results/ngsLD/hap2/all_popln/all_poplns_{hap2scaffold_prefix}_ld_output.ld"
  params:
    pop_size=lambda wildcards: sum(1 for _ in open("data/lists/hap2/all_populations_clipped_hap2.txt"))
  log:
    "results/logs/ngsLD/hap2/all_popln_{hap2scaffold_prefix}_ld_estimation.log"
  threads: 4
  shell:
    """
    # Filter beagle file to keep only markers from this scaffold
    zcat {input.beagle_file} | awk -v s="{wildcards.hap2scaffold_prefix}" 'NR==1 || $1 ~ "^"s"__"' > temp_{wildcards.hap2scaffold_prefix}.beagle.txt

    # Strip first three columns and keep GLs only
    awk 'NR==1 {{next}} {{for (i=4; i<=NF; i++) printf "%s%s", $i, (i==NF?"\\n":"\\t")}}' temp_{wildcards.hap2scaffold_prefix}.beagle.txt | gzip > {output.GL_file}
    rm temp_{wildcards.hap2scaffold_prefix}.beagle.txt

    # Filter MAF file and extract chrom+pos
    zcat {input.mafs_file} | awk -v s="{wildcards.hap2scaffold_prefix}" 'NR==1 || $1 ~ "^"s"__"' > temp_{wildcards.hap2scaffold_prefix}.mafs
    awk 'NR>1 {{print $1"\\t"$2}}' temp_{wildcards.hap2scaffold_prefix}.mafs | gzip > {output.pos_file}
    rm temp_{wildcards.hap2scaffold_prefix}.mafs

    # Count number of sites dynamically
    n_sites=$(zcat {output.pos_file} | wc -l)

    # Run ngsLD
    ngsLD --geno {output.GL_file} \
          --probs 1 \
          --n_ind {params.pop_size} \
          --n_sites $n_sites \
          --pos {output.pos_file} \
          --max_kb_dist 50 \
          --out {output.ld_output} \
          --n_threads {threads} \
          --extend_out \
          --verbose 1 \
          &> {log}
    """


# can add --weight_filter to filter edge weights, but need 4th column on weights
# Filter LD SNPs and use info for PCA below
rule prune_graph_all_poplns:
  input:
    ld_input="results/ngsLD/hap2/all_popln/all_poplns_{hap2scaffold_prefix}_ld_output.ld"
  output:
    tsv_file="results/ngsLD/hap2/all_popln/{hap2scaffold_prefix}_ld.tsv",
    pruned_list="results/ngsLD/hap2/all_popln/{hap2scaffold_prefix}_pruned_list.txt",
    pruned_snps="results/ngsLD/hap2/all_popln/{hap2scaffold_prefix}_pruned_snps.txt"
  params:
    prune_graph="/home/socamero/prune_graph/target/release/prune_graph"
  log:
    "results/logs/ngsLD/hap2/all_popln/{hap2scaffold_prefix}_prune_graph.log"
  threads: 4
  shell:
    """
    # Extract snp1, snp2, r2 from ngsLD output (assumes tab-delimited with header)
    awk 'NR>1 {{print $1"\\t"$2"\\t"$7}}' {input.ld_input} > {output.tsv_file}

    # Run prune_graph to get list of pruned SNPs
    {params.prune_graph} --in {output.tsv_file} \
    --out {output.pruned_list} \
    --out-excl {output.pruned_snps} \
    --weight-filter 'column_3 > 0.2' \
    2> {log}
    """

# !!! Need to check that output of Beagle is compatible with removing markers with pruned_list
# Combine ANGSD beagle outputs files 
rule filter_beagle_with_prune_graph:
  input:
    beagle="results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.beagle.gz",
    pruned_list=expand("results/ngsLD/hap2/all_popln/{hap2scaffold_prefix}_pruned_list.txt", hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    pruned_snps="results/angsd/hap2/canonical/pcangsd_input/all_poplns_pruned_LD_SNPs.beagle.gz"
  log:
    "results/logs/ngsLD/hap2/all_popln/filter_beagle_with_prune_graph.log"
  params:
    script="scripts/filter_beagle_with_prune_graph.py"
  threads: 4
  envmodules:
    "python/3.13.2"
  shell:
    """
    python {params.script} {input.beagle} {input.pruned_list} {output.pruned_snps} 2> {log}
    """


# Calculate PCA on genotype likelihoods using PCAngsd
# .Q output = admixture proportions

rule PCAngsd_all_populations:
  input:
    beagle="results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.beagle.gz"
  output:
    admix_Q="results/pcangsd/hap2/canonical/pcangsd_input/all_popln_canonical_SNP_pcangsd.admix.{K}.Q",
    admix_P="results/pcangsd/hap2/canonical/pcangsd_input/all_popln_canonical_SNP_pcangsd.admix.{K}.P"
  params:
    file_name="results/pcangsd/hap2/canonical/pcangsd_input/all_popln_canonical_SNP_pcangsd"
  log:
    "results/logs/pcangsd/hap2/canonical/pcangsd_input/all_popln_canonical_SNP_pcangsd.{K}.log"
  threads: 8
  envmodules:
    "python/3.12.4"
  shell:
    """
    pcangsd -b {input.beagle}\
    -o {params.file_name}\
    -t {threads}\
    --iter 1000\
    --admix\
    --admix_K {wildcards.K}\
    2> {log}
    """


rule PCAngsd_all_populations_LD_pruned:
  input:
    beagle="results/angsd/hap2/canonical/pcangsd_input/all_poplns_pruned_LD_SNPs.beagle.gz"
  output:
    admix_Q="results/pcangsd/hap2/canonical/pcangsd_input/all_popln_pruned_LD_SNP_pcangsd.admix.{K}.Q",
    admix_P="results/pcangsd/hap2/canonical/pcangsd_input/all_popln_pruned_LD_SNP_pcangsd.admix.{K}.P"
  params:
    file_name="results/pcangsd/hap2/canonical/pcangsd_input/all_popln_pruned_LD_SNP_pcangsd"
  log:
    "results/logs/pcangsd/hap2/canonical/pcangsd_input/all_popln_pruned_LD_SNP_pcangsd.{K}.log"
  threads: 8
  envmodules:
    "python/3.12.4"
  shell:
    """
    pcangsd -b {input.beagle}\
    -o {params.file_name}\
    -t {threads}\
    --iter 1000\
    --admix\
    --admix_K {wildcards.K}\
    2> {log}
    """


# Plot PCAs and admixture plots
# Remove LD pruned datasets for now until that has been resolved
rule PCAngsd_all_populations_plots:
  input:
    pop_info           = "data/lists/hap2/all_samples_to_popln_bam_order.csv",
    admixture_prop     = "results/pcangsd/hap2/canonical/pcangsd_input/all_popln_canonical_SNP_pcangsd.admix.{K}.Q"
  output:
    pca_plot        = "results/plots/hap2/PCAngsd/all_popln_SNP_K{K}_PCA.png",
    admix_plot      = "results/plots/hap2/PCAngsd/all_popln_SNP_K{K}_admix.png",
    #ld_pca_plot     = "results/plots/hap2/PCAngsd/all_popln_LD_pruned_SNP_K{K}_PCA.png",
    #ld_admix_plot   = "results/plots/hap2/PCAngsd/all_popln_LD_pruned_SNP_K{K}_admix.png"
  params:
    cov      = "results/pcangsd/hap2/canonical/pcangsd_input/all_popln_canonical_SNP_pcangsd.cov",
    #ld_cov   = "results/pcangsd/hap2/canonical/pcangsd_input/all_popln_pruned_LD_SNP_pcangsd.cov"
  threads: 2
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_PCA.R  {params.cov}      {input.pop_info} {output.pca_plot}
    Rscript scripts/plot_admixture.R {input.admixture_prop}    {input.pop_info} {output.admix_plot}
    """
    #input ld_admixture_prop  = "results/pcangsd/hap2/canonical/pcangsd_input/all_popln_pruned_LD_SNP_pcangsd.admix.{K}.Q"
    #Rscript scripts/plot_PCA.R  {params.ld_cov}   {input.pop_info} {output.ld_pca_plot}
    #Rscript scripts/plot_admixture.R {input.ld_admixture_prop} {input.pop_info} {output.ld_admix_plot}


# Calculate geographical distances between pairwise individuals between populations
rule PCA_calc_geo_distances:
  input:
    csv="data/lists/hap2/all_samples_geo_coord.csv"
  output:
    dist="results/pcangsd/hap2/canonical/pcangsd_input/pairwise_individ_geodist.csv"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/calculate_popln_geodist.R {input.csv} {output.dist}"


# PROGRAM ngsAdmix - estimate admixture using ngsAdmix
# You can run admixture using PCAngsd OR ngsAdmix.
rule ngsAdmix_analysis:
  input:
    beagle="results/angsd/hap2/canonical/pcangsd_input/all_poplns_canonical_SNPs.beagle.gz"
  output:
    qopt="results/ngsAdmix/hap2/canonical_K{K}.qopt",
    fopt="results/ngsAdmix/hap2/canonical_K{K}.fopt.gz",
  params:
    out_prefix="results/ngsAdmix/hap2/canonical_K{K}"
  log:
    "results/logs/ngsAdmix/hap2/canonical_K{K}.log"
  threads: 8
  envmodules:
    "angsd/0.940"
  shell:
    """
    NGSadmix -likes {input.beagle} \
      -K {wildcards.K} \
      -P {threads} \
      -o {params.out_prefix} &> {log}
    """

# sample_metadata.csv must be in same format as the bamlist from ANGSD
rule plot_ngsAdmix_pie_map:
  input:
    Q="results/ngsAdmix/hap2/canonical_K18.qopt",
    sample_metadata="data/lists/hap2/all_samples_to_popln_bam_order.csv",
    coords="data/lists/hap2/all_popln_geo_coord.csv",
    script="scripts/plot_ngsAdmix_pie_map.R"
  output:
    barplot="results/plots/hap2/ngsAdmix/K18_admixture_barplot.png",
    piemix="results/plots/hap2/ngsAdmix/K18_pie_charts_map.png",
    combined="results/plots/hap2/ngsAdmix/K18_combined_admixture_map.png"
  threads: 4
  envmodules:
    "r/4.4.0",
    "gdal/3.9.1"
  shell:
    """
    Rscript {input.script} \
      {input.Q} \
      {input.sample_metadata} \
      {input.coords} \
      {output.barplot} \
      {output.piemix} \
      {output.combined}
    """

# similar to pie_map except using the R package mapmixture
rule plot_ngsAdmix_mapmixture:
  input:
    Q             = "results/ngsAdmix/hap2/canonical_K18.qopt",
    sample_metadata = "data/lists/hap2/all_samples_to_popln_bam_order.csv",
    coords        = "data/lists/hap2/all_popln_geo_coord.csv",
    script        = "scripts/plot_ngsAdmix_mapmixture.R"
  output:
    barplot       = "results/plots/hap2/ngsAdmix/K18_admixture_barplot_2.png",
    piemix        = "results/plots/hap2/ngsAdmix/K18_pie_charts_map_mapmixture.png",
    combined      = "results/plots/hap2/ngsAdmix/K18_combined_mapmixture.png"
  threads: 4
  envmodules:
    "r/4.4.0",
    "gdal/3.9.1"
  shell:
    """
    Rscript {input.script} \
      {input.Q} \
      {input.sample_metadata} \
      {input.coords} \
      {output.barplot} \
      {output.piemix} \
      {output.combined}
    """