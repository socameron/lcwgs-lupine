## ANALYSIS 5: Thetas (nucleotide diversity, etc)
# PROGRAM: ANGSD | https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests

# We do not filter for MAFs
rule angsd_for_thetas:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.arg",
    mafs_file="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.mafs.gz",
    saf_1="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.idx",
    saf_2="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.gz",
    glf_file="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.glf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites",
    scaffolds="results/scaffolds/hap{hap}_24_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap{hap}/canonical/thetas/by_popln/angsd_canonical_sites_hap{hap}_{population}.log"
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
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -doMajorMinor 1\
    -doMaf 1\
    -doGlf 1\
    -doSaf 2\
    -doCounts 1\
    -anc {params.ref}\
    -doPost 2\
    &> {log}
    """

# NOTE: Did not filter with -minMAF 0.001
rule global_SFS_theta_by_population:
  input:
    saf_idx="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.idx",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    sfs="results/realSFS/hap{hap}/globalSFS/{population}_globalSFS_folded_theta.sfs"
  log:
    "results/logs/realSFS/hap{hap}/globalSFS/{population}_globalSFS_folded_theta.log"
  envmodules:
    "angsd/0.940"
  threads: 40
  shell:
    """
    realSFS {input.saf_idx}\
    -P {threads}\
    -seed 1\
    -fold 1\
    -anc {input.ref}\
    > {output.sfs}\
    2> {log}
    """


rule theta_prep_by_population:
  input:
    sfs="results/realSFS/hap{hap}/globalSFS/{population}_globalSFS_folded_theta.sfs",
    saf_idx="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.idx"
  output:
    theta="results/theta/hap{hap}/{population}_out.thetas.gz",
    index="results/theta/hap{hap}/{population}_out.thetas.idx"
  log:
    "results/logs/theta/hap{hap}/{population}_estimate_theta.log"
  params:
    file="results/theta/hap{hap}/{population}_out"
  envmodules:
    "angsd/0.940"
  threads: 4
  shell:
    """
    realSFS saf2theta {input.saf_idx}\
    -sfs {input.sfs}\
    -outname {params.file}\
    -fold 1\
    -P {threads}\
    &> {log}
    """


# Extract log scaled estimates of 
rule estimate_theta_by_sites:
  input:
    theta_index="results/theta/hap{hap}/{population}_out.thetas.idx"
  output:
    file_name="results/theta/hap{hap}/{population}_log_scale.out"
  envmodules:
    "angsd/0.940"
  threads: 4
  shell:
    """
    thetaStat print {input.theta_index}\
    > {output.file_name}
    """

rule estimate_theta_sliding_window:
  input:
    theta_index="results/theta/hap{hap}/{population}_out.thetas.idx"
  output:
    window="results/theta/hap{hap}/{population}_out.thetasWindow.gz.pestPG"
  params:
    file_name="results/theta/hap{hap}/{population}_out.thetasWindow.gz"
  log:
    "results/logs/theta/hap{hap}/{population}_estimate_theta_sliding_window.log"
  envmodules:
    "angsd/0.940"
  threads: 4
  shell:
    """
    thetaStat do_stat {input.theta_index}\
    -win 10000\
    -step 1000\
    -outnames {params.file_name}\
    &> {log}
    """

rule plot_thetas_output:
  input:
    theta_files=expand("results/theta/hap2/{population}_out.thetasWindow.gz.pestPG", population=POPULATIONS)
  output:
    tajimasD_plot="results/plots/hap2/theta/tajimasD_boxplot.tiff",
    nucleotide_div_plot="results/plots/hap2/theta/nucleotide_diversity_violinplot.tiff"
  log:
    "results/logs/theta/hap2/plot_theta_metrics.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_theta.R
    #Rscript scripts/plot_theta_by_scaffold.R
    """


rule plot_theta_vs_Nc:
  input:
    theta_files = expand("results/theta/hap{hap}/{population}_out.thetasWindow.gz.pestPG",
                         hap=[2], population=POPULATIONS),
    nc_csv      = "data/field_data/Nc_by_site.csv"       # <-- your CSV: Site_code, Site_code2, Range_position, Nc_2023, Nc_2024, Nc_2025
  output:
    scatter = "results/plots/hap{hap}/theta/pi_vs_Nc_{nc_year}.tiff",
    summary = "results/theta/hap{hap}/pi_summary_by_population_{nc_year}.csv"
  params:
    hap = 2,                 # switch to 1 if needed
    nc_year = "2025",        # pick "2023", "2024", or "2025"
    nsites_q = 0.15          # drop bottom 15% nSites within each pop
  log:
    "results/logs/theta/hap{hap}/pi_vs_Nc_{nc_year}.log"
  envmodules:
    "r/4.4.0"
  shell:
    r"""
    Rscript scripts/plot_theta_vs_Nc.R \
      --hap {params.hap} \
      --nc_csv {input.nc_csv} \
      --nsites_q {params.nsites_q} \
      --nc_year {params.nc_year} \
      --scatter {output.scatter} \
      --summary {output.summary} \
      &> {log}
    """