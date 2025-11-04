# ANALYSIS 6: Pairwise Fst Between Populations
# We re-use estimated .saf files from the 05_ANGSD_diversity_stats theta analysis


# For identifying populatio pairs for Fst analysis
POP_COMBINATIONS = list(itertools.combinations(POPULATIONS, 2))

# Varying levels of K for admixture analyses
K_values = [16,17,18,19,20]



# Use the following in rule all:
#expand("results/realSFS/hap2/fst/{pop1}_{pop2}_fst_global.txt", pop1=[x[0] for x in POP_COMBINATIONS], pop2=[x[1] for x in POP_COMBINATIONS])

rule fst_prep_by_population:
  input:
    pop1_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop1}_canonical_sites.saf.idx",
    pop2_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop2}_canonical_sites.saf.idx"
  output:
    sfs_prior="results/realSFS/hap2/fst/{pop1}_{pop2}_prior.ml"
  log:
    "results/logs/realSFS/hap2/globalSFS_Fst/{pop1}_{pop2}_prior.log"
  envmodules:
    "angsd/0.940"
  threads: 40
  shell:
    """
    realSFS {input.pop1_saf_idx} {input.pop2_saf_idx} \
      -P {threads} \
      -fold 1 \
      > {output.sfs_prior} \
      2> {log}
    """

rule fst_analysis:
  input:
    pop1_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop1}_canonical_sites.saf.idx",
    pop2_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop2}_canonical_sites.saf.idx",
    sfs_prior="results/realSFS/hap2/fst/{pop1}_{pop2}_prior.ml"
  output:
    fst_idx="results/realSFS/hap2/fst/{pop1}_{pop2}.fst.idx"
  params:
    fst_out_prefix="results/realSFS/hap2/fst/{pop1}_{pop2}"
  log:
    "results/logs/realSFS/hap2/globalSFS_Fst/{pop1}_{pop2}_Fst_estimation.log"
  envmodules:
    "angsd/0.940"
  threads: 40
  shell:
    """
    realSFS fst index {input.pop1_saf_idx} {input.pop2_saf_idx} \
      -sfs {input.sfs_prior} \
      -fold 1 \
      -fstout {params.fst_out_prefix} \
      -P {threads} \
      2> {log}
    """

# Review sliding window doesn't include regions where site is low
rule estimate_fst_stats:
  input:
    fst_idx="results/realSFS/hap2/fst/{pop1}_{pop2}.fst.idx"
  output:
    global_fst="results/realSFS/hap2/fst/{pop1}_{pop2}_fst_global.txt"
  params:
    window_size=50000,  # Window size for sliding window Fst
    step_size=10000  # Step size for sliding window Fst
  log:
    "results/logs/realSFS/hap2/globalSFS_Fst/{pop1}_{pop2}_Fst_extract.log"
  envmodules:
    "angsd/0.940"
  shell:
    """
    # Global Fst estimate
    realSFS fst stats {input.fst_idx} > {output.global_fst}
    """

    #window_fst="results/realSFS/hap2/fst/{pop1}_{pop2}_fst_windows.txt"
    # Fst in sliding windows
    # realSFS fst stats2 {input.fst_idx} -win {params.window_size} -step {params.step_size} > {output.window_fst}

# Isolation by distance
rule plot_fst_by_distance:
  input:
    global_fst=expand("results/realSFS/hap2/fst/{pop1}_{pop2}_fst_global.txt", pop1=[x[0] for x in POP_COMBINATIONS], pop2=[x[1] for x in POP_COMBINATIONS]),
    coords="data/lists/hap2/all_popln_geo_coord.csv"
  output:
    fst_plot="results/plots/hap2/fst/fst_by_distance.png"
  log:
    "results/logs/Fst/fst_by_distance.log"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/plot_Fst_by_distance.R"