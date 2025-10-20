## ANALYSIS 9: Contemporary Effective Population Size

# This analysis uses CurrentNe2, which builds upon the theory presented by commonly used NeEstimator
# CurrentNe accounts for full-siblings contributing to LD, while CurrentNe accounts for 
# Confidence intervals is also better estimated in CurrentNe than NeEstimator

# Some criteria for NeEstimator:
  # Since this LD method requires genotype calls, we filter out individuals with coverage lower than 8x
  # We then filter out rare alleles and set a stringent call of removing alleles with <2% freq (-minMaf)
  # We also only allow SNPs where 90% of all individuals are represented (-minInd)
  # Only call genotypes if it reaches beyond 95% cutoff (-postCutoff)

rule angsd_for_currentNe2:
  input:
    bam_list="data/lists/hap2/{population}_clipped_hap2.txt",
    canonical_sites="results/bed/hap2/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap2/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap2/lupinehap2.fasta.fai"
  output:
    arg_file="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.arg",
    mafs_file="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.mafs.gz",
    depth_sample="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.depthSample",
    depth_global="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.depthGlobal",
    bcf_file="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.bcf"
  params:
    ref="data/reference/hap2/lupinehap2.fasta",
    file_name="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs",
    scaffolds="results/scaffolds/hap2_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.9 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap2/canonical/currentNe2_input/{population}_angsd_canonical_SNPs.log"
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
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 50\
    -setMaxDepth 4000\
    -minMapQ 30\
    -minQ 20\
    -minInd {params.minInd}\
    -minMaf 0.01\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -GL 1\
    -doPost 1\
    -doMajorMinor 1\
    -doMaf 1\
    -dobcf 1\
    -doGeno 2\
    -doCounts 1\
    -doDepth 1\
    -postCutoff 0.90\
    &> {log}
    """

rule convert_bcf_to_vcf_currentNe2:
  input:
    bcf="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.bcf"
  output:
    vcf="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.vcf.gz",
    csi="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.vcf.gz.csi"
  threads: 2
  envmodules:
    "bcftools/1.19"
  shell:
    """
    # just convert the BCF (with GT) to compressed VCF - not necessary to use bcftools since not inferring GT from GL
    bcftools view -Oz -o {output.vcf} {input.bcf}
    bcftools index -f {output.vcf}
    """

rule convert_vcf_GP_to_GT_currentNe2:
  input:
    vcf="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs.vcf.gz"
  output:
    vcf_GT="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs_GT.vcf"
  threads: 2
  envmodules:
    "bcftools/1.19"
  shell:
    """
    # 1) Convert GP -> GT (calls GT if max GP >= 0.95, else sets ./.)
    bcftools +tag2tag {input.vcf} -- --GP-to-GT -t 0.95 \
      | bcftools annotate -x FORMAT/AD,FORMAT/DP,FORMAT/PL,FORMAT/GL,FORMAT/GQ,FORMAT/PGT,FORMAT/PS \
      -Ov -o {output.vcf_GT}

    # 2) Index the *VCF* IF wanted a zipped version (with -Oz above), which produces the .csi next to it
    #bcftools index -f {output.vcf_GT}
    """

# Use RLPLV .vcf to determine genome size (# of chromosomes)
rule genome_size_morgans_global:
  input:
    vcf = "results/angsd/hap2/canonical/currentNe2_input/RLPLV/RLPLV_canonical_SNPs_GT.vcf.gz"
  output:
    genome_size = "results/currentNe2/hap2/genome_size_morgans.txt"
  log:
    "results/logs/currentNe2/hap2/genome_size_morgans.log"
  threads: 4
  envmodules:
    "bcftools/1.19"
  shell:
    r"""
    mkdir -p results/currentNe2/hap2 results/logs/currentNe2/hap2
    if [ -n "{GENOME_MORGANS}" ] && [ "{GENOME_MORGANS}" != "None" ]; then
        echo {GENOME_MORGANS} > {output.genome_size}
        echo "Override genome_morgans={GENOME_MORGANS}" > {log}
    else
        echo "Counting CHROM from RLPLV: {input.vcf}" > {log}
        bcftools query -f '%CHROM\n' {input.vcf} \
          | sed 's/^[[:space:]]*//' \
          | awk 'NF' \
          | sort -u \
          | wc -l \
          | awk '{{print $1}}' > {output.genome_size}
        echo "Derived GS=$(cat {output.genome_size})" >> {log}
    fi
    """


# Run currentNe2 but estimate K parameter (breeding system)
rule currentne2_estK:
  input:
    vcf="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs_GT.vcf",
    gs="results/currentNe2/hap2/genome_size_morgans.txt"
  output:
    ne="results/currentNe2/hap2/{population}/{population}_estK_currentNe2_OUTPUT.txt",
    mix="results/currentNe2/hap2/{population}/{population}_estK_currentNe2_mix_OUTPUT.txt"
  params:
    filename="results/currentNe2/hap2/{population}/{population}_estK"
  log:
   "results/logs/currentNe2/hap2/{population}/{population}_estK.log"
  threads: 4
  shell:
    """
    set -euo pipefail
    mkdir -p results/currentNe2/hap2/{wildcards.population} results/logs/currentNe2/hap2/{wildcards.population}
    GS=$(cat {input.gs})
    echo "[currentNe2 estK] GS=$GS; VCF={input.vcf}" > {log}
    currentne2 -t {threads} -o {params.filename} {input.vcf} $GS >> {log} 2>&1
    """


# One job per k value (k may be 0 or a decimal like 0.2, 0.33, 1.0)
rule currentne2_k:
  input:
    vcf = "results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs_GT.vcf",
    gs  = "results/currentNe2/hap2/genome_size_morgans.txt"
  output:
    ne = "results/currentNe2/hap2/{population}/{population}_k{k}_currentNe2_OUTPUT.txt",
    mix = "results/currentNe2/hap2/{population}/{population}_k{k}_currentNe2_mix_OUTPUT.txt"
  params:
    filename="results/currentNe2/hap2/{population}/{population}_k{k}"
  log:
    "results/logs/currentNe2/hap2/{population}/{population}_k{k}.log"
  threads: 4
  shell:
    """
    set -euo pipefail
    mkdir -p results/currentNe2/hap2/{wildcards.population} results/logs/currentNe2/hap2/{wildcards.population}
    GS=$(cat {input.gs})
    echo "[currentNe2 k] pop={wildcards.population} k={wildcards.k} GS=$GS; VCF={input.vcf}" > {log}
    currentne2 -t {threads} -k {wildcards.k} -o {params.filename} {input.vcf} $GS >> {log} 2>&1
    """


# Run currentNe2 and invoke -x allow for population structure within subpoplns
rule currentne2_x:
  input:
    vcf="results/angsd/hap2/canonical/currentNe2_input/{population}/{population}_canonical_SNPs_GT.vcf",
    gs="results/currentNe2/hap2/genome_size_morgans.txt"
  output:
    ne="results/currentNe2/hap2/{population}/{population}_x_currentNe2_OUTPUT.txt"
  params:
    filename="results/currentNe2/hap2/{population}/{population}_x"
  log:
    "results/logs/currentNe2/hap2/{population}/{population}_x.log"
  threads: 4
  shell:
    """
    set -euo pipefail
    mkdir -p results/currentNe2/hap2/{wildcards.population} results/logs/currentNe2/hap2/{wildcards.population}
    GS=$(cat {input.gs})
    echo "[currentNe2 -x] GS=$GS; VCF={input.vcf}" > {log}
    currentne2 -t {threads} -x -o {params.filename} {input.vcf} $GS >> {log} 2>&1
    """
