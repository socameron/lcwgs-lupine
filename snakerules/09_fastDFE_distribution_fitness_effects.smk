## ANALYSIS 9: Distribution of Fitness Effect (DFE)

# Identify 4fold degenerate sites for inputs into ANGSD (which determines the SFS)
rule find_4fold_degenerate_sites:
  input:
    gff="data/annotation/Lupinus_perennis_genes.gff3",
    fasta="data/reference/hap2/lupinehap2.fasta"
  output:
    fourfold="results/DFE/sites_4D.txt",
    other="results/DFE/sites_CDSother.txt"
  log:
    "results/logs/DFE/find_4D_sites.log"
  shell:
    """
    python scripts/find_4fold_degen.py {input.gff} {input.fasta} > {log} 2>&1
    mv sites_4D.txt {output.fourfold}
    mv sites_CDSother.txt {output.other}
    """

rule convert_DFE_sites_txt_to_bed:
  input:
    fourfold="results/DFE/sites_4D.txt",
    other="results/DFE/sites_CDSother.txt"
  output:
    bed_4fold="results/DFE/sites_4D.BED",
    bed_other="results/DFE/sites_CDSother.BED"
  log:
    "results/logs/bedtools/convert_4fold_and_CDS.log"
  shell:
    """
    awk -F '[:]' '{{OFS="\\t"; print $1, $2-1, $2}}' {input.fourfold} > {output.bed_4fold} 2> {log}
    awk -F '[:]' '{{OFS="\\t"; print $1, $2-1, $2}}' {input.other} > {output.bed_other} 2> {log}
    """


# Filter each population's site list with the 4D and CDS sites
rule filter_by_popln_4D_and_CDS:
  input:
    fourfold_bed="results/DFE/sites_4D.BED",
    cds_bed="results/DFE/sites_CDSother.BED",
    filtered_dupHMM_bed="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.BED",
    filtered_dupHMM_txt="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt"
  output:
    filtered_4D_bed="results/bed/hap{hap}/DFE_sites/{population}_4D_filtered.BED",
    filtered_CDS_bed="results/bed/hap{hap}/DFE_sites/{population}_CDS_other_filtered.BED",
    filtered_4D_txt="results/bed/hap{hap}/DFE_sites/{population}_4D_filtered.txt",
    filtered_CDS_txt="results/bed/hap{hap}/DFE_sites/{population}_CDS_other_filtered.txt"
  log:
    dupHMM_log="results/logs/bedtools/hap{hap}/DFE_sites/{population}_filtered_DFE.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Intersect dupHMM-filtered regions with 4D and other CDS sites
    bedtools intersect -a {input.filtered_dupHMM_bed} -b {input.fourfold_bed} > {output.filtered_4D_bed}
    bedtools intersect -a {input.filtered_dupHMM_bed} -b {input.cds_bed} > {output.filtered_CDS_bed}

    # Convert intersected BED files to ANGSD-compatible -sites txt files
    awk '{{print $1, $3}}' {output.filtered_4D_bed} > {output.filtered_4D_txt}
    awk '{{print $1, $3}}' {output.filtered_CDS_bed} > {output.filtered_CDS_txt}
    """


rule index_DFE_4fold_sites:
  input: 
    filtered_4D_txt="results/bed/hap2/DFE_sites/{population}_4D_filtered.txt"
  output: 
    bin_index="results/bed/hap2/DFE_sites/{population}_4D_filtered.txt.bin",
    idx_index="results/bed/hap2/DFE_sites/{population}_4D_filtered.txt.idx"
  envmodules:
    "angsd/0.940"
  shell: 
    """
    angsd sites index {input.filtered_4D_txt}
    """


rule index_DFE_CDS_other_sites:
  input: 
    filtered_CDS_txt="results/bed/hap2/DFE_sites/{population}_CDS_other_filtered.txt"
  output: 
    bin_index="results/bed/hap2/DFE_sites/{population}_CDS_other_filtered.txt.bin",
    idx_index="results/bed/hap2/DFE_sites/{population}_CDS_other_filtered.txt.idx"
  envmodules:
    "angsd/0.940"
  shell: 
    """
    angsd sites index {input.filtered_CDS_txt}
    """


rule angsd_for_DFE_4fold:
  input:
    bam_list="data/lists/hap2/{population}_clipped_hap2.txt",
    canonical_sites="results/bed/hap2/DFE_sites/{population}_4D_filtered.txt",
    bin_index="results/bed/hap2/DFE_sites/{population}_4D_filtered.txt.bin",
    fasta_fai="data/reference/hap2/lupinehap2.fasta.fai"
  output:
    arg_file="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites.arg",
    mafs_file="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites.mafs.gz",
    saf_1="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites.saf.idx",
    saf_2="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites.saf.pos.gz",
    saf_3="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites.saf.gz",
    glf_file="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites.glf.gz"
  params:
    ref="data/reference/hap2/lupinehap2.fasta",
    file_name="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites",
    scaffolds="results/scaffolds/hap2_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap2/canonical/DFE/by_popln/angsd_4fold_sites_hap2_{population}.log"
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
    -nThreads {threads}\
    -uniqueOnly 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doGlf 1\
    -doSaf 1\
    -doCounts 1\
    -anc {params.ref}\
    &> {log}
    """


# We do not filter for MAFs
rule angsd_for_DFE_CDS:
  input:
    bam_list="data/lists/hap2/{population}_clipped_hap2.txt",
    canonical_sites="results/bed/hap2/DFE_sites/{population}_CDS_other_filtered.txt",
    bin_index="results/bed/hap2/DFE_sites/{population}_CDS_other_filtered.txt.bin",
    fasta_fai="data/reference/hap2/lupinehap2.fasta.fai"
  output:
    arg_file="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites.arg",
    mafs_file="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites.mafs.gz",
    saf_1="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites.saf.idx",
    saf_2="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites.saf.pos.gz",
    saf_3="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites.saf.gz",
    glf_file="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites.glf.gz"
  params:
    ref="data/reference/hap2/lupinehap2.fasta",
    file_name="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites",
    scaffolds="results/scaffolds/hap2_scaffolds.txt",
    minInd=lambda wildcards, input: max(1, int(0.8 * sum(1 for _ in open(input.bam_list))))
  log:
    "results/logs/angsd/hap2/canonical/DFE/by_popln/angsd_CDS_other_hap2_{population}.log"
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
    -nThreads {threads}\
    -uniqueOnly 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doGlf 1\
    -doSaf 1\
    -doCounts 1\
    -anc {params.ref}\
    &> {log}
    """
    

rule SFS_4fold_by_population:
  input:
    saf_idx="results/angsd/hap2/canonical/DFE/by_popln/{population}_4fold_sites.saf.idx",
    canonical_sites="results/bed/hap2/DFE_sites/{population}_4D_filtered.txt",
    ref="data/reference/hap2/lupinehap2.fasta"
  output:
    sfs="results/realSFS/hap2/DFE/by_popln/{population}_DFE_folded_4fold.sfs"
  log:
    "results/logs/realSFS/hap2/DFE/by_popln/{population}_DFE_folded_4fold.log"
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



rule SFS_CDSother_by_population:
  input:
    saf_idx="results/angsd/hap2/canonical/DFE/by_popln/{population}_CDS_other_sites.saf.idx",
    canonical_sites="results/bed/hap2/DFE_sites/{population}_CDS_other_filtered.txt",
    ref="data/reference/hap2/lupinehap2.fasta"
  output:
    sfs="results/realSFS/hap2/DFE/by_popln/{population}_DFE_folded_CDS.sfs"
  log:
    "results/logs/realSFS/hap2/DFE/by_popln/{population}_DFE_folded_CDS.log"
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


# Run DFE - note requires Python 3.10 or 3.12
# NOTE: MUST RUN SNAKEMAKE RULE USING VIRTUAL ENVIRONMENT envs/fastdfe-py312
rule fastdfe_by_population:
  input:
    neut = "results/realSFS/hap2/DFE/by_popln/{population}_DFE_folded_4fold.sfs",
    sel  = "results/realSFS/hap2/DFE/by_popln/{population}_DFE_folded_CDS.sfs"
  output:
    png   = "results/plots/hap2/fastdfe/{population}_fastdfe_plot.png",
    csv   = "results/fastdfe/hap2/by_popln/{population}_fastdfe_discretized.csv",
    json  = "results/fastdfe/hap2/by_popln/{population}_fastdfe_inference.json",
    txt   = "results/fastdfe/hap2/by_popln/{population}_fastdfe_summary.txt"
  log:
    "results/logs/fastdfe/hap2/by_popln/{population}.log"
  envmodules:
    "r/4.5.0",
    "python/3.12.4"
  threads: 2
  params:
    outdir    = "results/fastdfe/hap2/by_popln",
    png       = "results/fastdfe/hap2/by_popln/{population}_fastdfe_plot.png",
    png_location = "results/plots/hap2/fastdfe",
    n_runs    = 100,
    py        = "~/envs/fastdfe-py312/bin/python",
    Rscript     = "scripts/estimate_fastdfe_folded.R"
  shell:
    """
    set -euo pipefail
    mkdir -p $(dirname {log}) {params.outdir}
    MPLBACKEND=Agg \
    PYTHONNOUSERSITE=1 \
    RETICULATE_USE_UV=FALSE \
    RETICULATE_PYTHON="{params.py}" \
    Rscript --vanilla {params.Rscript} \
      --neut "{input.neut}" \
      --sel "{input.sel}" \
      --outdir "{params.outdir}" \
      --pop "{wildcards.population}" \
      --n_runs "{params.n_runs}" \
      > "{log}" 2>&1
    mv {params.png} {params.png_location}
    """





## ANALYSIS 10: Genotype Environment Association Test

# Data requirements include:
  # (1) Allele frequencies by population
  # (2) Historical and recent climate data
  # (3) Soil data collected from the field
  # (4)






## ANALYSIS 12 : SIFT analysis for mutational load

rule make_ploidy_from_bamlist:
  input:
    bamlist = "data/lists/hap2/{population}_clipped_hap2.txt"
  output:
    ploidy = "data/lists/hap2/ploidy/{population}_ploidy.txt"
  log:
    "results/logs/bcftools/hap2/ploidy/{population}.log"
  envmodules:
    "samtools/1.20"
  shell:
    r"""
    mkdir -p $(dirname {output.ploidy}) $(dirname {log})

    # Pull SM: from @RG lines; \K keeps only the part after SM:
    xargs -a {input.bamlist} -I{{}} bash -lc '
      samtools view -H "{{}}" |
      grep -oP "\tSM:\K[^\t]+" || true
    ' \
    | LC_ALL=C sort -u \
    | awk '{{print $1 "\t2"}}' > {output.ploidy} 2> {log}
    """



# NOTE: we restrict the analysis to the first 24 scaffolds
# SNP‐calling rule, one job per population:
rule bcftools_snp_call_by_population:
  input:
    ref     = "data/reference/hap2/lupinehap2.fasta",
    bamlist = "data/lists/hap2/{population}_clipped_hap2.txt",
    scaffolds = "results/scaffolds/hap2_scaffolds.txt"
  output:
    vcf = "results/bcftools/hap2/vcf_by_popln/raw/{population}.vcf.gz",
    tbi = "results/bcftools/hap2/vcf_by_popln/raw/{population}.vcf.gz.tbi"
  log:
    "results/logs/bcftools/SNP_call_{population}.log"
  threads: 12
  envmodules:
    "bcftools/1.22"
  shell:
    """
    bcftools mpileup -Ou \
      -f {input.ref} \
      --bam-list {input.bamlist} \
      -q 5 \
      --regions-file {input.scaffolds} \
      -I \
      -a FMT/AD,FMT/DP \
    | bcftools call \
      -V indels -f GQ -mv -Oz \
      --threads {threads} \
      -o {output.vcf} 2> {log}
    tabix -f {output.vcf}
    """

# NOTE: we restrict the analysis to non-paralogous regions and the first 24 scaffolds
rule bcftools_snp_call_canonical_by_population:
  input:
    ref     = "data/reference/hap2/lupinehap2.fasta",
    bamlist = "data/lists/hap2/{population}_clipped_hap2.txt",
    canonical_sites = "results/bed/hap2/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.BED"
  output:
    vcf = "results/bcftools/hap2/vcf_by_popln/canonical/{population}_canonical.vcf.gz",
    tbi = "results/bcftools/hap2/vcf_by_popln/canonical/{population}_canonical.vcf.gz.tbi"
  log:
    "results/logs/bcftools/SNP_call_{population}_canonical.log"
  threads: 12
  envmodules:
    "bcftools/1.22"
  shell:
    """
    bcftools mpileup -Ou \
      -f {input.ref} \
      --bam-list {input.bamlist} \
      -q 5 \
      --targets-file {input.canonical_sites} \
      --max-depth 3500 \
      -I \
      -a FMT/AD,FMT/DP \
      --threads {threads} \
    | bcftools call \
      -V indels -f GQ -mv -Oz \
      --threads {threads} \
      -o {output.vcf} 2> {log}
    tabix -f {output.vcf}
    """


## ANALYSIS 13: Phylogenetic tree based on ML

# samples: TSV mapping sample -> BAM path, pop -> population code
# ref: reference fasta


rule angsd_for_phylotree:
  input:
    bam_list="data/lists/hap2/all_populations_clipped_hap2.txt",
    canonical_sites="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt",
    bin_index="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.bin",
    idx_index="results/bed/hap2/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.idx"
  output:
    arg_file="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.arg",
    mafs_file="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.mafs.gz",
    hwe_file="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.hwe.gz",
    depth_sample="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.depthSample",
    depth_global="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.depthGlobal",
    saf_1="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.saf.idx",
    saf_2="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.saf.pos.gz",
    saf_3="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.saf.gz",
    beagle="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs.beagle.gz"
  params:
    ref="data/reference/hap2/lupinehap2.fasta",
    file_name="results/angsd/hap2/canonical/phylotree_input/all_poplns_{hap2scaffold}_canonical_SNPs",
    scaffolds="results/scaffolds/hap2_scaffolds.txt"
  log:
    "results/logs/angsd/hap2/canonical/phylotree_input/angsd_{hap2scaffold}_canonical_SNPs_all_poplns.log"
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
    -uniqueOnly 1\
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
    &> {log}
    """


rule ngsdist_pairwise:
  input:
    beagle = "results/angsd/hap2/gl_tree/all_samples.beagle.gz",
    mafs   = "results/angsd/hap2/gl_tree/all_samples.mafs.gz",
    labels = "data/sample_labels.txt"   # one label per line, in same order as beagle individuals
  output:
    dist = "results/tree/ngsdist/distances.mat"
  log:
    "results/logs/tree/hap2/ngsdist.log"
  threads: 16
  envmodules:
    "ngstools/1.0"  # adjust to your module naming
  shell:
    r"""
    mkdir -p $(dirname {output.dist}) $(dirname {log})
    ngsDist --geno {input.beagle} --probs \
            --labels {input.labels} \
            --out {output.dist} --n_threads {threads} \
            > {log} 2>&1
    """

rule fastme_tree:
  input:
    dist = "results/tree/ngsdist/distances.mat"
  output:
    tree = "results/tree/ngsdist/tree_fastme.nwk"
  log:
    "results/logs/tree/hap2/fastme.log"
  envmodules:
    "fastme/2.1.6"
  shell:
    r"""
    fastme -i {input.dist} -o {output.tree} -n > {log} 2>&1
    """

















## ARCHIVE

# Previously, the @SM sample metadata was not the same for re-sequenced libraries. 
# This mean't that .bam files from the same sample were mistakenly taken as seaprate indiiduals.
# This problem is now fixed!
# Below is the rule to fix the .bam files. 

rule fix_sm_lb_per_bam:
  input:
    bam="{bam}"
  output:
    done="{bam}.smfix.DONE"
  log:
    "{bam}.smfix.log"
  envmodules:
    "samtools/1.20"
  threads: 2
  shell:
    r"""
    set -euo pipefail
    BAM="{input.bam}"
    base="$(basename "$BAM")"
    sample="${{base%%_hap2_realign.bam}}"
    prefix="${{sample%%_*}}"
    lb="${{prefix}}_LB"

    tmp_hdr="$(mktemp)"
    samtools view -H "$BAM" \
    | awk -v sm="$prefix" -v lb="$lb" -F'\t' 'BEGIN{{OFS="\t"}}
        $1=="@RG"{{smset=0; lbset=0;
          for(i=1;i<=NF;i++){{ if($i ~ /^SM:/){{$i="SM:" sm; smset=1}}
                                else if($i ~ /^LB:/){{$i="LB:" lb; lbset=1}} }}
          if(!lbset){{$0=$0 OFS "LB:" lb}}
          if(!smset){{$0=$0 OFS "SM:" sm}}
          print; next}}
        {{print}}
      ' > "$tmp_hdr" 2>> {log}

    if samtools reheader -i "$tmp_hdr" "$BAM" 2>> {log}; then
      samtools index -@4 "$BAM" 2>> {log}
      rm -f "$tmp_hdr"
      echo "In-place OK: $BAM (SM=$prefix LB=$lb)" >> {log}
    else
      out="${{BAM}}.tmp.smfix.bam"
      echo "In-place failed; rewriting to $out" >> {log}
      samtools reheader "$tmp_hdr" "$BAM" > "$out" 2>> {log}
      samtools index -@4 "$out" 2>> {log}
      mv -f "$out" "$BAM"
      mv -f "$out.bai" "$BAM.bai"
      rm -f "$tmp_hdr"
    fi

    touch {output.done}
    """
