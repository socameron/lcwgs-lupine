# Population Genetics - _Lupinus perennis_ (lcWGS)

This is a Snakemake workflow for running the population genetic analyses on low to medium coverage data (ranges from 2x to 20x, with mean 10x). It primarily uses genotype likelihoods (estimated from ANGSD) rather than genotype calls (e.g PLINK or bcftools). It includes the following analyses:

1. Ortholog finder between _L. perennis_ and _L. mutabilis_* plus WG alignment (software [OrthoFinder](https://github.com/davidemms/OrthoFinder) and [MUMer](https://mummer4.github.io/)
2. Relatedness structure for identifying clones (software [ngsRelate](https://github.com/ANGSD/NgsRelate))
3. Individual-level inbreeding coefficients (software [ngsF-HMM](https://github.com/fgvieira/ngsF-HMM))
4. Runs of Homozygosity (ROH) for inbreeding (software [RZooRoH](https://cran.r-project.org/web/packages/RZooRoH/index.html))
5. Population structure and admixture (software [PCangsd](https://github.com/Rosemeis/pcangsd) and [ngsAdmix](https://github.com/aalbrechtsen/NGSadmix?tab=readme-ov-file))
6. Thetas for Taijama's D, nucleotide diversity, and Watterson's theta (software [ANGSD](https://github.com/ANGSD/angsd) with command `thetaStat`) 
7. Population-level pairwise Fst and Fst by Distance (IBD) (software [ANGSD](https://github.com/ANGSD/angsd) with command `realSFS`)
8. Linkage disequilibrium (software [ngsLD](https://github.com/fgvieira/ngsLD))
9. Contemporary effective population size (software [currentNe2](https://github.com/esrud/currentNe2))
10. Distribution of fitness effects (software [fastdfe](https://fastdfe.readthedocs.io/en/latest/))
11. Genotype-environment associations (software [Redundancy Analysis](https://github.com/Capblancq/RDA-landscape-genomics))

*(to include _L. polyphyllus_ when available)


Some analyses will be implemented include:

1. Mutation screens (dN/dS) (custom software, email dan.schoen@mcgill.ca for info)
2. Genetic off-sets (software [Redundancy Analysis](https://github.com/Capblancq/RDA-landscape-genomics))
3. Demographic histories (software [dadi](https://dadi.readthedocs.io/en/latest/))

NOTE: since there are so many analyses, one day I'll break up the Snakefile into separate components per each analysis.

## About the _Lupinus perennis_ dataset

The _Lupinus perennis_ genomic dataset represents 18 populations sampled across southern Ontario in Canada and the US Midwest. I collected leaf samples from 22+ Wild lupine (_Lupinus perennis_) populations in 2022, 2023, and 2024 and sequenced 18 populations using low coverage whole genome sequencing (lcWGS; about 10x coverage). Some populations were sampled by colllaborators from US state agencies. 

For 7 populations, the [University of Wisconsin Biotechnology Center DNA Sequencing Facility](https://dnaseq.biotech.wisc.edu/) completed our DNA extractions and the [McGill Genome Centre](https://www.mcgillgenomecentre.ca/) performed the sequencing using Illumina NovaSeq6000. For the remaining 11 populations, [Génome Québec](https://genomequebec.com/en/) completed our DNA extractions and sequencing using Illumina NovaSeq X+. 

Populations are distinguished based on their geographical position relative to the overall species' distribution. 8 populations are considered 'core' populations, and the other 10 are northern 'edge' populations. No populations in the southern range were sampled (i.e towards Florida). 

![Sequence Map](https://github.com/socameron/lcwgs-lupine/blob/f8d73f7020cf4419618df976f589a47877c65d3b/GEA_sampling_figure.svg)

## Quality control prior to analyses

Before running any analyses, I mapped my reads to a scaffold-based reference genome. Since I'm using ANGSD for a majority of my analyses, I also clipped overlapping reads with `bamutil/1.0.14`. Previously I realigned reads around indels using `gatk/3.8` but I no longer include this into my workflow. Realignment around indels is dealt with during ANGSD calls with `-baq 2`. I then filtered out potential paralogs using [ngsParalog](https://github.com/tplinderoth/ngsParalog) and its internal program called `dupHMM`. Please see the file 'Snakefile_dataprep'. Any downstream analyses then underwent its own filtering criteria, but was required to constrain the analysis to the first 24 scaffolds of our reference assembly, as these in theory should represent the 24 chromosomes of _Lupinus perennis_. In the future, we might utilize our Hi-C data to construct a chromosome-based reference assembly. 

## Snakemake

My analyses uses the Snakemake workflow and is adapted for reproducibility on [Digital Research Alliance of Canada](alliancecan.ca/en) clusters. I specifically used the #Beluga cluster because it's located in Montreal. Setting up this workflow would not have been possible without [Eric Anderson](https://github.com/eriqande/mega-lcwgs-pw-fst-snakeflow) (Research Geneticist, NOAA) whom I met at [ConGen2023](https://www.umt.edu/ces/conferences/congen/). For more information about Snakemake's amazing abilities, see Eric Anderson's lecture slides [here](https://eriqande.github.io/con-gen-2023/slides/snake-slides.html#/section). Please review my Snakemake tutorial if you'd like to setup the workflow manager on DRAC clusters.

## Python and R

If you take a look into my scripts folder, you'll notice a lot of other python and R scripts. I mostly use Python to run analyses that do NOT require graphing for publication/presentations. On the other hand, I use R largely for small-data processing and making pretty graphs. R's 'ggplot2' package has some pretty settings compared to Python. 

Some Python and R packages that I require include:

```
# Python packages
matplotlib
biopython
numpy
multiqc
snakemake

# R packages
tidyverse
purrr
Rcpp
geosphere
RColorBrewer

```


## Snakemake organization

Each analysis is written in separate `.smk` files (e.g `01_ngsRelate_relatedness.smk`), which are then compiled into the working `Snakefile`. This reduces the length of code in a single script. I also use an `all` rule to expand the number of files that I'd like to create.

## Installing PCAngsd

To use PCAngsd, you'll need to download it from github and install it with a compiler. 

```
module load gcc python/3.10 angsd
source ENV/bin/activate
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
pip install -r requirements.txt
python setup.py build_ext --inplace
pip install .
```

There are also a lot of other packages that need custom installation. If you need help, please contact your university/cluster services or email me at cameron.so@mail.mcgill.ca.

## Additional resources and references

[A beginner's guide to low-coverage whole genome sequencing for population genomics](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.16077)

[My Website](https://www.cameronso.ca)

[Snakemake Tutorial for DRAC clusters](https://github.com/socameron/snakemake-tutorial)

