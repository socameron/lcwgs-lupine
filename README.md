# Population Genetics - _Lupinus perennis_ (lcWGS)

This is a Snakemake workflow for running the population genetic analyses. It includes the following analyses:

1. Relatedness structure for identifying clones (software `ngsRelate`)
2. Individual-level inbreeding coeffcieints (software `ngsF-HMM`)
3. Runs of Homozygosity (ROH) for inbreeding (software `bcftools roh`
4. Population structure and admixture (software `PCangsd`)
5. Thetas for Taijama's D, nucleotide diversity, and Watterson's theta (software `angsd` with command `thetaStat`) 
6. Population-level pairwise Fst and Fst by Distance (IBD) (software `realSFS`)
7. Contemporary effective population size (software `GONE`)
8. Linkage disequilibrium (software `ngsLD`)
9. Mutation screens (dN/dS) (custom software, email dan.schoen@mcgill.ca for info)

## About the _Lupinus perennis_ dataset

The _Lupinus perennis_ genomic dataset represents 19 populations sampled across the southern Ontario and the US Midwest. I collected leaf samples from 22+ Wild lupine (_Lupinus perennis_) populations in 2022 and 2023 and sequenced 19 populations using low coverage whole genome sequencing (lcWGS; about 10x coverage). Some populations were sampled by colllaborators in American state agencies. 

The [University of Wisconsin Biotechnology Center DNA Sequencing Facility](https://dnaseq.biotech.wisc.edu/) completed our DNA extractions and the [McGill Genome Centre](https://www.mcgillgenomecentre.ca/) performed the sequencing using Illumina NovaSeq6000. 

Populations are distinguished based on their geographical position relative to the overall species' distribution. 8 populations are considered 'core' populations, and the other 11 are northern 'edge' populations. No populations in the southern range were sampled (i.e towards Florida). 

![Sequence Map](https://github.com/socameron/lcwgs-lupine/blob/main/19%20populations%20to%20sequence.png?raw=true)

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


## Some other notes

I use an `all` rule in Snakemake so I don't have to specify a file to request. Rather, I expand all of the files that I'm requesting in rule `all`. It also makes it easier to type in the command line.

In case you're also new to the world of bioinformatics, I also use Cyberduck to access my data plus VS Code to edit code directly on DRAC servers.  

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

[PCAngsd Version 1.2](https://github.com/Rosemeis/pcangsd)

[My Website](https://www.cameronso.ca)

[Snakemake Tutorial for DRAC clusters](https://github.com/socameron/snakemake-tutorial)

