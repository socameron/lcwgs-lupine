# Population Genetics - _Lupinus perennis_ (lcWGS)

This is a Snakemake workflow for creating PCA plots (so far) from raw fastq.gz files using PCAngsd. 

## About the _Lupinus perennis_ dataset

I collected dried leaf samples from 22 Wild lupine (_Lupinus perennis_) populations in 2022 and 2023. 16 populations, all of which are remnant, were selected for low coverage whole genome sequencing (lcWGS).
Leaf samples were sent to [University of Wisconsin Biotechnology Center DNA Sequencing Facility](https://dnaseq.biotech.wisc.edu/) for DNA extractions. After which, sequencing was completed at the [McGill Genome Centre](https://www.mcgillgenomecentre.ca/). 

The distinction between populations is their geographical position relative to the species' distribution. 8 populations are considered 'core' populations, and the other 8 are northern 'edge' populations. No populations in the southern range were sampled (i.e towards Florida). 

## Snakemake workflow

This Snakemake workflow is adapted for use on [Digital Research Alliance of Canada](alliancecan.ca/en) clusters. I specifically used the #Beluga cluster because it's located in Montreal. Setting up this workflow would not have been possible without [Eric Anderson](https://github.com/eriqande/mega-lcwgs-pw-fst-snakeflow) (Research Geneticist, NOAA) whom I met at [ConGen2023](https://www.umt.edu/ces/conferences/congen/). For more information about Snakemake's amazing abilities, see Eric Anderson's lecture slides [here](https://eriqande.github.io/con-gen-2023/slides/snake-slides.html#/section). 

**TL;DR** Snakemake allows you to parallelize multiple jobs at once, while submitting jobs once tasks are complete. Snakemake rules are immensely easy to read and makes research more reproducible.

Setting up this workflow on DRAC clusters is simplier than you think:

1. First load Python and create and a virtual environment.

```
# We create a virtual environment in the scratch folder. I still haven't figured out if it's a better idea to do this in your home directory.
cd scratch 
module load python/3.10
virtual --no-download ENV

# Just labelled the virtual environment folder as 'ENV'
# Activate virtual environment

source ENV/bin/activate
```

2. Install and upgrade pip (python installer) plus Snakemake

```
pip install --no-index --upgrade pip
pip install --no-index snakemake
pip install --no-index multiqc
pip install --no-index matplotlib

# Also installing multiqc and matplotlib for further down analyses

# Double check version of Snakemake
snakemake --version
```
3. Download Snakemake files from my repository

Using Snakemake on DRAC clusters requires 4 items. (1) The `Snakemake` file; (2) a snakeprofile folder containing (3) a `config.yaml` and (4) `status.sacct-robust.sh` file. The Snakemake file is where all rules are written for your workflow. The `config.yaml` specifies the resources allocated for each rule (although this can also be written into your Snakemake file) plus other system settings. The `status-sacct-robust.sh` file interacts Snakemake with the SLURM scheduler to send jobs via `sbatch`. 

```
git clone https://github.com/socameron/lcwgs-lupine.git
```

4. Check that data files have readable permissions

Wherever your storing your data, say in `scratch/data`, then use `chmod ug+x *.fastq.gz` within that folder (logged onto your cluster, of course) to change permission rules. 

5. Call Snakemake using the profile set up

```
snakemake --profile snakeprofile all
# Calls snakemake using --profile snakeprofile. Tells Snakemake to specific run the rule 'all'.
# Add -np to run a practice run
```

## Some other notes

I use an `all` rule in Snakemake so I don't have to specify a file to request. Rather, I expand all of the files that I'm requesting in rule `all`. It also makes it easier to type in the command line.

In case you're also new to the world of bioinformatics, I also use Cyberduck to access my data plus VS Code to edit code directly on DRAC servers.  

## Installing PCAngsd

To use PCAngsd, you'll need to install it from github and install it with a compiler. 

```
module load gcc python/3.10 angsd
source ENV/bin/activate
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
pip install -r requirements.txt
python setup.py build_ext --inplace
pip install .
```


