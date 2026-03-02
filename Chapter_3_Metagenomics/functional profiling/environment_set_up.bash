###############################################
# 1) FUNCTIONAL PROFILING SETUP
#
# Description:
#   Installs Miniforge and Mamba, and creates:
#     - shotgun_meta environment (fastqc, multiqc, kraken, etc.)
#     - biobakery environment (HUMAnN, MetaPhlAn, etc.)
#
# Notes:
#   This script is based on NEOF's Microbial Shotgun Genomics course:
#   https://neof-workshops.github.io/Shotgun_ld5ug2/Course/01-Shotgun_metagenomics.html
#
#   HUMAnN databases must be downloaded separately prior to running
#   this script (see section A.1.3, biobakery3.9, in the NEOF bookdown).
###############################################

### Set base project directory
PROJECT_DIR=~/sonip-mtg
DB_DIR=$PROJECT_DIR/databases
MINIFORGE_DIR=$PROJECT_DIR/miniforge3

mkdir -p "$PROJECT_DIR" "$DB_DIR"
cd "$PROJECT_DIR"

### Install Miniforge 
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    -O Miniforge3.sh

bash Miniforge3.sh -b -p "$MINIFORGE_DIR"

# Activate conda in shell
"$MINIFORGE_DIR/bin/conda" init

source ~/.bashrc

conda activate base
conda update -n base -c defaults conda -y

### Create shotgun_meta environment
mamba create -n shotgun_meta -y
mamba activate shotgun_meta

# Installing select tools from those listed in tutorial
mamba install -c bioconda \
    fastqc trim-galore multiqc bowtie2 kraken2 \
    krona bracken lefse flash \
    -y

# Update Krona taxonomy
ktUpdateTaxonomy.sh

# Deactivate
mamba deactivate

### Create biobakery environment
mamba create -n biobakery python=3.7 -y
mamba activate biobakery

# Install HUMAnN
mamba install -c biobakery humann=3.9 -y

# Install supporting tools
mamba install bioconda::hclust2=1.0.0 -y
mamba install -c conda-forge matplotlib -y

# Install MetaPhlAn
metaphlan --install \
    --index mpa_vOct22_CHOCOPhlAnSGB_202403

### HUMAnN database configuration (point to your existing DBs)
humann_config --update database_folders nucleotide "$DB_DIR/humann/chocophlan"
humann_config --update database_folders protein "$DB_DIR/humann/uniref"
humann_config --update database_folders utility_mapping "$DB_DIR/humann/utility_mapping"

# Test HUMAnN installation
humann_test

# Export environment file
conda env export > biobakery.yml

mamba deactivate

### Cleanup
conda clean --all -y
