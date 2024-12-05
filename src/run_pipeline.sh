#!/usr/bin/env sh
set -e
module load /project/itmatlab/sharedmodules/use.shared
module load python/3.11
module load sratoolkit-2.11.0
module load edirect-15.3
module load STAR-v2.7.10b
module load bedtools2/2.30.0
module load bedGraphToBigWig/4.0
module load R/4.0.2
source ../.venv/bin/activate

### SETUP NOTE:
# In order to use the following --profile lsf command
# you need to follow the instructions at https://github.com/Snakemake-Profiles/lsf
# and set up LSF support for snakemake

bsub -e logs/snakemake.err \
    -o logs/snakemake.out \
    snakemake --profile lsf -j 100 -c 100 \
    --resources ncbi_download=3 \
    --use-singularity --singularity-args "-B /project/itmatlab/index/" \
    "$@"
