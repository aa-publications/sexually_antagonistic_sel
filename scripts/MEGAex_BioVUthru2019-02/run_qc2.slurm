#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=100G
#SBATCH --time=00:20:00
#SBATCH --job-name=qc2
#SBATCH --output=qc2_%A_%a.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL,END

### wrapper script for qc_per_batch.py
### This script will perform GWAS QC for one plink dataset (i.e. batch)
### REQUIRES: bfile_names.txt (one line per plink file prefix)
module load Anaconda3
source activate py36
module load PLINK/1.9b_5.2

python qc2_MEGAex_BioVUthru2019-02.py
