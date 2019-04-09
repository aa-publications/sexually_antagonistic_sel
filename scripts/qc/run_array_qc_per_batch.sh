#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --job-name=qc
#SBATCH --array=1-15
#SBATCH --output=qc_%A_task_%a.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL,END

### wrapper script for qc_per_batch.py 
### This script will perform GWAS QC for one plink dataset (i.e. batch)
### REQUIRES: bfile_names.txt (one line per plink file prefix)

source activate verBio
module load PLINK/1.9b_5.2

data_dir="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/preQC"
output_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
file=$( awk -v line=${SLURM_ARRAY_TASK_ID} 'NR==line {print}' ${data_dir}/cox_bfile_prefix.txt )

python qc_per_batch.py $file $data_dir $output_dir