#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --job-name=qc
#SBATCH --array=1-15
#SBATCH --output=qc_%A_task_%a.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL,END

### wrapper script for qc_per_batch.py 
### This script will perform GWAS QC for one plink dataset (i.e. batch)

iter=SLURM_ARRAY_TASK_ID

source activate verBio
module load PLINK/1.9b_5.2


data_dir="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/postQC"
# file="MEGA_ex_Array_Ancestry_MEGA_postQC_GRID"
output_dir="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/plink_qc_v2"


file=$( awk -v line=${SLURM_ARRAY_TASK_ID} 'NR==line {print}' ${data_dir}/bfile_names.txt )

python qc_per_batch.py $file $data_dir $output_dir