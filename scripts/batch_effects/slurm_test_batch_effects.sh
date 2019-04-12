#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=25G
#SBATCH --time=0-06:00:00
#SBATCH --job-name=batch_ef
#SBATCH --array=2,3,5,6,8,9,10,12
#SBATCH --output=batch_ef_%A_%a.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL,END

### wrapper script test_batch_effects.py
###         for each batch, test each snp for batch effects 
### REQUIRES: bfile_names.txt (one line per plink file prefix)

source activate py36_r_ml
module load PLINK/1.9b_5.2

list_of_batches_file="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/raw_preQC/cox_bfile_prefix.txt"
frq_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
shar_snps_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_shared_snps.txt" 
target_batch_prefix=$( awk -v line=${SLURM_ARRAY_TASK_ID} 'NR==line {print}' ${list_of_batches_file} )
output_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/batch_effects"

python test_batch_effects_parallel.py $frq_dir $shar_snps_file $target_batch_prefix $output_dir
