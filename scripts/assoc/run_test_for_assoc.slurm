#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=1-23:10:00
#SBATCH --mem=18G
#SBATCH --output=assoc_%J
#SBATCH --job-name=assoc

module load Anaconda3
module load PLINK/2.00-alpha1
source activate py36_r_ml


plink2 --bfile /dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/merge_batch_without_batch4/final_maf_filtered__merged_MEGA_no_batch4_2019-04-15 --covar /dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/covariates/no_batch4_merged_pca/covar_CEU_no_batch4_merged_mega_w_age_5PC_batch_2019-04-18.tsv --glm no-x-sex  hide-covar firth-fallback --vif 1500 --out /dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_18_glm_no_batch4/no_batch4_glm  >> glm.out

# Pass your Python script to the Anaconda3 python intepreter for execution
# python test_for_assoc.py
