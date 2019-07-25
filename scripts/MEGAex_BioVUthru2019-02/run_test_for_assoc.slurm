#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=4:50:00
#SBATCH --mem=12G
#SBATCH --output=assoc_%J.out
#SBATCH --job-name=assoc

module load Anaconda3
module load PLINK/2.00-alpha1
source activate py36_r_ml



# Pass your Python script to the Anaconda3 python intepreter for execution
# python test_for_assoc.py


PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/maf_filtered/final_maf_filtered__final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"
# COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/covariates/covar_white_YOB_5PC__2019-07-05.tsv"
COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/covariates/covar_white_YOB_12PC__2019-07-21.tsv"
LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic"



# -----------
# MAIN
# -----------
#

# run using plink2 w/ glm
# sex is not included for ANY of the variants with ' no-x-sex '
# covar-variance-standardize mean 0, variance 1

echo "plink2 --bfile ${PLINK_FILE} --covar ${COVAR_FILE} --glm firth-fallback intercept --covar-variance-standardize YOB --vif 1500 --out ${LOGISTIC_RESUTLS_OUTPUT}"

#plink2 --bfile ${PLINK_FILE} --covar ${COVAR_FILE} --glm firth-fallback intercept --covar-variance-standardize YOB --vif 1500 --out ${LOGISTIC_RESUTLS_OUTPUT} >> glm.out



plink --bfile ${PLINK_FILE} --covar ${COVAR_FILE} --logistic --out ${LOGISTIC_RESUTLS_OUTPUT} >> glm.out
