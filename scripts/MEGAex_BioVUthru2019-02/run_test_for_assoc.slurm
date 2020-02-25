#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=4:50:00
#SBATCH --mem=12G
#SBATCH --output=assoc_%J.out
#SBATCH --job-name=assoc 

module load PLINK/1.9b_5.2



# Pass your Python script to the Anaconda3 python intepreter for execution
# python test_for_assoc.py


# PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02_keep_missingness/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"
# COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/covariates/covar_white_YOB_5PC__2019-07-05.tsv"
COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/covariates/covar_white_YOB_12PC__2019-07-21.tsv"
LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic"
bv_sig_snps_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/MEGAex_BioVUthru2019-02/biovu_sig_snps.txt"

missingness_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/intermediate"
missingness_file=${missingness_dir}"/inter__cctrl_miss_snps_to_rm_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3.txt"

# -----------
# MAIN
# -----------
#

PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02_keep_missingness/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/temp_white_only_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"
LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2020_01_27_rerun_missingness_biovu_snps"

plink --bfile ${PLINK_FILE} --extract $missingness_file --covar $COVAR_FILE --logistic  --out ${LOGISTIC_RESUTLS_OUTPUT}"/rerun_missingness_assoc"
# 
# plink --bfile ${PLINK_FILE} --extract ${bv_sig_snps_file} --freqx --hardy midp --missing --test-mishap --mendel --het --ibc --fst case-control --out ${LOGISTIC_RESUTLS_OUTPUT}"/bv_sig_hits_stats"
# plink --bfile ${PLINK_FILE} --extract ${bv_sig_snps_file} --freq case-control --out ${LOGISTIC_RESUTLS_OUTPUT}"/bv_sig_hits_stats"
