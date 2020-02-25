#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=4:50:00
#SBATCH --mem=12G
#SBATCH --output=assoc_%J.out
#SBATCH --job-name=assoc remote



# RUN STATS ON SIGNIFICANT SNPS
module load Anaconda3
# module load PLINK/2.00-alpha1
source activate py36_r_ml



# Pass your Python script to the Anaconda3 python intepreter for execution
# python test_for_assoc.py


PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/maf_filtered/final_maf_filtered__final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"
LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic"
bv_sig_snps_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/MEGAex_BioVUthru2019-02"

cases_ids="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/maf_filtered/cases_fam_ids_for_plink.txt"
control_ids="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/maf_filtered/controls_fam_ids_for_plink.txt"

# -----------
# MAIN
# -----------
#





# plink --bfile ${PLINK_FILE} --covar ${COVAR_FILE} --extract ${bv_sig_snps_file} --freqx --hardy midp --missing --test-mishap --mendel --het --ibc --fst case-control >> bv_sig_snps.out

plink --bfile ${PLINK_FILE} --keep ${cases_ids} --freq --hardy midp --missing  --out ${LOGISTIC_RESUTLS_OUTPUT}"/2019_07_21_logistic_cases_frq_hwe.out"
plink --bfile ${PLINK_FILE} --keep ${control_ids} --freq --hardy midp --missing  --out ${LOGISTIC_RESUTLS_OUTPUT}"/2019_07_21_logistic_controls_frq_hwe.out"
