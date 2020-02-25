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


bfile="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/impute/final_imputed/test/mega_test"
outfile="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/impute/final_imputed/test/pl2_ffallback"

plink2 --bfile $bfile --glm no-x-sex firth-fallback --vif 1500 --out $outfile

# Pass your Python script to the Anaconda3 python intepreter for execution
# python test_for_assoc.py
