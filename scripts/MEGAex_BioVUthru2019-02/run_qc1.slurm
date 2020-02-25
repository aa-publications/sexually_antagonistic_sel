#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=200G
#SBATCH --time=1-02:00:00
#SBATCH --job-name=qc
#SBATCH --output=qc_%A_%a.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL,END

### wrapper script for qc_per_batch.py
### This script will perform GWAS QC for one plink dataset (i.e. batch)
### REQUIRES: bfile_names.txt (one line per plink file prefix)

source activate py36_r_ml
module load PLINK/1.9b_5.2

data_dir="/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/MEGAex_BioVUthru2019-02_BestOfMultipleCalls"
bfile="MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"
# bfile="temp"
output_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02"

# python qc1_MEGAex_BioVUthru2019-02.py $bfile $data_dir $output_dir

# data_dir=""
# bfile="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/maf_filtered/final_maf_filtered__final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"



data_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02_keep_missingness/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc"
bfile="final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"
output_w_missing_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02_keep_missingness"
python qc1_MEGAex_BioVUthru2019-02.py $bfile $data_dir $output_w_missing_dir
