#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=08:10:00
#SBATCH --mem=25G
#SBATCH --output=blat_ukk_wcsf
#SBATCH --job-name=blat_ukk_wcsf
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL




# database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
# query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/seqList_hg19_sequences.fa"
# output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/blat_results_$(date +%F).txt"

database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/fa/UKBB_WCSF_probes.fa"
output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/blat_output/UKBB_WCSF_probes_blat.txt"


echo " "
echo "Running:"
echo "blat $database_file $query_file $output_file -out=blast"
/dors/capra_lab/opt/kent-tools/blat/blat $database_file $query_file $output_file -out=blast8


echo "Done! Check ${output_file}"
