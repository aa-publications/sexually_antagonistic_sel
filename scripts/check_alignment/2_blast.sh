#!/bin/bash



database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/seqList_hg19_sequences.fa"

output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/blat_results_$(date +%F).txt"



echo "Running..."

echo "blat $database_file $query_file $output_file -out=blast"
/dors/capra_lab/opt/kent-tools/blat/blat $database_file $query_file $output_file -out=blast8


echo "Done!"
