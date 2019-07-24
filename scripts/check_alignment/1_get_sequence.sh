#!/bin/bash



twobit_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
seqList_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/seqList_for_twoBitToFa_2019-07-24.txt"
output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/seqList_hg19_sequences.fa"


echo "twoBitToFa $twobit_file  -seqList=$seqList_file $output_file"
twoBitToFa $twobit_file  -seqList=$seqList_file $output_file

echo "DONE!"
