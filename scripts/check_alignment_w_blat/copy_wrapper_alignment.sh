#!/bin/bash

now=$(date +"%Y_%m_%d")

SLURM_ARRAY_TASK_ID=1
my_array=(dummy 150)
flank_bp=${my_array[${SLURM_ARRAY_TASK_ID}]}


output_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_08_21_blat_pipeline/5e-6_gwas_threshold/${flank_bp}bp_flanking_region"
echo "writing output to: ${output_dir}"
echo " "
echo " "

#0
# ASSOC_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/20190722_gwas_pca12_age.csv"
ASSOC_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL.csv"
GWAS_PTHRESH=0.000005
PLINK_OUTPUT_DIR="${output_dir}/seqList_for_twoBitToFa_${now}.tsv"
PLINK_OUTPUT_FILE_KEY="${output_dir}/key_for_seqList_for_twoBitToFa_${now}.tsv"
FLANK_BP=${flank_bp%bp}
# echo "python 0_seqlist_for_2bit.py $ASSOC_FILE $GWAS_PTHRESH $FLANK_BP $PLINK_OUTPUT_DIR $PLINK_OUTPUT_FILE_KEY"
python 0_seqlist_for_2bit.py $ASSOC_FILE $GWAS_PTHRESH $FLANK_BP $PLINK_OUTPUT_DIR $PLINK_OUTPUT_FILE_KEY



#1
SEQ_OUTPUT="${output_dir}/seqList_hg19_sequences.fa"
./1_get_sequence.sh $PLINK_OUTPUT_DIR $SEQ_OUTPUT


#2
BLAT_OUTPUT="${output_dir}/blat_results_${now}.txt"
./2_blat.sh $SEQ_OUTPUT $BLAT_OUTPUT


#3
FORMATTED_BLAT_OUTPUT="${output_dir}/blat_results_w_assoc_${now}.tsv"
python 3_format_blat_results.py $BLAT_OUTPUT $PLINK_OUTPUT_FILE_KEY $FORMATTED_BLAT_OUTPUT



