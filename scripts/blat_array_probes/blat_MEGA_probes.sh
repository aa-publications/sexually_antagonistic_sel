#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --mem=25G
#SBATCH --output=blat_mega.out
#SBATCH --job-name=blat_mega
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=debug
# database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
# query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/seqList_hg19_sequences.fa"
# output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/blat_results_$(date +%F).txt"

database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"



# BLAT OUTPUT
# query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/check_probe_seqs/MEGAEX_probes_v1.fa"
# output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/check_probe_seqs/MEGAEx_probes_v1_blat.txt"
#
# echo "Running:"
# echo "blat $database_file $query_file $output_file -out=blast"
# # /dors/capra_lab/opt/kent-tools/blat/blat $database_file $query_file $output_file -out=blast8



# PSL OUTPUT

query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/sig_hits_biovu_probes.fa"
psl_output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/MEGAEx_probes_v1_blat_wublast.txt"

echo "blat $database_file $query_file $psl_output_file -out=pslx"
/dors/capra_lab/opt/kent-tools/blat/blat $database_file $query_file $psl_output_file -out=wublast
echo "Done! Check ${psl_output_file}"
