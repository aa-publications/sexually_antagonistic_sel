#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --mem=25G
#SBATCH --output=blat_mega
#SBATCH --job-name=blat_mega
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=debug
# database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
# query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/seqList_hg19_sequences.fa"
# output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/blat_results_$(date +%F).txt"

database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
# query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/MEGAEX_probes_v1.fa"
# query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/test.fa"
# output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/blat_psl/MEGAEx_probes_v1_blat.txt"


query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/sig_hits_biovu_probes.fa"
output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/webblat_wublast_format/MEGAEx_probes_v1_blat_wublast.txt"

echo "Running:"
echo "blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 $database_file $query_file $output_file -out=wublast"
/dors/capra_lab/opt/kent-tools/blat/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 $database_file $query_file $output_file -out=wublast



echo "Done! Check ${output_file}"
