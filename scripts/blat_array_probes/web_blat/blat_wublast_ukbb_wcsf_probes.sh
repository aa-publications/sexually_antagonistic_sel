#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --mem=25G
#SBATCH --output=blat_mega
#SBATCH --job-name=blat_mega
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=debug

database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
query_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/fa/split_to_n_files/ukbb_wcsf"


query_file=$( ls ${query_dir} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' )



query_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/fa/sig_hits_ukbb_wcsf_probes.fa"
output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/webblat_wublast_format/ukbb_wcsf_probes_v1_blat_wublast"


echo "/dors/capra_lab/opt/kent-tools/blat/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 $database_file $query_file $output_file -out=wublast"

/dors/capra_lab/opt/kent-tools/blat/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 $database_file $query_file $output_file -out=wublast
echo "Done! Check ${output_file}"

