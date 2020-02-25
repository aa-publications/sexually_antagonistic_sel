#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=04:10:00
#SBATCH --mem=25G
#SBATCH --job-name=blat_uk_bil
#SBATCH --output=slurm_out/blat_uk_bil_%A_%a.out
#SBATCH --array=1-26
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL

database_file="/dors/capra_lab/data/dna/human/hg19/hg19.2bit"
query_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/fa/split_to_n_files/ukbb_UKBiL"


query_file=$( ls ${query_dir} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' )

output_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/web_blat/uk_bil/webblat_${query_file}.out"


echo "/dors/capra_lab/opt/kent-tools/blat/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 $database_file $query_file $output_file"
/dors/capra_lab/opt/kent-tools/blat/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 $database_file ${query_dir}/$query_file $output_file
echo "Done! Check ${output_file}"

