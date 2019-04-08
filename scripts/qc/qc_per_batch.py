#!/bin/python
# This script will perform QC on each batch obtained from BioVU.
#
#
#
# Abin Abraham
# created on: 2019-04-04 08:30:42


import os
import sys
import numpy as np
import pandas as pd
import subprocess

import glob
import argparse
import time
from datetime import datetime

sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
from func_run_shell_cmd import run_shell_cmd
from func_rm_dups_individuals import rm_dup_individuals
from func_rm_discordant_sex import rm_discordant_sex
from func_rm_high_sample_missing_rate import rm_high_sample_missing_rate
from func_rm_high_het_rate import rm_high_het_rate
from func_rm_relatedness import rm_relatedness
from func_rm_high_variant_missing_rate import rm_high_variant_missing_rate
from func_rm_dups_varID import rm_dup_var_id
from func_add_pheno import add_pheno
from func_test_missing import test_missing
from func_track import get_num_lines

DATE = datetime.now().strftime('%Y-%m-%d')
start = time.time()

# # retrieve command line args 
# parser = argparse.ArgumentParser(description='Example with nonoptional arguments')
# parser.add_argument('bfile', action='store', type=str, help="plink bfile prefix")
# parser.add_argument('data_dir', action='store', type=str, help="dir to plink input files")
# parser.add_argument('output_dir', action='store')
  
# results = parser.parse_args()
 
# plink_prefix = results.bfile
# data_dir = results.data_dir
# output_dir = results.output_dir

# -----------
# FUNCTIONS
# -----------

def update_counts(label, fam_ct, bim_ct):

    TRACK_INDVID_DICT[label] = fam_ct
    TRACK_SNP_DICT[label] = bim_ct




# -----------
# MAIN
# -----------

plink_prefix = "MEGA_ex_Array_Ancestry_MEGA_preQC_GRID" # prefix before .bed/.bim/.ped
data_dir = "/Users/abin-personal/Documents/katja_biobank/katja_biobank/data" 
output_dir = "/Users/abin-personal/Documents/katja_biobank/katja_biobank/data"

print("Running qc_per_batch.py on {}...\
        \n\tfile: {}\n\tdata_dir: {}\n\toutput_dir: {}\n\n".format(DATE, plink_prefix, data_dir, output_dir))

raw_plink_prefix = os.path.join(data_dir, plink_prefix)
base_prefix = plink_prefix

TRACK_INDVID_DICT = {}
TRACK_SNP_DICT = {}
TRACK_INDVID_DICT['raw_data'] = get_num_lines(raw_plink_prefix+".fam")
TRACK_SNP_DICT['raw_data'] = get_num_lines(raw_plink_prefix+".bim")


#
#
#


# read .fam file
# write a new three column file with the third column being male or female 


#  QC ON SAMPLES
(dups_removed_plink_prefix, plink_stdout), fam_ct, bim_ct = rm_dup_individuals(raw_plink_prefix, output_dir, base_prefix)
update_counts('rm_dup_indvid', fam_ct, bim_ct)
(cleaned_sex_plink_prefix, all_stdout), fam_ct, bim_ct = rm_discordant_sex(dups_removed_plink_prefix, output_dir, base_prefix) # note: --set-hh-missing is run along with checking for sex
update_counts('rm_dscrd_sex_individ', fam_ct, bim_ct)
(ind_miss_removed_plink_prefix, rm_miss_indiv_stdout), fam_ct, bim_ct = rm_high_sample_missing_rate(cleaned_sex_plink_prefix, output_dir, base_prefix)
update_counts('rm_miss_indvid', fam_ct, bim_ct)
(het_removed_plink_file, ind_snp_output_file, all_stdout), fam_ct, bim_ct = rm_high_het_rate(ind_miss_removed_plink_prefix, output_dir, base_prefix)
update_counts('rm_het_indvid', fam_ct, bim_ct)
(related_indiv_removed_plink_prefix, all_stdout), fam_ct, bim_ct = rm_relatedness(het_removed_plink_file, output_dir,ind_snp_output_file, base_prefix)
update_counts('rm_related_indvid', fam_ct, bim_ct)

# add phenotype code 
(phenotyped_plink_prefix, make_pheno_stdout), fam_ct, bim_ct = add_pheno(related_indiv_removed_plink_prefix, output_dir, base_prefix, cases='female')
update_counts('add_pheno', fam_ct, bim_ct)
#  QC ON SNPS
(vars_miss_removed_plink_prefix, rm_miss_var_stdout),fam_ct, bim_ct = rm_high_variant_missing_rate(phenotyped_plink_prefix, output_dir, base_prefix)
update_counts('rm_miss_snps', fam_ct, bim_ct)
(no_dups_vars_plink_prefix, all_stdout),fam_ct, bim_ct = rm_dup_var_id(vars_miss_removed_plink_prefix, output_dir, base_prefix)
update_counts('rm_dup_snps', fam_ct, bim_ct)
(vars_miss_removed_plink_prefix, rm_miss_var_stdout),fam_ct, bim_ct = test_missing(no_dups_vars_plink_prefix, output_dir, base_prefix)
update_counts('rm_test_miss_snps', fam_ct, bim_ct)

# count lines in beween each step 

# mv files 1) .log, 2) tmp 3) inter 


#
# SUMMARIZE
#

columns = ["start", "removed_after_geno_mind_hwe", "removed_after_heterozygozity", "removed_after_relatedness", "final", "filename"]
indvid_line = [num_individ, 
                num_individ-num_individ_filter1, 
                num_individ_filter1-num_individ_after_het,
                num_individ_after_het-num_individ_after_relatives, 
                num_individ_after_relatives, basefile]

var_line = [num_variants,
                num_variants-num_variants_filter1,
                num_variants_filter1-num_variants_after_het,
                num_variants_after_het-num_variants_after_relatives,
                num_variants_after_relatives, basefile]

summary_file = os.path.join(output_dir, basefile+"_summary.tsv")
df = pd.DataFrame([indvid_line, var_line], index=['individual_count', 'variant_count'], columns=columns)
df.to_csv(summary_file, sep="\t")

#
# CLEAN UP TEMP FILES 
#

## make log director 
log_path=os.path.join(output_dir,'log')
summary_path=os.path.join(output_dir,'summary')

if not os.path.isdir(log_path):
    os.mkdir(log_path)

if not os.path.isdir(summary_path):
    os.mkdir(summary_path)

# move logs to log directory
for logfile in glob.glob(os.path.join(output_dir, "*.log")):
    os.rename(logfile, os.path.join(log_path, os.path.basename(logfile)))

# move summary files  to summary directory
for raw_plink_prefix in glob.glob(os.path.join(output_dir, "*_summary.tsv")):
    os.rename(raw_plink_prefix, os.path.join(summary_path, os.path.basename(raw_plink_prefix)))


# remove all temp files
[os.remove(x) for x in glob.glob(filter1_out_file + ".*")]
[os.remove(x) for x in glob.glob(ind_snp_out_file + ".*")]
[os.remove(x) for x in glob.glob(het_out_file+ ".*") ] 
[os.remove(x) for x in glob.glob(het_out_file+".het_exclude.tsv") ] 
[os.remove(x) for x in glob.glob(het_removed_out_file+".*") ] 
[os.remove(x) for x in glob.glob(related_genome_out_file+".*") ] 
[os.remove(x) for x in glob.glob(relatives_out_file) ] 
[os.remove(x) for x in glob.glob(rm_relatives_out_file+".nosex") ] 


print("Done processing {}. Took {} minutes.".format(plink_prefix, (time.time()-start)/60))