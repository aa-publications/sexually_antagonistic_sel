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
from collections import OrderedDict

# sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')

from func_calc_stats import calc_stats
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

# retrieve command line args
parser = argparse.ArgumentParser(description='Example with nonoptional arguments')
parser.add_argument('bfile', action='store', type=str, help="plink bfile prefix")
parser.add_argument('data_dir', action='store', type=str, help="dir to plink input files")
parser.add_argument('output_dir', action='store')
results = parser.parse_args()

plink_prefix = results.bfile
data_dir = results.data_dir
output_dir = results.output_dir

# -----------
# FUNCTIONS
# -----------

def update_counts(label, fam_ct, bim_ct):

    TRACK_INDVID_DICT[label] = fam_ct
    TRACK_SNP_DICT[label] = bim_ct


def safe_mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


# -----------
# MAIN
# -----------

#### !!!!! ##### !!!!! #### !!!!! ##### !!!!! #### !!!!! ##### !!!!! ##### !!!!! ##### !!!!!
#               ONLY FOR DEV/TESTING SCRIPTS... 
# plink_prefix = "MEGA_ex_Array_Ancestry_MEGA_preQC_GRID"  # prefix before .bed/.bim/.ped
# data_dir = "/Users/abin-personal/Documents/katja_biobank/katja_biobank/data"
# data_dir = "/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/raw_preQC"
# output_dir = "/Users/abin-personal/Documents/katja_biobank/katja_biobank/data"
# output_dir = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank"
#### !!!!! ##### !!!!! #### !!!!! ##### !!!!! #### !!!!! ##### !!!!! ##### !!!!! ##### !!!!!

# make a folder to hold all the analysis & outputs
os.mkdir(os.path.join(output_dir, plink_prefix+"_qc"))
output_dir = os.path.join(output_dir, plink_prefix+"_qc")

# =============  OUTPUT FILES =============
track_snp_df_out_file = os.path.join(output_dir, 'snp_count_{}.tsv'.format(plink_prefix))
track_indiv_df_out_file = os.path.join(output_dir, 'indvid_count_{}.tsv'.format(plink_prefix))


# =============  SET UP =============
print("Running qc_per_batch.py on {}...\
        \n\tfile: {}\n\tdata_dir: {}\n\toutput_dir: {}\n\n".format(DATE, plink_prefix, data_dir, output_dir))

raw_plink_prefix = os.path.join(data_dir, plink_prefix)
base_prefix = plink_prefix

TRACK_INDVID_DICT = OrderedDict()
TRACK_SNP_DICT = OrderedDict()


# =============  CREATE A COPY AND CALC BASIC STATS =============
(frq_and_miss_plink_prefix, plink_stdout), fam_ct, bim_ct = calc_stats(raw_plink_prefix, output_dir, base_prefix)
update_counts('raw_data', fam_ct, bim_ct)

#
#  QC ON SAMPLES
#

(dups_removed_plink_prefix, plink_stdout), fam_ct, bim_ct = rm_dup_individuals(raw_plink_prefix, output_dir, base_prefix)
update_counts('rm_dup_indvid', fam_ct, bim_ct)

# note: --set-hh-missing is run along with checking for sex
(cleaned_sex_plink_prefix, all_stdout), fam_ct, bim_ct = rm_discordant_sex(
    dups_removed_plink_prefix, output_dir, base_prefix)
update_counts('rm_dscrd_sex_individ', fam_ct, bim_ct)

(ind_miss_removed_plink_prefix, rm_miss_indiv_stdout), fam_ct, bim_ct = rm_high_sample_missing_rate(
    cleaned_sex_plink_prefix, output_dir, base_prefix)
update_counts('rm_miss_indvid', fam_ct, bim_ct)

(het_removed_plink_file, ind_snp_output_file, all_stdout), fam_ct, bim_ct = rm_high_het_rate(
    ind_miss_removed_plink_prefix, output_dir, base_prefix)
update_counts('rm_het_indvid', fam_ct, bim_ct)

(related_indiv_removed_plink_prefix, all_stdout), fam_ct, bim_ct = rm_relatedness(
    het_removed_plink_file, output_dir, ind_snp_output_file, base_prefix)
update_counts('rm_related_indvid', fam_ct, bim_ct)

# add phenotype code to plink files
(phenotyped_plink_prefix, make_pheno_stdout), fam_ct, bim_ct = add_pheno(
    related_indiv_removed_plink_prefix, output_dir, base_prefix, cases='female')
update_counts('add_pheno', fam_ct, bim_ct)

#
#  QC ON SNPS
#

(vars_miss_removed_plink_prefix, rm_miss_var_stdout), fam_ct, bim_ct = rm_high_variant_missing_rate(
    phenotyped_plink_prefix, output_dir, base_prefix)
update_counts('rm_miss_snps', fam_ct, bim_ct)

(no_dups_vars_plink_prefix, all_stdout), fam_ct, bim_ct = rm_dup_var_id(
    vars_miss_removed_plink_prefix, output_dir, base_prefix)
update_counts('rm_dup_snps', fam_ct, bim_ct)

(vars_miss_removed_plink_prefix, rm_miss_var_stdout), fam_ct, bim_ct = test_missing(
    no_dups_vars_plink_prefix, output_dir, base_prefix)
update_counts('rm_test_miss_snps', fam_ct, bim_ct)


# =============  CALC FINAL STATS =============
(frq_and_miss_plink_prefix, plink_stdout), fam_ct, bim_ct = calc_stats(raw_plink_prefix, output_dir, base_prefix, prefix='inter_final_stats')
update_counts('final_stats', fam_ct, bim_ct)

#
#   CONVERT COUNTS TO DF
#

ind_ct_df = pd.DataFrame(TRACK_INDVID_DICT, index=[base_prefix])
snp_ct_df = pd.DataFrame(TRACK_SNP_DICT, index=[base_prefix])

ind_ct_df = ind_ct_df.reset_index().rename(columns={'index':'batch'})
snp_ct_df = snp_ct_df.reset_index().rename(columns={'index':'batch'})

ind_count_file = os.path.join(output_dir, 'counts_individuals_{}.tsv'.format(base_prefix))
snps_count_file = os.path.join(output_dir, 'counts_snps_{}.tsv'.format(base_prefix))
ind_ct_df.to_csv(ind_count_file, sep="\t", header=True, index=False)
snp_ct_df.to_csv(snps_count_file, sep="\t", header=True, index=False)


#
#   CLEAN UP FILES
#

log_files_path = os.path.join(output_dir, 'log')
intermediate_files_path = os.path.join(output_dir, 'intermediate')
temp_files_path = os.path.join(output_dir, 'temp')

safe_mkdir(log_files_path)
safe_mkdir(intermediate_files_path)
safe_mkdir(temp_files_path)

for file in glob.glob(output_dir+"/*.log"):
    os.rename(file, os.path.join(log_files_path, os.path.basename(file)))

for file in glob.glob(output_dir+"/inter_*"):
    os.rename(file, os.path.join(intermediate_files_path, os.path.basename(file)))

for file in glob.glob(output_dir+"/temp_*"):
    os.rename(file, os.path.join(temp_files_path, os.path.basename(file)))


print("\n\nFINISHED! processing {}. Took {:.3f} minutes.".format(plink_prefix, (time.time()-start)/60))
