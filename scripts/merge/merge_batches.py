#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2019-04-12 09:53:27



import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')

sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')
from func_rm_fids import rm_fids
from func_test_missing import test_missing
from func_run_shell_cmd import run_shell_cmd
from func_maf_filter import maf_filter
from func_hwe import hwe
from func_pca import pca

start = time.time()
# =============  PATHS =============

# INPUTS:
PLINK_TO_MERGE_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/batch_effects/plink_batch_effect_snps_removed"
plink_file_prefix = os.path.join(PLINK_TO_MERGE_DIR, "no_batch_snps_{}")

DUP_FID_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/dups_fids_per_study"
dup_fid_file_prefix = os.path.join(DUP_FID_DIR, "fid_dups_{}.tsv")

# OUTPUTS
OUTPUT_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches"
PLINK_FILES_TO_MERGE = os.path.join(OUTPUT_DIR, 'merge_list_plink_prefixes.txt')
MERGED_PLINK_FILE = os.path.join(OUTPUT_DIR, 'merged_MEGA_{}'.format(DATE))

# batches
batches = ["MEGA_ex_Array_Batch10_Cox_14_GenderDysphoria9_preQC_GRID",
           "MEGA_ex_Array_Batch11_Cox_12_Batch1_preQC_GRID",
           "MEGA_ex_Array_Batch12_Cox_12_Batch2_preQC_GRID",
           "MEGA_ex_Array_Batch13_Cox_17_preQC_GRID",
           "MEGA_ex_Array_Batch1_Cox_01_06_preQC_GRID",
           "MEGA_ex_Array_Batch2_Cox_06_08_preQC_GRID",
           "MEGA_ex_Array_Batch3_Cox_05_preQC_GRID",
           "MEGA_ex_Array_Batch4_Cox07_Cox07Shadow_preQC_GRID",
           "MEGA_ex_Array_Batch5_Cox_09_preQC_GRID",
           "MEGA_ex_Array_Batch6_Cox_03_15_preQC_GRID",
           "MEGA_ex_Array_Batch7_Cox_11_preQC_GRID",
           "MEGA_ex_Array_Batch8_Cox_13_02_preQC_GRID",
           "MEGA_ex_Array_Batch9_Cox_04_10_preQC_GRID"]

# -----------
# FUNCTION
# -----------


def safe_mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


# -----------
# MAIN
# -----------

#
#   REMOVE DUPLICATE INDVIDUALS
# # 
# for batch in batches:

#     input_prefix = plink_file_prefix.format(batch)
#     dups_ids_file = dup_fid_file_prefix.format(batch)
#     pprefix = "temp_fid_dedups"
#     (no_dups_plink_prefix, plink_stdout), fam_ct, bim_ct = rm_fids(
#         dups_ids_file, input_prefix, OUTPUT_DIR, batch, prefix=pprefix)


# #
# #   MERGE IN PLINK
# #

# # write a txt file with plink prefix fiels to merge
# with open(PLINK_FILES_TO_MERGE, 'w') as fout:
#     for batch in batches:
#         fout.write(os.path.join(OUTPUT_DIR, pprefix+"_"+batch+"\n"))

# merge_cmd = "plink --merge-list {} --out {}".format(PLINK_FILES_TO_MERGE, MERGED_PLINK_FILE)
# shell_stdout = run_shell_cmd(merge_cmd)

_, merged_plink_prefix = os.path.split(MERGED_PLINK_FILE)
# #
# #   DOWNSTREAM ANALYSIS
# #

# # -test-missing
# (vars_miss_removed_plink_prefix, rm_miss_var_stdout), fam_ct, bim_ct = test_missing(
#     MERGED_PLINK_FILE, OUTPUT_DIR, merged_plink_prefix, prefix='temp_merged')

# # -maf filter
# (maf_filtered_plink_file, plink_stdout), fam_ct, bim_ct = maf_filter(
#     vars_miss_removed_plink_prefix, OUTPUT_DIR, merged_plink_prefix, prefix='temp_maf_filtered_')


# # HWE SNPs all
# hwe_outputs, plink_stdout = hwe(maf_filtered_plink_file, OUTPUT_DIR, merged_plink_prefix, prefix='inter_hwe_')
# # HWE SNPs males
# hwe_outputs, plink_stdout = hwe(maf_filtered_plink_file, OUTPUT_DIR,
#                                      merged_plink_prefix, prefix='inter_hwe_males', plink_opts="--filter-males")
# # HWE SNPs females
# hwe_outputs, plink_stdout = hwe(maf_filtered_plink_file, OUTPUT_DIR,
#                                      merged_plink_prefix, prefix='inter_hwe_females', plink_opts="--filter-females")

# calc PCA across all individuals
maf_filtered_plink_file='/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/temp_maf_filtered__merged_MEGA_2019-04-12'
pca_outputs, plink_stdout = pca(maf_filtered_plink_file, OUTPUT_DIR, merged_plink_prefix, prefix='inter_pca_')
print(plink_stdout.splitlines())
sys.stdout.flush()

#
#   CLEAN UP FILES
#

log_files_path = os.path.join(OUTPUT_DIR, 'log')
temp_files_path = os.path.join(OUTPUT_DIR, 'temp')
inter_files_path = os.path.join(OUTPUT_DIR, 'intermediate')

safe_mkdir(log_files_path)
safe_mkdir(temp_files_path)
safe_mkdir(inter_files_path)

# for file in glob.glob(OUTPUT_DIR+"/*.log"):
#     os.rename(file, os.path.join(log_files_path, os.path.basename(file)))

# for file in glob.glob(OUTPUT_DIR+"/temp_*"):
#     os.rename(file, os.path.join(temp_files_path, os.path.basename(file)))

# for file in glob.glob(OUTPUT_DIR+"/inter_*"):
#     os.rename(file, os.path.join(inter_files_path, os.path.basename(file)))


print("\n\nFINISHED! processing {}. Took {:.3f} minutes.".format(merged_plink_prefix, (time.time()-start)/60))
