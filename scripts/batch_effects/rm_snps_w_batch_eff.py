#!/bin/python
# This script will aggregate snp with batch effects and then remove them each batch file.
#
#
#
# Abin Abraham
# created on: 2019-04-12 08:55:13



import os
import sys
import glob
import pandas as pd
from datetime import datetime

sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')
from func_run_shell_cmd import run_shell_cmd
from func_track import get_num_lines

DATE = datetime.now().strftime('%Y-%m-%d')


# INPUT FILES
SNPS_TO_RM_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/batch_effects"
PLINK_FILES_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch/{}_qc"
plink_prefix = "final_qc_{}"

# OUTPUT FILES
OUTPUT_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/batch_effects"
BATCH_EF_SNP_TO_RM_FILE = os.path.join(OUTPUT_DIR, 'uniq_batch_eff_snps_to_rm_file_{}.txt'.format(DATE))
OUTPUT_BEFFECTS_RM_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/batch_effects/plink_batch_effect_snps_removed"

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
# MAIN
# -----------

# load all the batch effect snps from each batch, dedups, and then write to plink compatible text file 
all_df = pd.DataFrame()
for batch_file in glob.glob(SNPS_TO_RM_DIR + "/batch_eff_snps_*.txt"):

    df = pd.read_csv(batch_file, sep="\t", header=None, names=['snps'])
    all_df = all_df.append(df)


uniq_snps_to_rm = df.snps.unique()

with open(BATCH_EF_SNP_TO_RM_FILE, 'w') as fout:
    for line in uniq_snps_to_rm:
        fout.writelines(line+"\n")


# run plink to remove snps 
if os.stat(BATCH_EF_SNP_TO_RM_FILE).st_size != 0:
    for batch in batches:
        rm_dups_cmd = "plink --bfile {} --exclude {} --make-bed --out {}".format(
            os.path.join(PLINK_FILES_DIR.format(batch), plink_prefix.format(batch)),
            BATCH_EF_SNP_TO_RM_FILE,
            os.path.join(OUTPUT_BEFFECTS_RM_FILE, "no_batch_snps_"+batch))

        plink_stdout = run_shell_cmd(rm_dups_cmd)

        print("Done removing batch snps for {}".format(batch))

print("Plink output written to {}".format(OUTPUT_BEFFECTS_RM_FILE))