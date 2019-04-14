#!/bin/python
# 
# 
#           This script will remove SNPs with same CHR:POS from '.assoc.fisher' plink output. (Will retain one of the duplicates randomly)
#           Also REMOVE C_A and C_U columns from plink output so its is ready for FUMA.
#
#
# Abin Abraham
# created on: 2019-04-13 13:17:28


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
 
DATE = datetime.now().strftime('%Y-%m-%d')




# -----------
# PATHS
# ----------- 
ASSOC_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_19_fisher_exact/merged_MEGA.assoc.fisher"
CLEAN_ASSOC_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_19_fisher_exact/dedup_snps_merged_MEGA.assoc.fisher"
SNPS_TO_REMOVE_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/var_w_same_chr_pos_to_rm.txt"

# -----------
# MAIN
# ----------- 

assoc_df = pd.read_csv(ASSOC_FILE, sep="\s+")
var_to_rm = pd.read_csv(SNPS_TO_REMOVE_FILE,header=None)
snp_set = set(var_to_rm.values.flatten())

filtered_df = assoc_df[~assoc_df.SNP.isin(snp_set)].copy()
filtered_df.drop(['C_A','C_U'], axis=1, inplace=True)

print("Removed (n={:,}) all but one SNP with duplicated CHR and POS.".format(assoc_df.shape[0] - filtered_df.shape[0]))
print("Final number of SNPs: {:,}".format(filtered_df.shape[0]))

filtered_df.to_csv(CLEAN_ASSOC_OUTPUT, sep="\t", index=False)
print("output written to:\n{}".format(CLEAN_ASSOC_OUTPUT))