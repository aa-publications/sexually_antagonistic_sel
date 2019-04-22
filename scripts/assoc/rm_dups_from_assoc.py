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
ASSOC_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_18_glm_no_batch4/no_batch4_glm.PHENO1.glm.logistic.hybrid"
CLEAN_ASSOC_OUTPUT=os.path.join(os.path.split(ASSOC_FILE)[0], "dedup_clean_"+os.path.split(ASSOC_FILE)[1])
SNPS_TO_REMOVE_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/merge_batch_without_batch4/var_w_same_chr_pos_to_rm.txt"

# -----------
# MAIN
# ----------- 

assoc_df = pd.read_csv(ASSOC_FILE, sep="\s+")
assoc_df.rename(columns={'ID':'SNP', '#CHROM':'CHR'}, inplace=True)
var_to_rm = pd.read_csv(SNPS_TO_REMOVE_FILE,header=None)
snp_set = set(var_to_rm.values.flatten())

filtered_df = assoc_df[~assoc_df.SNP.isin(snp_set)].copy()

# note: this is for plink2 output for --glm
final_df = filtered_df.loc[:, ['CHR','POS','SNP','REF','ALT', 'A1', 'OR','SE', 'P' ]].copy()
no_na_final_df = final_df[~final_df.P.isnull()].copy()


print("{:,} SNPs removed due to NA for p-values". format(final_df.P.isnull().sum()))
print("Removed (n={:,}) all but one SNP with duplicated CHR and POS.".format(assoc_df.shape[0] - final_df.shape[0]))
print("Final number of SNPs: {:,}".format(final_df.shape[0]))

final_df.to_csv(CLEAN_ASSOC_OUTPUT, sep="\t", index=False)
print("output written to:\n{}".format(CLEAN_ASSOC_OUTPUT))