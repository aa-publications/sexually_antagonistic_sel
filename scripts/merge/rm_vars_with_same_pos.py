#!/bin/python
# This scripts works on an already merged (combined batches) plink MEGA file. 
#           * identifty variants with the same CHR and POS 
#           * pick one to remove and write 
# 
# Abin Abraham
# created on: 2019-04-13 12:42:42



import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')

sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')
from func_run_shell_cmd import run_shell_cmd


start = time.time()



# =============  PATHS =============

# INPUTS:
MERGED_PLINK_PREFIX = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/merged_MEGA_2019-04-12"
plink_file_prefix = os.path.basename(MERGED_PLINK_PREFIX)

OUTPUT_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches"
LIST_DUPS_OUTPUT=os.path.join(OUTPUT_DIR, 'dup_vars_merged_MEGA_2019-04-12')
IDS_TO_RM_FILE=os.path.join(OUTPUT_DIR, 'var_w_same_chr_pos_to_rm.txt')

# =============  MAIN =============

shell_cmd = "plink --bfile ${bfile_path} --list-duplicate-vars --out id_dup_vars".format(MERGED_PLINK_PREFIX, LIST_DUPS_OUTPUT)
plink_stdout = run_shell_cmd(shell_cmd)

# pick all but one snp to remove 
dups_df = pd.read_csv(LIST_DUPS_OUTPUT+".dupvar", sep="\t")
dups_df['ID_to_rm']  = dups_df.IDS.apply(lambda x: ' '.join(x.split()[1:]))
dups_to_rm_df = pd.DataFrame(dups_df['ID_to_rm'].str.split(' ').tolist()).stack()

df_to_write = dups_to_rm_df.to_frame(name='snps_to_rm')
df_to_write.to_csv(IDS_TO_RM_FILE, index=False, header=False )
