#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2019-04-12 20:26:10



import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
import time 

sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')
from func_run_shell_cmd import run_shell_cmd

DATE = datetime.now().strftime('%Y-%m-%d')
start = time.time()


# =============  PATHS =============

MERGED_PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/merge_batch_without_batch4/final_maf_filtered__merged_MEGA_no_batch4_2019-04-15"
FISHER_OUTPUT_PLINK_PREFIX="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_18_glm_no_batch4/no_batch4_fisher"

LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_18_glm_no_batch4/no_batch4_glm"

COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/covariates/no_batch4_merged_pca/covar_CEU_no_batch4_merged_mega_w_age_5PC_batch_2019-04-18.tsv"


# -----------
# MAIN
# ----------- 

# run plink assoc 
assoc_cmd = "plink --bfile {} --assoc fisher counts --out {}".format(MERGED_PLINK_FILE, FISHER_OUTPUT_PLINK_PREFIX)
# shell_stdout = run_shell_cmd(assoc_cmd)

# run plink with --model  
assoc_cmd = "plink --bfile {} --model fisher --out {}".format(MERGED_PLINK_FILE, MODEL_OUTPUT_PLINK_PREFIX)
# shell_stdout_model = run_shell_cmd(assoc_cmd)

# run using plink2 w/ glm
lr_cmd = "plink2 --bfile {} --covar {} --glm no-x-sex  hide-covar firth-fallback  --vif 1500 --out {}".format(MERGED_PLINK_FILE, COVAR_FILE, LOGISTIC_RESUTLS_OUTPUT)

print("Running:\n{}".format(lr_cmd))
shell_stdout_model = run_shell_cmd(lr_cmd)





# logistic regression
# lr_cmd = "plink --bfile {} --covar {} --logistic no-x-sex --out {}".format(MERGED_PLINK_FILE, COVAR_FILE, LOGISTIC_RESUTLS_OUTPUT)
