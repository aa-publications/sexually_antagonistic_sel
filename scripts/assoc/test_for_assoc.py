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

MERGED_PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/merged_MEGA_2019-04-12"
FISHER_OUTPUT_PLINK_PREFIX="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_13_fisher_exact/merged_MEGA"
MODEL_OUTPUT_PLINK_PREFIX="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_13_geno_assoc/merged_MEGA_assoc"
LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_04_14_logistic_assoc/merged_MEGA_w_covar"


COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/covariates/covar_CEU_merged_mega_age_5PC_2019-04-14.tsv"


# -----------
# MAIN
# ----------- 

# run plink assoc 
assoc_cmd = "plink --bfile {} --assoc fisher counts --out {}".format(MERGED_PLINK_FILE, FISHER_OUTPUT_PLINK_PREFIX)
# shell_stdout = run_shell_cmd(assoc_cmd)

# run plink with --model  
assoc_cmd = "plink --bfile {} --model fisher --out {}".format(MERGED_PLINK_FILE, MODEL_OUTPUT_PLINK_PREFIX)
# shell_stdout_model = run_shell_cmd(assoc_cmd)

# logistic regression
lr_cmd = "plink --bfile {} --covar {} --covar-name YOB, ARRAY, PC1, PC2, PC3, PC4, PC5 --adjust gc --logistic no-x-sex --ci 0.95 --out {}".format(MERGED_PLINK_FILE, COVAR_FILE, LOGISTIC_RESUTLS_OUTPUT)

print("Running:\n{}".format(lr_cmd))
shell_stdout_model = run_shell_cmd(lr_cmd)
print(lr_cmd.splitlines())