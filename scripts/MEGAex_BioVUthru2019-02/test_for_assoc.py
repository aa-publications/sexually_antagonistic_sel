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

MERGED_PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/maf_filtered/final_maf_filtered__final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3"
COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/covariates/covar_white_YOB_5PC__2019-07-05.tsv"

LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/{}_logistic_w_sex_correction".format(DATE)


# -----------
# MAIN
# -----------
#

# run using plink2 w/ glm
# sex is not included for ANY of the variants with ' no-x-sex '
# covar-variance-standardize mean 0, variance 1
# lr_cmd = "plink2 --bfile {} --covar {} --glm no-x-sex firth-fallback intercept --covar-variance-standardize YOB --vif 1500 --out {}".format(MERGED_PLINK_FILE, COVAR_FILE, LOGISTIC_RESUTLS_OUTPUT)

lr_cmd = "plink2 --bfile {} --covar {} --glm firth-fallback intercept --covar-variance-standardize YOB --vif 1500 --out {}".format(MERGED_PLINK_FILE, COVAR_FILE, LOGISTIC_RESUTLS_OUTPUT)

print("Running:\n{}".format(lr_cmd))
shell_stdout_model = run_shell_cmd(lr_cmd)

print("Finished: {} mins.".format((time.time() - start)/60)
