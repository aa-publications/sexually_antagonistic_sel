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

MERGED_PLINK_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/impute/final_imputed/merged_filt_imp"

LOGISTIC_RESUTLS_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_11_21_imputed_mega"

COVAR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/impute/final_imputed/covar_white_YOB_12PC__2019-07-21.tsv"


# -----------
# MAIN
# ----------- 

# logistic regression
lr_cmd = "plink --bfile {} --covar {} --logistic no-x-sex --out {}".format(MERGED_PLINK_FILE, COVAR_FILE, LOGISTIC_RESUTLS_OUTPUT)
plink --bfile $MERGED_PLINK_FILE --covar $COVAR_FILE --logistic no-x-sex --out $LOGISTIC_RESUTLS_OUTPUT