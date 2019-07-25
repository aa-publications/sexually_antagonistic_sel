#!/bin/python
# This script will perform QC on MEGA data "BioVUthru2019-02"
#
#
#
# Abin Abraham
# created on: 2019-04-12 09:53:27
# modified on: 2019-07-03 08:58:00


import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')

sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')
from func_test_missing import test_missing
from func_run_shell_cmd import run_shell_cmd

sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/merge')
from func_maf_filter import maf_filter
from func_hwe import hwe


start = time.time()


# =============  PATHS =============
# INPUTS:
PLINK_INPUT_PATH ="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3" # w/ plink prefix , no extension
input_prefix = os.path.basename(PLINK_INPUT_PATH)

# OUTPUTS
OUTPUT_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02"

print("Running {} on {}...\
        \n\tplink_input: {}\n\tdata_dir: {}\n\toutput_dir: {}\n\n".format(os.path.basename(__file__), DATE, input_prefix, os.path.dirname(PLINK_INPUT_PATH), OUTPUT_DIR))


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

# -maf filter
(maf_filtered_plink_file, plink_stdout), fam_ct, bim_ct = maf_filter(
        PLINK_INPUT_PATH, OUTPUT_DIR, input_prefix, prefix='final_maf_filtered_')

# -hwe
# hwe SNPs all
hwe_outputs, plink_stdout = hwe(maf_filtered_plink_file, OUTPUT_DIR, input_prefix, prefix='inter_hwe_')

# hwe SNPs males
hwe_outputs, plink_stdout = hwe(maf_filtered_plink_file, OUTPUT_DIR,
                                     input_prefix, prefix='inter_hwe_males', plink_opts="--filter-males")
# hwe SNPs females
hwe_outputs, plink_stdout = hwe(maf_filtered_plink_file, OUTPUT_DIR,
                                     input_prefix, prefix='inter_hwe_females', plink_opts="--filter-females")



#
#   CLEAN UP FILES
#

log_files_path = os.path.join(OUTPUT_DIR, 'log')
temp_files_path = os.path.join(OUTPUT_DIR, 'temp')
inter_files_path = os.path.join(OUTPUT_DIR, 'intermediate')

safe_mkdir(log_files_path)
safe_mkdir(temp_files_path)
safe_mkdir(inter_files_path)

for file in glob.glob(OUTPUT_DIR+"/*.log"):
    os.rename(file, os.path.join(log_files_path, os.path.basename(file)))

for file in glob.glob(OUTPUT_DIR+"/temp_*"):
    os.rename(file, os.path.join(temp_files_path, os.path.basename(file)))

for file in glob.glob(OUTPUT_DIR+"/inter_*"):
    os.rename(file, os.path.join(inter_files_path, os.path.basename(file)))


print("\n\nFINISHED! processing {}. Took {:.3f} minutes.".format(input_prefix, (time.time()-start)/60))
