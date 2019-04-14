#!/bin/python
# This script will will output a tsv with individual (GRID) and batch that it was sequenced in. (no header)
#
#
#
# Abin Abraham
# created on: 2019-04-14 10:04:38

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
 
DATE = datetime.now().strftime('%Y-%m-%d')


PLINK_PREFIX="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/merge_batches/temp/temp_fid_dedups_{}"
OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/covariates/create_individ_array_table_{}.py".format(DATE)

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


# MAIN 


# =============  MAIN =============

final_df = pd.DataFrame(columns=['GRID','batch'])

for batch in batches: 
    df = pd.read_csv(PLINK_PREFIX.format(batch)+".fam", sep="\s+", header=None, names=['GRID'], usecols=[0] ) 
    df['batch'] = batch

    final_df = final_df.append(df)


final_df.to_csv(OUTPUT_FILE, sep="\t", index=False, header=False)