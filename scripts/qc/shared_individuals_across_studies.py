#!/bin/python
# This script will create a list of FID to remove from each study/batch due to duplicates. 
#
#   OUTPUTS: a file per batch with IDs to remove 
#
# Abin Abraham
# created on: 2018-10-05 13:45:20


import os
import sys
import numpy as np
import pandas as pd
import time 

start = time.time()

# ******************************************
#          USER MODIFY THESE FILES 
# ******************************************
data_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/batch_effects/plink_batch_effect_snps_removed"
output_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data"
plink_prefix = os.path.join(data_dir,"no_batch_snps_{}.fam")

# -----------
# MAIN
# ----------- 

# create directory to store id's to remove 
fid_to_remove_dir = os.path.join(output_dir, 'dups_fids_per_study')
if not os.path.isdir(fid_to_remove_dir): 
    os.mkdir(fid_to_remove_dir)




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



# concat all FID with study to one df
fid_df = pd.DataFrame()
for batch in batches: 
    print("Appending {}.".format(batch))

    df = pd.read_csv(plink_prefix.format(batch), sep="\s+", usecols=[1],names=['FID'])
    df['batch'] = batch

    fid_df = fid_df.append(df)
print("{:,} FIDs across all {} studies/batches".format(fid_df.shape[0], fid_df.batch.nunique()))

# identify duplicate FID 
# then keep only the first instance, and then record the batch also.... 
dups_id = fid_df[fid_df.duplicated(subset='FID')].copy()
dups_id['IID'] = dups_id['FID']
print("{:,} FIDs are duplicates.".format(dups_id.FID.nunique()))

# write duplicate FID to file (in plink-ready format)
# for each batch, write a file with FIDS to remove 
for batch in dups_id.batch.unique(): 
    fid_to_remove_file = os.path.join(fid_to_remove_dir, batch+"_fid_dups.tsv")
    dups_id.loc[dups_id['batch']==batch, ['FID', 'IID']].to_csv(fid_to_remove_file, sep="\t", index=False, header=False)
    
    print("Wrote FIDs to remove for {}".format(batch))


print("Done. Took {:.2f} minutes.".format(time.time()-start))    