#!/bin/python
# This script will create a list of FID to remove from each study/batch due to duplicates. 
#
#
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
study_prefix_file = "/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/plink_qc_v2/bfile_names.txt"
data_dir="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/plink_qc_v2/"

# data_dir="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/postQC"
output_dir="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/plink_qc_v2/"


# -----------
# MAIN
# ----------- 

# create directory to store id's to remove 
fid_to_remove_dir = os.path.join(output_dir, 'dups_fids_per_study')

if not os.path.isdir(fid_to_remove_dir): 
    os.mkdir(fid_to_remove_dir)

# create lsit of file prefixies for studies 
study_prefix_list = []
with open(study_prefix_file, 'r') as fread: 
    for line in fread: 
        study_prefix_list.append(line.splitlines()[0]+"_clean")
        # study_prefix_list.append(line.splitlines()[0])

# concat all FID with study to one df
fid_df = pd.DataFrame()

for study in study_prefix_list: 
    print("Appending {}.".format(study))

    df = pd.read_csv(os.path.join(data_dir,study+".fam"), sep="\s+", usecols=[1],names=['FID'])
    df['study'] = study

    fid_df = fid_df.append(df)

print("{:,} FIDs across all {} studies/batches".format(fid_df.shape[0], fid_df.study.nunique()))

# duplicate FIDs
dups_id = fid_df[fid_df.duplicated(subset='FID')].copy()
dups_id['IID'] = dups_id['FID']
print("{:,} FIDs are duplicates.".format(dups_id.FID.nunique()))

# write duplicate FID to file (in plink-ready format)
for study in dups_id.study.unique(): 
    fid_to_remove_file = os.path.join(fid_to_remove_dir, study+"_fid_dups.tsv")
    dups_id.loc[dups_id['study']==study, ['FID', 'IID']].to_csv(fid_to_remove_file, sep="\t", index=False, header=False)
    
    print("Wrote FIDs to remove for {}".format(study))

print("Done. Took {:.2f} minutes.".format(time.time()-start))    