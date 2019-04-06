#!/bin/python
# This script will look through all studies/batches and create a list of snps that are shared across all studies.
#
#
#
# Abin Abraham
# created on: 2018-10-05 12:37:33

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
final_snps_file="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/plink_qc_v2/shared_snps.txt"


# create lsit of prefixies for studies 
study_prefix_list = []
with open(study_prefix_file, 'r') as fread: 
    for line in fread: 
        study_prefix_list.append(line.splitlines()[0])


# create a set of snps (rsID) that are shared across all studies/batches
snp_set = set()
for index, study in enumerate(study_prefix_list):
    
    print("Getting snps from {}".format(study))
    df = pd.read_csv(os.path.join(data_dir, study+"_clean.bim"), sep="\s+", usecols=[1], names=['rsID'])
    
    if index == 0: 
        snp_set.update(df.rsID.values)
    else: 
        this_set = set(df.rsID.values)

        snp_set.intersection_update(this_set)
        

final_snps = list(snp_set)


with open(final_snps_file,'w') as fout: 

    fout.writelines("\n".join(final_snps))


print("Done! Check {}".format(final_snps_file))
print("It took {:.2f} seconds.".format(time.time()-start))