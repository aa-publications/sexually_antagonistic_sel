#!/bin/python
# This script will look through all studies/batches and create a list of snps that are shared across all studies.
# 
#
#   USER MODIFIED PATHS:
#       - study_prefix_file: full path to txt file with one line per plink prefix for each batch 
#       - data_dir: full path to dir with plink bfiles 
#       - output_snps_file: full path to output .txt file with common snps
#
#
#
#
#
#
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
study_prefix_file = "/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/raw_preQC/cox_bfile_prefix.txt"
data_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
output_snps_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_shared_snps.txt"


# create list of prefixies for studies 
study_prefix_list = []
with open(study_prefix_file, 'r') as fread: 
    for line in fread: 
        study_prefix_list.append(line.splitlines()[0])


# create a set of snps (rsID) that are shared across all studies/batches
snp_set = set()
for index, study in enumerate(study_prefix_list):
    
    print("Getting snps from {}".format(study))
    df = pd.read_csv(os.path.join(data_dir,study+"_qc", "final_qc_"+study+".bim"), sep="\s+", usecols=[1], names=['rsID'])
    
    if index == 0: 
        snp_set.update(df.rsID.values)
    else: 
        this_set = set(df.rsID.values)

        snp_set.intersection_update(this_set)
        

final_snps = list(snp_set)


with open(output_snps_file,'w') as fout: 

    fout.writelines("\n".join(final_snps))


print("Done! Check {}".format(output_snps_file))
print("It took {:.2f} seconds.".format(time.time()-start))