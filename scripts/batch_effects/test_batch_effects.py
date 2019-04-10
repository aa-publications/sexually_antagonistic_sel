#!/bin/python
# Test for batch effects in MEGA data by testing for differences in genotype count between target batch and all other batches. 
#
#
#       INPUTS:
#           * see argparse help  
#
#
#       OUTPUTS: 
#           * 'batch_eff_{target_batch_plink_prefix}.tsv : snp and pvalue from fisher's exact test 
#           * 'rm_snps_batch_eff_{target_batch_plink_prefix}.txt: snps to remove that pass genome wide significance for batch effects
#
# 
#        DEPENDENCIES: 
#           * the list of batches plink prefix is hard-coded
#           * frqx_file_path: the path to look for the *.frqx file is HARD CODED! 
#
#
# Abin Abraham
# created on: 2019-04-10 14:08:38

import os
import sys
import time 
import argparse
import logging 
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')
sstart = time.time()

from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

# argparse
parser = argparse.ArgumentParser(description='Example with nonoptional arguments')
parser.add_argument('frq_dir', action='store', type=str, help="full path to directory with .frqx plink output files")
parser.add_argument('shared_snp_file', action='store', type=str, help="full path to txt file with shared snps across all batches")
parser.add_argument('target_batch_prefix', action='store', type=str, help="target batch plink prefix")
parser.add_argument('output_dir', action='store', type=str, help="full path to output dir")

results = parser.parse_args()
 
FRQ_DIR = results.frq_dir
SHARED_SNPS_FILE = results.shared_snp_file
TARGET_BATCH = results.target_batch_prefix
OUTPUT_DIR = results.output_dir


# logger 
logging_file = os.path.join(OUTPUT_DIR, 'rm_batch_eff_{}.log'.format(TARGET_BATCH))
logging.basicConfig(level=logging.INFO, filename=logging_file, filemode='w')
logger = logging.getLogger()
logger.addHandler(logging.FileHandler(logging_file, 'w'))
print = logger.info

# -----------
# FUNCTION
# -----------

def load_shared_snp_set(shared_snps_file):

    shared_snps = pd.read_csv(shared_snps_file, header=None)
    shared_snps_array = shared_snps.iloc[:, 0].values
    shared_snps_set = set(shared_snps_array)
    assert len(shared_snps_array) == len(shared_snps_set), "There are duplicate snpsIDs in the shared snp set."

    return shared_snps_set


# -----------
# MAIN
# -----------

PVAL_THRESHOLD = 5*(10**-8) # pval threshold to determien snps w/ batch effects 

# ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! #
#  FOR TESTING/DEV
# FRQ_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
# SHARED_SNPS_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_shared_snps.txt"
# TARGET_BATCH = "MEGA_ex_Array_Batch13_Cox_17_preQC_GRID"
# OUTPUT_DIR ="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
# ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! #

# all batches plink prefix
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

frqx_file_path = "{}_qc/temp/temp_final_stats_{}.frqx"
FRQ_FILE = frqx_file_path.format(TARGET_BATCH, TARGET_BATCH)
batches.remove(TARGET_BATCH)

# OUTPUT PATHS
snp_pval_output = os.path.join(OUTPUT_DIR, 'batch_eff_{}.tsv'.format(TARGET_BATCH))
snps_to_rm_output = os.path.join(OUTPUT_DIR, 'batch_eff_snps_{}.txt'.format(TARGET_BATCH))

print("Running {} on {}....".format(sys.argv[0], DATE))
print("\t target batch file:\n\t\t{}".format(os.path.join(FRQ_DIR, FRQ_FILE)))
print("\t outputs written to:\n\t\t{}".format(OUTPUT_DIR))
print("\t p-value threshold:\n\t\t{}".format(PVAL_THRESHOLD))
print("\t using shared snps from:\n\t\t{}".format(SHARED_SNPS_FILE))


#
#   LOAD TARGET BATCH
#

shared_snps_set = load_shared_snp_set(SHARED_SNPS_FILE)

full_target_df = pd.read_csv(os.path.join(FRQ_DIR, FRQ_FILE), sep="\t")
target_df = full_target_df[full_target_df.SNP.isin(shared_snps_set)].copy()
target_df.set_index('SNP', inplace=True)
target_df = target_df.loc[:, ['C(HOM A1)',  'C(HET)',  'C(HOM A2)', 'C(MISSING)']].copy()


#
#   SUM COUNTS FOR GENOTYEPS OVER ALL BATCHES
#

# concatenate each batch ...
prealloc_counts_tensor = np.ones((1, len(shared_snps_set), 4))*-1
for this_batch in batches:

    print("Collapsing batch: {}".format(this_batch))
    
    full_batch_df = pd.read_csv(os.path.join(
        FRQ_DIR, frqx_file_path.format(this_batch, this_batch)), sep="\t")
    
    batch_df = full_batch_df[full_batch_df.SNP.isin(shared_snps_set)].copy()

    # convert to numpy array; ensure the order of snps match the target_df
    batch_df.set_index('SNP', inplace=True)
    reindex_batch_df = batch_df.reindex(index=target_df.index)

    this_array = reindex_batch_df.loc[:, ['C(HOM A1)',  'C(HET)',  'C(HOM A2)',
                                          'C(MISSING)']].values.reshape((1, len(shared_snps_set), 4))

    prealloc_counts_tensor = np.concatenate((prealloc_counts_tensor, this_array), axis=0)

# drop the preallocated frame
counts_tensor = prealloc_counts_tensor[1:, :, :]

# sum over all batches
summed_counts = np.sum(counts_tensor, axis=0)
all_batch_df = pd.DataFrame(summed_counts, index=target_df.index, columns=[
                            'C(HOM A1)',  'C(HET)',  'C(HOM A2)', 'C(MISSING)'])

#
#   FISHER EXACT TEST PER SNP
#

stats = importr('stats')
pval_df = pd.DataFrame({'snp': target_df.index,
                        'pval': np.nan*np.ones(target_df.shape[0])})

for counter, snp in enumerate(target_df.index):
    print("Testing {:,} out of {:,}".format(counter, target_df.shape[0])) if (counter % 100000 == 0) else None

    table = np.array([target_df.loc[snp, :].values, all_batch_df.loc[snp, :].values])
    res = stats.fisher_test(table,workspace=2e8)

    pval_df.loc[pval_df['snp'] == snp, 'pval'] = res[0][0]


pval_df['bonferroni_thresh'] = PVAL_THRESHOLD
pval_df['pass_mult_test'] = pval_df.pval < PVAL_THRESHOLD

snps_to_rm_df = pval_df.loc[pval_df['pass_mult_test'] == True, 'snp'].copy()

print("\n*** {:,} out of {:,} SNPs show significant evidence for batch effects.\n".format(snps_to_rm_df.shape[0],pval_df.shape[0] ))
# write 
pval_df.to_csv(snp_pval_output, sep="\t", index=False, header=True)
snps_to_rm_df.to_csv(snps_to_rm_output, header=False, index=False)
print("Output written to: {}".format(OUTPUT_DIR))
print("Done! Took {:.2f} minutes".format( (time.time() - sstart)/60 ))