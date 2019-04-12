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


from functools import partial
from multiprocessing import Pool, cpu_count


from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

DATE = datetime.now().strftime('%Y-%m-%d')
sstart = time.time()


# ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! #
#  FOR TESTING/DEV
#FRQ_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
#SHARED_SNPS_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_shared_snps.txt"
#TARGET_BATCH = "MEGA_ex_Array_Batch13_Cox_17_preQC_GRID"
#OUTPUT_DIR ="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
# ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! # ! #



# # argparse
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
fmt='%(asctime)s:%(message)s'
logging.basicConfig(level=logging.INFO, format=fmt, filename='{}_{}.log'.format(TARGET_BATCH, DATE), filemode='w')
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)


# -----------
# FUNCTION
# -----------

def load_shared_snp_set(shared_snps_file):

    shared_snps = pd.read_csv(shared_snps_file, header=None)
    shared_snps_array = shared_snps.iloc[:, 0].values
    shared_snps_set = set(shared_snps_array)
    assert len(shared_snps_array) == len(shared_snps_set), "There are duplicate snpsIDs in the shared snp set."

    return shared_snps_set


def calc_fishers(this_snp):

    table = np.array([target_df.loc[this_snp, :].values, all_batch_df.loc[this_snp, :].values])
    res = stats.fisher_test(table,workspace=2e8)
    
    pval = res[0][0]
    # pval_df.loc[pval_df['snp'] == snp, 'pval'] = res[0][0]

    return (this_snp, pval)

# -----------
# MAIN
# -----------

PVAL_THRESHOLD = 5*(10**-8) # pval threshold to determien snps w/ batch effects 

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

frqx_file_path = "{}_qc/final_qc_{}.frqx"
FRQ_FILE = frqx_file_path.format(TARGET_BATCH, TARGET_BATCH)
batches.remove(TARGET_BATCH)

# OUTPUT PATHS
snp_pval_output = os.path.join(OUTPUT_DIR, 'batch_eff_{}.tsv'.format(TARGET_BATCH))
snps_to_rm_output = os.path.join(OUTPUT_DIR, 'batch_eff_snps_{}.txt'.format(TARGET_BATCH))

logging.info("Running {} on {}....".format(sys.argv[0], DATE))
logging.info("target batch file:\n{}".format(os.path.join(FRQ_DIR, FRQ_FILE)))
logging.info("outputs written to:\n{}".format(OUTPUT_DIR))
logging.info("p-value threshold:\n{}".format(PVAL_THRESHOLD))
logging.info("using shared snps from:\n{}".format(SHARED_SNPS_FILE))


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
logging.info("concatenating all non-target batch genotype counts...")
prealloc_counts_tensor = np.ones((1, len(shared_snps_set), 4))*-1
for this_batch in batches:

    logging.info("Collapsing batch: {}".format(this_batch))
    
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
all_batch_df = pd.DataFrame(summed_counts, index=target_df.index, columns=['C(HOM A1)',  'C(HET)',  'C(HOM A2)', 'C(MISSING)'])


#
#   FISHER EXACT TEST PER SNP
#

#benchmark 1000 snps in 2.16 mins 

bts = time.time()
logging.info("Calculating fisher's test with {} core".format(os.cpu_count()))

stats = importr('stats')

partial_calc_fish = partial(calc_fishers)
pool= Pool()
results = pool.map(partial_calc_fish, shared_snps_set, 100)

pool.close()
pool.join()

logging.info("Done in {:.2f} mins for fisher's test.".format((time.time() - bts)/60))


pval_df = pd.DataFrame(results, columns=['snp','pval'])
pval_df['bonferroni_thresh'] = PVAL_THRESHOLD
pval_df['pass_mult_test'] = pval_df.pval < PVAL_THRESHOLD


snps_to_rm_df = pval_df.loc[pval_df['pass_mult_test'] == True, 'snp'].copy()

logging.info("*** {:,} out of {:,} SNPs show significant\
         evidence for batch effects.\n".format(snps_to_rm_df.shape[0],pval_df.shape[0]))

# write 
pval_df.to_csv(snp_pval_output, sep="\t", index=False, header=True)
snps_to_rm_df.to_csv(snps_to_rm_output, header=False, index=False)
logging.info("Output written to: {}".format(OUTPUT_DIR))
logging.info("Done! Took {:.2f} minutes".format( (time.time() - sstart)/60 ))
