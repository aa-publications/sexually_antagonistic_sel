#!/bin/python
# This script will take plink assoc output and create a list of sequence to pull using twoBItToFa
#
#
#
# Abin Abraham
# created on: 2019-07-24 10:36:57

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')


OUTPUT_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'seqList_for_twoBitToFa_{}.txt'.format(DATE))
OUTPUT_FILE_KEY = os.path.join(OUTPUT_DIR, 'key_for_seqList_for_twoBitToFa_{}.tsv'.format(DATE))


pval_thresh = 0.5*10**-5
SEPERATOR=","

###
#   INPUT FILES
###

# plink output from assoc
# must include the following headers:
#       >  BP
#       >  CHR
#       >  P
#       >  TEST

assoc_file ="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/20190722_gwas_pca12_age.csv"




###
#   MAIN
###

# select SNPS to create seqList
as_df = pd.read_csv(assoc_file, sep=SEPERATOR)
sig_snps_df = as_df.loc[ (as_df['P'] < pval_thresh) & (as_df['P'] !=0) & (as_df['TEST'] == "ADD" ) ].copy()


# create list for twoBItToFa
# e.g.
# chr1:0-189
# [start,end)

sig_snps_df['seq_start'] = sig_snps_df.BP - 26 # take 25 bases to the left of index position
sig_snps_df['seq_end'] = sig_snps_df.BP + 27 # take 25 bases to the left of index position inclusive
sig_snps_df['seqList'] = 'chr' + sig_snps_df.CHR.map(str) + ":" + sig_snps_df.seq_start.map(str) + "-" +sig_snps_df.seq_end.map(str)

# write seqlist
sig_snps_df['seqList'].to_csv(OUTPUT_FILE, index=False,header=False)
sig_snps_df.to_csv(OUTPUT_FILE_KEY, index=False,header=True, sep="\t")

print("wrote seqlist to {}".format(OUTPUT_DIR))
