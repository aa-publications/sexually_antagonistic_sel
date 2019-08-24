#!/bin/python
# This script will take plink assoc output and create a list of sequence to pull using twoBItToFa
#
#
#
#
#
# assoc_file: full path to plink output from assoc analysis
# this file must include the following headers:
#       >  BP
#       >  CHR
#       >  P
#       >  TEST
#
# Abin Abraham
# created on: 2019-07-24 10:36:57

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

import argparse

DATE = datetime.now().strftime('%Y-%m-%d')

###
#   INPUT FILES
###


### REQUIRED ARGUMENTS IN ORDER
parser = argparse.ArgumentParser(description='')
parser.add_argument('assoc_file', action='store', type=str)
parser.add_argument('gwas_pval_tresh', action='store', type=int)
parser.add_argument('flank_bp', action='store', type=int)
parser.add_argument('output_file', action='store', type=str)
parser.add_argument('output_file_key', action='store', type=str)
results = parser.parse_args()


ASSOC_FILE = results.assoc_file
PVAL_THRESH = results.gwas_pval_tresh
FLANK_BP = results.flank_bp
OUTPUT_FILE = results.output_file
OUTPUT_FILE_KEY = results.output_file_key


# pval_thresh = 0.5*10**-5
# seperator used by pd.read_csv (changed depeneding which plink version was used to genreate the ASSOC_FILE)
SEPERATOR=","




# ASSOC_FILE = '/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL.csv'
# PVAL_THRESH=1
# FLANK_BP=50
# OUTPUT_DIR = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results"
# OUTPUT_FILE = '/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_08_21_blat_pipeline/all_gwas_hits/50bp/seqList_for_twoBitToFa_2019_08_24.tsv'
# OUTPUT_FILE_KEY = '/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_08_21_blat_pipeline/all_gwas_hits/50bp/key_for_seqList_for_twoBitToFa_2019_08_24.tsv'


###
#   MAIN
###
print("Running:\n{}".format('\n\t'.join(sys.argv)))
print(" ")

# select SNPS to create seqList
as_df = pd.read_csv(ASSOC_FILE, sep=SEPERATOR)
sig_snps_df = as_df.loc[ (as_df['P'] < PVAL_THRESH) & (as_df['P'] !=0) & (as_df['TEST'] == "ADD" ) ].copy()


# create list for twoBItToFa
# e.g.
# chr1:0-189
# [start,end)

left_flank=int(FLANK_BP/2)
right_flank=int(FLANK_BP/2)
print("num bases on left and rigth: {}, {}".format(left_flank, right_flank))
sig_snps_df['seq_start'] = sig_snps_df.BP - left_flank # take 25 bases to the left of index position
sig_snps_df['seq_end'] = sig_snps_df.BP + right_flank # take 25 bases to the left of index position inclusive
sig_snps_df['seqList'] = 'chr' + sig_snps_df.CHR.map(str) + ":" + sig_snps_df.seq_start.map(str) + "-" +sig_snps_df.seq_end.map(str)

# write seqlist
sig_snps_df['seqList'].to_csv(OUTPUT_FILE, index=False,header=False)
sig_snps_df.to_csv(OUTPUT_FILE_KEY, index=False,header=True, sep="\t")

print("wrote seqlist to {}".format(os.path.basename(OUTPUT_FILE)))
