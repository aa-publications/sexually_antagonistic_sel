#!/bin/python
# This script identifies individuals that are > 3 S.D from mean for obersved heterozygosity rate.
#               * Requires output of plink --het as input.
#               * calcualtes obserevd heterozygosity rate for each individuals
#               * outputs individuals who are > 3 S.D. away (to be removed)
#
# 
#           WARNING: if there are no idnividuals with outlier heterozygosity values, output file is not written
# 
#
# Abin Abraham
# created on: 2018-10-04 11:44:01


import os
import sys
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='calc individuals with outlier heterozygozity')

# REQUIRED ARGUMENTS IN ORDER
parser.add_argument('het_file', action='store', type=str)
parser.add_argument('output_file', action='store', type=str)

# retrieve passed arguments
args = parser.parse_args()
het_file_path = args.het_file
output_file_path = args.output_file

print("Running: {}".format(sys.argv[:]))

df = pd.read_csv(het_file_path, sep="\s+")
df['obs_het_rate'] = (df['N(NM)'] - df['O(HOM)'])/df['N(NM)']

obs_het_rate_mean = df['obs_het_rate'].mean()
obs_het_rate_std_threshold = df['obs_het_rate'].std()*3

high_thresh = obs_het_rate_mean + obs_het_rate_std_threshold
low_thresh = obs_het_rate_mean - obs_het_rate_std_threshold


removed_df = df.loc[(df['obs_het_rate'] > high_thresh) | (df['obs_het_rate'] < low_thresh)].copy()
print("Number of individuals with het rate > 3 S.D from mean: {}".format(removed_df.shape[0]))


## output tsv with two columsn FID and IID 
if (removed_df.shape[0] > 0):
    removed_df.to_csv(output_file_path, sep="\t", index=False, header=False, columns = ['FID','IID'])
else: 
    print("No individuals with outlier heterozygosity values...")

print("Done! Check: {}".format(output_file_path))