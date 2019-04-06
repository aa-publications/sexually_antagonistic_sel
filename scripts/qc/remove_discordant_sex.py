#!/bin/python
# This script will parse *.check_sex plink output and return 1) distribution plot of F-estimate and 2) individuals w/ discordant sex to remove. 
#   - output are written to hte same directory as input file 
#
#
# Abin Abraham
# created on: 2019-04-05 11:02:56


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns 

DATE = datetime.now().strftime('%Y-%m-%d')


import argparse
parser = argparse.ArgumentParser(description='Example with nonoptional arguments')
 
### REQUIRED ARGUMENTS IN ORDER
parser.add_argument('sex_check_file', action='store', type=str)
parser.add_argument('prefix', action='store', help="output_prefix", type=str)
parser.add_argument('fig_output_file', action='store', help="full path to figure file", type=str)
parser.add_argument('ids_output_file', action='store', help="full path to ids_to_remove file", type=str)
 
results = parser.parse_args()
sex_chk_file = results.sex_check_file
output_prefix= results.prefix
output_fig_file = results.fig_output_file
output_file = results.ids_output_file


# OUTPUT FILES
fig_out_file = output_fig_file
grids_to_remove_file = output_file


# !!!! MODIFY BEFORE DEPLOYMENT !!! 
# testing
# sex_chk_file = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGA_ex_Array_Batch13_Cox_17_preQC_GRID_check_sex.sexcheck"
# output_prefix = "demo"

# -----------
# MAIN
# ----------- 
print("Identifying individuals with discordant sex based on F coef calculated on X Chromosome...")
df = pd.read_csv(sex_chk_file, sep="\s+")
failed_df = df.loc[df['STATUS'] != 'OK'].copy()

sns.distplot(df.F, rug = True,  color='blue', label="OK(n={:,})".format(df.shape[0]-failed_df.shape[0]), 
                      kde_kws={"alpha": 0.7, "lw": 2, })

sns.distplot(failed_df.F, rug = True, color='y',label="PROBLEM(n={:,})".format(failed_df.shape[0]))    
plt.title("F_coef for X Chromosome\nDataset:{}".format(output_prefix))
plt.axvline(x=0.2, color='r', alpha = 0.6, linestyle='dashed')
plt.axvline(x=0.8, color='r', alpha = 0.6, linestyle='dashed')
plt.legend(loc='upper left')
  
plt.savefig(fig_out_file)
print("chrX_Fcoef distribution saved to:\n\t{}".format(fig_out_file))


# calc individuals to remove 
dedup_df = failed_df[~failed_df.duplicated(subset=['FID','IID'], keep='first')]
dedup_df.loc[:, ['FID','IID']].to_csv(grids_to_remove_file, index=False, sep="\t", header=False)
print("Discordant IDs written to:\n\t{}".format(grids_to_remove_file))
