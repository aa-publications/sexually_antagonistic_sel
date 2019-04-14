#!/bin/python
# This script will create the covariates table for PheWAS run on all individuals with MEGAx array data in BioVU.
#           * keeps 5PCs from plink pca 
#           * add age as YOB
#           * add batch
#           * geneder is removed! 
#
# Abin Abraham
# created on: 2019-04-14 10:21:24


import os
import sys
import numpy as np
import pandas as pd
import datetime


DATE = datetime.date.today()
# -----------
# PATHS
# -----------

# 
# **************** USER MUST MODIFY ****************
# 

DEMOGRAPHICS_FILE = "/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/demographics/Capra_Preterm_A2_Phenotype.xlsx"
ARRAY_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/covariates/mega_fid_batch_list_2019-04-14.tsv"
PCA_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/covariates/inter_pca__merged_MEGA_2019-04-13.eigenvec"

OUTPUT_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/covariates/covar_CEU_merged_mega_age_5PC_{}.tsv".format(DATE)


use_plink_pca=True

# -----------
# MAIN
# -----------

# LOAD DATA
demo_df = pd.read_excel(DEMOGRAPHICS_FILE, skiprows=0)
demo_df.drop(['Unnamed: 1', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15'], axis=1, inplace=True)

array_df = pd.read_csv(ARRAY_FILE, sep="\t", header=None, names=['GRID', 'ARRAY'])
array_df.ARRAY = array_df.ARRAY.apply(lambda x: x[:-4])

pca_df = pd.read_csv(PCA_FILE, sep="\t")
pca_df.rename(columns={'FID': 'GRID'}, inplace=True)


# FORMAT DOB
demo_df['DOB'] = pd.to_datetime(demo_df['DOB'], format='%Y')
demo_df['YOB'] = demo_df.DOB.map(lambda x: x.year)


# INCLUDE COVARAITES OF INTEREST

if use_plink_pca:
    # - use pca from plink run
    pca_keep_df = pca_df.loc[:, ['GRID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']]

    demo_keep_df = demo_df.loc[:, ['GRID', 'GENDER', 'YOB']]
    temp_merge_df = pd.merge(demo_keep_df, pca_keep_df, how='inner', on='GRID')

    # add array information
    final_df = pd.merge(temp_merge_df, array_df, how='inner', on="GRID")


else:
    # - use pca provided in demographics file from BioVU
    # also add array information
    demo_keep_df = demo_df.loc[:, ['GRID', 'GENDER', 'YOB', 'CEU','YRI','JPT_CHB']]
    final_df = pd.merge(demo_keep_df, array_df, how='inner', on='GRID')


dups = final_df[final_df.duplicated(keep='first')]
duplicate_index = dups.index
print("number of duplciate GRIDS removed: {}".format(dups.GRID.nunique()))

final_df.drop(duplicate_index, axis=0, inplace=True)

# TABULATE MISSING DATA
final_grids = set(final_df.GRID.values)
set1 = set(demo_keep_df.GRID.unique())
try:
    set2 = set(pca_keep_df.GRID.unique())
except NameError:
    set2 = {}

set3 = set(array_df.GRID.unique())

all_sets = set1.union(set2).union(set3)

grids_dropped_num = len(all_sets.difference(final_grids))

print("Number of GRIDS lost in merges: {}".format(grids_dropped_num))

final_df.drop(['GENDER'], axis=1, inplace=True)

# final clean up
final_no_dups = final_df[~final_df.duplicated(subset=['GRID'])].copy()
final_no_dups["IID"] = final_no_dups['GRID']
final_no_dups.rename(columns={'GRID':'FID'}, inplace=True)
final_no_dups_reorder = final_no_dups.loc[:, ['FID','IID', 'YOB', 'PC1', 'PC2', 'PC3','PC4', 'PC5', 'ARRAY']].copy()

# one hot encode
concat_df = pd.concat([final_no_dups_reorder, pd.get_dummies(final_no_dups_reorder.ARRAY)], axis=1)
concat_df.drop('ARRAY', axis=1, inplace=True)

concat_df.to_csv(OUTPUT_FILE, index=False, sep="\t")
print("Done. check: {}".format(OUTPUT_FILE))
