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

# needs ['GRID', 'GENDER', 'YOB']
DEMOGRAPHICS_FILE = "/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/BV227_A3_Capra_Preterm_A3_Phenotype_Demographics_date_formatted.csv"
PCA_FILE = "/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/pca/15pca_race_is_w__biovu_thru_2019_02.eigenvec"
WHITE_GRIDS_FILE="/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/grid_lists/race_is_w_grids.txt"

# TO DO: rename
OUTPUT_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/covariates/covar_white_YOB_12PC__{}.tsv".format(DATE)




# -----------
# MAIN
# -----------

# LOAD DATA
demo_df = pd.read_csv(DEMOGRAPHICS_FILE, sep="," )
pca_df = pd.read_csv(PCA_FILE, sep="\t")
pca_df.rename(columns={'FID': 'GRID'}, inplace=True)
white_grids_df = pd.read_csv(WHITE_GRIDS_FILE, sep="\t", names=['GRID'])


# FORMAT DOB
demo_df['DOB'] = pd.to_datetime(demo_df['DOB'], format='%Y-%m-%d')
demo_df['YOB'] = demo_df.DOB.map(lambda x: x.year)


# INCLUDE COVARAITES OF INTEREST

# - use pca from plink run
pca_keep_df = pca_df.loc[:, ['GRID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12' ]]
demo_keep_df = demo_df.loc[:, ['GRID', 'GENDER', 'YOB']]
temp_merge_df = pd.merge(demo_keep_df, pca_keep_df, how='inner', on='GRID')

# keep only 'white' grids
final_df = temp_merge_df[temp_merge_df.GRID.isin(white_grids_df.GRID)].copy()


# rm dups
dups = final_df[final_df.duplicated(keep='first')]
duplicate_index = dups.index
print("number of duplciate GRIDS removed: {}".format(dups.GRID.nunique()))
final_df.drop(duplicate_index, axis=0, inplace=True)

# drop GENDER
final_df.drop(['GENDER'], axis=1, inplace=True)



# final clean up
final_no_dups = final_df[~final_df.duplicated(subset=['GRID'])].copy()
final_no_dups["IID"] = final_no_dups['GRID']
final_no_dups.rename(columns={'GRID':'FID'}, inplace=True)
final_no_dups_reorder = final_no_dups.loc[:, ['FID','IID', 'YOB', 'PC1', 'PC2', 'PC3','PC4', 'PC5','PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12']].copy()

# one hot encode
# concat_df = pd.concat([final_no_dups_reorder, pd.get_dummies(final_no_dups_reorder.ARRAY)], axis=1)
# concat_df.drop('ARRAY', axis=1, inplace=True)

final_no_dups_reorder.to_csv(OUTPUT_FILE, index=False, sep="\t")
print("Done. check: {}".format(OUTPUT_FILE))
