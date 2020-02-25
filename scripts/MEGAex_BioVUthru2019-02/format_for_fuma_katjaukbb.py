#!/bin/python
# This script will conver the output of plink2 --glm to be ready for FUMA.
#
#
#
# Abin Abraham
# created on: 2019-07-06 09:43:36


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')


# -----------
# MAIN
# -----------

FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL.csv"
OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL_for_fuma_missP_removed.csv"

#
raw_df = pd.read_csv(FILE, sep=",")
df = raw_df.loc[raw_df['TEST'] == "ADD", :].copy()

df['effect_allle'] = df.A1
# df['non_effect_allele'] = df.apply(lambda x: x.REF if x.A1 == x.ALT else x.ALT, axis=1)

df.rename(columns={'#CHROM':'CHR', 'POS':'BP', 'ID':'SNP'}, inplace=True)


final_df = df.loc[:, ['SNP','CHR','BP','effect_allle','P','OR','missing_p']].copy()
final_df.sort_values('P', inplace=True)


# remove XY  to X region
final_df.loc[final_df['CHR'] =='XY', 'CHR'] = "X"


# remove MT
temp_to_write_df = final_df.loc[final_df['CHR'] !='MT'].copy()

# remove missing SNSPs
to_write_df = temp_to_write_df.loc[temp_to_write_df['missing_p'] > 0.00001, :]


to_write_df.to_csv(OUTPUT_FILE, sep="\t", index=False, float_format="%e")

#
#
#
# SNP | snpid | markername | rsID: rsID
# CHR | chromosome | chrom: chromosome
# BP | pos | position: genomic position (hg19)
# A1 | effect_allele | allele1 | alleleB: affected allele
# A2 | non_effect_allele | allele2 | alleleA: another allele
# P | pvalue | p-value | p_value | frequentist_add_pvalue | pval: P-value (Mandatory)
# OR: Odds Ratio
# Beta | be: Beta
# SE: Standard error
