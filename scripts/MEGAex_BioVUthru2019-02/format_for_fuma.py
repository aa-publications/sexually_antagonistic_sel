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

FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/2019_07_21_logistic.assoc.logistic"
OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/2019_07_21_logistic.assoc.logistic.for_fuma.tsv1"

# FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_11_21_imputed_mega/snps_only_pl2_covar_standar.PHENO1.glm.logistic.hybrid"
# OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_11_21_imputed_mega/snps_only_pl2_covar_standar.PHENO1.glm.logistic.hybrid.for.fuma.tsv"


is_plink2=False

#

if is_plink2:
    col_names=["#CHROM","BP","SNP","REF", "ALT", "A1", "FIRTH?", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "P"]
    raw_df = pd.read_csv(FILE, sep="\t", names=col_names)
    df = raw_df.loc[raw_df['TEST'] == "ADD", :].copy()
else:
    col_names=["CHR","SNP", "BP", "A1", "TEST", "NMISS", "OR", "STAT", "P"]
    raw_df = pd.read_csv(FILE, sep="\s+", names=col_names)
    df = raw_df.loc[raw_df['TEST'] == "ADD", :].copy()

df['effect_allle'] = df.A1
# df['non_effect_allele'] = df.apply(lambda x: x.REF if x.A1 == x.ALT else x.ALT, axis=1)

df.rename(columns={'#CHROM':'CHR', 'POS':'BP', 'ID':'SNP'}, inplace=True)


filtered_df = df.loc[:, ['SNP','CHR','BP','effect_allle','P','OR']].copy()
filtered_df['CHR'] = filtered_df['CHR'].map(int)
# keep only autosoems
final_df = filtered_df.loc[filtered_df['CHR'] <23].copy()


final_df['float_P'] = final_df.P.apply(lambda x: float(x))
final_df.sort_values('float_P', inplace=True)
final_df.drop(columns={'float_P'}, inplace=True)


final_df.to_csv(OUTPUT_FILE, sep="\t", index=False, float_format="%e")

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
