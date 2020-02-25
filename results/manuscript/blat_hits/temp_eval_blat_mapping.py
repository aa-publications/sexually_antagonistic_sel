#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'



import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
fpath='/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
prop = fm.FontProperties(fname=fpath, size=16)
bigprop = fm.FontProperties(fname=fpath, size=20)






# CONSTANTS:
GWAS_P_THRESH = 5*10**-8
SUGG_GWAS_P_THRESH = 1*10**-6
MISS_P_THRESHOLD = 0.00001


# -----------
# PATHS
# -----------

# WEBBLAT
WEBBLAT_UKBIL="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/web_blat/ukbil_w_score/chrxy_UKBil_probes_webscores.txt"
WEBBLAT_UKWCSF="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/web_blat/ukwcsf_w_score/chrxy_ukwcsf_probes_webscores.txt"
WEBBLAT_BV="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/blat_psl/chrxy_MEGAx_probes_v1_blat_webscores.txt"

# GWAS SUMMARY STATS
GWAS_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/"
UK_GWAS_FILE = os.path.join(GWAS_DIR, "data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL.csv")
BV_GWAS_FILE = os.path.join(GWAS_DIR, "results/2019_07_21_logistic/2019_07_21_logistic.assoc.logistic")

# OUTPUT
OUTPUT_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits"

MEGA_PROBES_FILE="/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/MEGAEx_BioVU_15075710_A1_name_snp_probe_chr_pos.csv"

# -----------
# MAIN
# -----------
sys.path.append("/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits")
from plot_blat_analysis import load_and_format_gwas, load_and_format_blat


# load gwas
uk_gwas_df, bv_gwas_df = load_and_format_gwas(UK_GWAS_FILE, BV_GWAS_FILE)
ukbil_webblat_df, ukwcsf_webblat_df, bv_webblat_df = load_and_format_blat(WEBBLAT_UKBIL, WEBBLAT_UKWCSF, WEBBLAT_BV)


keep = pd.DataFrame()
for snp in bv_gwas_df.SNP.unique():

    keep = keep.append(bv_webblat_df[bv_webblat_df.queryName.apply(lambda x: x.find(snp)) > -1])

keep.shape


bv_gwas_df.head()
bv_webblat_df.head(15)

bv_gwas_df.SNP.nunique()
bv_gwas_df.SNP.isin(bv_webblat_df.bim_mapped_rsID).sum()


bv_webblat_df.bim_mapped_rsID.unique()[-300:-200]


bv_webblat_df.head()

bv_gwas_df[~bv_gwas_df.SNP.isin(bv_webblat_df.bim_mapped_rsID)]
