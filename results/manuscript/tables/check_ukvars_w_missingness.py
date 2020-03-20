#!/bin/python
# This script will check uk variants that show differentital missingness if they have any BLAT hits with X or Y.
#
#
#
# Abin Abraham
# created on: 2020-03-20 08:27:04



import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')



uk_bil_file = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/web_blat/ukbil_w_score/chrxy_UKBil_probes_webscores.txt"
uk_wcsf_file = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/web_blat/ukwcsf_w_score/chrxy_ukwcsf_probes_webscores.txt"
var_file = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/tables/ukbb_sig_snps_missing_bysex.csv"



###
###    FUNCTIONS
###


def load_and_format_blat(uk_bil_file, uk_wcsf_file ):
    # LOAD WEBBLAT HITS
    ukbil_webblat_df = pd.read_csv(uk_bil_file, sep="\t", names=['targetChr','targetStart','targetEnd','queryName','score','per_identity'])
    ukwcsf_webblat_df = pd.read_csv(uk_wcsf_file, sep="\t", names=['targetChr','targetStart','targetEnd','queryName','score','per_identity'])


    # format webblat resutls
    get_start = lambda x: int(x.split(":")[-1].split("-")[0])
    get_end = lambda x: int(x.split(":")[-1].split("-")[1])


    # add query start and end
    ukbil_webblat_df['q.start'] = ukbil_webblat_df.queryName.apply(get_start)
    ukbil_webblat_df['q.end'] = ukbil_webblat_df.queryName.apply(get_end)

    ukwcsf_webblat_df['q.start'] = ukwcsf_webblat_df.queryName.apply(get_start)
    ukwcsf_webblat_df['q.end'] = ukwcsf_webblat_df.queryName.apply(get_end)

    ukbil_webblat_df['q.length'] = ukbil_webblat_df['q.end'] - ukbil_webblat_df['q.start']
    ukwcsf_webblat_df['q.length'] = ukwcsf_webblat_df['q.end'] -ukwcsf_webblat_df['q.start']

    # add bim_mapped_rsID
    ukbil_webblat_df['bim_mapped_rsID'] = ukbil_webblat_df.queryName.apply(lambda x: x.split(',')[-1].split(":")[0])
    ukwcsf_webblat_df['bim_mapped_rsID'] = ukwcsf_webblat_df.queryName.apply(lambda x: x.split(',')[-1].split(":")[0])


    return ukbil_webblat_df, ukwcsf_webblat_df



def uk_array_filter(og_blat_df, filter_perID=False):
    """ required each blat hit to be a) ≥ 40 bp in length, b) overlap the middle base (where the variant is), and c) ≥ 90% id"""

    blat_df = og_blat_df.copy()

    uk_qlength_bool = blat_df['q.length'] >= 40 # webblat is 0 start
    incl_center_bool = (blat_df['q.end'] > 35) & (blat_df['q.start'] < 35)
    uk_perID_bool = blat_df['per_identity'] >= 90

    if filter_perID:
        blat_df['pass_blat_filters'] = uk_qlength_bool  & incl_center_bool & uk_perID_bool
    else:
        blat_df['pass_blat_filters'] = uk_qlength_bool  & incl_center_bool


    return blat_df, og_blat_df


def pick_best_blat_by_score(blat_df, blat_filter=True):
    # pick BLAT hit with the highest blat score as teh 'best' blat hit (imposes a 1:1 relationship between variant and blat hit)
    # if blat_filter == True, then filter out bad blat hits and then choose teh 'best' blat hit

    temp_blat_df = blat_df.copy()

    if blat_filter:
        good_blat_hits = temp_blat_df.loc[temp_blat_df['pass_blat_filters']==True].reset_index(drop=True)
    else:
        good_blat_hits = temp_blat_df.reset_index(drop=True)

    best_df = good_blat_hits.loc[good_blat_hits.groupby('bim_mapped_rsID')['score'].idxmax()].reset_index(drop=True)

    return best_df


# %%
# -----------
# MAIN
# -----------


###
###    READ
###

ukbil_webblat_df, ukwcsf_webblat_df = load_and_format_blat(uk_bil_file, uk_wcsf_file)

# load variants that show differential missingness between males and females
miss_vars_df = pd.read_csv(var_file, sep=",")
miss_vars = miss_vars_df.SNP.unique()


# keep only variants with missingness between fmeales and males
miss_vars_ukbil_df = ukbil_webblat_df[ukbil_webblat_df.bim_mapped_rsID.isin(miss_vars)].copy()
miss_vars_ukwcsf_df = ukwcsf_webblat_df[ukwcsf_webblat_df.bim_mapped_rsID.isin(miss_vars)].copy()


###
###    FILTER BASED on 'HOMOLOGY' criteria
###

filt_miss_vars_ukbil_df,_ = uk_array_filter(miss_vars_ukbil_df, filter_perID=True)
filt_miss_vars_ukwcsf_df,_ = uk_array_filter(miss_vars_ukwcsf_df, filter_perID=True)


###
###    PICK THE BEST BLAT HIT (based on highest blat score)
###


# pool the two arrays used in UK Biobank
pooled_uk_df = pd.concat([filt_miss_vars_ukbil_df, filt_miss_vars_ukwcsf_df], axis=0)
best_blat_hit_uk_df = pick_best_blat_by_score(pooled_uk_df, blat_filter=True)


best_blat_hit_uk_df.to_csv('uk_var_w_missingness_best_blatscore_xy.tsv', sep="\t",index=False)

best_blat_hit_uk_df.bim_mapped_rsID.nunique()
best_blat_hit_uk_df.shape


# temp_file ="best_blatscore_xy_hit_length_filtered_uk_gwas.tsv.gz"
# temp_df = pd.read_csv(temp_file,sep='\t')
# temp_df[temp_df.SNP.isin(miss_vars)]
# temp_df.head()
