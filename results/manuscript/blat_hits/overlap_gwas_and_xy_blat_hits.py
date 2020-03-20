
#!/bin/python
# This script will
#
#
#
# Abin Abraham
# created on: 2020-02-08 14:55:47


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




# %%
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
OUTPUT_TABLE_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/tables"


# %%
# -----------
# FUNCTIONS
# -----------

def load_and_format_gwas(uk_file, bv_file):

    # LOAD GWAS HITS
    raw_uk_df = pd.read_csv( uk_file, sep=",")
    raw_bv_df = pd.read_csv( bv_file, sep="\s+")

    # format results
    uk_gwas_df = raw_uk_df.loc[raw_uk_df.missing_p > MISS_P_THRESHOLD].copy()
    print("UKBB: Removed {:,} out of {:,} snps with singificant missingness b/w cases and controls.".format(raw_uk_df.shape[0] - uk_gwas_df.shape[0], raw_uk_df.shape[0]))
    uk_gwas_df['chr_pos'] = uk_gwas_df.CHR.map(str) + ":" + uk_gwas_df.BP.map(str)

    bv_gwas_df = raw_bv_df.loc[(raw_bv_df.TEST == "ADD") & (raw_bv_df.CHR < 23),].copy()
    bv_gwas_df['chr_pos'] = bv_gwas_df.CHR.map(str) + ":" + bv_gwas_df.BP.map(str)

    # add significance column
    uk_gwas_df['stat_sig'] = False
    uk_gwas_df.loc[uk_gwas_df['P']< GWAS_P_THRESH, 'stat_sig'] = True

    bv_gwas_df['stat_sig'] = False
    bv_gwas_df.loc[bv_gwas_df['P']< GWAS_P_THRESH, 'stat_sig'] = True


    return uk_gwas_df, bv_gwas_df

def load_and_format_blat(uk_bil_file, uk_wcsf_file, bv_file):
    # LOAD WEBBLAT HITS
    ukbil_webblat_df = pd.read_csv(uk_bil_file, sep="\t", names=['targetChr','targetStart','targetEnd','queryName','score','per_identity'])
    ukwcsf_webblat_df = pd.read_csv(uk_wcsf_file, sep="\t", names=['targetChr','targetStart','targetEnd','queryName','score','per_identity'])
    bv_webblat_df = pd.read_csv(bv_file, sep="\t", names=['targetChr','targetStart','targetEnd','queryName','score','per_identity'])

    # format webblat resutls
    get_start = lambda x: int(x.split(":")[-1].split("-")[0])
    get_end = lambda x: int(x.split(":")[-1].split("-")[1])


    # add query start and end
    ukbil_webblat_df['q.start'] = ukbil_webblat_df.queryName.apply(get_start)
    ukbil_webblat_df['q.end'] = ukbil_webblat_df.queryName.apply(get_end)

    ukwcsf_webblat_df['q.start'] = ukwcsf_webblat_df.queryName.apply(get_start)
    ukwcsf_webblat_df['q.end'] = ukwcsf_webblat_df.queryName.apply(get_end)

    bv_webblat_df['q.start'] = bv_webblat_df.queryName.apply(get_start)
    bv_webblat_df['q.end'] = bv_webblat_df.queryName.apply(get_end)

    ukbil_webblat_df['q.length'] = ukbil_webblat_df['q.end'] - ukbil_webblat_df['q.start']
    ukwcsf_webblat_df['q.length'] = ukwcsf_webblat_df['q.end'] -ukwcsf_webblat_df['q.start']
    bv_webblat_df['q.length'] = bv_webblat_df['q.end'] -bv_webblat_df['q.start']

    # add bim_mapped_rsID
    ukbil_webblat_df['bim_mapped_rsID'] = ukbil_webblat_df.queryName.apply(lambda x: x.split(',')[-1].split(":")[0])
    ukwcsf_webblat_df['bim_mapped_rsID'] = ukwcsf_webblat_df.queryName.apply(lambda x: x.split(',')[-1].split(":")[0])
    bv_webblat_df['bim_mapped_rsID'] = bv_webblat_df.queryName.apply(lambda x: x.split(',')[1])


    return ukbil_webblat_df, ukwcsf_webblat_df, bv_webblat_df

def bv_blat_filter(og_blat_df):
    """ required each blat hit to be a) ≥ 40 bp in length, b) overlap the last base, and c) ≥ 90% id"""


    blat_df = og_blat_df.copy()

    qlength_bool = blat_df['q.length'] >= 40
    end_match_bool = blat_df['q.end'] == 50
    perID_bool = blat_df['per_identity'] >= 90

    blat_df['pass_blat_filters'] = qlength_bool  & end_match_bool & perID_bool

    return blat_df, og_blat_df

def uk_array_filter(og_blat_df):
    """ required each blat hit to be a) ≥ 40 bp in length, b) overlap the middle base (where the variant is), and c) ≥ 90% id"""

    blat_df = og_blat_df.copy()

    uk_qlength_bool = blat_df['q.length'] >= 40 # webblat is 0 start
    incl_center_bool = (blat_df['q.end'] > 35) & (blat_df['q.start'] < 35)
    uk_perID_bool = blat_df['per_identity'] >= 90

    blat_df['pass_blat_filters'] = uk_qlength_bool  & incl_center_bool & uk_perID_bool

    return blat_df, og_blat_df

def pick_best_blat(blat_df, blat_filter=True):
    # pick BLAT hit with the highest %ID as teh 'best' blat hit (imposes a 1:1 relationship between variant and blat hit)
    # if blat_filter == True, then filter out bad blat hits and then choose teh 'best' blat hit

    temp_blat_df = blat_df.copy()

    if blat_filter:
        good_blat_hits = temp_blat_df.loc[temp_blat_df['pass_blat_filters']==True].reset_index(drop=True)
    else:
        good_blat_hits = temp_blat_df.reset_index(drop=True)

    best_df = good_blat_hits.loc[good_blat_hits.groupby('bim_mapped_rsID')['per_identity'].idxmax()].reset_index(drop=True)

    return best_df

def merge_blat_and_gwas(blat_df , gwas_df):

    # merge  gwas with BLAT hits
    to_merge_df = blat_df.copy()
    to_merge_df['target_Start-End'] = to_merge_df.targetStart.map(str) + "-" + to_merge_df.targetEnd.map(str)
    to_merge_df['query_Start-End'] = to_merge_df['q.start'].map(str) + "-" + to_merge_df['q.end'].map(str)
    keep_cols = ['targetChr','target_Start-End', 'score', 'per_identity','query_Start-End', 'q.length','bim_mapped_rsID','pass_blat_filters']

    merged_df = pd.merge(gwas_df, to_merge_df.loc[:, keep_cols], left_on = "SNP", right_on="bim_mapped_rsID", how='left')

    print(f"Out of {blat_df.bim_mapped_rsID.nunique():,} SNPs with X or Y hits,")
    print(f"{(~merged_df.bim_mapped_rsID.isna()).sum():,}, SNPs in the GWAS summary stats were mapped.")

    return merged_df

def patch_violinplot(color_list ):
    from matplotlib.collections import PolyCollection
    ax = plt.gca()

    violins = [art for art in ax.get_children() if isinstance(art, PolyCollection)]
    colors = color_list*len(violins)
    for i in range(len(violins)):
        violins[i].set_edgecolor(colors[i])
        violins[i].set_lw(0.7)

def format_to_write(gwas_blat_df):
    keep_cols1 = ['chr_pos','SNP','A1','NMISS','OR','STAT','P', 'targetChr','target_Start-End', 'score','per_identity','query_Start-End','q.length','pass_blat_filters']
    temp_df = gwas_blat_df.loc[~gwas_blat_df.isna().any(1), keep_cols1]
    temp_df['target_seq'] = temp_df['targetChr'] + ":" + temp_df['target_Start-End']
    keep_cols2 = ['chr_pos','SNP','A1','NMISS','OR','STAT','P','target_seq', 'score','per_identity','query_Start-End','q.length','pass_blat_filters']
    new_names = ['chr_pos', 'snp','a1','n_samples','or','stat','p','target_seq', 'score','per_identity','query_Start-End','q.length','pass_blat_filters']

    to_write_df = temp_df.loc[:, keep_cols2].copy()
    to_write_df.columns = new_names

    return to_write_df

def format_to_write_gwas_sig_vars(og_gwas_blat_df):
    gwas_blat_df = og_gwas_blat_df.copy()

    gwas_blat_df['target_chr_pos'] = gwas_blat_df['targetChr'] + "_" +  gwas_blat_df['target_Start-End']
    keep_cols = ['chr_pos', 'SNP', 'A1', 'NMISS', 'OR', 'STAT', 'P',
        'target_chr_pos', 'score', 'per_identity', 'query_Start-End', 'q.length']

    to_write_df = gwas_blat_df.loc[gwas_blat_df['stat_sig']==True, keep_cols].copy()
    new_cols = ['chr_pos', 'snp', 'a1', 'n_samples', 'or', 'stat', 'gwas_p', 'blat_hit_chr_pos', 'blat_score', 'per_identity', 'query_start-end', 'query_length']
    to_write_df.columns = new_cols

    return to_write_df



# %%
if __name__ == '__main__':

    # -----------
    # MAIN
    # -----------


    # load data
    uk_gwas_df, bv_gwas_df = load_and_format_gwas(UK_GWAS_FILE, BV_GWAS_FILE)
    og_ukbil_webblat_df, og_ukwcsf_webblat_df, og_bv_webblat_df = load_and_format_blat(WEBBLAT_UKBIL, WEBBLAT_UKWCSF, WEBBLAT_BV)



    # -----------
    # APPLY BLAT FILTER
    # -----------
    # note: no variatns are removed; a column for filtering is added


    # biovu
    bv_webblat_df, all_bv_webblat_df =  bv_blat_filter(og_bv_webblat_df)

    # ukbb
    ukbil_webblat_df, all_ukbil_webblat_df = uk_array_filter(og_ukbil_webblat_df)
    ukwcsf_webblat_df, all_ukwcsf_webblat_df = uk_array_filter(og_ukwcsf_webblat_df)

    pd.value_counts(bv_webblat_df.pass_blat_filters)
    pd.value_counts(ukbil_webblat_df.pass_blat_filters)
    pd.value_counts(ukwcsf_webblat_df.pass_blat_filters)


    # -----------
    # pick the 'best' blat hit
    # -----------

    # pool both uk arrays
    both_uk_webblat_df = pd.concat([ukbil_webblat_df, ukwcsf_webblat_df],axis=0)

    # after filtering out bad blat hits
    uk_best_webblast_df = pick_best_blat(both_uk_webblat_df)
    bv_best_webblast_df = pick_best_blat(bv_webblat_df)

    # w/o filtering out bad blat hits
    uk_no_filt_best_webblast_df = pick_best_blat(both_uk_webblat_df, blat_filter=False)
    bv_no_filt_best_webblast_df = pick_best_blat(bv_webblat_df, blat_filter=False)


    # -----------
    # counts and overlap of variants and blat hits
    # -----------

    bv_gwas_df.shape
    uk_gwas_df.shape

    # num gwas variants that overlap at least one  blat xy hit - unfiltered
    bv_gwas_df[bv_gwas_df.SNP.isin(bv_webblat_df.bim_mapped_rsID)].shape
    uk_gwas_df[uk_gwas_df.SNP.isin(ukwcsf_webblat_df.bim_mapped_rsID)].shape

    (83083/798051)*100
    (128090/620040)*100

    blatxy_bv_one_to_many_df = bv_webblat_df[bv_webblat_df.bim_mapped_rsID.isin(bv_gwas_df.SNP)].copy()
    blatxy_uk_one_to_many_df = ukwcsf_webblat_df[ukwcsf_webblat_df.bim_mapped_rsID.isin(uk_gwas_df.SNP)].copy()

    blatxy_bv_one_to_many_df.groupby('bim_mapped_rsID').size().describe()
    blatxy_uk_one_to_many_df.groupby('bim_mapped_rsID').size().describe()


    # repeat above but filter to keep only high-quality hits
    blat_filt = lambda df: df.loc[df['pass_blat_filters']==True].copy()

    filt_bv_blat_df = blat_filt(blatxy_bv_one_to_many_df)
    filt_uk_blat_df = blat_filt(blatxy_uk_one_to_many_df)


    filt_uk_blat_df.sort_values( ['q.length', 'per_identity'], ascending=False)

    filt_bv_blat_df.bim_mapped_rsID.nunique()
    filt_uk_blat_df.bim_mapped_rsID.nunique()

    (4994/798051)*100
    (28555/620040)*100



    # -----------
    # merge blat hits with gwas summary stats
    # -----------

    uk_gwas_blat_df = merge_blat_and_gwas(uk_best_webblast_df, uk_gwas_df)
    bv_gwas_blat_df = merge_blat_and_gwas(bv_best_webblast_df, bv_gwas_df)


    # unfiltred ukbb
    uk_gwas_unfilt_blat_df = merge_blat_and_gwas(uk_no_filt_best_webblast_df, uk_gwas_df)
    uk_gwas_unfilt_blat_df.loc[uk_gwas_unfilt_blat_df['stat_sig']==True]

    # not significant variants
    suspicious_nonsig_bv_df = bv_gwas_blat_df[(~bv_gwas_blat_df.bim_mapped_rsID.isna()) & (bv_gwas_blat_df['stat_sig']==False)].sort_values(['q.length','per_identity'], ascending=False)
    suspicious_nonsig_uk_df =uk_gwas_blat_df[(~uk_gwas_blat_df.bim_mapped_rsID.isna()) & (uk_gwas_blat_df['stat_sig']==False)].sort_values(['q.length','per_identity'], ascending=False)

    high_xy_hits_for_non_sig_vars_bv.shape
    high_xy_hits_for_non_sig_vars_bv = suspicious_nonsig_bv_df.loc[(suspicious_nonsig_bv_df['per_identity']==100) & (suspicious_nonsig_bv_df['q.length']==50) ].copy()
    high_xy_hits_for_non_sig_vars_uk = suspicious_nonsig_uk_df.loc[(suspicious_nonsig_uk_df['per_identity']==100) & (suspicious_nonsig_uk_df['q.length']==71) ].copy()

    high_xy_hits_for_non_sig_vars_bv.to_csv(os.path.join(OUTPUT_DIR, 'high_xy_hits_for_non_sig_vars_bv.tsv'), sep="\t", index=False)
    high_xy_hits_for_non_sig_vars_uk.to_csv(os.path.join(OUTPUT_DIR,'high_xy_hits_for_non_sig_vars_uk.tsv'), sep="\t", index=False)


    uk_gwas_blat_clean_df = format_to_write(uk_gwas_blat_df)
    bv_gwas_blat_clean_df = format_to_write(bv_gwas_blat_df)

    uk_sig_gwas_blat_clean_df = format_to_write_gwas_sig_vars(uk_gwas_blat_df)
    bv_sig_gwas_blat_clean_df = format_to_write_gwas_sig_vars(bv_gwas_blat_df)




    # -----------
    # write tables ...
    # -----------

    # 1) raw blat hits to x or y

    og_ukbil_webblat_df.head()




    og_ukbil_webblat_df.rename(columns={'bim_mapped_rsID':'UK_BilEVE_rsID'}).to_csv(os.path.join(OUTPUT_DIR, 'raw_blat_xy_hits_ukbil.tsv'), sep="\t", index=False)
    og_ukwcsf_webblat_df.rename(columns={'bim_mapped_rsID':'UK_Axiom_rsID'}).to_csv(os.path.join(OUTPUT_DIR, 'raw_blat_xy_hits_ukwcsf.tsv'), sep="\t", index=False)
    og_bv_webblat_df.rename(columns={'bim_mapped_rsID':'BioVU_MEGAEx_rsID'}).to_csv(os.path.join(OUTPUT_DIR, 'raw_blat_xy_hits_bv.tsv'), sep="\t", index=False)


    uk_gwas_blat_clean_df.to_csv(os.path.join(OUTPUT_DIR, 'uk_gwas_best_blat_xy_hits.tsv'), sep="\t", index=False)
    bv_gwas_blat_clean_df.to_csv(os.path.join(OUTPUT_DIR, 'bv_gwas_best_blat_xy_hits.tsv'), sep="\t", index=False)

    uk_sig_gwas_blat_clean_df.to_csv(os.path.join(OUTPUT_DIR, 'uk_sig_gwas_best_blat_xy_hits.tsv'), sep="\t", index=False)
    bv_sig_gwas_blat_clean_df.to_csv(os.path.join(OUTPUT_DIR, 'bv_sig_gwas_best_blat_xy_hits.tsv'), sep="\t", index=False)

    # -- dictionary map SNP to whether it is GWAS SIGNIFICANT
    uksig_dict = dict(zip(uk_gwas_df.SNP, uk_gwas_df.stat_sig))
    bvsig_dict = dict(zip(bv_gwas_df.SNP, bv_gwas_df.stat_sig))


    all_bv_webblat_df['stat_sig'] = all_bv_webblat_df.bim_mapped_rsID.map(bvsig_dict)
    all_ukwcsf_webblat_df['stat_sig'] = all_ukwcsf_webblat_df.bim_mapped_rsID.map(uksig_dict)
    all_ukbil_webblat_df['stat_sig'] = all_ukbil_webblat_df.bim_mapped_rsID.map(uksig_dict)

    bv_best_webblast_df['stat_sig'] = bv_best_webblast_df.bim_mapped_rsID.map(bvsig_dict)
    uk_best_webblast_df['stat_sig'] = uk_best_webblast_df.bim_mapped_rsID.map(uksig_dict)



