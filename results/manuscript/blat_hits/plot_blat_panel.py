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
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator,MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
fpath='/dors/capra_lab/users/abraha1/conda/envs/py36_r_ml/lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf'
prop = fm.FontProperties(fname=fpath, size=16)
smallprop = fm.FontProperties(fname=fpath, size=12)
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
OUTPUT_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits/bivariate_blat_distributions"
OUTPUT_TABLE_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/tables"

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

def bv_blat_filter(og_blat_df, filter_perID=False):
    """ required each blat hit to be a) ≥ 40 bp in length, b) overlap the last base, and c) ≥ 90% id"""


    blat_df = og_blat_df.copy()

    qlength_bool = blat_df['q.length'] >= 40
    end_match_bool = blat_df['q.end'] == 50
    perID_bool = blat_df['per_identity'] >= 90

    if filter_perID:
        blat_df['pass_blat_filters'] = qlength_bool  & end_match_bool & perID_bool
    else:
        blat_df['pass_blat_filters'] = qlength_bool  & end_match_bool



    return blat_df, og_blat_df

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

def pick_best_blat_by_score(blat_df, blat_filter=True):
    # pick BLAT hit with the highest %ID as teh 'best' blat hit (imposes a 1:1 relationship between variant and blat hit)
    # if blat_filter == True, then filter out bad blat hits and then choose teh 'best' blat hit

    temp_blat_df = blat_df.copy()

    if blat_filter:
        good_blat_hits = temp_blat_df.loc[temp_blat_df['pass_blat_filters']==True].reset_index(drop=True)
    else:
        good_blat_hits = temp_blat_df.reset_index(drop=True)

    best_df = good_blat_hits.loc[good_blat_hits.groupby('bim_mapped_rsID')['score'].idxmax()].reset_index(drop=True)

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

def plot_bivar_best_seqid_vs_length(best_plt_df, dataset, ax):

    temp_sig_stat= best_plt_df.loc[best_plt_df['stat_sig']==True].copy()



    clip_min, clip_max =  best_plt_df['prop_of_probe_length'].min(), best_plt_df['prop_of_probe_length'].max()
    yclip_min, yclip_max =  best_plt_df['per_identity'].min(), best_plt_df['per_identity'].max()


    sns.kdeplot(best_plt_df['prop_of_probe_length'], best_plt_df['per_identity'], ax=ax, shade=True, shade_lowest=True,  cbar=True, clip=((clip_min, clip_max),(90,100)))

    sns.scatterplot(data=temp_sig_stat, x='prop_of_probe_length', y='per_identity', alpha=1, s=100, marker="o", color='indianred', ax=ax)


    ax.set_xlim(clip_min, clip_max)
    ax.set_title(f"{dataset} Array")
    ax.set_xlabel("Match length\n(% probe length)")
    ax.set_ylabel("Match % identity\nto chrX or chrY")

    sns.despine( ax=ax, top=True, right=True)

    return fig

def plot_bivar_best_score_vs_length(best_plt_df, dataset, ax):

    temp_sig_stat= best_plt_df.loc[best_plt_df['stat_sig']==True].copy()



    clip_min, clip_max =  best_plt_df['score'].min(), best_plt_df['score'].max()
    yclip_min, yclip_max =  best_plt_df['per_identity'].min(), best_plt_df['per_identity'].max()


    sns.kdeplot(best_plt_df['score'], best_plt_df['per_identity'], ax=ax, shade=True, shade_lowest=True,  cbar=True, clip=((clip_min, clip_max),(90,100)))

    sns.scatterplot(data=temp_sig_stat, x='score', y='per_identity', alpha=1, s=100, marker="o", color='indianred', ax=ax)


    ax.set_xlim(clip_min, clip_max)
    ax.set_title(f"{dataset} Array")
    ax.set_xlabel("Match BLAT score")
    ax.set_ylabel("Match % identity\nto chrX or chrY")

    sns.despine( ax=ax, top=True, right=True)

    return fig


def plot_bivar_many___vs_length(best_plt_df, dataset, ax, xvar):
    temp_sig_stat= best_plt_df.loc[best_plt_df['stat_sig']==True].copy()

    clip_min, clip_max =  best_plt_df[xvar].min(), best_plt_df[xvar].max()
    yclip_min, yclip_max =  best_plt_df['per_identity'].min(), best_plt_df['per_identity'].max()


    kdeax = sns.kdeplot(best_plt_df[xvar], best_plt_df['per_identity'], ax=ax, shade=True, shade_lowest=True,  cbar=True, clip=((clip_min, clip_max),(90,100)))
    sns.scatterplot(data=temp_sig_stat, x=xvar, y='per_identity', alpha=1, s=120,  color='indianred', style='SNP', ax=ax)


    ax.set_xlim(clip_min, clip_max)
    ax.set_ylim(yclip_min, yclip_max)
    ax.set_title(f"{dataset}", fontproperties=prop)

    if xvar=="prop_of_probe_length":
        ax.set_xlabel("Match length (% of probe length)",fontproperties=prop, labelpad=0)
    else:
        ax.set_xlabel("BLAT score", fontproperties=prop, labelpad=0)


    ax.set_ylabel("Sequence similarity (%)", fontproperties=prop, labelpad=0)
    ax.legend(loc='upper left', bbox_to_anchor=(1.28, 0.5), prop=fm.FontProperties(fname=fpath, size=11))

    ax.tick_params(length=3, width=1, axis='both', which='major', labelsize=12)
    sns.despine( ax=ax, top=True, right=True)

    return ax, kdeax





# %%
if __name__ == '__main__':

    # -----------
    # MAIN
    # -----------


    # load data
    uk_gwas_df, bv_gwas_df = load_and_format_gwas(UK_GWAS_FILE, BV_GWAS_FILE)
    og_ukbil_webblat_df, og_ukwcsf_webblat_df, og_bv_webblat_df = load_and_format_blat(WEBBLAT_UKBIL, WEBBLAT_UKWCSF, WEBBLAT_BV)


    # pool uk biobank arrays
    og_uk_webblat_df = pd.concat( [og_ukbil_webblat_df, og_ukwcsf_webblat_df], axis=0)


    # write inputs used for analysis:
    og_ukbil_webblat_df.to_csv(os.path.join(OUTPUT_TABLE_DIR, 'raw_ukbil_webblat_xy_df.tsv'), sep="\t", index=False)
    og_ukwcsf_webblat_df.to_csv(os.path.join(OUTPUT_TABLE_DIR, 'raw_ukwcsf_webblat_xy_df.tsv'), sep="\t", index=False)
    og_bv_webblat_df.to_csv(os.path.join(OUTPUT_TABLE_DIR, 'raw_bv_webblat_xy_df.tsv'), sep="\t", index=False)


    # -----------
    # BLAT Score -- no filtering on seq. length
    # -----------


    # for ukbb, pick the best blat hit using blat score, and merge with gwas data
    uk_best_by_score_no_filter_df = pick_best_blat_by_score(og_uk_webblat_df, blat_filter=False)
    uk_gwas_best_score_df = merge_blat_and_gwas(uk_best_by_score_no_filter_df, uk_gwas_df)

    # for biovu, do the same thing
    bv_best_by_score_no_filter_df = pick_best_blat_by_score(og_bv_webblat_df, blat_filter=False)
    bv_gwas_best_score_df = merge_blat_and_gwas(bv_best_by_score_no_filter_df, bv_gwas_df)

    # remove nas
    uk_gwas_best_score_no_na_df = uk_gwas_best_score_df[~uk_gwas_best_score_df.score.isna()].copy()
    bv_gwas_best_score_no_na_df = bv_gwas_best_score_df[~bv_gwas_best_score_df.score.isna()].copy()

    # for probes without a x or y hit, convert that to 0
    bv_gwas_best_score_df.score = bv_gwas_best_score_df.score.fillna(0)
    uk_gwas_best_score_df.score = uk_gwas_best_score_df.score.fillna(0)

    bv_gwas_best_score_df.query('score !=0').shape[0]
    uk_gwas_best_score_df.query('score !=0').shape[0]
    uk_gwas_best_score_df.query('score ==0').shape[0] - uk_gwas_best_score_df.shape[0]

    bv_gwas_best_score_df.shape[0]
    uk_gwas_best_score_df.shape[0]



    # %%
    ###
    ###    Plot - BLAT Score distribution
    ###


    bv_gwas_best_score_df.query('score == 0').shape[0]/bv_gwas_best_score_df.shape[0]
    uk_gwas_best_score_df.query('score == 0').shape[0]/uk_gwas_best_score_df.shape[0]


    # %%

    # set up canvas
    sns.set(style="ticks", context='talk', font_scale=0.8, rc={"figure.figsize": (6, 8)})
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=True)
    bv_axins = inset_axes(axs[0], width="73%", height="65%", loc=1)
    uk_axins = inset_axes(axs[1], width="73%", height="65%", loc=1)

    # bv
    sns.distplot(bv_gwas_best_score_df.score, hist=True, kde=False, color='blue', norm_hist=True, ax=axs[0])
    sns.distplot(bv_gwas_best_score_no_na_df.score, hist=False,  kde=True, color='gray', norm_hist=True, ax=bv_axins, kde_kws=dict(linewidth=1.5,  clip=(20,49)))
    bv_axins.set_xlim(bv_gwas_best_score_no_na_df.score.min(),bv_gwas_best_score_no_na_df.score.max())

    # uk
    sns.distplot(uk_gwas_best_score_df.score, hist=True, kde=False, color='blue', norm_hist=True, ax=axs[1])
    sns.distplot(uk_gwas_best_score_no_na_df.score, hist=False, kde=True, color='gray', norm_hist=True, ax=uk_axins, kde_kws=dict(linewidth=1.5))
    uk_axins.set_xlim(uk_gwas_best_score_no_na_df.score.min(),uk_gwas_best_score_no_na_df.score.max())


    # mark 90 and 99th percentile
    uk_90_score, uk_95_score, uk_99_score = uk_gwas_best_score_df.score.quantile([0.9,0.95, 0.99]).values
    bv_90_score, bv_95_score, bv_99_score = bv_gwas_best_score_df.score.quantile([0.9,0.95, 0.99]).values

    # bv_axins.axvline(bv_90_score, lw=1, color='indianred')
    # uk_axins.axvline(uk_90_score, lw=1, color='indianred')

    bv_axins.axvline(bv_95_score, lw=1, color='indianred')
    uk_axins.axvline(uk_95_score, lw=1, color='indianred')

    bv_axins.axvline(bv_99_score, lw=1, color='indianred')
    uk_axins.axvline(uk_99_score, lw=1, color='indianred')


    bv_axins.annotate('99th', xy=(bv_99_score, 0.12), xytext=(bv_99_score+0.3, 0.12) , fontproperties=smallprop)
    uk_axins.annotate('99th', xy=(uk_99_score, 0.09), xytext=(uk_99_score+0.3, 0.09) , fontproperties=smallprop)
    bv_axins.annotate('95th', xy=(bv_95_score, 0.12), xytext=(bv_95_score+0.3, 0.12) , fontproperties=smallprop)
    uk_axins.annotate('95th', xy=(uk_95_score, 0.09), xytext=(uk_95_score+0.3, 0.09) , fontproperties=smallprop)


    # overlap significant variants
    plot_uk_sig_vars_score = uk_gwas_best_score_no_na_df.loc[uk_gwas_best_score_no_na_df['stat_sig']==True, ['SNP','score']].values.tolist()
    plot_bv_sig_vars_score = bv_gwas_best_score_no_na_df.loc[bv_gwas_best_score_no_na_df['stat_sig']==True, ['SNP','score']].values.tolist()

    for sigsnp, sigscore in plot_uk_sig_vars_score:
        # uk_axins.axvline(sigscore, ymax=0.25, color = 'blue', lw=1)
        uk_axins.plot(sigscore, 0,fillstyle='none',markersize=15, marker="v",color='blue')
    for sigsnp, sigscore in plot_bv_sig_vars_score:
        # bv_axins.axvline(sigscore, ymax=0.25, color = 'blue', lw=1)
        bv_axins.plot(sigscore, 0,fillstyle='none',markersize=15, marker="v",color='blue')



    # modify ticks
    axs[0].tick_params(length=3, width=1)
    axs[1].tick_params(length=3, width=1)
    bv_axins.tick_params(length=2, width=1, axis='both', which='major', labelsize=10)
    uk_axins.tick_params(length=2, width=1, axis='both', which='major', labelsize=10)

    axs[0].tick_params(axis='both', which='major', labelsize=10)
    axs[1].tick_params(axis='both', which='major', labelsize=10)

    uk_axins.xaxis.set_major_locator(MultipleLocator(10))
    bv_axins.xaxis.set_major_locator(MultipleLocator(10))
    bv_axins.yaxis.set_major_locator(MaxNLocator(4))
    uk_axins.yaxis.set_major_locator(MaxNLocator(4))

    # labels
    bv_axins.set_xlabel("")
    uk_axins.set_xlabel("")
    axs[0].set_xlabel("")

    bv_axins.set_xticklabels([int(x) for x in bv_axins.get_xticks()], fontproperties=fm.FontProperties(fname=fpath, size=12))
    uk_axins.set_xticklabels([int(x) for x in uk_axins.get_xticks()] , fontproperties=fm.FontProperties(fname=fpath, size=12))
    bv_axins.set_yticklabels(["{:.2f}".format(x) for x in bv_axins.get_yticks()], fontproperties=fm.FontProperties(fname=fpath, size=12))
    uk_axins.set_yticklabels(["{:.2f}".format(x) for x in uk_axins.get_yticks()] , fontproperties=fm.FontProperties(fname=fpath, size=12))

    axs[0].set_yticklabels(["{:.2f}".format(x) for x in axs[0].get_yticks()], fontproperties=fm.FontProperties(fname=fpath, size=12))
    axs[1].set_yticklabels(["{:.2f}".format(x) for x in axs[1].get_yticks()] , fontproperties=fm.FontProperties(fname=fpath, size=12))
    axs[0].set_xticklabels([int(x) for x in axs[0].get_xticks()], fontproperties=fm.FontProperties(fname=fpath, size=12))
    axs[1].set_xticklabels([int(x) for x in axs[1].get_xticks()] , fontproperties=fm.FontProperties(fname=fpath, size=12))

    fig.subplots_adjust(hspace=0.35)
    axs[0].set_title('BioVU',fontproperties=prop)
    axs[1].set_title('UK Biobank',fontproperties=prop)
    axs[1].set_xlabel('BLAT Score', fontproperties=prop, labelpad=10)

    # spines
    sns.despine(ax=axs[0], top=True, right=True)
    sns.despine(ax=axs[1], top=True, right=True)
    sns.despine(ax=bv_axins, top=True, right=True)
    sns.despine(ax=uk_axins, top=True, right=True)

    for spine_pos in ['bottom','left']:
        axs[0].spines[spine_pos].set_linewidth(0.5)
        axs[1].spines[spine_pos].set_linewidth(0.5)
        uk_axins.spines[spine_pos].set_linewidth(1)
        bv_axins.spines[spine_pos].set_linewidth(1)

    # plt.show()
    plt.savefig(os.path.join(OUTPUT_DIR, '{}_blat_score_w_inset.pdf'.format(DATE)))



    # %%
    # -----------
    # Filter on sequence length and overlap -- plot % seq id
    # -----------


    # filter on match length but not on sequence identity
    bv_filtered_blat_df, _ =  bv_blat_filter(og_bv_webblat_df, filter_perID=False)

    ukbil_filtered_blat_df, _ = uk_array_filter(og_ukbil_webblat_df, filter_perID=False)
    ukwcsf_filtered_blat_df, _ = uk_array_filter(og_ukwcsf_webblat_df, filter_perID=False)
    uk_filtered_blat_df = pd.concat([ukbil_filtered_blat_df,ukwcsf_filtered_blat_df], axis=0)



    # pick the best blat hit using blat score, and merge with gwas data
    # uk
    uk_best_by_score_filter_df = pick_best_blat_by_score(uk_filtered_blat_df, blat_filter=True)
    uk_gwas_best_score_filtered_df = merge_blat_and_gwas(uk_best_by_score_filter_df, uk_gwas_df)

    # biovu
    bv_best_by_score_filter_df = pick_best_blat_by_score(bv_filtered_blat_df, blat_filter=True)
    bv_gwas_best_score_filtered_df = merge_blat_and_gwas(bv_best_by_score_filter_df, bv_gwas_df)

    # remove nas
    uk_gwas_best_score_filtered_no_na_df = uk_gwas_best_score_filtered_df[~uk_gwas_best_score_filtered_df.score.isna()].copy()
    bv_gwas_best_score_filtered_no_na_df = bv_gwas_best_score_filtered_df[~bv_gwas_best_score_filtered_df.score.isna()].copy()


    assert np.all(uk_gwas_best_score_filtered_df.groupby('SNP').size() == 1 ), "more than one row per SNP"


    uk_xy_hit_prop = 1 - uk_gwas_best_score_filtered_df.score.isna().sum()/uk_gwas_best_score_filtered_df.SNP.nunique()
    bv_xy_hit_prop = 1 - bv_gwas_best_score_filtered_df.score.isna().sum()/bv_gwas_best_score_filtered_df.SNP.nunique()

    # now filter on %sequence similarity ≥ 90%
    # proportion of probes with ≥90% seq simliarity and length restrictions out of all probes...
    uk_high_qual_hits_df = uk_gwas_best_score_filtered_df.loc[ (~uk_gwas_best_score_filtered_df['score'].isna()) & (uk_gwas_best_score_filtered_df['per_identity'] >= 90)].copy()
    uk_high_qual_hits = uk_gwas_best_score_filtered_df.loc[ (~uk_gwas_best_score_filtered_df['score'].isna()) & (uk_gwas_best_score_filtered_df['per_identity'] >= 90)].SNP.nunique()
    uk_high_qual_hits/uk_gwas_best_score_filtered_df.SNP.nunique() * 100

    bv_high_qual_hits_df = bv_gwas_best_score_filtered_df.loc[ (~bv_gwas_best_score_filtered_df['score'].isna()) & (bv_gwas_best_score_filtered_df['per_identity'] >= 90)].copy()
    bv_high_qual_hits = bv_gwas_best_score_filtered_df.loc[ (~bv_gwas_best_score_filtered_df['score'].isna()) & (bv_gwas_best_score_filtered_df['per_identity'] >= 90)].SNP.nunique()
    bv_high_qual_hits/bv_gwas_best_score_filtered_df.SNP.nunique() * 100


    # %%
    ###
    ###    write these dfs
    ###

    #
    # ---- write best blat XY hit based on blat score for GWAS significant variants ----
    #

    keep_cols = ['SNP','chr_pos','A1','OR','P','targetChr','target_Start-End','score','per_identity','q.length']
    uk_best_blat_score_sig_hits_table = uk_gwas_best_score_filtered_df.query('stat_sig == True').loc[:, keep_cols].reset_index(drop=True).copy()
    bv_best_blat_score_sig_hits_table = bv_gwas_best_score_filtered_df.query('stat_sig == True').loc[:, keep_cols].reset_index(drop=True).copy()


    rename_bv_snps = {'JHU_3.16652239' :'rs9870157',
    'JHU_13.20119335' :'rs9508454',
    '14:35761675-C-G' :'rs1048990'}

    bv_best_blat_score_sig_hits_table.loc[bv_best_blat_score_sig_hits_table.SNP.isin(rename_bv_snps.keys()), 'SNP'] = bv_best_blat_score_sig_hits_table.loc[bv_best_blat_score_sig_hits_table.SNP.isin(rename_bv_snps.keys()), 'SNP'].map(rename_bv_snps).values

    bv_best_blat_score_sig_hits_table['dataset'] = 'BioVU'
    uk_best_blat_score_sig_hits_table['dataset'] = 'UKBiobank'


    combined_df = pd.concat([bv_best_blat_score_sig_hits_table, uk_best_blat_score_sig_hits_table])
    combined_df.to_csv(os.path.join(OUTPUT_TABLE_DIR, 'uk_bv_gwas_sig_hits_w_best_blat_score_xy_match.tsv'), sep="\t", index=False)

    #
    # ---- write for GWAS varaints with an X or Y BLAT hit that passes the filtering, write the best BLAT hit based on score ----
    #

    uk_high_qual_hits_df.to_csv(os.path.join(OUTPUT_TABLE_DIR, 'uk_best_blatscore_xy_hit_and_filtered.tsv'), sep="\t", index=False)
    bv_high_qual_hits_df.to_csv(os.path.join(OUTPUT_TABLE_DIR, 'bv_best_blatscore_xy_hit_and_filtered.tsv'), sep="\t", index=False)


    # %%
    ###
    ###    Plot - % Sequence ID historgram
    ###


    # num probes w/ ≥90% seq similarity out of all autosomal probes
    bv_gwas_best_score_filtered_no_na_df.query("per_identity >=90").shape[0]/bv_gwas_best_score_df.shape[0]*100
    uk_gwas_best_score_filtered_no_na_df.query("per_identity >=90").shape[0]/uk_gwas_best_score_df.shape[0]*100

    bv_gwas_best_score_filtered_no_na_df.query("per_identity >=90").shape[0]
    uk_gwas_best_score_filtered_no_na_df.query("per_identity >=90").shape[0]

    # %%
    sns.set(style="ticks", context='talk', font_scale=0.8, rc={"figure.figsize": (6, 8)})
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=False)
    sns.distplot(bv_gwas_best_score_filtered_no_na_df.per_identity, ax=axs[0], hist=True, kde=True, color='gray', norm_hist=True, kde_kws=dict(linewidth=1.5, bw=4))
    sns.distplot(uk_gwas_best_score_filtered_no_na_df.per_identity, ax=axs[1], hist=True, kde=True, color='gray', norm_hist=True, kde_kws=dict(linewidth=1.5, bw=4))

    axs[0].axvline(90, c='indianred', lw=1)
    axs[1].axvline(90, c='indianred', lw=1)

    sns.despine(ax=axs[0], top=True, right=True)
    sns.despine(ax=axs[1], top=True, right=True)

    for spine_pos in ['bottom','left']:
        axs[0].spines[spine_pos].set_linewidth(1)
        axs[1].spines[spine_pos].set_linewidth(1)

    axs[0].tick_params(length=3, width=1, axis='both', which='major', labelsize=10)
    axs[1].tick_params(length=3, width=1, axis='both', which='major', labelsize=10)

    axs[0].set_xlabel('', fontproperties=prop, labelpad=10)
    axs[1].set_xlabel('Sequence identity(%)', fontproperties=prop, labelpad=10)

    axs[0].set_yticklabels(["{:.2f}".format(x) for x in axs[0].get_yticks()], fontproperties=fm.FontProperties(fname=fpath, size=12))
    axs[1].set_yticklabels(["{:.2f}".format(x) for x in axs[1].get_yticks()] , fontproperties=fm.FontProperties(fname=fpath, size=12))
    axs[0].set_xticklabels([int(x) for x in axs[0].get_xticks()], fontproperties=fm.FontProperties(fname=fpath, size=12))
    axs[1].set_xticklabels([int(x) for x in axs[1].get_xticks()] , fontproperties=fm.FontProperties(fname=fpath, size=12))

    axs[0].set_xlim(right=100)
    axs[1].set_xlim(right=100)

    plot_uk_sig_vars_perid = uk_gwas_best_score_no_na_df.loc[uk_gwas_best_score_no_na_df['stat_sig']==True, ['SNP','per_identity']].values.tolist()
    plot_bv_sig_vars_perid = bv_gwas_best_score_no_na_df.loc[bv_gwas_best_score_no_na_df['stat_sig']==True, ['SNP','per_identity']].values.tolist()

    for sigsnp, sigscore in plot_uk_sig_vars_perid:
        # uk_axins.axvline(sigscore, ymax=0.25, color = 'blue', lw=1)
        axs[1].plot(sigscore, 0.004,fillstyle='none',markersize=15, marker="v",color='blue')
    for sigsnp, sigscore in plot_bv_sig_vars_perid:
        # bv_axins.axvline(sigscore, ymax=0.25, color = 'blue', lw=1)
        axs[0].plot(sigscore, 0.004,fillstyle='none',markersize=15, marker="v",color='blue')


    fig.subplots_adjust(hspace=0.35)
    axs[0].set_title('BioVU',fontproperties=prop)
    axs[1].set_title('UK Biobank',fontproperties=prop)

    # plt.savefig(os.path.join(OUTPUT_DIR, '{}_seq_identity.pdf'.format(DATE)))
    # plt.show()


    # %%
    #rename
    rename_bv_snps = {'JHU_3.16652239' :'rs9870157',
    'JHU_13.20119335' :'rs9508454',
    '14:35761675-C-G' :'rs1048990'}



    # %%
    ###
    ###    Bivariate Plots of Score vs. Sequence ID
    ###

    prop_probe_length = lambda df,probe_length: df['q.length']/probe_length
    filt_perid = lambda df: df.loc[df['per_identity']>= 90].copy()

    bv_filt_df = filt_perid(bv_gwas_best_score_filtered_no_na_df)
    uk_filt_df = filt_perid(uk_gwas_best_score_filtered_no_na_df)
    bv_filt_df.loc[bv_filt_df['SNP'].isin(rename_bv_snps.keys()), 'SNP'] = bv_filt_df.loc[bv_filt_df['SNP'].isin(rename_bv_snps.keys())].SNP.map(rename_bv_snps).values

    bv_filt_df['prop_of_probe_length'] = prop_probe_length(bv_filt_df, 50)
    uk_filt_df['prop_of_probe_length'] = prop_probe_length(uk_filt_df, 71)



    # %% prop length vs. % id
    sns.set(context='paper', font_scale=1.5, style='ticks', rc={'figure.figsize':(15,5)})
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)

    bvmanyfig, bvkdeax = plot_bivar_many___vs_length(bv_filt_df, 'BioVU', axs[0], 'prop_of_probe_length')
    ukmanyfig, ukkdeax = plot_bivar_many___vs_length(uk_filt_df, 'UK Biobank', axs[1], 'prop_of_probe_length')


    plt.tight_layout()
    savefile = os.path.join(OUTPUT_DIR, f'{DATE}_uk_bv_bivar_prop_length_vs_seqidentity.pdf')
    plt.savefig(savefile)


    # %%

    sns.set(context='paper', font_scale=1.5, style='ticks', rc={'figure.figsize':(15,5)})
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)

    ukmanyfig, ukkdeax = plot_bivar_many___vs_length(uk_filt_df, 'UK Biobank', axs[1], 'score')
    bvmanyfig, bvkdeax = plot_bivar_many___vs_length(bv_filt_df, 'BioVU', axs[0], 'score')


    plt.tight_layout()


    savefile = os.path.join(OUTPUT_DIR, f'{DATE}_uk_bv_bivar_score_vs_seqidentity.pdf')
    plt.savefig(savefile)


# %%
