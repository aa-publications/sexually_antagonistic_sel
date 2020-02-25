
#!/bin/python
# This script will compare e-values distributions from BLAT on variants w/ a certain flankign region (from GWAS assoc resutls provided by Katja.)
#
#
#
# Abin Abraham
# created on: 2019-08-22 13:11:52

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime


import glob
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.graphics.mosaicplot import mosaic
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("/dors/capra_lab/users/abraha1/bin/python_modules/assocplots")
from assocplots.manhattan import *

DATE = datetime.now().strftime('%Y-%m-%d')



# %%
# -----------  PATHS  -----------

KATJA_ASSOC_RESULTS="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL.csv"

BLAT_250_OUTPUTS = '/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_08_21_blat_pipeline/all_gwas_hits/250bp/blat_results_w_assoc_2019_08_24.tsv'
BLAT_150_OUTPUTS = '/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_08_21_blat_pipeline/all_gwas_hits/150bp/blat_results_w_assoc_2019_08_25.tsv'
BLAT_50_OUTPUTS = '/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_08_21_blat_pipeline/all_gwas_hits/50bp/blat_results_w_assoc_2019_08_24.tsv'

OUTPUT_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_08_21_blat_pipeline/all_gwas_hits/"

# -----------  CONSTANTS  -----------

GWAS_P_THRESH = 5*10**-8
MISS_P_THRESHOLD = 0.00001

# -----------
# FUNCTIONS
# -----------


def plot_e_distribution(blat_output_file, title="", ax=None, keep_min_evalue=False):

    blat_df = pd.read_csv(blat_output_file, sep="\t")

    # 1) remove hits with high missing
    # 2) keep only variants with BLAT hits to X or Y
    # assuming 'p_missing' refers to high missing rate between cases and controls
    print("*** FILTER OUT MISSING VALUES ***")
    xy_blat_df = blat_df.loc[ ( (blat_df['missing_p'] > 0.00001) & ((blat_df['Subject'] == "chrX") | (blat_df['Subject'] == "chrY"))) , :].copy()
    xy_blat_df['gwas_signif'] = xy_blat_df.P.apply(lambda x: 'signif' if x < GWAS_P_THRESH else 'no_signif')

    plot_df = xy_blat_df.loc[:, ['SNP','P','%id','e-value','gwas_signif']].copy()
    plot_df['log_evalue'] = plot_df['e-value'].apply(lambda x: np.log10(x))

    if keep_min_evalue:
        plot_df.sort_values(['SNP','e-value'], inplace=True, ascending=True)
        plot_df = plot_df[plot_df.duplicated('SNP', keep='first')].copy()
        print("keeping only the best hit by min. e-value")


    # 3) Plot e-value distributions

    if ax is None:
        fig, ax = plt.subplots()

    sns.set(style="whitegrid",  font_scale=1.5, rc={"figure.figsize": (12, 4)})
    sns.violinplot(data=plot_df, y='gwas_signif', x='log_evalue', scale='width',inner='quartile',  ax=ax)
    sns.stripplot(data=plot_df, y='gwas_signif', x='log_evalue', color='black', edgecolor="black", ax=ax, size=3, alpha=0.8,jitter=False)
    ax.axvline(-50, color='r')
    ax.axvline(-2, color='b')
    plt.title(title)
    plt.xlabel('log10(E-value)')

    return ax



def plot_e_distribution_by_x_or_y(blat_output_file, title="", ax=None, keep_min_evalue=False):

    blat_df = pd.read_csv(blat_output_file, sep="\t")

    # 1) remove hits with high missing
    # 2) keep only variants with BLAT hits to X or Y
    # assuming 'p_missing' refers to high missing rate between cases and controls
    print("*** FILTER OUT MISSING VALUES ***")
    xy_blat_df = blat_df.loc[ ( (blat_df['missing_p'] > 0.00001) & ((blat_df['Subject'] == "chrX") | (blat_df['Subject'] == "chrY"))) , :].copy()
    xy_blat_df['gwas_signif'] = xy_blat_df.P.apply(lambda x: 'signif' if x < GWAS_P_THRESH else 'no_signif')


    plot_df = xy_blat_df.loc[:, ['SNP','P','%id','e-value','gwas_signif','Subject']].copy()
    plot_df['log_evalue'] = plot_df['e-value'].apply(lambda x: np.log10(x))

    if keep_min_evalue:
        plot_df.sort_values(['SNP','e-value'], inplace=True, ascending=True)
        plot_df = plot_df[plot_df.duplicated('SNP', keep='first')].copy()
        print("keeping only the best hit by min. e-value")


    # 3) Plot e-value distributions

    if ax is None:
        fig, ax = plt.subplots()

    sns.set(style="whitegrid",  font_scale=1.5, rc={"figure.figsize": (12, 4)})
    sns.violinplot(data=plot_df, y='gwas_signif', x='log_evalue', scale='width',inner='quartile', hue="Subject", ax=ax)
    sp = sns.stripplot(data=plot_df, y='gwas_signif', x='log_evalue', color='black', edgecolor="black",ax=ax,
                        hue="Subject", size=3, alpha=0.8,jitter=False, dodge=True)

    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


    ax.axvline(-50, color='r')
    ax.axvline(-2, color='b')
    plt.title(title)
    plt.xlabel('log10(E-value)')

    return ax




def plot_perID_distribution(blat_output_file, title="", ax=None):

    blat_df = pd.read_csv(blat_output_file, sep="\t")

    # 1) remove hits with high missing
    # 2) keep only variants with BLAT hits to X or Y
    # assuming 'p_missing' refers to high missing rate between cases and controls

    xy_blat_df = blat_df.loc[ ( (blat_df['missing_p'] > MISS_P_THRESHOLD) & ((blat_df['Subject'] == "chrX") | (blat_df['Subject'] == "chrY"))) , :].copy()
    xy_blat_df['gwas_signif'] = xy_blat_df.P.apply(lambda x: 'signif' if x < GWAS_P_THRESH else 'no_signif')

    plot_df = xy_blat_df.loc[:, ['SNP','P','%id','e-value','gwas_signif']].copy()

    # 3) Plot e-value distributions


    if ax is None:
        fig, ax = plt.subplots()

    sns.set(style="whitegrid",  font_scale=1.5, rc={"figure.figsize": (12, 4)})
    sns.violinplot(data=plot_df, y='gwas_signif', x='%id', scale='width',inner='quartile',  ax=ax)
    sns.stripplot(data=plot_df, y='gwas_signif', x='%id', color='black', edgecolor="black", ax=ax, size=3, alpha=0.8,jitter=False)
    plt.title(title)
    plt.xlabel('%id')

    return ax

def calc_or(BLAT_OUTPUT_FILE, GWAS_P_THRESH, BLAT_P_THRESH, label):

    blat_df = pd.read_csv(BLAT_OUTPUT_FILE, sep="\t")

    print("*** FILTER OUT MISSING VALUES ***")
    print("*** KEEPING BLAT HITS w/ E-VALUE <= {}".format(BLAT_P_THRESH))
    no_miss_blat_df = blat_df.loc[ (blat_df['missing_p'] > MISS_P_THRESHOLD) & (blat_df['e-value'] <= BLAT_P_THRESH), :].copy()

    no_miss_blat_df['gwas_signif'] = no_miss_blat_df.P.apply(lambda x: 'gwas_signif' if x < GWAS_P_THRESH else 'gwas_no_signif')
    no_miss_blat_df['blat_XY_hit'] = no_miss_blat_df['Subject'].apply(lambda x: 'blat_XY' if ( (x == 'chrX') | (x == 'chrY'))  else 'blat_no_XY')

    contig_table = sm.stats.Table(pd.crosstab(no_miss_blat_df.gwas_signif, no_miss_blat_df.blat_XY_hit).loc[['gwas_signif','gwas_no_signif'], ['blat_XY', 'blat_no_XY']])
    counts = contig_table.table_orig
    odds_ratio, pval = stats.fisher_exact(counts)


    return pd.DataFrame({'label':[label] , 'gwas_p_thresh':[GWAS_P_THRESH], 'blat_p_thresh': [BLAT_P_THRESH],
                  'gwas_signif_blat_XY_signif'   : [counts.loc['gwas_signif', 'blat_XY']],
                  'gwas_signif_blat_XY_no_signif': [counts.loc['gwas_signif', 'blat_no_XY']],
                  'gwas_no_signif_blat_XY_signif': [counts.loc['gwas_no_signif', 'blat_XY']],
                  'gwas_no_signif_blat_XY_no_signif': [counts.loc['gwas_no_signif', 'blat_no_XY']],
                  'odds_ratio':[odds_ratio], 'pval':[pval]})




# %%
# -----------
# MAIN
# -----------

# -----------  Load association stats  -----------

# load gwas assoc results
assoc_df = pd.read_csv(KATJA_ASSOC_RESULTS, sep=",")

# remove hits with high missingness b/w case and controls
clean_assoc_df = assoc_df.loc[assoc_df['missing_p'] > 0.00001].copy()
print("{:,} variants removed due to missingnesss..".format(assoc_df.shape[0] - clean_assoc_df.shape[0]))



# %ll
###
#   Manhattan Plots of GWAS Summary
###

mpl.rcParams['figure.figsize']=12, 6
manhattan(assoc_df.P , assoc_df.BP, assoc_df.CHR, '', title="katja_assoc_pca12_centers_age",
               lines= [6, 8],
               lines_colors=['b', 'r'])

# remove hits with significant missing rate between cases and controls

manhattan(clean_assoc_df.P , clean_assoc_df.BP, clean_assoc_df.CHR, '', title="katja_assoc_pca12_centers_age\n(rm variatns w/ high missing rate b/w cases and controls)",
           lines= [6, 8],
           lines_colors=['b', 'r'],)


# %%
###
#   Calc OR for BLAT HIT TO SEX CHROMOSOMES to GWAS SIGNIFICANCE
###

or_df  = pd.DataFrame()
for label, blat_file in {'50bp': BLAT_50_OUTPUTS, '150bp':BLAT_150_OUTPUTS, '250bp':BLAT_250_OUTPUTS}.items():

    for BLAT_P_THRESH in [0.01]:
        df  = calc_or(blat_file, GWAS_P_THRESH, BLAT_P_THRESH, label)
        or_df = or_df.append(df)


or_df




# %%
###
#   BLAT results - Distribution of  E-Value between signif and non- signic
###

#
# PLOT X AND Y
#

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12,8), sharex=True)
plot_e_distribution_by_x_or_y(BLAT_50_OUTPUTS, title="",ax=axs[0])
plot_e_distribution_by_x_or_y(BLAT_150_OUTPUTS, title="",ax=axs[1])
plot_e_distribution_by_x_or_y(BLAT_250_OUTPUTS, title="",ax=axs[2])

axs[0].set_ylabel("50 BP")
axs[1].set_ylabel("250 BP")
axs[2].set_ylabel("500 BP")

axs[0].set_xlabel("")
axs[1].set_xlabel("")

#
# PLOT X OR Y vs. none
#



fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12,8), sharex=True)
plot_e_distribution(BLAT_50_OUTPUTS, title="",ax=axs[0])
plot_e_distribution(BLAT_150_OUTPUTS, title="", ax=axs[1])
plot_e_distribution(BLAT_250_OUTPUTS, title="", ax=axs[2])

axs[0].set_ylabel("50 BP")
axs[1].set_ylabel("250 BP")
axs[2].set_ylabel("500 BP")

axs[0].set_xlabel("")
axs[1].set_xlabel("")
plt.savefig(os.path.join(OUTPUT_DIR, 'evalue_distribution_by_gwas_significance_{}.png'.format(DATE)))

# %%
#
# PLOT MIN E-VALUE
#


fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12,8), sharex=True)
plot_e_distribution(BLAT_50_OUTPUTS, title="", ax=axs[1], keep_min_evalue=True)
plot_e_distribution(BLAT_150_OUTPUTS, title="", ax=axs[1], keep_min_evalue=True)
plot_e_distribution(BLAT_250_OUTPUTS, title="", ax=axs[2], keep_min_evalue=True)


axs[0].set_ylabel("50 BP")
axs[1].set_ylabel("250 BP")
axs[2].set_ylabel("500 BP")

axs[0].set_xlabel("")
axs[1].set_xlabel("")
plt.savefig(os.path.join(OUTPUT_DIR, 'best_hit_by_min_evalue_distribution_by_gwas_significance_{}.png'.format(DATE)))


# %%
###
#   BLAT results - %id
###

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12,8), sharex=True)
plot_perID_distribution(BLAT_50_OUTPUTS, title="", ax=axs[0])
plot_perID_distribution(BLAT_150_OUTPUTS, title="", ax=axs[1])
plot_perID_distribution(BLAT_250_OUTPUTS, title="", ax=axs[2])

axs[0].set_ylabel("50 BP")
axs[1].set_ylabel("250 BP")
axs[2].set_ylabel("500 BP")

axs[0].set_xlabel("")
axs[1].set_xlabel("")

plt.savefig(os.path.join(OUTPUT_DIR, 'percent_id_distribution_by_gwas_significance_{}.png'.format(DATE)))
#



###
#   plot only top hit
###

hits_df = pd.read_csv(BLAT_50_OUTPUTS, sep="\t")
keep_cols = ['CHR','SNP','BP', 'A1','OR', 'STAT','P','missing_p', 'Query','Subject','q.start','q.end','e-value','bit score']
hits_df = hits_df.loc[:, keep_cols].copy()

hits_df.groupby('SNP').apply(lambda x: x[x['e-value'] == x['e-value'].min()])

hits_df.head()


hits_df.columns