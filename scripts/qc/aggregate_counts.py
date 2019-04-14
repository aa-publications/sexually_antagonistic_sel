#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2019-04-13 14:16:10   

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
 
DATE = datetime.now().strftime('%Y-%m-%d')


import re
import glob


#PATH 

DATA_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
OUTPUT_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/mega_data/post_qc_per_batch"
snps_summary_output = os.path.join(OUTPUT_DIR, 'snps_qc_count.tsv')
individual_summary_output = os.path.join(OUTPUT_DIR, 'individual_qc_count.tsv')
# -----------
# FUNCTION
# ----------- 


def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)



def format_count(df):

    df['batch_num'] = df.batch.apply(lambda x: x.split("MEGA_ex_Array_")[1].split("_Cox")[0] )
    df.drop(['batch'], axis=1, inplace=True)
    df.rename(columns={'batch_num':'batch'}, inplace=True)

    df = df.T
    df.columns = df.loc['batch',:]
    df.drop(['batch'], axis=0, inplace=True)

    final_df = df.loc[:, sorted_nicely(df.columns)].copy()
    final_df.columns.name=''
    final_df.index.name='qc_steps'

    diff_df = final_df.apply(pd.DataFrame.diff)
    diff_df.fillna(0)
    for index, col in final_df.iteritems(): 
        for row_index, row in col.iteritems(): 
            
            final_df.loc[row_index, index] = "{:,}({:,})".format(final_df.loc[row_index, index], diff_df.loc[row_index, index])

    return final_df

# =============  MAIN =============

snp_df = pd.DataFrame()
for this_batch in glob.glob(DATA_DIR+"/*/counts_snps*"):
    df = pd.read_csv(this_batch, sep="\t")
    snp_df =snp_df.append(df)


individ_df = pd.DataFrame()
for this_batch in glob.glob(DATA_DIR+"/*/counts_individuals*"):
    df = pd.read_csv(this_batch, sep="\t")
    individ_df =individ_df.append(df)


final_snp_df = format_count(snp_df)
final_individ_df = format_count(individ_df)


final_snp_df.to_csv(snps_summary_output, sep="\t",index=True)
final_individ_df.to_csv(individual_summary_output, sep="\t",index=True)

