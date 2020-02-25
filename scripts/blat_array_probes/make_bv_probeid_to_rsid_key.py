#!/bin/python
# This script will 
#
#
#
# Abin Abraham
# created on: 2020-01-25 14:13:22

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime


from glob import glob
DATE = datetime.now().strftime('%Y-%m-%d')


# PATHS
bv_bim_file= "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3.bim"
bv_annot_file="/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/MEGAEx_BioVU_15075710_A1_no_header.csv"


KEY_FILE_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/biovu.annot.key.to.bim.tsv"

# -----------
# FUNCTIONS
# -----------


def intersect_uk(uk_bim, anno_df, uk_col, anno_col): 
        

    n_all = uk_bim[uk_col].nunique()
    both = set(uk_bim[uk_col].unique()).intersection(set(anno_df[anno_col].unique()))
    bim_only = set(uk_bim[uk_col].unique()).difference(set(anno_df[anno_col].unique()))
    anno_only = set(anno_df[anno_col].unique()).difference(set(uk_bim[uk_col].unique()))

    print("{:,} out of {:,} SNPs mapped by {}.".format(len(both), n_all, anno_col))
    
    return both, bim_only, anno_only


def intersect(a, b): 
    
    both = a.intersection(b)
    a_only = a.difference(b)
    b_only = b.difference(a)
    
    print("{:,} shared, {:,} only in first set,  {:,} only in second set".format(len(both), len(a_only), len(b_only)))
    
    return both, a_only, b_only

# %%

# -----------
# MAIN
# -----------

# load annotation file 
anno_df = pd.read_csv(bv_annot_file, low_memory=False)
keep_cols = ['IlmnID', 'Name', 'IlmnStrand', 'SNP', 'AlleleA_ProbeSeq', 'AlleleB_ProbeSeq', 'Chr', 'MapInfo']
anno_df = anno_df.loc[~anno_df.MapInfo.isna(), keep_cols].copy()
anno_df.MapInfo = anno_df.MapInfo.map(int)

# load bim file 
bv_bim_df = pd.read_csv(bv_bim_file, sep="\t", names=['chr','SNP','dummy','pos','A1','A2'])



# %%
###
#   COMPARE rsIDS
###

anno_col='Name'
rsid_both, rsid_bim_only, rsid_anno_only = intersect_uk(bv_bim_df, anno_df, 'SNP', anno_col)
# NAME columns captures SNP in BIM compleemte 



# -----------
# write a key 
# -----------

# for each probe in the annotation file, provide the rsID in the bim file 
anno_df.head(2)


key_anno_df = anno_df.copy()
key_anno_df['bim_SNP'] = "None"
key_anno_df['bim_SNP'] = key_anno_df['Name']
key_anno_df.to_csv(KEY_FILE_OUTPUT, sep="\t", index=False)

