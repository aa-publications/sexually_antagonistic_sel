#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2019-07-24 11:10:22

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')

###
#   PATHS
###

# BLAT OUTPUT MUST BE IN 'BLAST OUTPUT FORMAT #8'
BLAT_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/blat_results_2019-07-24.txt"
PLINK_ASSOC_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/key_for_seqList_for_twoBitToFa_2019-07-24.tsv"

OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/2019_07_21_logistic/katja_results/blat_results_w_assoc.tsv"

###
#   MAIN
###
pl_output=pd.read_csv(PLINK_ASSOC_OUTPUT, sep="\s+")
bl_output=pd.read_csv(BLAT_OUTPUT, sep="\t",header=None, names=['Query', 'Subject', '%id', 'alignment_length',
                                                                 'mistmatches', 'gap_openings', 'q.start', 'q.end', 's.start',
                                                                  's.end', 'e-value', 'bit score'])


keep_pl_output = pl_output.loc[:, ['CHR', 'SNP', 'BP',  'OR', 'P', 'seqList']].copy()

# merge
final_df = pd.merge(pl_output, bl_output, left_on='seqList', right_on='Query', how='inner')


# write
final_df.to_csv(OUTPUT_FILE, index=False, sep="\t", header=True)
