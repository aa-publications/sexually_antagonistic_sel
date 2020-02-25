#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-01-23 09:44:14



import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')



# PATHS

ukbb_bil_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/UKBiLEVE_probes_blat.txt"
ukbb_wcsf_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/UKBB_WCSF_probes_blat.txt"
# mega_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/check_probe_seqs/MEGAEx_probes_v1_blat.txt"

# 
output_ukbb_bil_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/chrx_y_hits_only_UKBiLEVE_probes_blat.txt"
output_ukbb_wcsf_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/chrx_y_hits_only_UKBB_WCSF_probes_blat.txt"
# output_mega_file="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/check_probe_seqs/chrx_y_hits_only_MEGAEx_probes_v1_blat.txt"
# 
# 


# MAIN 
col_names= ["QueryId","SubjectId","per_identity","alignment_length","mismatches","gap_openings","q.start","q.end","s.start","s.end","e-value","bit_score"]

bil_df = pd.read_csv(ukbb_bil_file, sep="\t", names=col_names)
ukbb_df = pd.read_csv(ukbb_wcsf_file, sep="\t", names=col_names)
# mega_df = pd.read_csv(mega_file, sep="\t", names=col_names)


# keep only hits to chrX and chrY
keep_subjects = ['chrX','chrY']
filtered_bil_df = bil_df.loc[bil_df['SubjectId'].isin(keep_subjects) ].copy()
filtered_ukbb_df = ukbb_df.loc[ukbb_df['SubjectId'].isin(keep_subjects) ].copy()
# filtered_mega_df = mega_df.loc[mega_df['SubjectId'].isin(keep_subjects) ].copy()


# write 
filtered_bil_df.to_csv(output_ukbb_bil_file, sep="\t", index=False)
filtered_ukbb_df.to_csv(output_ukbb_wcsf_file, sep="\t", index=False)
# filtered_mega_df.to_csv(output_mega_file, sep="\t", index=False)
