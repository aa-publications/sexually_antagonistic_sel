#!/bin/python
# See function doc string...
#
#
#
# Abin Abraham
# created on: 2019-04-07 13:01:13

import os
import sys
import subprocess
import pandas as pd 

sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd

@track_fx
def rm_dup_individuals(input_prefix, output_dir, base_prefix):
    """
    Plink wrapper to remove duplicate individuals.
        - The first instance is retained, all other are removed. 
        - See OUTPUT FILES for all outputs files created

    WARNING: 
        - The original .fam file will be overwritten with dupllicated. 
        - A copy of the original is made before overwritting.
        - If there are no duplicates, then --make-bed is still run (to fit with the pipeline)

    Parameters
    ----------
    input_prefix : str
        full path with plink prefix of file 
    ouput_dir : str
        full path to directory to write outputs
    base_prefix : str
        the original plink prefix to be modified for output plink prefix 

    Returns
    -------
    dups_removed_plink_prefix : str
        - full path with plink prefix with duplicate samples removed
    plink_stdout : str 
        - STDOUT from running plink command

    
    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============

    copy_of_orig_fam_file = os.path.join(output_dir, "original_"+base_prefix+".fam")
    duplicated_samples_file = os.path.join(output_dir, 'tmp_dups_samples_to_rm{}.csv'.format(base_prefix))
    no_dups_plink_prefix = os.path.join(output_dir, "tmp_dedups_{}".format(base_prefix))


    # =============  REMOVE DUPLICATE SAMPLES =============
    
    # read and write a copy of fam file 
    fam_df = pd.read_csv(input_prefix+".fam", sep="\s+", names=['FID','IID','c3','c4','c5','c6'])
    fam_df.to_csv(copy_of_orig_fam_file, sep=" ", header=None, index=None)

    
    assert fam_df[~(fam_df.FID == fam_df.IID)].shape[0] == 0,\
         "FID and IID are *not* the same in this file:\n{}".format(input_prefix+".fam")

    # identify duplicated FID&IID
    dup_index = fam_df[fam_df.duplicated(subset=['FID','IID'], keep='first')].index
    dup_fids = fam_df.iloc[dup_index, :].FID.unique()

    # each duplicate FID & IID, except for the first instance, will be have "_[counter]" appened
    for this_fid in dup_fids: 
        for counter, index_row in enumerate(fam_df.loc[fam_df['FID']== this_fid].iterrows()):
            index, this_row = index_row
            if counter == 0: 
                continue 
            else: 
                fam_df.loc[index, ['FID', 'IID']] = fam_df.loc[index, ['FID', 'IID']].apply(lambda x: x+"_{}".format(counter))

    # write duplicated FID and IID to file 
    if (fam_df.loc[dup_index, ['FID','IID']].shape[0] > 0): 
        fam_df.loc[dup_index, ['FID','IID']].to_csv(duplicated_samples_file, sep=" ", header=None, index=None)

    # OVERWRITE existing .fam to tagging duplicates
    fam_df.to_csv(input_prefix+".fam", sep=" ", header=None, index=None)

    # plink to rm duplicates
    if (fam_df.loc[dup_index, ['FID','IID']].shape[0] > 0):
        rm_dups_cmd = "plink --bfile {} --remove {} --make-bed --out {}".format( input_prefix, duplicated_samples_file, no_dups_plink_prefix)
    else: 
        rm_dups_cmd = "plink --bfile {} --make-bed --out {}".format( input_prefix, no_dups_plink_prefix)

    plink_stdout = run_shell_cmd(rm_dups_cmd)

    return no_dups_plink_prefix, plink_stdout
