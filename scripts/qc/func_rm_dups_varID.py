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
def rm_dup_var_id(input_prefix, output_dir, base_prefix):
    """
    Plink wrapper to remove duplicate individuals.
        - The first instance is retained, all other are removed. 
        - See OUTPUT FILES for all outputs files created

    WARNING: 
        - The original .fam file will be overwritten with dupllicated. 
        - A copy of the original is made before overwritting.
        - If plink IDs no dup variants, then this will still run --make-bed 
        to be compatiable w/ the pipeline

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
    copy_of_orig_bim_file = os.path.join(output_dir, "original_"+base_prefix+".bim")
    dup_vars_output = os.path.join(output_dir, "inter_dup_vars_"+base_prefix)
    no_dups_vars_plink_prefix = os.path.join(output_dir, "temp_dedups_vars_{}".format(base_prefix))


    # =============  REMOVE MONOALLELIC & DUPLICATE SAMPLES =============
    
    # read and write a copy of fam file 
    bim_df = pd.read_csv(input_prefix+".bim", sep="\s+", names=['chrm','varID','morgans','basepair','A1','A2'])
    bim_df.to_csv(copy_of_orig_bim_file, sep=" ", header=None, index=None)

    # identify & write duplicated varIDs
    dup_index = bim_df[bim_df.duplicated(subset=['varID'], keep='first')].index
    dup_varID = bim_df.iloc[dup_index, :].varID.unique()

    # each duplicate varID, except for the first instance, will be have "_[counter]" appened
    for this_varID in dup_varID: 
        for counter, index_row in enumerate(bim_df.loc[bim_df['varID']== this_varID].iterrows()):
            index, this_row = index_row
            if counter == 0: 
                continue 
            else: 
                bim_df.loc[index, ['varID']] = bim_df.loc[index, ['varID']].apply(lambda x: x+"_{}".format(counter))

    # OVERWRITE existing .fam to taggin duplicates
    bim_df.to_csv(input_prefix+".bim", sep=" ", header=None, index=None)
    
    # call plink to identify dups 
    #   if there are no dups, then plink writes an empty file
    id_dup_vars_cmd = "plink --bfile {} --list-duplicate-vars suppress-first ids-only  --out {}".format(input_prefix,  dup_vars_output)
    ls_dups_vars_stdout = run_shell_cmd(id_dup_vars_cmd)

    
    # plink to rm duplicates
    if os.stat(dup_vars_output+".dupvar").st_size != 0:
        rm_dups_cmd = "plink --bfile {} --exclude {} --make-bed --out {}".format( input_prefix, dup_vars_output+".dupvar", no_dups_vars_plink_prefix)
    else: 
        rm_dups_cmd = "plink --bfile {} --make-bed --out {}".format( input_prefix,  no_dups_vars_plink_prefix)
    
    plink_stdout = run_shell_cmd(rm_dups_cmd)

    all_stdout = (ls_dups_vars_stdout, plink_stdout)
    return no_dups_vars_plink_prefix, all_stdout
