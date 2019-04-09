#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2019-04-07 14:03:47



import os
import sys
import subprocess
import pandas as pd

sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd

@track_fx
def rm_high_sample_missing_rate(input_prefix, output_dir, base_prefix):
    """
    Uses plink to remove samples with high missing rate.

    Parameters
    ----------
    input_prefix : str
        full path with plink prefix of file
    ouput_dir : str
        full path to directory to write outputs

    Returns
    -------
    ind_miss_removed_plink_prefix : str
        - full path with plink prefix removing samples with high missing rate 
    rm_miss_indiv_stdout : str
        - STDOUT from running plink command
    base_prefix : str
        the original plink prefix to be modified for output plink prefix 

    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    ind_miss_removed_plink_prefix = os.path.join(output_dir, "temp_sample_missing_{}".format(base_prefix))
    

    # ============= REMOVE SAMPLES WITH HIGH MISSING RATE =============

    shell_cmd = ("plink --bfile {}"
                " --mind 0.05"
                " --make-bed"
                " --out {}").format(input_prefix, ind_miss_removed_plink_prefix)

    rm_miss_indiv_stdout = run_shell_cmd(shell_cmd)


    return ind_miss_removed_plink_prefix, rm_miss_indiv_stdout