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
def calc_stats(input_prefix, output_dir, base_prefix, prefix='temp_frq_and_miss'):
    """
    Uses plink to calc frequency and missingness and outputs new plink bfiles. 
        - See OUTPUT FILES for all outputs files created

    Parameters
    ----------
    input_prefix : str
        full path with plink prefix of file 
    ouput_dir : str
        full path to directory to write outputs
    base_prefix : str
        the original plink prefix to be modified for output plink prefix 
    prefix : str
        prefix to append to plink output

    Returns
    -------
    dups_removed_plink_prefix : str
        - full path with plink prefix with duplicate samples removed
    plink_stdout : str 
        - STDOUT from running plink command


    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    frq_and_miss_plink_prefix = os.path.join(output_dir, "{}_{}".format(prefix, base_prefix))

    # =============  RUN BASIC STATS =============
    shell_cmd = "plink --bfile {} --missing --freqx --make-bed --out {}".format(input_prefix, frq_and_miss_plink_prefix) 

    plink_stdout = run_shell_cmd(shell_cmd)

    return frq_and_miss_plink_prefix, plink_stdout
