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

# sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')

from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd


@track_fx
def rm_hwe_snps(input_prefix, output_dir, base_prefix, prefix='temp_hwe_filtered_'):
    """
    Plink wrapper to remove variants not in Hardy Weinberg Equilibirum.
        - See OUTPUT FILES for all outputs files created


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
    hwe_cmd : str
        - full path with plink prefix with hwe variants removed 
    hwe_stdout : str 
        - STDOUT from running plink command


    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    no_hwe_output = os.path.join(output_dir, "{}_{}".format(prefix, base_prefix))
    
    # =============  REMOVE MONOALLELIC & DUPLICATE SAMPLES =============


    # call plink to identify dups
    #   if there are no dups, then plink writes an empty file
    hwe_cmd = "plink --bfile {} --hwe 1e-50 midp include-nonctrl --make-bed --out {}".format(input_prefix, no_hwe_output)
    hwe_stdout = run_shell_cmd(hwe_cmd)


    return no_hwe_output, hwe_stdout
