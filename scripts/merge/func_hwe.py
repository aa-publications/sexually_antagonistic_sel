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

from func_run_shell_cmd import run_shell_cmd


def hwe(input_prefix, output_dir, base_prefix, prefix='temp_hwe_', plink_opts=None):
    """
    Run HWE in plink. 

    Parameters
    ----------
    input_prefix : str
        full path with plink prefix of file 
    ouput_dir : str
        full path to directory to write outputs
    base_prefix : str
        the original plink prefix to be modified for output plink prefix 
    plink_opts : str
        plink options to add on. no extra space. 
            e.g.: '--filter-cases


    Returns
    -------
    hwe_outputs : str
        - prefix for plink hwe outputs
    plink_stdout : str 
        - STDOUT from running plink command


    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    hwe_outputs = os.path.join(output_dir, "{}_{}".format(prefix, base_prefix))

    # =============  RM FIDS  =============
    if plink_opts:
        hardy_cmd = "plink --bfile {} {} --hardy midp --out {}".format(input_prefix, plink_opts, hwe_outputs)
    else: 
        hardy_cmd = "plink --bfile {} --hardy midp --out {}".format(input_prefix, hwe_outputs)
    plink_stdout = run_shell_cmd(hardy_cmd)

    return hwe_outputs, plink_stdout
