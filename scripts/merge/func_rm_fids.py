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
def rm_fids(dups_file, input_prefix, output_dir, base_prefix, prefix='temp_dedups_fids'):
    """
    Plink wrapper to remove duplicate individuals GIVEN A FILE W/ IDS TO REMOVE

    WARNING:
        - If there are no duplicates, then --make-bed is still run (to fit with the pipeline)


    Parameters
    ----------
    dups_file : str
        full path to tsv file with one line per FID and IID to remove
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
    no_dups_plink_prefix = os.path.join(output_dir, "{}_{}".format(prefix, base_prefix))

    # =============  RM FIDS  =============

    if os.path.isfile(dups_file):
        if  os.path.getsize(dups_file) > 0:
            rm_dups_cmd = "plink --bfile {} --remove {} --make-bed --out {}".format(
                input_prefix, dups_file, no_dups_plink_prefix)
    else:
        rm_dups_cmd = "plink --bfile {} --make-bed --out {}".format(input_prefix, no_dups_plink_prefix)

    plink_stdout = run_shell_cmd(rm_dups_cmd)

    return no_dups_plink_prefix, plink_stdout
