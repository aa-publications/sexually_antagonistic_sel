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



def pca(input_prefix, output_dir, base_prefix, prefix='temp_pca_'):
    """
    Run pca from plink

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
    hwe_outputs : str
        - prefix for plink hwe outputs
    plink_stdout : str
        - STDOUT from running plink command


    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    pca_outputs = os.path.join(output_dir, "{}_{}".format(prefix, base_prefix))

    # =============  RM FIDS  =============
    pca_cmd = "plink --bfile {} --pca 10 header tabs --out {}".format(input_prefix, pca_outputs)
    plink_stdout = run_shell_cmd(pca_cmd)

    return pca_outputs, plink_stdout
