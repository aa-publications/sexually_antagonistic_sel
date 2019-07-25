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
def extract_race_white(input_prefix, output_dir, base_prefix, prefix="temp_white_only"):
    """
    Filter and keep on 'white' individuals.
        - LIST OF 'WHITE' GRIDS ARE HARD CODED!


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
    WHITE_ONLY_PLINK_GRIDS="/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/grid_lists/race_is_w_grids_for_plink.txt"
    # CEU_ONLY_PLINK_GRIDS="/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/demographics/CEU_at_least_90percent_GRIDS_only_for_plink.txt"

    full_path, pprefix = os.path.split(input_prefix)
    # =============  OUTPUT FILES =============
    white_only_plink_prefix = os.path.join(output_dir, "temp_white_only_{}".format(base_prefix))

    # =============  FILTER AND INCLUDE ONLY EUROPEAN SAMPLES =============
    keep_white_cmd = "plink --bfile {} --keep {} --make-bed --out {}".format(
            input_prefix, WHITE_ONLY_PLINK_GRIDS, white_only_plink_prefix)

    plink_stdout = run_shell_cmd(keep_white_cmd)

    return white_only_plink_prefix, plink_stdout
