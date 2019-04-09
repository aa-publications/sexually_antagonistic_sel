#!/bin/python
# See function doc string
#
#
#
# Abin Abraham
# created on: 2019-04-07 14:52:34

import os
import sys
import subprocess
import pandas as pd 


# sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')

from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd

@track_fx
def rm_discordant_sex(input_prefix, output_dir, base_prefix):
    """
        Remove individuals with discordant sex using plink. 
            - See OUTPUT FILES for all outputs files created
                - writes a list of individuals that are removed (ids_to_remove_file)
                - plots a figure of F-measure of X chromosome (chrx_Fmeasure_fig_file)

        WARNING: 
            * after removing individuals w/  discordant sex; the  --set-hh-missing is run! 
                - set-hh-missing converts all het. haploid and non-missing male calls to missing.
            * if there are no discordant IDs to remove, then --make-bed is still run (to be compatible with pipeline.)

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
        cleaned_sex_plink_prefix : str
            - full path with plink prefix with discordant sex values removed
        all_stdout : tuple of strings 
            - each element is STDOUT from running all plink commands

        
    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    split_x_pprefix = os.path.join(output_dir, "temp_split_X_{}".format(base_prefix))
    chk_sex_out_file = os.path.join(output_dir, "inter_check_sex_{}".format(base_prefix))
    chrx_Fmeasure_fig_file = os.path.join(output_dir, 'chrX_F_measure_{}.png'.format(base_prefix))
    ids_to_remove_file = os.path.join(output_dir,  'inter_ids_to_remove_discord_sex_{}.txt'.format(base_prefix))
    cleaned_sex_plink_prefix = os.path.join(output_dir, "temp_concdordant_sex_{}".format(base_prefix))

    # ============= REMOVE DISCORDANT SEX =============

    #   1. split PAR of X chromosome (must do this first and seperately before checking sex)
    split_x_cmd = ("plink --bfile {}"
                    " --split-x b37 no-fail"
                    " --make-bed"
                    " --out {}").format(input_prefix, split_x_pprefix)
    splitx_stdout = run_shell_cmd(split_x_cmd)

    #   2. check sex 
    chk_sex_cmd = ("plink --bfile {}"
                    " --check-sex"
                    " --out {}").format(split_x_pprefix, chk_sex_out_file)
    chk_sex_stdout = run_shell_cmd(chk_sex_cmd)

    #   3. identify individuals to remove 
    #       - if there are no individuals to remove, then no ids_to_remove_file is created...
    discordant_sex_cmd = ("python helper_remove_discordant_sex.py"
                    " {} {} {} {}").format(chk_sex_out_file+".sexcheck",
                                           pprefix, 
                                           chrx_Fmeasure_fig_file, 
                                           ids_to_remove_file) 
    py_id_individ_stdout = run_shell_cmd(discordant_sex_cmd)

    
    #   4. remove individuals with discordant IDs  & set het. haploid and non-missing haploids as missing
    if os.path.isfile(ids_to_remove_file):
        rm_discord_sex_cmd = ("plink --bfile {}"
                            " --remove {} --set-hh-missing --make-bed"
                            " --out {}").format( split_x_pprefix, 
                                                ids_to_remove_file,
                                                cleaned_sex_plink_prefix)
        
    else: 
        rm_discord_sex_cmd = ("plink --bfile {}"
                            " --set-hh-missing --make-bed"
                            " --out {}").format(split_x_pprefix, cleaned_sex_plink_prefix)
        
    rm_discord_sex_stdout = run_shell_cmd(rm_discord_sex_cmd)

    all_stdout = (splitx_stdout, chk_sex_stdout, py_id_individ_stdout, rm_discord_sex_stdout)

    return cleaned_sex_plink_prefix, all_stdout