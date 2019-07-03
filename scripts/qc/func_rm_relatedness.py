#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2019-04-07 14:28:03

import os
import sys
import subprocess
import pandas as pd


# sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')

from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd

@track_fx
def rm_relatedness(input_prefix, output_dir, ind_snp_output_file, base_prefix):
    """
        Remove individuals with high relatedness.
            - See OUTPUT FILES for all outputs files created


        WARINING: if no highly related individuals present,
        this will still run --make-bed (so that this is  compatible with pipeline)



        Parameters
        ----------
        input_prefix : str
            full path with plink prefix of file
        ouput_dir : str
            full path to directory to write outputs
        ind_snp_output_file : str
            full path plink prefix obtaiend from calculatign indep pairwise SNPs
            * note: should not include the ".prune.in" or "prune.out" suffix
            * see func_rm_high_het_rate.py
        base_prefix : str
            the original plink prefix to be modified for output plink prefix

        Returns
        -------
        related_indiv_removed_plink_prefix : str
            - full path with plink prefix with related individuals removed
        all_stdout : tuple of strings
            - each element is STDOUT from running all plink commands


    """

    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    related_genome_output = os.path.join(output_dir, "inter_genome_{}".format(base_prefix))
    ids_to_remove_relatedness = os.path.join(output_dir, "inter_ids_to_remove_relatedness_{}.txt".format(base_prefix))
    related_indiv_removed_plink_prefix = os.path.join(output_dir, "temp_no_relatedness_{}".format(base_prefix))

    # ============= REMOVE HIGHLY RELATED INDIVIDUALS =============


    # 1) calculate relatedness
    related_cmd = ("plink --bfile {}"
                " --extract {}"
                " --genome --min 0.2"
                " --out {}").format(input_prefix, ind_snp_output_file+".prune.in", related_genome_output)

    related_stdout = run_shell_cmd(related_cmd)

    # 2) identify related individuals to remove
    relatives_cmd = ("python /dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc/helper_remove_related_individuals.py"
                    " {}"
                    " {}").format(related_genome_output+".genome", ids_to_remove_relatedness)

    id_relatives_stdout = run_shell_cmd(relatives_cmd)

    # 3) remove related individuals
    if os.path.isfile(ids_to_remove_relatedness):
        rm_relatives_cmd = ("plink --bfile {}"
                            " --remove {} --make-bed"
                            " --out {}").format(input_prefix, ids_to_remove_relatedness, related_indiv_removed_plink_prefix)
    else:
        rm_relatives_cmd = ("plink --bfile {}"
                            " --make-bed"
                            " --out {}").format(input_prefix, related_indiv_removed_plink_prefix)
    rm_relatives_output = run_shell_cmd(rm_relatives_cmd)


    all_stdout = (related_stdout, id_relatives_stdout, rm_relatives_output)
    return related_indiv_removed_plink_prefix, all_stdout
