#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2019-04-07 14:15:11

import os
import sys
import subprocess
import pandas as pd 


# sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')

from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd

@track_fx
def rm_high_het_rate(input_prefix, output_dir, base_prefix): 
    """
        Remove individuals with high heterozygosity rate on autosomes
            - individuals with outlier heterozygosity on the autosomes may relfect contamination or poor DNA quality
            - note: expected heterozygosity will vary between populations

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
        het_removed_plink_file : str
            - full path with plink prefix high autosomal heterozygosity values removed.
        ind_snp_output_file : str
            - full path w/ prefix to a set of independent snps calculated with plink 
            - * note:  this is only the prefix. must append "prune.in" or "prune.out"
        all_stdout : tuple of strings 
            - each element is STDOUT from running all plink commands

        
    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    ind_snp_output_file = os.path.join(output_dir, "temp_indepSNP-{}".format( base_prefix))
    het_out_file = os.path.join(output_dir, "inter_hetro_rate_{}".format(base_prefix))
    ids_to_remove_by_het_file = het_out_file+".het_exclude.tsv"
    het_removed_plink_file = os.path.join(output_dir,"temp_clean_het_rate_{}".format(base_prefix))
    

    # =============  REMOVE HIGH HET INDVIDUALS =============

    # 1) identify indepdnent snps
    ind_snp_cmd = ("plink --bfile {}"
                " --autosome"
                " --indep-pairwise 50 5 0.2"
                " --out {}").format(input_prefix, ind_snp_output_file)

    ind_snp_stdout = run_shell_cmd(ind_snp_cmd)

    # 2) calcualte heterozygosity using independent snps 
    het_cmd = ("plink --bfile {}"
            " --extract {}"
            " --het"
            " --out {}").format(input_prefix, ind_snp_output_file+".prune.in", het_out_file)
    het_stdout = run_shell_cmd(het_cmd)

    # 3) identify smaples to remove based on hetrozygozity
    het_individ_cmd = ("python helper_calc_het_outliers.py"
                    " {} {}").format(het_out_file+".het", ids_to_remove_by_het_file)

    het_ids_stdout = run_shell_cmd(het_individ_cmd)

    # 4) remove individuals based on heterozygosity
    if os.path.isfile(ids_to_remove_by_het_file):
        rm_het_cmd = ("plink --bfile {}"
                    " --remove {}"
                    " --make-bed"
                    " --out {}").format(input_prefix, ids_to_remove_by_het_file, het_removed_plink_file)
    else: 
        rm_het_cmd = ("plink --bfile {}"
                    " --make-bed"
                    " --out {}").format(input_prefix, het_removed_plink_file)

    clean_het_stdout = run_shell_cmd(rm_het_cmd)

    all_stdout = (ind_snp_stdout, het_stdout, het_ids_stdout, clean_het_stdout)

    return het_removed_plink_file, ind_snp_output_file,  all_stdout