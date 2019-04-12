#!/bin/python
# See function doc string.
#
#
#
# Abin Abraham
# created on: 2019-04-07 14:03:47



import os
import sys
import subprocess
import pandas as pd

# sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
sys.path.append('/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/scripts/qc')

from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd

@track_fx
def test_missing(input_prefix, output_dir, base_prefix, pval_threshold=0.00001, prefix="inter_"):
    """
    Tests for differences in missing call counts for each variant between cases and controls.
    Vairiants with significant p-value after mult. testing is removed. 
    
    Parameters
    ----------
    input_prefix : str
        full path with plink prefix of file
    ouput_dir : str
        full path to directory to write outputs
    base_prefix : str
        the original plink prefix to be modified for output plink prefix 
    pval_threshold : float
        variatns w/ p values less than this will be reported for --test-missing
    prefix : str
        prefix to append to plink output

    WARNING: 
        IF no snps pass pval threshold to be removed, --make-bed is still run,
        so that it is compatiable with pipeline...

    Returns
    -------
    clean_case_ctrl_snps_plink_prefix : str
        - full path with plink prefix after removing variants with singificant difference between cases and controls
    all_stdout : str
        - STDOUT from running plink command


    """
    full_path, pprefix = os.path.split(input_prefix)

    # =============  OUTPUT FILES =============
    test_missing_output = os.path.join(output_dir, "{}_cctrl_miss_snps_{}".format(prefix, base_prefix))
    snps_to_rm  = os.path.join(output_dir, "{}_cctrl_miss_snps_to_rm_{}.txt".format(prefix, base_prefix))
    clean_case_ctrl_snps_plink_prefix = os.path.join(output_dir, "{}_cctrl_miss_snps_{}".format(prefix, base_prefix))

    
    # ============= TEST FOR MISSING RATE BETWEEN CASES AND CONTROL =============
    #   note: if no snps pass the pval threshol, plink still writes a file with headers but no other rows...
    shell_cmd = ("plink --bfile {}"
                " --test-missing midp"
                " --pfilter {}"
                " --out {}").format(input_prefix, pval_threshold, test_missing_output)

    test_miss_stdout = run_shell_cmd(shell_cmd)


    miss_df = pd.read_csv(test_missing_output+".missing", sep="\s+")

    if (miss_df.shape[0] > 0): 
        miss_df.SNP.to_csv(snps_to_rm, sep=" ", header=False, index=False)
        rm_snps_cmd = ("plink --bfile {}"
                " --exclude {}"
                " --make-bed"
                " --out {}").format(input_prefix, snps_to_rm, clean_case_ctrl_snps_plink_prefix)
    else:
        print("No SNPs pass pval threshold for different missingness b/w cases and control...")
        rm_snps_cmd = ("plink --bfile {}"
                " --make-bed"
                " --out {}").format(input_prefix, clean_case_ctrl_snps_plink_prefix)

    

    rm_snps_stdout = run_shell_cmd(rm_snps_cmd)
    all_stdout = (test_miss_stdout, rm_snps_stdout)

    return clean_case_ctrl_snps_plink_prefix, all_stdout