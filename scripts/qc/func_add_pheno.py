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

sys.path.append('/Users/abin-personal/Documents/katja_biobank/katja_biobank/scripts/qc')
from func_track import track_fx
from func_run_shell_cmd import run_shell_cmd

@track_fx
def add_pheno(input_prefix, output_dir, base_prefix, cases='female'):
    """
    Tests for differences in missing call counts for each variant between cases and controls.
    Vairiants with significant p-value after mult. testing is removed. 
    See output files for intermediate files created.
    
    Parameters
    ----------
    input_prefix : str
        full path with plink prefix of file
    ouput_dir : str
        full path to directory to write outputs
    base_prefix : str
        the original plink prefix to be modified for output plink prefix 
    cases : str (default='female')
        'male' or 'female' to set as cases


    Returns
    -------
    phenotyped_plink_prefix : str
        - full path with plink prefix adding phenotype codes to fam file 
    rm_miss_indiv_stdout : str
        - STDOUT from running plink command


    """
    full_path, pprefix = os.path.split(input_prefix)
  


    # =============  OUTPUT FILES =============
    pheno_labels_file = os.path.join(full_path, "inter_pheno_labels_{}.txt".format(base_prefix))
    phenotyped_plink_prefix =os.path.join(full_path,  "temp_w_pheno_{}__{}_as_cases".format(base_prefix, cases))

    # ============= REMOVE SAMPLES WITH HIGH MISSING RATE =============

    fam_df = pd.read_csv(input_prefix+".fam", sep="\s+", header=None, names=['FID','IID','col3','col4','sex','pheno'])

    fam_df['pheno'] = 0
    fam_df.loc[fam_df['sex']==1, 'pheno'] = 'male'
    fam_df.loc[fam_df['sex']==2, 'pheno'] = 'female'
    fam_df.to_csv(pheno_labels_file, sep=" ", header=None, index=False, columns=['FID','IID','pheno'])

    shell_cmd = ("plink --bfile {}"
                " --make-pheno {} {}" # have females be listed as cases
                " --make-bed --out {}").format(input_prefix, pheno_labels_file,  cases, phenotyped_plink_prefix)

    make_pheno_stdout = run_shell_cmd(shell_cmd)


    return phenotyped_plink_prefix, make_pheno_stdout