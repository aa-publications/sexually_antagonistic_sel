#!/bin/python
# This script will perform QC on each batch obtained from BioVU.
#
#
#
# Abin Abraham
# created on: 2019-04-04 08:30:42


import os
import sys
import numpy as np
import pandas as pd
import subprocess

import glob
import argparse
import time
from datetime import datetime
 
DATE = datetime.now().strftime('%Y-%m-%d')
start = time.time()

# retrieve command line args 
parser = argparse.ArgumentParser(description='Example with nonoptional arguments')
parser.add_argument('bfile', action='store', type=str, help="plink bfile prefix")
parser.add_argument('data_dir', action='store', type=str, help="dir to plink input files")
parser.add_argument('output_dir', action='store')
  
results = parser.parse_args()
 
plink_prefix = results.bfile
data_dir = results.data_dir
output_dir = results.output_dir

# -----------
# FUNCTIONS
# -----------

def get_num_lines(file):

    cmd = "wc -l {}".format(os.path.join(data_dir, file))

    run_output = subprocess.run(cmd.split(), stderr=subprocess.PIPE,  stdout=subprocess.PIPE)
    exit_status = run_output.returncode
    output = run_output.stdout.decode('utf-8')
    assert exit_status == 0, "wc -l command had {} exit status.\nOUTPUT:\n{}".format(exit_status, output)

    num_lines = int(output.split()[0])

    return num_lines

def run_shell_cmd(cmd):
    run_output = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
    exit_status = run_output.returncode
    output = run_output.stdout.decode('utf-8')

    assert exit_status == 0, "COMMAND:\n{}\n\nEXIT STATUS:\n{}\n\nOUTPUT: {}".format(cmd, exit_status, output)

    return output

# -----------
# MAIN
# -----------

plink_prefix = "MEGA_ex_Array_Ancestry_MEGA_preQC_GRID" # prefix before .bed/.bim/.ped
data_dir = "/Users/abin-personal/Documents/katja_biobank/katja_biobank/data" 
output_dir = "/Users/abin-personal/Documents/katja_biobank/katja_biobank/data"

print("Running qc_per_batch.py on {}...\
        \n\tfile: {}\n\tdata_dir: {}\n\toutput_dir: {}".format(DATE, plink_prefix, data_dir, output_dir))


tracking_indvids = {}
tracking_snps = {}

#
#   RAW DATA 
#

raw_file_path = os.path.join(data_dir, plink_prefix)

tracking_indvids['raw_data'] = get_num_lines(raw_file_path+".fam")
tracking_snps['raw_data'] = get_num_lines(raw_file_path+".bim")


#
#   REMOVE DUPLICATE INDIVIDUALS 
#           - area for improvment: don't just randomly remove duplicate, see if there are diff in missing rate/sex errors...

dup_samples_out_file = os.path.join(data_dir, 'dups_samples_{}.csv'.format(plink_prefix))
fam_df = pd.read_csv(raw_file_path+".fam", sep="\s+", names=['FID','IID','c3','c4','c5','c6'])


# create a copy of the original .fam file 
fam_df.to_csv(raw_file_path+"_orignal.fam", sep=" ", header=None, index=None)

# check that FID and IID are equal 
assert fam_df[~(fam_df.FID == fam_df.IID)].shape[0] == 0, "FID and IID are *not* equal in this file:\n{}".format(raw_file_path+".fam")

# identify & write duplicated FID & IID rows
dup_index = fam_df[fam_df.duplicated(subset=['FID','IID'], keep='first')].index
dup_fids = fam_df.iloc[dup_index, :].FID.unique()

# each duplicate, except for the first instance, will be append "_[counter]" to FID and IID 
for this_fid in dup_fids: 
    for counter, index_row in enumerate(fam_df.loc[fam_df['FID']== this_fid].iterrows()):
        index, this_row = index_row
        if counter == 0: 
            continue 
        else: 
            fam_df.loc[index, ['FID', 'IID']] = fam_df.loc[index, ['FID', 'IID']].apply(lambda x: x+"_{}".format(counter))

# write duplicated FID and IID to file 
fam_df.loc[dup_index, ['FID','IID']].to_csv( dup_samples_out_file, sep=" ", header=None, index=None)

# modify existing .fam to include dups
fam_df.to_csv(raw_file_path+".fam", sep=" ", header=None, index=None)

# remove duplciates
rm_dups_out_file = os.path.join(output_dir, "{}_dedups".format(plink_prefix))
rm_dups_cmd = "plink --bfile {} --remove {} --make-bed --out {}".format( raw_file_path, dup_samples_out_file, rm_dups_out_file)
rm_dups_out = run_shell_cmd(rm_dups_cmd)


#
#   REMOVE INDIVIDUALS W/ DISCORDANT SEX
#

#   1. split PAR of X chromosome (must do this first and seperately before checking sex)
split_x_out_file = "{}_splitX".format(raw_file_path)
split_x_cmd = ("plink --bfile {}"
                   " --split-x b37 no-fail"
                   " --make-bed"
                   " --out {}").format(rm_dups_out_file, split_x_out_file)
split_x_out = run_shell_cmd(split_x_cmd)

#   2. check sex 
chk_sex_out_file = "{}_check_sex".format(raw_file_path)
chk_sex_cmd = ("plink --bfile {}"
                   " --check-sex"
                   " --out {}").format(split_x_out_file, chk_sex_out_file)
chk_sex_out = run_shell_cmd(chk_sex_cmd)

#   3. identify individual to remove 
fig_file = os.path.join(output_dir, '{}_chrX_F_measure.png'.format(plink_prefix))
ids_to_remove_file = os.path.join(output_dir,  '{}_ids_to_remove.txt'.format(plink_prefix))
discordant_sex_cmd = ("python remove_discordant_sex.py"
                   " {}"
                   " {}"
                   " {}"
                   " {}").format(chk_sex_out_file+".sexcheck", plink_prefix, fig_file, ids_to_remove_file) 
discordant_out = run_shell_cmd(discordant_sex_cmd)

#   4. remove individuals with discordant IDs  & set het. haploid and non-missing haploids as missing
rm_discord_sex_out_file = os.path.join(output_dir, "{}_concdordant_sex".format(plink_prefix))
rm_discord_sex_cmd = "plink --bfile {} --remove {} --set-hh-missing --make-bed --out {}".format( split_x_out_file, ids_to_remove_file, rm_discord_sex_out_file)
rm_discordant_out = run_shell_cmd(rm_discord_sex_cmd)

tracking_indvids['concord_sex'] = get_num_lines(rm_discord_sex_out_file+".fam")
tracking_snps['concord_sex'] = get_num_lines(rm_discord_sex_out_file+".fam")


#
#   EXCLUDE ON HIGH MISSING GENOTYPE AND SAMPLE CALL RATE 
#

missing_samples_file = "{}_sample_missing".format(os.path.join(output_dir, plink_prefix))
filter1_cmd = ("plink --bfile {}"
               " --mind 0.05"
               " --make-bed"
               " --out {}").format(rm_discord_sex_out_file, missing_samples_file)

filter1_output = run_shell_cmd(filter1_cmd)

# get number of individuals and variants
tracking_indvids['sample_missing'] = get_num_lines(missing_samples_file+".fam")
tracking_snps['sample_missing'] = get_num_lines(missing_samples_file+".fam")

#
#  REMOVE INDIVIDUALS W/ HIGH HETEROZYGOSITY RATE 
#
#   individuals with outlier heterozygosity on the autosome may relfect contamination or poor DNA quality 
#   note: expected heterozygosity will vary between populations

# 1) identify indepdnent snps
ind_snp_out_file = "{}_indepSNP".format(os.path.join(output_dir, plink_prefix))
ind_snp_cmd = ("plink --bfile {}"
               " --autosome"
               " --indep-pairwise 50 5 0.2"
               " --out {}").format(missing_samples_file, ind_snp_out_file)

ind_snp_output = run_shell_cmd(ind_snp_cmd)

# 2) calcualte heterozygosity using independent snps 
het_out_file = "{}_het".format(os.path.join(output_dir, plink_prefix))
het_cmd = ("plink --bfile {}"
           " --extract {}"
           " --het"
           " --out {}").format(missing_samples_file, ind_snp_out_file+".prune.in", het_out_file)

het_cmd_output = run_shell_cmd(het_cmd)

# 3) identify smaples to remove based on hetrozygozity
het_individ_cmd = ("python calc_het_outliers.py"
                   " {}"
                   " {}").format(het_out_file+".het", het_out_file+".het_exclude.tsv")

het_individ_output = run_shell_cmd(het_individ_cmd)

# 4) remove individuals based on heterozygosity
het_removed_out_file = "{}_clean_het_rate".format(os.path.join(output_dir, plink_prefix))
rm_het_cmd = ("plink --bfile {}"
              " --remove {}"
              " --make-bed"
              " --out {}").format(missing_samples_file, het_out_file+".het_exclude.tsv", het_removed_out_file)

rm_het_output = run_shell_cmd(rm_het_cmd)

# get number of individuals and variants
tracking_indvids['heterozygosity'] = get_num_lines(het_removed_out_file+".fam")
tracking_snps['heterozygosity'] = get_num_lines(het_removed_out_file+".fam")


#
#   REMOVE INDIVIDUALS WITH HIGH RELATEDNESS
#

# 1) calculate relatedness
related_genome_out_file = "{}_genome".format(os.path.join(output_dir, plink_prefix))
related_cmd = ("plink --bfile {}"
               " --extract {}"
               " --genome --min 0.2"
               " --out {}").format(het_removed_out_file, ind_snp_out_file+".prune.in", related_genome_out_file)

related_output = run_shell_cmd(related_cmd)

# 2) identify related individuals to remove
relatives_out_file = "{}_related_fid_exclude.tsv".format(os.path.join(output_dir, plink_prefix))
relatives_cmd = ("python remove_related_individuals.py"
                 " {}"
                 " {}").format(related_genome_out_file+".genome", relatives_out_file)

relatives_output = run_shell_cmd(relatives_cmd)

# 3) remove related individuals
rm_relatives_out_file = "{}_clean".format(os.path.join(output_dir, plink_prefix))
rm_relatives_cmd = ("plink --bfile {}"
                    " --remove {}"
                    " --make-bed"
                    " --out {}").format(het_removed_out_file, relatives_out_file, rm_relatives_out_file)

rm_relatives_output = run_shell_cmd(rm_relatives_cmd)

# 4) get number of individuals and variants

tracking_indvids['heterozygosity'] = get_num_lines(rm_relatives_out_file+".fam")
tracking_snps['heterozygosity'] = get_num_lines(rm_relatives_out_file+".fam")

#
# SUMMARIZE
#

columns = ["start", "removed_after_geno_mind_hwe", "removed_after_heterozygozity", "removed_after_relatedness", "final", "filename"]
indvid_line = [num_individ, 
                num_individ-num_individ_filter1, 
                num_individ_filter1-num_individ_after_het,
                num_individ_after_het-num_individ_after_relatives, 
                num_individ_after_relatives, basefile]

var_line = [num_variants,
                num_variants-num_variants_filter1,
                num_variants_filter1-num_variants_after_het,
                num_variants_after_het-num_variants_after_relatives,
                num_variants_after_relatives, basefile]

summary_file = os.path.join(output_dir, basefile+"_summary.tsv")
df = pd.DataFrame([indvid_line, var_line], index=['individual_count', 'variant_count'], columns=columns)
df.to_csv(summary_file, sep="\t")

#
# CLEAN UP TEMP FILES 
#

## make log director 
log_path=os.path.join(output_dir,'log')
summary_path=os.path.join(output_dir,'summary')

if not os.path.isdir(log_path):
    os.mkdir(log_path)

if not os.path.isdir(summary_path):
    os.mkdir(summary_path)

# move logs to log directory
for logfile in glob.glob(os.path.join(output_dir, "*.log")):
    os.rename(logfile, os.path.join(log_path, os.path.basename(logfile)))

# move summary files  to summary directory
for raw_file_path in glob.glob(os.path.join(output_dir, "*_summary.tsv")):
    os.rename(raw_file_path, os.path.join(summary_path, os.path.basename(raw_file_path)))


# remove all temp files
[os.remove(x) for x in glob.glob(filter1_out_file + ".*")]
[os.remove(x) for x in glob.glob(ind_snp_out_file + ".*")]
[os.remove(x) for x in glob.glob(het_out_file+ ".*") ] 
[os.remove(x) for x in glob.glob(het_out_file+".het_exclude.tsv") ] 
[os.remove(x) for x in glob.glob(het_removed_out_file+".*") ] 
[os.remove(x) for x in glob.glob(related_genome_out_file+".*") ] 
[os.remove(x) for x in glob.glob(relatives_out_file) ] 
[os.remove(x) for x in glob.glob(rm_relatives_out_file+".nosex") ] 


print("Done processing {}. Took {} minutes.".format(plink_prefix, (time.time()-start)/60))