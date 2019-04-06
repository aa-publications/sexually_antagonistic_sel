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
 
file_base = results.bfile
data_dir = results.data_dir
output_dir = results.output_dir

# -----------
# FUNCTIONS
# -----------

def get_num_lines(file):

    cmd = "wc -l {}".format(os.path.join(data_dir, file))

    run_output = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
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

file_base = "MEGA_ex_Array_Batch13_Cox_17_preQC_GRID" # prefix before .bed/.bim/.ped
data_dir = "/dors/capra_lab/users/abraha1/projects/PTB_phewas/data/biovu_samples_MEGAx_phewas/raw_preQC" 
output_dir = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data"

print("Running qc_per_batch.py on {}...\
        \n\tfile: {}\n\tdata_dir: {}\n\toutput_dir: {}".format(DATE, file_base, data_dir, output_dir))


tracking_indvids = {}
tracking_snps = {}

#
#   RAW DATA 
#

raw_file_path = os.path.join(data_dir, file_base)
basefile = os.path.basename(raw_file_path)

tracking_individs['raw_data'] = get_num_lines(raw_file_path+".fam")
tracking_snps['raw_data'] = get_num_lines(raw_file_path+".fam")




#
#   REMOVE DUPLICATE INDIVIDUALS 
#

# 1) identify duplicates in .fam file 
# 2) rename each duplicate with _1, _2, _3 etc. 



df = pd.read_csv(raw_file_path+".fam", sep="\s+", columns=['FID','IID','col3','col4','col5','col6'])

#
#   REMOVE INDIVIDUALS W/ DISCORDANT SEX
#

#   1. split PAR of X chromosome 
split_x_out_file = "{}_splitX".format(os.path.join(output_dir, basefile))
split_x_cmd = ("plink --bfile {}"
                   " --split-x b37 no-fail"
                   " --make-bed"
                   " --out {}").format(raw_file_path, split_x_out_file)
split_x_out = run_shell_cmd(split_x_cmd)

#   2. check sex 
chk_sex_out_file = "{}_check_sex".format(os.path.join(output_dir, basefile))
chk_sex_cmd = ("plink --bfile {}"
                   " --check-sex"
                   " --out {}").format(split_x_out_file, chk_sex_out_file)
chk_sex_out = run_shell_cmd(chk_sex_cmd)

#   3. identify individual to remove 
fig_file = os.path.join(output_dir, '{}_chrX_F_measure.png'.format(basefile))
ids_to_remove_file = os.path.join(output_dir,  '{}_ids_to_remove.txt'.format(basefile))
discordant_sex_cmd = ("python remove_discordant_sex.py"
                   " {}"
                   " {}"
                   " {}"
                   " {}").format(chk_sex_out_file+".sexcheck", basefile, fig_file, ids_to_remove_file) 
discordant_out = run_shell_cmd(discordant_sex_cmd)

#   4. remove individuals with discordant IDs 
rm_discord_sex_out_file = "{}_concdordant_sex".format(os.path.join(output_dir, basefile))
rm_discord_sex_cmd = "plink --bfile {} --remove {} --make-bed --out {}".format( raw_file_path, ids_to_remove_file, rm_discord_sex_out_file)
rm_discordant_out = run_shell_cmd(rm_discord_sex_cmd)

tracking_individs['concord_sex'] = get_num_lines(rm_discord_sex_out_file+".fam")
tracking_snps['concord_sex'] = get_num_lines(rm_discord_sex_out_file+".fam")


#
#   EXCLUDE ON HIGH MISSING GENOTYPE AND SAMPLE CALL RATE 
#

filter1_out_file = "{}_temp1".format(os.path.join(output_dir, basefile))
filter1_cmd = ("plink --bfile {}"
               " --geno 0.05"
               " --mind 0.05"
               " --make-bed"
               " --out {}").format(raw_file_path, filter1_out_file)

filter1_output = run_shell_cmd(filter1_cmd)

# get number of individuals and variants
num_individ_filter1 = get_num_lines(filter1_out_file+".fam")
num_variants_filter1 = get_num_lines(filter1_out_file+".bim")

#
#  REMOVE INDIVIDUALS W/ HIGH HETEROZYGOSITY RATE 
#

# 1) identify indepdnent snps
ind_snp_out_file = "{}_indepSNP".format(os.path.join(output_dir, basefile))
ind_snp_cmd = ("plink --bfile {}"
               " --autosome"
               " --indep-pairwise 50 5 0.2"
               " --out {}").format(filter1_out_file, ind_snp_out_file)

ind_snp_output = run_shell_cmd(ind_snp_cmd)

# 2) calcualte heterozygosity using independent snps 
het_out_file = "{}_het".format(os.path.join(output_dir, basefile))
het_cmd = ("plink --bfile {}"
           " --extract {}"
           " --het"
           " --out {}").format(filter1_out_file, ind_snp_out_file+".prune.in", het_out_file)

het_cmd_output = run_shell_cmd(het_cmd)

# 3) identify smaples to remove based on hetrozygozity
het_individ_cmd = ("python calc_het_outliers.py"
                   " {}"
                   " {}").format(het_out_file+".het", het_out_file+".het_exclude.tsv")

het_individ_output = run_shell_cmd(het_individ_cmd)

# 4) remove individuals based on heterozygosity
het_removed_out_file = "{}_temp2".format(os.path.join(output_dir, basefile))
rm_het_cmd = ("plink --bfile {}"
              " --remove {}"
              " --make-bed"
              " --out {}").format(filter1_out_file, het_out_file+".het_exclude.tsv", het_removed_out_file)

rm_het_output = run_shell_cmd(rm_het_cmd)

# get number of individuals and variants
num_individ_after_het = get_num_lines(het_removed_out_file+".fam")
num_variants_after_het = get_num_lines(het_removed_out_file+".bim")

#
#   REMOVE DUPLICATE INDIVIDUALS
#



#
#   REMOVE INDIVIDUALS WITH HIGH RELATEDNESS
#

# 1) calculate relatedness
related_genome_out_file = "{}_genome".format(os.path.join(output_dir, basefile))
related_cmd = ("plink --bfile {}"
               " --extract {}"
               " --genome --min 0.2"
               " --out {}").format(het_removed_out_file, ind_snp_out_file+".prune.in", related_genome_out_file)

related_output = run_shell_cmd(related_cmd)

# 2) identify related individuals to remove
relatives_out_file = "{}_related_fid_exclude.tsv".format(os.path.join(output_dir, basefile))
relatives_cmd = ("python remove_related_individuals.py"
                 " {}"
                 " {}").format(related_genome_out_file+".genome", relatives_out_file)

relatives_output = run_shell_cmd(relatives_cmd)

# 3) remove related individuals
rm_relatives_out_file = "{}_clean".format(os.path.join(output_dir, basefile))
rm_relatives_cmd = ("plink --bfile {}"
                    " --remove {}"
                    " --make-bed"
                    " --out {}").format(het_removed_out_file, relatives_out_file, rm_relatives_out_file)

rm_relatives_output = run_shell_cmd(rm_relatives_cmd)

# 4) get number of individuals and variants
num_individ_after_relatives = get_num_lines(rm_relatives_out_file+".fam")
num_variants_after_relatives = get_num_lines(rm_relatives_out_file+".bim")

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


print("Done processing {}. Took {} minutes.".format(file_base, (time.time()-start)/60))