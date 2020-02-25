

#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-01-22 20:14:30


import os
import sys


#
MEGA_PROBES_FILE="/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/MEGAEx_BioVU_15075710_A1_name_snp_probe_chr_pos.csv"
OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/check_probe_seqs/temp_MEGAEX_probes_v1.fa"

# %%

fw = open(OUTPUT_FILE, 'w')

with open(MEGA_PROBES_FILE, 'r') as fo:

    for num, line in enumerate(fo):
        if num ==0:
            continue


        split_line = line.split(',')

        split_line


        a1=split_line[1].split('/')[0][1:]
        a2=split_line[1].split('/')[1][:-1]
        seqA = f">seqA,{split_line[0]},{split_line[4]},{split_line[5].splitlines()[0]},{a1},{a2}\n{split_line[2]}\n"
        seqB = f">seqB,{split_line[0]},{split_line[4]},{split_line[5].splitlines()[0]},{a1},{a2}\n{split_line[3]}\n"


        if split_line[4] in ["MT", "X", "XY","Y"]:
            continue


        if (len(split_line[3]) > 1) and (len(split_line[2]) > 1):
            line_to_write = seqA+seqB
        else:
            line_to_write = seqA


        fw.write(line_to_write)



fw.close()
