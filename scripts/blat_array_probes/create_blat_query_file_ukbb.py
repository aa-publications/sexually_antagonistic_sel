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
BiL_PROBES_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/UKBiLEVE.annot.key.to.bim.tsv"
UKBB_PROBES_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/UKB_WCSG.annot.key.to.bim.tsv"

BIL_OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/ukbb_UKBiLEVE_probes.fa"
UKBB_OUTPUT_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/UKBB_WCSF_probes.fa"

BIL_ERROR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/ukbb_UKBiLEVE_probes_error.fa"
UKBB_ERROR_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/UKBB_WCSF_probes_error.fa"


def seq_to_write(seq): 
    
    

    
    try:
        probID,affyID,dbSNPID=split_line[0:3]
        
        first = split_line[6].split('[')[0]
        end = split_line[6].split('[')[1].split(']')[1]
        a1 = split_line[6].split('/')[0][-1]
        a2 = split_line[6].split('/')[1][0]
    
        return "{}{}{}".format(first, a1, end), "{}{}{}".format(first, a2, end)
    
    except IndexError as error:
        print("ERROR!")
        ferr.write(line)   
        
        return "", ""
        

# %%

fw = open(BIL_OUTPUT_FILE, 'w') 
ferr = open(BIL_ERROR_FILE, 'w') 
with open(BiL_PROBES_FILE, 'r') as fo: 
    
    for num, line in enumerate(fo): 
        
        if num ==0: 
            continue

        split_line = line.split('\t')
        
        
        
        seq_a, seq_b = seq_to_write(split_line)
           
        if len(seq_a) == 0: 
            continue

        probID,affyID,dbSNPID,chrom,pos = split_line[0:5]
        a1,a2 = split_line[7:9]
        mappedID = split_line[-1].splitlines()[0]

        line_to_write_seq_a = f">seq_a,{probID},{affyID},{dbSNPID},{chrom}:{pos},{a1},{a2},{mappedID}\n{seq_a}\n"
        line_to_write_seq_b = f">seq_b,{probID},{affyID},{dbSNPID},{chrom}:{pos},{a1},{a2},{mappedID}\n{seq_b}\n"

    
        fw.write(line_to_write_seq_a)
        fw.write(line_to_write_seq_b)
        


fw.close()
ferr.close()

# UKBB
fw = open(UKBB_OUTPUT_FILE, 'w') 
ferr = open(UKBB_ERROR_FILE, 'w')
with open(UKBB_PROBES_FILE, 'r') as fo: 
    
    for num, line in enumerate(fo): 
        
        if num ==0: 
            continue
        
        
        split_line = line.split('\t')    
        
            
        
        seq_a, seq_b = seq_to_write(split_line)
        
        if len(seq_a) == 0: 
            continue
        
        probID,affyID,dbSNPID,chrom,pos = split_line[0:5]
        a1,a2 = split_line[7:9]
        mappedID = split_line[-1].splitlines()[0]

        line_to_write_seq_a = f">seq_a,{probID},{affyID},{dbSNPID},{chrom}:{pos},{a1},{a2},{mappedID}\n{seq_a}\n"
        line_to_write_seq_b = f">seq_b,{probID},{affyID},{dbSNPID},{chrom}:{pos},{a1},{a2},{mappedID}\n{seq_b}\n"

    
        fw.write(line_to_write_seq_a)
        fw.write(line_to_write_seq_b)
        


fw.close()
ferr.close()