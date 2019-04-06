#!/bin/python
# This script will look throug pairs of related individuals and remove one randomly. 
#          * input: plink --genome output file with *.genome
#
#          * output: 2 column tab seperated file with IID and FID to remove
#
#
# Abin Abraham
# created on: 2018-10-04 14:19:39



import os
import sys
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='output a list of individuals to remove due to relatedness')

# REQUIRED ARGUMENTS IN ORDER
parser.add_argument('genome_file', action='store', type=str)
parser.add_argument('output_file', action='store', type=str)

# retrieve passed arguments
args = parser.parse_args()
genome_file = args.genome_file
output_file = args.output_file

print("Running: {}".format(sys.argv[:]))

# read *.genome
df = pd.read_csv(genome_file, sep="\s+", usecols=[0,2])

# a running list of FIDs removed
fid_removed = set()


for _, row in df.iterrows():

    fid1, fid2 = row[0], row[1]

    if (fid1 not in fid_removed) and (fid2 not in fid_removed): 
        # neither fids have already been removed 
        
        # so remove 1 of them
        fid_removed.add(fid2)

final_df = pd.DataFrame( {'FID': list(fid_removed),'IID':list(fid_removed)})
final_df.to_csv(output_file, sep="\t", index=False, header=False)

print("Total FIDs removed {}".format(len(fid_removed)))
print("Done. Removed one FID for related FID pairs. Check:\n{}".format(output_file))
