#!/bin/python
# This script will 
#
#
#
# Abin Abraham
# created on: 2020-01-25 14:13:22

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime


from glob import glob
DATE = datetime.now().strftime('%Y-%m-%d')


# PATHS
uk_ano_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes"
uk_bil_annot_file = os.path.join(uk_ano_dir, "Axiom_UKBiLEVE.na34.annot_clean.csv")
uk_wcsg_annot_file = os.path.join(uk_ano_dir, "Axiom_UKB_WCSG.na34.annot_clean.csv")
uk_bim_dir = "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/ukbb_bim_files"


bv_bim_file= "/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/MEGAex_BioVUthru2019-02/MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3_qc/final_qc_MEGAex_BioVUthru2019-02_BestOfMultipleCalls_Capra_Preterm_A3.bim"
bv_annot_file="/dors/capra_lab/users/abraha1/data/biovu_mega_ex_2019_02_capra_preterm_a3/MEGAEx_BioVU_15075710_A1_name_snp_probe_chr_pos.csv"


KEY_FILE_OUTPUT="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/ukbb_probes/UKBiLEVE.annot.key.to.bim.tsv"

# -----------
# FUNCTIONS
# -----------


def intersect_uk(uk_bim, anno_df, uk_col, anno_col): 
        

    n_all = uk_bim[uk_col].nunique()
    both = set(uk_bim[uk_col].unique()).intersection(set(anno_df[anno_col].unique()))
    bim_only = set(uk_bim[uk_col].unique()).difference(set(anno_df[anno_col].unique()))
    anno_only = set(anno_df[anno_col].unique()).difference(set(uk_bim[uk_col].unique()))

    print("{:,} out of {:,} SNPs mapped by {}.".format(len(both), n_all, anno_col))
    
    return both, bim_only, anno_only


def intersect(a, b): 
    
    both = a.intersection(b)
    a_only = a.difference(b)
    b_only = b.difference(a)
    
    print("{:,} shared, {:,} only in first set,  {:,} only in second set".format(len(both), len(a_only), len(b_only)))
    
    return both, a_only, b_only

# %%

# -----------
# MAIN
# -----------

# load uk bim file 
uk_bim = pd.DataFrame()
for bfile in glob(uk_bim_dir+"/*.bim"): 
    
    bdf = pd.read_csv(bfile, sep="\s+", names=['chr','SNP','dummy','pos','A1','A2'])
    
    uk_bim = uk_bim.append(bdf)
    

keep_cols = ['Probe Set ID', 'Affy SNP ID', 'dbSNP RS ID', 'Chromosome', 'Physical Position', 'Flank', 'Allele A', 'Allele B', 'Strand', 'Annotation Notes']
anno_df = pd.read_csv(uk_bil_annot_file, sep=",", low_memory=False, usecols=keep_cols)

uk_bim['chr_pos_A1_A2'] = uk_bim['chr'].map(str) +"_"+  uk_bim['pos'].map(str) +"_"+ uk_bim['A1'].map(str) +"_"+ uk_bim['A2'].map(str) 
anno_df['chr_pos_A1_A2'] = anno_df['Chromosome'].map(str) +"_"+  anno_df['Physical Position'].map(str) +"_"+ anno_df['Allele A'].map(str) +"_"+ anno_df['Allele B'].map(str) 
anno_df['chr_pos_A2_A1'] = anno_df['Chromosome'].map(str) +"_"+  anno_df['Physical Position'].map(str) +"_"+ anno_df['Allele B'].map(str) +"_"+ anno_df['Allele A'].map(str) 


uk_bim_pos_dict = dict(zip(uk_bim.chr_pos_A1_A2, uk_bim.SNP))


# %%
###
#   COMPARE rsIDS
###

uk_bim.head()
all_bim_snps = set(uk_bim['SNP'].unique())
all_bim_chpos = set(uk_bim['chr_pos_A1_A2'].unique())


anno_col='dbSNP RS ID'
rsid_both, rsid_bim_only, rsid_anno_only = intersect_uk(uk_bim, anno_df, 'SNP', anno_col)

anno_col='Affy SNP ID' # not as good 
affy_both, affy_bim_only, affy_anno_only = intersect_uk(uk_bim, anno_df, 'SNP', anno_col)


# any shared between rsid and affy
both, a_only, b_only = intersect(rsid_both, affy_both)

# combine affy and dBSNP IDs 
rsid_affy_ids = rsid_both.union(affy_both)
both, a_only, b_only = intersect(rsid_affy_ids, all_bim_snps)


# map by chr_pos
chpos1_both, chpos1_bim_only, chpos1_anno_only = intersect_uk(uk_bim, anno_df, 'chr_pos_A1_A2', 'chr_pos_A1_A2')
chpos2_both, chpos2_bim_only, chpos2_anno_only = intersect_uk(uk_bim, anno_df, 'chr_pos_A1_A2', 'chr_pos_A2_A1')
chpos_snps = chpos1_both.union(chpos2_both)

both, a_only, b_only = intersect(chpos_snps, all_bim_chpos)



# USE both rsID and affy and chr_pos 

rsid_both # rsIDs in BIM and ANNO File 
affy_both  # affyID in BIM and ANNO 
chpos_snps # chpos in both bim and anno 

matched_anno_df = anno_df.loc[anno_df['dbSNP RS ID'].isin(rsid_both)| 
            anno_df['Affy SNP ID'].isin(affy_both)| 
            anno_df['chr_pos_A1_A2'].isin(chpos_snps)|
            anno_df['chr_pos_A2_A1'].isin(chpos_snps)].copy()

anno_df.shape
matched_anno_df.shape            


print("{:,} out of {:,} ({:,} remaining) probes in the annotations files were matched.".format(matched_anno_df.shape[0], anno_df.shape[0], anno_df.shape[0]-matched_anno_df.shape[0]))

matched_bim_df = uk_bim.loc[ uk_bim['SNP'].isin(rsid_both)| 
                             uk_bim['SNP'].isin(affy_both)| 
                             uk_bim['chr_pos_A1_A2'].isin(chpos_snps)].copy()


uk_bim.shape
matched_bim_df.shape                         
print("{:,} out of {:,} ({:,} remaining) SNPs in the BIM file matched.".format(matched_bim_df.shape[0], uk_bim.shape[0], uk_bim.shape[0]-matched_bim_df.shape[0]))
# WRITE A KEY 
# for each probe in the annotation file, provide the rsID in the bim file 




key_anno_df = matched_anno_df.copy()
key_anno_df['bim_SNP'] = "None"



matched_anno_df.loc[matched_anno_df['dbSNP RS ID'].isin(rsid_both), 'bim_SNP'] = matched_anno_df['dbSNP RS ID']
matched_anno_df.loc[matched_anno_df['Affy SNP ID'].isin(rsid_both), 'bim_SNP'] = matched_anno_df['Affy SNP ID']
            

matched_anno_df.loc[matched_anno_df['chr_pos_A1_A2'].isin(chpos_snps), 'bim_SNP'] = matched_anno_df.loc[matched_anno_df['chr_pos_A1_A2'].isin(chpos_snps), 'chr_pos_A1_A2'].map(uk_bim_pos_dict)            
matched_anno_df.loc[matched_anno_df['chr_pos_A2_A1'].isin(chpos_snps), 'bim_SNP'] = matched_anno_df.loc[matched_anno_df['chr_pos_A2_A1'].isin(chpos_snps), 'chr_pos_A2_A1'].map(uk_bim_pos_dict)            
    
    

# check that all rsIDs got filled    
matched_anno_df.loc[matched_anno_df['bim_SNP']=="None"].shape[0]    


write_cols = ['Probe Set ID', 'Affy SNP ID', 'dbSNP RS ID', 'Chromosome',
       'Physical Position', 'Strand', 'Flank', 'Allele A', 'Allele B', 'bim_SNP']


final_df = matched_anno_df.loc[:, write_cols].copy()
final_df.to_csv(KEY_FILE_OUTPUT, sep="\t", index=False)

