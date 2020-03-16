#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-01-25 10:46:09

import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors
from datetime import datetime

from glob import glob

DATE = datetime.now().strftime('%Y-%m-%d')
import rpy2.rinterface
%load_ext rpy2.ipython


###
###    PATHS & THRESHOLDS
###


# %%
GWAS_P_THRESH = 5*10**-8
SUGG_GWAS_P_THRESH = 1*10**-6
MISS_P_THRESHOLD = 0.00001


#BLAT HITS mapped to GWAS and best hit picked by highest BLAT score
BEST_BLAT_DIR="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits/bivariate_blat_distributions"
UK_BLAT_FILE = os.path.join(BEST_BLAT_DIR, 'best_blatscore_xy_hit_length_filtered_uk_gwas.tsv')
BV_BLAT_FILE = os.path.join(BEST_BLAT_DIR, 'best_blatscore_xy_hit_length_filtered_bv_gwas.tsv')


# GWAS SUMMARY STATS
root="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/"
ukbb_path = os.path.join(root, "data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL.csv")
biovu_path= os.path.join(root, "results/2019_07_21_logistic/2019_07_21_logistic.assoc.logistic")

# OUTPUT
output_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/manhattan_plots/manhattan_w_best_blat_score"


# %%
###
###    python functions
###

def format_for_manh_plot(og_df):

    df = og_df.copy()
    df['Subject'] = df.targetChr.fillna('No X/Y Hit')
    df['cat_perID'] = df.per_identity.apply(lambda x: map_catperID(x))
    df['cat_perID'].unique()

    return df.loc[:, ['CHR','SNP','BP','P','Subject','cat_perID']].copy()




def map_catperID(x):
    if (np.isnan(x) or x <90):
        return "<90%"
    elif ((x >=90) and (x < 95)):
        return "≥90%"
    elif ( x >=95):
        return "≥95%"
    else:
        return np.nan



###
###    set up R env & functions
###

# %% - load env
%%R
library(ggplot2)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


# %% - set up manhattan plot params
%%R


# set plot parameters
theme_Publication <- function(base_size=14, base_family="Arial") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "plain", size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "plain",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.major.y = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="plain")
          ))

}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}


# %% - fx for plotting manhattan plt
%%R
gg.manhattan.blat.size <- function(df, threshold, hlight, ylims, title){
    threshold=5e-8 # for adding SNP rsID labels
    sig = 5e-8 # significant threshold line
    sugg = 1e-6 # suggestive threshold line

    # adjust x coordinate based on chromosome
    df.tmp = df %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% left_join(df, ., by=c("CHR"="CHR"))
    df.tmp = df.tmp %>% arrange(CHR, BP) %>% mutate(BPcum=BP+tot) %>% mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>% mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
    axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # create factors
    df.tmp$Subject = factor(df.tmp$Subject, levels = c("No X/Y Hit", "chrX", "chrY"))
    df.tmp$cat_perID = factor(df.tmp$cat_perID, levels = c("<90%", "≥90%", "≥95%"))
    subj_col = c("gray25", "maroon3", "palegreen3")
    size_vals = c(1.25, 3.25, 5)
    alpha_vals = c(0.3,1)

    # plot
    plt = ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=Subject, size=cat_perID, alpha=-log10(P))) +
            scale_color_manual(values = subj_col) +
            scale_size_manual(values= size_vals) +
            scale_alpha(range= alpha_vals ) +

             # custom X axis:
            scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
            scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis

            # add plot and axis titles
            ggtitle(paste0(title)) +
            labs(x = "Chromosome") +     # add genome-wide sig and sugg lines
            geom_hline(yintercept = -log10(sig), color='red2') +

            # Add label using ggrepel to avoid overlapping
            geom_text_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP)), size=5, color="gray41", force=1.3) +

            # set theme
            theme_bw(base_size = 22) +
            theme(
              plot.title = element_text(hjust = 0.5),
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())


        return(plt)
}


# %%
###
###    MAIN
###


# load gwas sum stats with best blat hit
uk_df = pd.read_csv(UK_BLAT_FILE, sep="\t")
bv_df = pd.read_csv(BV_BLAT_FILE, sep="\t")

bv_df.query('stat_sig == True')
uk_df.query('stat_sig == True')




# -----------
# format dfs for plotting
# -----------

# create column named "Subject" w/ values: c("No X/Y Hit", "chrX", "chrY")
# create column named "cat_perID" w/ values: "<90%", "≥90%", "≥95%"

bv_plt_df = format_for_manh_plot(bv_df)
uk_plt_df = format_for_manh_plot(uk_df)


# rename
rename_bv_snps = {'JHU_3.16652239' :'rs9870157',
'JHU_13.20119335' :'rs9508454',
'14:35761675-C-G' :'rs1048990'}

bv_plt_df.SNP.map(rename_bv_snps)



bv_plt_df.loc[bv_plt_df['SNP'].isin(rename_bv_snps.keys()), 'SNP'] = bv_plt_df.loc[bv_plt_df['SNP'].isin(rename_bv_snps.keys()), 'SNP'].map(rename_bv_snps).values


# %% -- UK BIOBANK
%%R -i output_dir -i uk_plt_df -h 500 -w 1200


plot_df = uk_plt_df
highlight_snps = (plot_df %>% filter(P < 5e-8) %>% select(SNP))$
SNP
hlight=highlight_snps

ylims=c(0,-1*log10(min(plot_df$P))+2)
title="UK Biobank"


plt = gg.manhattan.blat.size(plot_df, threshold, hlight, ylims, title) + theme_Publication(base_size=18)
plt = plt + guides(alpha = FALSE, color = guide_legend(override.aes = list(size=5))) +
                                    theme( legend.position=c(.9, .8),
                                    legend.margin = margin(t = 6, r = 3, b = 6, l = 3, unit = "pt"),
                                    legend.key.size= unit(0.2, "cm"),
                                    legend.background = element_rect(fill=alpha('gray', 0), color="black", size=0),
                                    legend.key = element_rect(fill = "transparent", colour = "transparent"),
                                    legend.box = "vertical", legend.direction="vertical", legend.box.just = "left",
                                    legend.title=element_text(face="bold", size=16)) +
                                    labs(size = "Sequence similarity", color="Top Blat Hit")



# png_file = file.path(output_dir, paste( Sys.Date(), "_ukbb_manhplot_best_blatscore_colorsize.pdf", sep=""))
# ggsave(png_file, plot=plt)

no_lgd_plt = plt + theme(legend.position = "none")
nolgd_png_file = file.path(output_dir, paste( Sys.Date(), "_ukbb_manhplot_best_blatscore_colorsize_nolegend.png", sep=""))
ggsave(nolgd_png_file, plot=no_lgd_plt)
# ggsave(nolgd_png_file, plot=no_lgd_plt,  device=cairo_pdf)



# %% -- BioVU BIOBANK
%%R -i output_dir -i bv_plt_df -h 500 -w 1200


plot_df = bv_plt_df
highlight_snps = (plot_df %>% filter(P < 5e-8) %>% select(SNP))$SNP
hlight=highlight_snps

ylims=c(0,-1*log10(min(plot_df$P))+2)
title="BioVU"


plt = gg.manhattan.blat.size(plot_df, threshold, hlight, ylims, title) + theme_Publication(base_size=18)
plt = plt + guides(alpha = FALSE, color = guide_legend(override.aes = list(size=5))) +
                                    theme( legend.position=c(.9, .8),
                                    legend.margin = margin(t = 6, r = 3, b = 6, l = 3, unit = "pt"),
                                    legend.key.size= unit(0.2, "cm"),
                                    legend.background = element_rect(fill=alpha('gray', 0), color="black", size=0),
                                    legend.key = element_rect(fill = "transparent", colour = "transparent"),
                                    legend.box = "vertical", legend.direction="vertical", legend.box.just = "left",
                                    legend.title=element_text(face="bold", size=16)) +
                                    labs(size = "Sequence similarity", color="Top Blat Hit")

png_file = file.path(output_dir, paste( Sys.Date(), "_biovu_manhplot_best_blatscore_colorsize.png", sep=""))
ggsave(png_file, plot=plt)


# no_lgd_plt = plt + theme(legend.title = element_blank())
# nolgd_png_file = file.path(output_dir, paste( Sys.Date(), "_biovu_manhplot_best_blatscore_colorsize_nolegend.png", sep=""))


# ggsave(nolgd_png_file, plot=no_lgd_plt, device=cairo_pdf)
# ggsave(nolgd_png_file, plot=no_lgd_plt)

