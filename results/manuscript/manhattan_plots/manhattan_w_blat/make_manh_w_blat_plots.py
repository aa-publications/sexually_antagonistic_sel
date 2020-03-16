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


# %%
GWAS_P_THRESH = 5*10**-8
SUGG_GWAS_P_THRESH = 1*10**-6
MISS_P_THRESHOLD = 0.00001
# %%
# FILE PATHS
BLAT_RSID_UK ="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/uk_bv_blat_w_assoc_rsID/uk_bil_blat_mapped_rsID.tsv"
BLAT_RSID_UKW ="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/uk_bv_blat_w_assoc_rsID/uk_wcsf_blat_mapped_rsID.tsv"
BLAT_RSID_BV ="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/uk_bv_blat_w_assoc_rsID/bv_blat_mapped_rsID.tsv"

# GWAS SUMMARY STATS
root="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/"
ukbb_path = os.path.join(root, "data/assoc_from_katja_2019_08_20/20190819_gwas_pca12_centers_age_FINAL.csv")
biovu_path= os.path.join(root, "results/2019_07_21_logistic/2019_07_21_logistic.assoc.logistic")

# OUTPUT
output_dir="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/manhattan_plots"

# %% markdown
# # load and filter
# %% markdown
# ## load gwas hits
# %%
# LOAD
raw_uk_df = pd.read_csv( ukbb_path, sep=",")
raw_bv_df = pd.read_csv( biovu_path, sep="\s+")

# %%
# format results
uk_gwas_df = raw_uk_df.loc[raw_uk_df.missing_p > MISS_P_THRESHOLD].copy()
print("Removed {:,} out of {:,} snps with singificant missingness b/w cases and controls.".format(raw_uk_df.shape[0] - uk_gwas_df.shape[0], raw_uk_df.shape[0]))
uk_gwas_df['chr_pos'] = uk_gwas_df.CHR.map(str) + ":" + uk_gwas_df.BP.map(str)


bv_gwas_df = raw_bv_df.loc[(raw_bv_df.TEST == "ADD") & (raw_bv_df.CHR < 23),].copy()
bv_gwas_df['chr_pos'] = bv_gwas_df.CHR.map(str) + ":" + bv_gwas_df.BP.map(str)
# %%
sug_sig_uk_df = uk_gwas_df.loc[uk_gwas_df.P < SUGG_GWAS_P_THRESH].copy()
sug_sig_bv_df = bv_gwas_df.loc[(bv_gwas_df.P < SUGG_GWAS_P_THRESH)].copy()
# %%
print("UKBB: Number of suggestive genome wide significant {:,}".format(sug_sig_uk_df.shape[0]))
print("BV: Number of suggestive genome wide significant {:,}".format(sug_sig_bv_df.shape[0]))
# %% markdown
# # merge blat hits with association hits
# %%
# load blat hits
ukbil_blat_df = pd.read_csv(BLAT_RSID_UK, sep="\t")
ukw_blat_df = pd.read_csv(BLAT_RSID_UKW, sep="\t")
bv_blat_df = pd.read_csv(BLAT_RSID_BV , sep="\t")
# %%



bv_blat_df.head()
pd.value_counts(bv_blat_df['q.length']).reset_index().sort_values('index')
# %% markdown
# # manhattan plots
# %% markdown
# ## set up R env
# %%
%%R
# manhattan function
library(ggplot2)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
# %%
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
# %% markdown
# ## plot
# %%
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
            geom_hline(yintercept = -log10(sugg), linetype="dashed", color='red4') +

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
%%R
gg.manhattan.blat <- function(df, threshold, hlight, ylims, title){
    threshold=5e-8 # for adding SNP rsID labels
    sig = 5e-8 # significant threshold line
    sugg = 1e-6 # suggestive threshold line

    # adjust x coordinate based on chromosome
    df.tmp = df %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% left_join(df, ., by=c("CHR"="CHR"))
    df.tmp = df.tmp %>% arrange(CHR, BP) %>% mutate(BPcum=BP+tot) %>% mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>% mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
    axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # create factors
    df.tmp$Subject = factor(df.tmp$Subject, levels = c("No Hit", "chrX", "chrY"))
    df.tmp$BlatHit = factor(df.tmp$BlatHit, levels = c('No X/Y Hit', 'chrX <90%', 'chrX ≥90%','chrX ≥95%',
                                                       'chrY <90%', 'chrY ≥90%', 'chrY ≥95%' ))

    hit_labels = c('chrX <90%', 'chrX ≥90%','chrX ≥95%','chrY <90%', 'chrY ≥90%', 'chrY ≥95%')
    chr.x.cols = rep(brewer.pal(n = 9, name = 'PiYG')[1],3)
    chr.y.cols = rep(brewer.pal(n = 9, name = 'PiYG')[8],3)
    subj_col = c(chr.x.cols, chr.y.cols)
    subj_points =c(20,18,8,20,18,8)


    # plot
    plt = ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
            geom_point(data=df.tmp[df.tmp$BlatHit=="No X/Y Hit",], color='gray21',alpha=0.7, size=2.5, shape=19) +
            geom_point(data=df.tmp[df.tmp$BlatHit!="No X/Y Hit",], aes(color=BlatHit, shape=BlatHit), alpha=1, size=5.5) +
            scale_alpha_discrete(range = c(0.3, 1,1,1,1,1,1,1)) +
            scale_color_manual(name="Sex Chromosome Similarity(%id)",labels=hit_labels, values = subj_col) +
            scale_shape_manual(name="Sex Chromosome Similarity(%id)",labels=hit_labels, values = subj_points) +


             # custom X axis:
            scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
            scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis

            # add plot and axis titles
            ggtitle(paste0(title)) +
            labs(x = "Chromosome") +     # add genome-wide sig and sugg lines
            geom_hline(yintercept = -log10(sig), color='firebrick') +
            geom_hline(yintercept = -log10(sugg), linetype="dashed", color='firebrick') +

            # Add label using ggrepel to avoid overlapping
            geom_text_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP)), size=5, color="gray41", force=1.3) +

            # set theme
            theme_bw(base_size = 22) +
            theme(
              plot.title = element_text(hjust = 0.5),
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
                legend.position=c(.9, .8),
                legend.margin = margin(t = 6, r = 3, b = 6, l = 3, unit = "pt"),
                legend.key.size= unit(0.2, "cm"),
                legend.background = element_rect(fill=alpha('gray', 0), color="black", size=0),
                legend.key = element_rect(fill = "transparent", colour = "transparent"),
                legend.box = "vertical", legend.direction="vertical", legend.box.just = "left",
                legend.title=element_text(face="bold", size=16)) + guides(alpha = FALSE)
        return(plt)
}
# %%
%%R -i output_dir -i uk_merged_df -h 500 -w 1200

##
## UKBB
##

plot_df = uk_merged_df
highlight_snps = (plot_df %>% filter(P < 5e-8) %>% select(SNP))$SNP
hlight=highlight_snps

ylims=c(0,-1*log10(min(plot_df$P))+2)
title="UKBiobank"


plt = gg.manhattan.blat(plot_df, threshold, hlight, ylims, title) + theme_Publication(base_size=18)
plt = plt + theme( legend.position=c(.9, .8),
        legend.margin = margin(t = 6, r = 3, b = 6, l = 3, unit = "pt"),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.spacing.x = unit(0.4, 'cm'),
        legend.key.size= unit(0.2, "cm"),
        legend.background = element_rect(fill=alpha('gray', 0), color="black", size=0),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.box = "vertical", legend.direction="vertical", legend.box.just = "left",
        legend.title=element_text(face="bold", size=16)) +
        labs(size = "Sequence similarity", color="Top Blat Hit") + guides(color = guide_legend(ncol=2, override.aes = list(size=5)))

plt
png_file = file.path(output_dir, paste( Sys.Date(), "_ukbb_manhplot_blat_colorshape.png", sep=""))
ggsave(png_file, plot=plt)

# %%
bv_best_hit_df.sort_values('P').head()
# %%
bv_merged_df.sort_values('P').head()
# %%
%%R -i output_dir -i bv_merged_df -h 500 -w 1200

##
## BIOVU
##

plot_df = bv_merged_df
highlight_snps = (plot_df %>% filter(P < 5e-8) %>% select(SNP))$SNP
hlight=highlight_snps

ylims=c(0,-1*log10(min(plot_df$P))+2)
title="BioVU"


plt = gg.manhattan.blat(plot_df, threshold, hlight, ylims, title) + theme_Publication(base_size=18)
plt = plt + theme( legend.position=c(.9, .8),
        legend.margin = margin(t = 6, r = 3, b = 6, l = 3, unit = "pt"),
        legend.key.size= unit(0.2, "cm"),
        legend.background = element_rect(fill=alpha('gray', 0), color="black", size=0),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.box = "vertical", legend.direction="vertical", legend.box.just = "left",
        legend.title=element_text(face="bold", size=16)) +
        labs(size = "Sequence similarity", color="Top Blat Hit") +  guides(color = guide_legend(override.aes = list(size=5)))

# plt
png_file = file.path(output_dir, paste( Sys.Date(), "_bv_manhplot_blat_colorshape.png", sep=""))
ggsave(png_file, plot=plt)
# %% markdown
# ### use color and size
# %%
%%R -i output_dir -i uk_merged_df -h 500 -w 1200


plot_df = uk_merged_df
highlight_snps = (plot_df %>% filter(P < 5e-8) %>% select(SNP))$SNP
hlight=highlight_snps

ylims=c(0,-1*log10(min(plot_df$P))+2)
title="UKBiobank"


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

png_file = file.path(output_dir, paste( Sys.Date(), "_ukbb_manhplot_blat_colorsize.png", sep=""))
ggsave(png_file, plot=plt)

# %%
%%R -i output_dir -i bv_merged_df -h 500 -w 1200


plot_df = bv_merged_df
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

png_file = file.path(output_dir, paste( Sys.Date(), "_bv_manhplot_blat_colorsize.png", sep=""))
ggsave(png_file, plot=plt)

# %%
