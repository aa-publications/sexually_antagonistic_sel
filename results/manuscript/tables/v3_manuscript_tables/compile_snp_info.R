nonsig_snps <- read.csv("ukbb_sig_snps_missing_bysex.csv", stringsAsFactors=FALSE)
nonsig_snps$SIGNIF <- FALSE
nonsig_snps$DATASET <- "ukb"
sig_snps <- read.table("sigSNP_HWE_w_NMISS.csv", stringsAsFactors=FALSE, sep="\t", header=TRUE)
sig_snps$SIGNIF <- TRUE
names(sig_snps)[names(sig_snps) == "P"] <- "HWE_PValue"
stopifnot(all(names(sig_snps) %in% names(nonsig_snps)))
for (xn in setdiff(names(nonsig_snps), names(sig_snps))) { sig_snps[[xn]] <- NA }
all_snps <- rbind(sig_snps, nonsig_snps)
all_snps$TEST[all_snps$TEST == "AFF"] <- "FEMALE"
all_snps$TEST[all_snps$TEST == "UNAFF"] <- "MALE"
names(all_snps)[names(all_snps) == "TEST"] <- "SEX"

all_snps$hom1 <- as.numeric(sapply(strsplit(all_snps$GENO, "/"), "[[", 1))
all_snps$het <- as.numeric(sapply(strsplit(all_snps$GENO, "/"), "[[", 2))
all_snps$hom2 <- as.numeric(sapply(strsplit(all_snps$GENO, "/"), "[[", 3))
all_snps$total <- (all_snps$hom1 + all_snps$het + all_snps$hom2)
all_snps$freq <- (0.5 * all_snps$het + all_snps$hom1)/all_snps$total
all_snps$exp_hom1 <- all_snps$freq^2 * all_snps$total

# get number missing 
sample_counts <- data.frame(DATASET=rep(c("ukb", "biovu"), each=3),
                            SEX=c("FEMALE", "MALE", "ALL"),
                            total_genotyped=c(264813, 223478, 264813+223478, 34269, 27491, 34269+27491))
all_snps$num_missing <- NA
for (j in 1:nrow(all_snps)) {
    k <- which(sample_counts$SEX == all_snps$SEX[j] & sample_counts$DATASET == all_snps$DATASET[j])
    all_snps$num_missing[j] <- sample_counts$total_genotyped[k] - all_snps$total[j]
}


# p-value for excess alternate homozygosity:
all_snps$hom_pval <- NA
for (j in 1:nrow(all_snps)) {
    all_snps$hom_pval[j] <- with(all_snps[j,], 
        pbinom(hom1, size=total, prob=freq^2))
}


# p-value for HWE:
# this p-value differs from HWE_PValue - why?
all_snps$chisq_pval <- NA
for (j in 1:nrow(all_snps)) {
    all_snps$chisq_pval[j] <- with(all_snps[j,], 
        chisq.test(c(hom1, het, hom2), p=c(freq^2, 2*freq*(1-freq), (1-freq)^2)))$p.value
}

plot(HWE_PValue ~ chisq_pval, data=all_snps)
abline(0,1)

plot(chisq_pval ~ hom_pval, data=all_snps)
abline(0,1)


# plot observed and expected hom1
plot(hom1 ~ exp_hom1, data=subset(all_snps, SEX != "ALL"),
     col=ifelse(SEX=="FEMALE", "red", "blue"), pch=20, cex=2+pmin(5, -log10(hom_pval)/3))
abline(0,1)
points(hom1 ~ exp_hom1, data=subset(all_snps, SEX != "ALL" & SIGNIF), cex=10)

with(all_snps, {
    plot(1+exp_hom1, 1+hom1, log='xy',
         col=ifelse(SEX=="FEMALE", "red", "blue"), pch=20, cex=2+pmin(5, -log10(hom_pval)/3));
    abline(0,1);
    points(1+exp_hom1[SIGNIF], 1+hom1[SIGNIF], cex=10)
})

plot(hom1 - exp_hom1 ~ exp_hom1, data=subset(all_snps, SEX != "ALL"),
     col=ifelse(SEX=="FEMALE", "red", "blue"), pch=20, cex=2+pmin(5, -log10(hom_pval)/3),
     log='x')
points(hom1 - exp_hom1 ~ exp_hom1, data=subset(all_snps, SEX != "ALL" & SIGNIF),
       pch="*", cex=3)
abline(h=0)

### make SNP summary

snps <- data.frame(SNP=unique(all_snps$SNP), stringsAsFactors=FALSE)
both <- subset(all_snps, SEX=="ALL")
for (xn in c("hom_pval", "DATASET", "CHR", "SIGNIF")) {
    snps[[xn]] <- both[[xn]][match(snps$SNP, both$SNP)]
}
fem <- subset(all_snps, SEX=="FEMALE")
male <- subset(all_snps, SEX=="MALE")
for (xn in c("hom1", "het", "hom2", "total", "freq", "exp_hom1", "num_missing")) {
    snps[[paste0("fem_", xn)]] <- fem[[xn]][match(snps$SNP, fem$SNP)]
    snps[[paste0("male", xn)]] <- male[[xn]][match(snps$SNP, male$SNP)]
}

snps$missing_pval <- NA
for (j in 1:nrow(snps)) {
    snps$missing_pval[j] <- with(snps[j,],
                                 fisher.test(rbind(c(fem_total, maletotal),
                                              c(fem_num_missing, malenum_missing)))$p.value)
}

plot(missing_pval ~ hom_pval, data=snps)

# how many are below 1e-4 for each test?
addmargins(with(snps, table(missing=missing_pval < 1e-4, homozygotes=hom_pval < 1e-4)))

##
# load homology

ukb_blat <- read.table("biovu_ukbb_blat_for_sighits_ukb.csv", header=TRUE, sep='\t')
biovu_blat <- read.table("biovu_ukbb_blat_for_sighits_biovu.csv", header=TRUE, sep='\t')
ukb_blat$DATASET <- "ukb"
biovu_blat$DATASET <- "biovu"
stopifnot(all(names(ukb_blat) == names(biovu_blat)))
blat <- rbind(ukb_blat, biovu_blat)
names(blat)[names(blat) == "X.id"] <- "percent_identity"

biovu_tab <- read.table("biovu_sex_sel_hits_100flank_blat.tab", skip=3, header=FALSE)
names(biovu_tab) <- c("ACTIONS", "MORE_ACTIONS", "QUERY", "SCORE", "START", "END", "QSIZE", "IDENTITY", "CHROM", "STRAND", "START", "END", "SPAN")
biovu_tab$SNP <- gsub("_hg19.*", "", biovu_tab$QUERY) 
biovu_tab$SNP[biovu_tab$SNP == "14.35761675-C-G"] <- "14:35761675-C-G"
biovu_tab$SNP[biovu_tab$SNP == "JHU_13.20119335_chr13"] <- "JHU_13.20119335"
biovu_tab$percent_identity <- as.numeric(gsub("%", "", biovu_tab$IDENTITY))

# max homology for 40bp match
max_40bp_homology <- with(subset(blat, Subject %in% c("chrX", "chrY") & abs(s.end-s.start)>=40), tapply(percent_identity, SNP, max))
snps$max_40bp_homology <- max_40bp_homology[match(snps$SNP, names(max_40bp_homology))]

# matches at 40bp, 85%?
matches_X <- with(subset(blat, Subject == "chrX" & percent_identity > 85), tapply(abs(s.end-s.start), SNP, max) >= 40)
matches_Y <- with(subset(blat, Subject == "chrY" & percent_identity > 85), tapply(abs(s.end-s.start), SNP, max) >= 40)
matches_X[is.na(matches_X)] <- FALSE
matches_Y[is.na(matches_Y)] <- FALSE
biovu_matches_X <- with(subset(biovu_tab, CHROM == "chrX" & percent_identity > 85), tapply(SPAN, SNP, max) >= 40)
biovu_matches_Y <- with(subset(biovu_tab, CHROM == "chrY" & percent_identity > 85), tapply(SPAN, SNP, max) >= 40)
snps$matches_X <- NA
snps$matches_Y <- NA
snps$matches_X[match(names(matches_X), snps$SNP)] <- matches_X
snps$matches_X[match(names(biovu_matches_X), snps$SNP)] <- (snps$matches_X[match(names(biovu_matches_X), snps$SNP)] | biovu_matches_X)
snps$matches_Y[match(names(matches_Y), snps$SNP)] <- matches_Y
snps$matches_Y[match(names(biovu_matches_Y), snps$SNP)] <- (snps$matches_Y[match(names(biovu_matches_Y), snps$SNP)] | biovu_matches_Y)


# Write this out

with(snps, 
    write.csv(data.frame(STUDY=DATASET,
                         CHR=CHR,
                         SNP=SNP,
                         FEM_HOM1=fem_hom1,
                         FEM_HET=fem_het,
                         FEM_HOM2=fem_hom2,
                         FEM_MISSING=fem_num_missing,
                         MALE_HOM1=malehom1,
                         MALE_HET=malehet,
                         MALE_HOM2=malehom2,
                         MALE_MISSING=malenum_missing,
                         FEM_EXP_HOM1=fem_exp_hom1,
                         MALE_EXP_HOM1=maleexp_hom1,
                         MISSING_PVAL=missing_pval,
                         HOM1_PVAL=hom_pval
                         # , MATCHES_X=matches_X, MATCHES_Y=matches_Y
                         ),
     file="snp_summaries.csv", row.names=FALSE))
