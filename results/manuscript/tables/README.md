Abin Abraham
2020-03-01 21:18:54
updated: 2020-03-20 09:59:02






Contains tables pertinent to the entire sexually antagonist selection in humans manuscript.



sup_table_ukbb_biovu_sig_snps_HWE_genotype_counts.xlsx
    * For GWAS significant SNPs in BIOVU and UKBB, provides HWE and genotype counts.

sup_table_x_gwas_signif_missingness_rate_ukbb_biovu.xlsx
    * For GWAS significant SNPs in BIOVU and UKBB, provides p-value for different missingness rate between males and females.


raw_bv_webblat_xy_df.tsv
raw_ukbil_webblat_xy_df.tsv
raw_ukwcsf_webblat_xy_df.tsv
    * raw BLAT HITS to X or Y in UK Biobank (ukbil and ukwcsf are the two arrays) and BioVU (MEGA array) before any processing.
    * script:/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits/plot_blat_panel.py

bv_best_blatscore_xy_hit_and_filtered.tsv
uk_best_blatscore_xy_hit_and_filtered.tsv
    * GWAS variants and their best BLAT hit to X or Y based on BLAT score and filtering (≥90% seqID, ≥40bp length, overlap variant)
    * script:/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits/plot_blat_panel.py
    * script:/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits/plot_blat_panel.py


uk_bv_gwas_sig_hits_w_best_blat_score_xy_match.tsv
    * a subset of *bv_best_blatscore_xy_hit_and_filtered.tsv to include only gwas significant varaints
    * script:/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/results/manuscript/blat_hits/plot_blat_panel.py



ukbb_sig_snps_missing_bysex.csv
    * contains UKBB 'significant variants' with a statistically significant difference in missingness between males and females


uk_var_w_missingness_best_blatscore_xy.tsv
    *best blat xy hit based on score for uk variants with stat. sig. differnece in missing rate between males and females
    * script: check_ukvars_w_missingness.py