---
title: "Differential expression analysis"
author: "Rylee K. Hackley"
output: 
  html_document:
    keep_md: true
---





load ggf data info for informative gene descriptions and counts files

```
## [1] 3863   29
```

### Differential Expression Analysis

```r
dds_mydata <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = counts_meta,
  design = ~ batch + group
)

# Estimate Size Factors to inspect the raw count data which is generated from various growth experiments, flowcells and bio reps.
dds_mydata <- estimateSizeFactors(dds_mydata)
ddsDE <- DESeq(dds_mydata, test = "Wald")

# Estimate of sample quality; a plot of per-gene dispersion estimates together with the fitted mean-dispersion relationship
plotDispEsts(ddsDE)
```

![](01_DeSeq2_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

#### PCA

```r
# The variance stabilizing transformation is obtained using the vst() function.
vsd <- vst(ddsDE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)

# Principal Components Analysis
plotPCA(vsd, intgroup = c("genotype", "glucose"), returnData = T) -> PCA
(cbind(PCA, "pyrF"= counts_meta$pyrF, "experiment"= str_sub(counts_meta$batch, start= -4)) -> PCA)
##                                 PC1         PC2     group genotype glucose
## 10_WT_B_S18              -10.438954 -18.0192675   WT:mGlu       WT    mGlu
## 11_WT_C_S19                6.828975 -15.7542289   WT:mGlu       WT    mGlu
## 12_WT_D_S20              -11.723334 -11.6963856   WT:mGlu       WT    mGlu
## 13_TRMB_A_S21            -29.293681  11.8539029 trmB:mGlu     trmB    mGlu
## 14_TRMB_B_S22            -28.624264  11.3852119 trmB:mGlu     trmB    mGlu
## 15_TRMB_C_S23            -30.795449  12.5220458 trmB:mGlu     trmB    mGlu
## 16_TRMB_D_S24            -27.544185  14.5314203 trmB:mGlu     trmB    mGlu
## 1_WT_A_S9                 18.981809  -0.6555241   WT:pGlu       WT    pGlu
## 2_WT_B_S10                15.434573  -0.1242029   WT:pGlu       WT    pGlu
## 3_WT_C_S11                19.127755   2.9230367   WT:pGlu       WT    pGlu
## 4_WT_D_S12                13.301701  -2.3744928   WT:pGlu       WT    pGlu
## 5_TRMB_A_S13              17.579184   3.3102823 trmB:pGlu     trmB    pGlu
## 6_TRMB_B_S14              14.948109   1.1031914 trmB:pGlu     trmB    pGlu
## 7_TRMB_C_S15              20.610970   4.7440739 trmB:pGlu     trmB    pGlu
## 8_TRMB_D_S16              19.685553   3.2474883 trmB:pGlu     trmB    pGlu
## 9_WT_A_S17                -8.078762 -16.9965516   WT:mGlu       WT    mGlu
## aBWT_C_noglu_24_S20_SP   -13.494908 -12.4118297   WT:mGlu       WT    mGlu
## aCWT_D_noglu_24_S21_SP   -25.033075 -16.8748503   WT:mGlu       WT    mGlu
## aFTrmB_C_noglu_24_S48_S1 -42.716396   8.6421278 trmB:mGlu     trmB    mGlu
## aGTrmB_D_noglu_24_S42_S1 -44.073413  10.4828360 trmB:mGlu     trmB    mGlu
## RWT_A_glu_24_S41_S1       21.233823   0.2493970   WT:pGlu       WT    pGlu
## SWT_B_glu_24_S16_SP       12.158238  10.0171725   WT:pGlu       WT    pGlu
## TWT_C_glu_24_S17_SP       13.392279   7.0281015   WT:pGlu       WT    pGlu
## UWT_D_glu_24_S46_S1       20.475351  -1.7718851   WT:pGlu       WT    pGlu
## VTrmB_A_glu_24_S25_S1     22.420761   2.9081360 trmB:pGlu     trmB    pGlu
## WTrmB_B_glu_24_S11_S1     20.607003   3.9376374 trmB:pGlu     trmB    pGlu
## XTrmB_C_glu_24_S26_S1     12.977465  12.2414060 trmB:pGlu     trmB    pGlu
## YTrmB_D_glu_24_S18_SP     21.044465   3.6095038 trmB:pGlu     trmB    pGlu
## ZWT_A_noglu_24_S47_S1    -18.991593 -28.0577529   WT:mGlu       WT    mGlu
##                                              name  pyrF experiment
## 10_WT_B_S18                           10_WT_B_S18 FALSE       2022
## 11_WT_C_S19                           11_WT_C_S19 FALSE       2022
## 12_WT_D_S20                           12_WT_D_S20 FALSE       2022
## 13_TRMB_A_S21                       13_TRMB_A_S21 FALSE       2022
## 14_TRMB_B_S22                       14_TRMB_B_S22  TRUE       2022
## 15_TRMB_C_S23                       15_TRMB_C_S23  TRUE       2022
## 16_TRMB_D_S24                       16_TRMB_D_S24 FALSE       2022
## 1_WT_A_S9                               1_WT_A_S9 FALSE       2022
## 2_WT_B_S10                             2_WT_B_S10 FALSE       2022
## 3_WT_C_S11                             3_WT_C_S11 FALSE       2022
## 4_WT_D_S12                             4_WT_D_S12 FALSE       2022
## 5_TRMB_A_S13                         5_TRMB_A_S13 FALSE       2022
## 6_TRMB_B_S14                         6_TRMB_B_S14  TRUE       2022
## 7_TRMB_C_S15                         7_TRMB_C_S15 FALSE       2022
## 8_TRMB_D_S16                         8_TRMB_D_S16  TRUE       2022
## 9_WT_A_S17                             9_WT_A_S17 FALSE       2022
## aBWT_C_noglu_24_S20_SP     aBWT_C_noglu_24_S20_SP FALSE       2021
## aCWT_D_noglu_24_S21_SP     aCWT_D_noglu_24_S21_SP FALSE       2021
## aFTrmB_C_noglu_24_S48_S1 aFTrmB_C_noglu_24_S48_S1  TRUE       2021
## aGTrmB_D_noglu_24_S42_S1 aGTrmB_D_noglu_24_S42_S1  TRUE       2021
## RWT_A_glu_24_S41_S1           RWT_A_glu_24_S41_S1 FALSE       2021
## SWT_B_glu_24_S16_SP           SWT_B_glu_24_S16_SP FALSE       2021
## TWT_C_glu_24_S17_SP           TWT_C_glu_24_S17_SP FALSE       2021
## UWT_D_glu_24_S46_S1           UWT_D_glu_24_S46_S1 FALSE       2021
## VTrmB_A_glu_24_S25_S1       VTrmB_A_glu_24_S25_S1  TRUE       2021
## WTrmB_B_glu_24_S11_S1       WTrmB_B_glu_24_S11_S1  TRUE       2021
## XTrmB_C_glu_24_S26_S1       XTrmB_C_glu_24_S26_S1  TRUE       2021
## YTrmB_D_glu_24_S18_SP       YTrmB_D_glu_24_S18_SP  TRUE       2021
## ZWT_A_noglu_24_S47_S1       ZWT_A_noglu_24_S47_S1 FALSE       2021
ggplot(data = PCA) +
  geom_point(aes(x = PC1, y = PC2, fill = group, color = pyrF), size = 7, shape = 21, stroke = 2.5)+
  scale_color_manual(values = c("grey", "black")) + scale_fill_viridis(discrete = T)+
  theme_bw()+ theme(text = element_text(size = 18)) -> pca.plot
pca.plot
```

![](01_DeSeq2_files/figure-html/clustering-1.png)<!-- -->

### Differential Expression Analysis (no pyrF revertants)

```r
counts_meta2 <- counts_meta[counts_meta$pyrF == F,]
counts_data2 <- counts_data[counts_meta2$sample]

dds.mydata2 <- DESeqDataSetFromMatrix(
  countData = counts_data2,
  colData = counts_meta2,
  design = ~ batch + group
)
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
## design formula are characters, converting to factors

## Estimate Size Factors to inspect the raw count data which is generated from various lanes and bio reps.
dds.mydata2 <- estimateSizeFactors(dds.mydata2)
ddsDE2 <- DESeq(dds.mydata2, test = "Wald")
## using pre-existing size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 61 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testing
```

### Save normalized counts

```r
# Convert the raw count to the normalized count
normalized_counts <- as.data.frame(counts(ddsDE, normalized = TRUE)) %>% rownames_to_column()
write_csv(normalized_counts, file = "01_deseq2_output/normalised_counts_all.csv")

normalized_counts <- as.data.frame(counts(ddsDE2, normalized = TRUE)) %>% rownames_to_column()
write_csv(normalized_counts, file = "01_deseq2_output/normalised_counts_filtered.csv")
```

### What are the DEGs that are unique to pyrF revertant samples? What genes are DE across both analyses?

```r
all <- results(ddsDE, name = "group_independent_vs_dependent", lfcThreshold = 1, alpha = 0.05)
all <- all %>% # Make a result table
  data.frame() %>%
  rownames_to_column(var = "locus_tag") %>%
  as_tibble()
all <- left_join(all, gff_df[c(3:4, 10)], by = "locus_tag")
all_sig <- all %>% filter(padj < 0.05)

no_revert <- results(ddsDE2, name = "group_independent_vs_dependent", lfcThreshold = 1, alpha = 0.05)
no_revert_all <- no_revert %>% # Make a result table
  data.frame() %>%
  rownames_to_column(var = "locus_tag") %>%
  as_tibble()
no_revert_all <- left_join(no_revert_all, gff_df[c(3:4, 10)], by = "locus_tag")
no_revert_sig <- no_revert_all %>% filter(padj < 0.05)

#export DEGs from both analyses
write_csv(all, "01_deseq2_output/group_all.csv")
write_csv(no_revert_all, "01_deseq2_output/group_filtered.csv")
```



```r
# Genes that are unique to pryF
(all_sig[!(all_sig$locus_tag %in% no_revert_sig$locus_tag), ] -> tmp) # 1 gene: pyrF
## # A tibble: 1 × 9
##   locus_tag   baseMean log2FoldChange lfcSE  stat   pvalue    padj old_locus_tag
##   <chr>          <dbl>          <dbl> <dbl> <dbl>    <dbl>   <dbl> <chr>        
## 1 HAH_RS10105    3182.           6.05 0.998  5.06  4.14e-7 7.85e-5 HAH_2085     
## # ℹ 1 more variable: annotation <chr>

# Genes that are unique to non-revertants
(no_revert_sig[!(no_revert_sig$locus_tag %in% all_sig$locus_tag), ] -> tmp)  # 13 genes: most in HAH_4316-21 (ABC transporter for branched chain amino acids?)
## # A tibble: 13 × 9
##    locus_tag   baseMean log2FoldChange lfcSE  stat  pvalue    padj old_locus_tag
##    <chr>          <dbl>          <dbl> <dbl> <dbl>   <dbl>   <dbl> <chr>        
##  1 HAH_RS04570     925.          -1.77 0.203 -3.81 1.38e-4 1.40e-2 HAH_0941     
##  2 HAH_RS06045    7195.          -2.81 0.424 -4.27 1.95e-5 2.28e-3 HAH_1242     
##  3 HAH_RS06505    1678.          -2.40 0.300 -4.68 2.84e-6 3.78e-4 HAH_1336     
##  4 HAH_RS12385    2229.           1.48 0.138  3.49 4.87e-4 4.36e-2 HAH_2550     
##  5 HAH_RS16590     289.          -2.71 0.407 -4.19 2.84e-5 3.12e-3 HAH_4316     
##  6 HAH_RS16595     325.          -3.02 0.453 -4.47 7.91e-6 1.01e-3 HAH_4317     
##  7 HAH_RS16600     517.          -3.02 0.482 -4.20 2.70e-5 3.05e-3 HAH_4318     
##  8 HAH_RS16605     588.          -2.78 0.476 -3.73 1.92e-4 1.80e-2 HAH_4319     
##  9 HAH_RS16610    5543.          -3.18 0.458 -4.77 1.84e-6 2.53e-4 HAH_4320     
## 10 HAH_RS16615    1527.          -3.14 0.404 -5.29 1.19e-7 2.87e-5 HAH_4321     
## 11 HAH_RS16730     255.          -3.19 0.584 -3.74 1.81e-4 1.74e-2 HAH_4344     
## 12 HAH_RS16750     253.          -2.73 0.335 -5.17 2.29e-7 4.91e-5 HAH_4348     
## 13 HAH_RS18730     284.          -2.14 0.291 -3.93 8.50e-5 9.09e-3 HAH_5293     
## # ℹ 1 more variable: annotation <chr>

# comparing all
tmp <- inner_join(all[c(1,3,8:9)], no_revert_all[c(1,3)], by = "locus_tag", suffix = c(".pyr", ".nopyr")) 

#plot correlations
ggscatter(tmp, 
  x = "log2FoldChange.pyr", y = "log2FoldChange.nopyr",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson", size = 2.5
) +
  labs(y = "all samples", x = "no pyrF samples") +
  theme_bw() +
  theme(text = element_text(size = 18))
## Warning: Removed 1 rows containing non-finite values (`stat_smooth()`).
## Warning: Removed 1 rows containing non-finite values (`stat_cor()`).
## Warning: Removed 1 rows containing missing values (`geom_point()`).
```

![](01_DeSeq2_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r

# of the genes that are shared, how well does expression correlate?
(tmp <- inner_join(all_sig[c(1:3,8:9)], no_revert_sig[c(1:3)], by = "locus_tag", suffix = c(".pyr", ".nopyr")))
## # A tibble: 32 × 7
##    locus_tag   baseMean.pyr log2FoldChange.pyr old_locus_tag annotation         
##    <chr>              <dbl>              <dbl> <chr>         <chr>              
##  1 HAH_RS19805        4565.              -4.63 HAH_0188      hypothetical prote…
##  2 HAH_RS02790         456.              -1.48 HAH_0576      hypothetical prote…
##  3 HAH_RS04315         367.              -1.36 HAH_0887      cysteine synthase A
##  4 HAH_RS04465       55413.               3.33 HAH_0919      glutamate synthase…
##  5 HAH_RS04485        1207.              -1.61 HAH_0923      transcriptional re…
##  6 HAH_RS05300         544.              -2.22 HAH_1090      TRAP transporter f…
##  7 HAH_RS05305         449.              -1.83 HAH_1091      TAXI family TRAP t…
##  8 HAH_RS06030       15966.              -3.05 HAH_1239      Glu/Leu/Phe/Val de…
##  9 HAH_RS06035        9977.              -2.28 HAH_1240      arginase           
## 10 HAH_RS06050        2486.              -2.87 HAH_1243      amino acid ABC tra…
## # ℹ 22 more rows
## # ℹ 2 more variables: baseMean.nopyr <dbl>, log2FoldChange.nopyr <dbl>

#save DEG list
write_csv(tmp, "01_deseq2_output/group_final.csv")

ggscatter(tmp, 
  x = "log2FoldChange.pyr", y = "log2FoldChange.nopyr",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson", size = 2.5
) +
  labs(y = "all samples", x = "no pyrF samples") +
  theme_bw() +
  theme(text = element_text(size = 18)) -> correlation
correlation
```

![](01_DeSeq2_files/figure-html/unnamed-chunk-7-2.png)<!-- -->
