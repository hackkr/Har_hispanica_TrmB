---
title: "Combine experiments for count matrix"
author: Rylee K. Hackley
output: 
  html_document:
    keep_md: true
---


```r
library(tidyverse)
library(magrittr)
library(purrr)
library(DESeq2)
```

#wrangling (combining) counts matrices
data from 2022

```r
# read in featureCounts output for new data
count_mat_new <- read_tsv("000_counts_2022/HCA_featureCounts_2022.txt", skip = 1)
count_mat2022 <- as.matrix(count_mat_new[-c(1:6)]) # remove extra columns
rownames(count_mat2022) <- count_mat_new$Geneid # set rownames

# read in meta
meta2022 <- read_csv("000_counts_2022/000_meta.csv")

# check colnames with metadata rownames
colnames(count_mat2022) <- str_remove(colnames(count_mat2022), "03_aligned_sorted/")
colnames(count_mat2022) <- str_remove(colnames(count_mat2022), "_sorted.bam")
colnames(count_mat2022) == meta2022$sample
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [16] TRUE TRUE

counts_col22 <- meta2022
counts_col22 <- as.data.frame(count_mat2022)

# add pyrF column to metadata
counts_col22["HAH_RS10105",] %>%
  t() -> tmp
tmp; as.vector(tmp > 67)
##               HAH_RS10105
## 10_WT_B_S18            43
## 11_WT_C_S19            34
## 12_WT_D_S20            32
## 13_TRMB_A_S21          17
## 14_TRMB_B_S22        5214
## 15_TRMB_C_S23        5077
## 16_TRMB_D_S24          36
## 17_WT_ref_S25        4589
## 1_WT_A_S9              67
## 2_WT_B_S10             60
## 3_WT_C_S11             46
## 4_WT_D_S12             38
## 5_TRMB_A_S13           35
## 6_TRMB_B_S14        11627
## 7_TRMB_C_S15           50
## 8_TRMB_D_S16        11180
## 9_WT_A_S17             62
##  [1] FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE
## [13] FALSE  TRUE FALSE  TRUE FALSE
meta2022$pyrF <- as.vector(tmp > 67) #max reads mapping to pyrF locus allowed
```

data from 2021

```r
# read in featureCounts output for 2021 data
count_mat2021 <- read_tsv("000_counts_2021/HCA_featureCounts_2021.txt", skip = 1)
colnames(count_mat2021) <- str_remove(colnames(count_mat2021), "bams/")
colnames(count_mat2021) <- str_remove(colnames(count_mat2021), "_sorted.bam")
counts2021 <- count_mat2021

# read in meta
meta2021 <- read_csv("000_counts_2021/000_meta.csv")

#remove rows from data that are not in metadata (48hr samples)
counts2021 <- dplyr::select(counts2021, c("Geneid", meta2021$sample))

# check colnames with metadata rownames
colnames(counts2021) %in% meta2021$sample
##  [1] FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [13]  TRUE  TRUE  TRUE  TRUE

# add pyrF column
counts2021[counts2021$Geneid == "HAH_RS10105",][-1] %>%
  t() -> tmp
tmp; as.vector(tmp > 153)
##                           [,1]
## aAWT_B_noglu_24_S19_SP      29
## aBWT_C_noglu_24_S20_SP      80
## aCWT_D_noglu_24_S21_SP      63
## aFTrmB_C_noglu_24_S48_S1  7088
## aGTrmB_D_noglu_24_S42_S1  1928
## HCA_ref_RNA_S7_S1         5294
## RWT_A_glu_24_S41_S1         95
## SWT_B_glu_24_S16_SP        116
## TWT_C_glu_24_S17_SP        153
## UWT_D_glu_24_S46_S1        122
## VTrmB_A_glu_24_S25_S1    11534
## WTrmB_B_glu_24_S11_S1    18022
## XTrmB_C_glu_24_S26_S1     5632
## YTrmB_D_glu_24_S18_SP    20685
## ZWT_A_noglu_24_S47_S1       58
##  [1] FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE
## [13]  TRUE  TRUE FALSE
meta2021$pyrF <- as.vector(tmp > 153) #max reads mapping to pyrF locus allowed
```

combine 2 experiments

```r
meta <- rbind(meta2022, meta2021)
counts <- cbind(count_mat2022, counts2021[-1])

# calc total reads per sample
data.frame("total counts" = colSums(as.matrix(counts))) -> cnt_sums
cnt_sums <- rownames_to_column(cnt_sums, "sample_name")

counts <- rownames_to_column(counts, "locus_tag")
```


```r
# export combined counts
write_csv(counts, "00_deseq2_input/00_combined_counts.csv")
write_csv(meta, "00_deseq2_input/00_combined_meta.csv")
write_csv(cnt_sums, "00_deseq2_input/00_counts_sums.csv")
```
