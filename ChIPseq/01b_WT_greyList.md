Grey_list regions
================
Rylee K. Hackley

Make greylist to mask peaks that reproducibly occur in the no-tag
controls.

``` r
# BiocManager::install("GreyListChIP")
library("GreyListChIP")

ctrls <- read.csv("01b_HCA_chip_meta.csv") # update with correct PATH-TO/BAM/FILES
(filenames <- ctrls$bamReads)
```

    ## [1] "00_sorted_bams/4774_H3_S63_L004_R1_001_trimmed_sorted.bam"  
    ## [2] "00_sorted_bams/4774_A8_S96_L004_R1_001_trimmed_sorted.bam"  
    ## [3] "00_sorted_bams/4774_C5_S74_L004_R1_001_trimmed_sorted.bam"  
    ## [4] "00_sorted_bams/4774_D2_S51_L004_R1_001_trimmed_sorted.bam"  
    ## [5] "00_sorted_bams/4774_B4_S65_L004_R1_001_trimmed_sorted.bam"  
    ## [6] "00_sorted_bams/4774_A11_S120_L004_R1_001_trimmed_sorted.bam"
    ## [7] "00_sorted_bams/4774_E7_S92_L004_R1_001_trimmed_sorted.bam"  
    ## [8] "00_sorted_bams/4774_E10_S116_L004_R1_001_trimmed_sorted.bam"

Load in data

``` r
set.seed(123)

# chromosome information:
gl <- new("GreyList", karyoFile = "01b_greyList/HCA_karotype.txt")

for (i in 1:2) {
  gl <- countReads(gl, bamFile = filenames[i])
}

gl <- calcThreshold(gl, sampleSize = 3000) # see below for results with varying sample sizes
(gl <- makeGreyList(gl, maxGap = 150)) # use same bounds as will be used for peak calling
```

    ## coverage: 44543 bp (1.15%)

    ## GreyList on karyotype file HCA_karotype.txt
    ##   tiles: 7600
    ##   files: 00_sorted_bams/4774_H3_S63_L004_R1_001_trimmed_sorted.bam, 00_sorted_bams/4774_A8_S96_L004_R1_001_trimmed_sorted.bam
    ##   size (mean): 83.0697211909628
    ##   mu (mean): 473.784506667988
    ##   params: reps=100, sample size=3000, p-value=0.99
    ##   threshold: 614
    ##   regions: 26
    ##   coverage: 1.15%

``` r
export(gl, con = "01b_greyList/WT_greyList.bed")
```
