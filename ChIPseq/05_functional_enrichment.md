---
title: "Functional enrichment of gene near peaks"
author: Rylee K. Hackley
output: 
  html_document:
    keep_md: true
---

```r
# BiocManager::install(c("GenomicRanges","rtracklayer", "ChIPseeker", "IRanges"))
```


```r
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.2     ✔ readr     2.1.4
## ✔ forcats   1.0.0     ✔ stringr   1.5.0
## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
## ✔ purrr     1.0.1     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

```r
library(GenomicRanges)
```

```
## Loading required package: stats4
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:lubridate':
## 
##     intersect, setdiff, union
## 
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
## 
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
## 
## Loading required package: S4Vectors
## 
## Attaching package: 'S4Vectors'
## 
## The following objects are masked from 'package:lubridate':
## 
##     second, second<-
## 
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
## 
## The following object is masked from 'package:tidyr':
## 
##     expand
## 
## The following object is masked from 'package:utils':
## 
##     findMatches
## 
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
## 
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## 
## The following object is masked from 'package:lubridate':
## 
##     %within%
## 
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
## 
## The following object is masked from 'package:purrr':
## 
##     reduce
## 
## The following object is masked from 'package:grDevices':
## 
##     windows
## 
## Loading required package: GenomeInfoDb
```

```r
library(IRanges)
library(readxl)
library(rtracklayer)
library(viridis)
```

```
## Loading required package: viridisLite
```

```r
# Note: can also do cogtest by different ID, just replace "cogfile$locus_tag" argument with other type of identifier, e.g. "cogfile$acc")
# Calculate statistical information about the cogs represented in the list of genes using the function below.
# Cogtest: three inputs to this function:
# 1) list of gene names from the cluster (namelist)
# 2) COG file
# 3) p-value cutoff
# 4) method for multiple testing correction

cogtest <- function(namelist, cogfile, pvalue, method = "fdr", cutoff = 5) {
  cogs <- subset(cogfile, is.element(cogfile$locus_tag, unique(namelist$locus_tag)) == TRUE)
  clust <- table(cogs$COG_category) %>%
    subset(. > cutoff) # apply cutoff
  res <- data.frame(matrix(0, length(clust), 4)) # create 0 matrix
  rownames(res) <- names(clust)
  colnames(res) <- c("probability", "expect", "count", "p_adjust")
  all <- table(cogfile$COG_category)
  for (i in 1:length(clust)) { # calc expected frequencies and pval by hypergeo and append to DF
    all2 <- all[names(all) == names(clust[i])][[1]]
    res[i, 1] <- phyper(clust[[i]], all2, sum(all) - all2, nrow(cogs), lower.tail = F)
    res[i, 2] <- all2 * (nrow(cogs) / nrow(cogfile))
    res[i, 3] <- clust[[i]]
  }
  # multiple testing correction:
  res$p_adjust <- p.adjust(res$probability, method = method)
  fin <- subset(res, p_adjust <= pvalue)
  fin <- rownames_to_column(fin, var = "COG")
  return(fin)
}

# Note: the COGcategory name must be entered in quotes and add a space to the end of the category name of interest, e.g. 'transport '
## Use the following function to look at the genes in your cluster associated with a particular COG
cogset <- function(namelist, cogfile, COGcategory) {
  subset(cogfile, is.element(cogfile$locus_tag, namelist$locus_tag) & is.element(cogfile$COG_category, COGcategory) == TRUE)
}
```



```r
gff_df <- read_csv("../00_genome_files/GCF_000223905.1_gff.key.csv")[c(1:4, 10)]

# load arcogs, from EggNOG mapper
hca_cogs <- read_csv("../00_genome_files/HCA_eggNOG_mapper.csv")[c(1, 6:8)]
colnames(hca_cogs)[1] <- "old_locus_tag"
tmp <- inner_join(gff_df, hca_cogs)
colnames(hca_cogs)[1] <- "locus_tag"
tmp2 <- inner_join(gff_df, hca_cogs)
hca_cogs <- rbind(tmp, tmp2) %>% arrange(locus_tag)

# peaks
peaks <- read_csv("04a_peak_annotation/04a_consensus_genelist.csv")

# motifs
all_motif <- read_csv("04d_motif_annotation/motifs_annotated.csv")
peak_motif <- all_motif[all_motif$locus_tag %in% peaks$locus_tag,]

# look at the cogs in peaks
hca_cogs[hca_cogs$locus_tag %in% peaks$locus_tag, ]
```

```
## # A tibble: 25 × 8
##    chr         acc   locus_tag old_locus_tag annotation COG_category Description
##    <chr>       <chr> <chr>     <chr>         <chr>      <chr>        <chr>      
##  1 NC_015948.1 WP_0… HAH_RS00… HAH_0189      hypotheti… -            -          
##  2 NC_015948.1 WP_0… HAH_RS04… HAH_0887      cysteine … E            Cysteine s…
##  3 NC_015948.1 WP_0… HAH_RS04… HAH_1011      CBS domai… S            Signal tra…
##  4 NC_015948.1 WP_0… HAH_RS04… HAH_1012      Gfo/Idh/M… S            dehydrogen…
##  5 NC_015948.1 WP_0… HAH_RS06… HAH_1264      universal… T            COG0589 Un…
##  6 NC_015948.1 WP_0… HAH_RS06… HAH_1265      phosphoen… C            phosphoeno…
##  7 NC_015948.1 WP_0… HAH_RS06… HAH_1365      FAD-depen… C            COG1018 Fl…
##  8 NC_015948.1 WP_0… HAH_RS06… HAH_1366      2-oxoacid… C            COG1014 Py…
##  9 NC_015948.1 WP_0… HAH_RS06… HAH_1418      3-hydroxy… I            3-hydroxya…
## 10 NC_015948.1 WP_0… HAH_RS06… HAH_1419      class 1 f… G            D-fructose…
## # ℹ 15 more rows
## # ℹ 1 more variable: Preferred_name <chr>
```

HYPERGEO TESTS

```r
cogtest(all_motif, hca_cogs, 0.05) -> motif
# carbohydrate is most significant, but doesn't pass threshold for multiple testing

# use only motifs in promoter regions
all_motif %>% filter(type == "promoter") -> motif_pro
motif_pro <- motif_pro[, c(12, 13, 14)]
cogtest(motif_pro, hca_cogs, 0.05) -> pro_motif

length(unique(motif_pro$locus_tag))
```

```
## [1] 59
```

```r
# peaks with motifs
cogtest(peak_motif, hca_cogs, 0.05, cutoff = 1) -> pk_motif

motif$set <- rep("wg", nrow(motif))
pro_motif$set <- rep("promoter", nrow(pro_motif))
pk_motif$set <- rep("peaks", nrow(pk_motif))

rbind(motif, pro_motif, pk_motif)
```

```
##   COG  probability    expect count     p_adjust      set
## 1   G 2.929575e-03 5.5036067    12 3.515490e-02       wg
## 2   E 1.882910e-02 4.7969543     9 3.765820e-02 promoter
## 3   G 1.499577e-05 1.5228426     8 5.998307e-05 promoter
## 4   G 8.956013e-06 0.5343308     5 5.373608e-05    peaks
```

Look at involved genes

```r
cogset(all_motif, hca_cogs, "G")
```

```
## # A tibble: 12 × 8
##    chr         acc   locus_tag old_locus_tag annotation COG_category Description
##    <chr>       <chr> <chr>     <chr>         <chr>      <chr>        <chr>      
##  1 NC_015948.1 WP_0… HAH_RS04… HAH_0943      MFS trans… G            COG0477 Pe…
##  2 NC_015948.1 WP_0… HAH_RS05… HAH_1188      phosphope… G            Belongs to…
##  3 NC_015948.1 WP_0… HAH_RS06… HAH_1255      pyruvate … G            Belongs to…
##  4 NC_015948.1 WP_0… HAH_RS06… HAH_1419      class 1 f… G            D-fructose…
##  5 NC_015948.1 WP_0… HAH_RS07… HAH_1549      ABC trans… G            COG1653 AB…
##  6 NC_015948.1 WP_0… HAH_RS11… HAH_2323      phosphoen… G            Catalyzes …
##  7 NC_015948.1 WP_0… HAH_RS11… HAH_2420      MFS trans… G            COG0477 Pe…
##  8 NC_015948.1 WP_0… HAH_RS13… HAH_2730      type II g… G            Belongs to…
##  9 NC_015948.1 WP_0… HAH_RS13… HAH_2804      phosphogl… G            Belongs to…
## 10 NC_015943.1 WP_0… HAH_RS16… HAH_4332      glucose 1… G            Catalyzes …
## 11 NC_015944.1 WP_0… HAH_RS17… HAH_5129      sugar por… G            COG0477 Pe…
## 12 NC_015944.1 WP_0… HAH_RS18… HAH_5177      PTS fruct… G            COG1299 Ph…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(motif_pro, hca_cogs, "G")
```

```
## # A tibble: 8 × 8
##   chr         acc    locus_tag old_locus_tag annotation COG_category Description
##   <chr>       <chr>  <chr>     <chr>         <chr>      <chr>        <chr>      
## 1 NC_015948.1 WP_01… HAH_RS06… HAH_1255      pyruvate … G            Belongs to…
## 2 NC_015948.1 WP_01… HAH_RS06… HAH_1419      class 1 f… G            D-fructose…
## 3 NC_015948.1 WP_02… HAH_RS07… HAH_1549      ABC trans… G            COG1653 AB…
## 4 NC_015948.1 WP_01… HAH_RS11… HAH_2323      phosphoen… G            Catalyzes …
## 5 NC_015948.1 WP_00… HAH_RS13… HAH_2730      type II g… G            Belongs to…
## 6 NC_015943.1 WP_01… HAH_RS16… HAH_4332      glucose 1… G            Catalyzes …
## 7 NC_015944.1 WP_01… HAH_RS17… HAH_5129      sugar por… G            COG0477 Pe…
## 8 NC_015944.1 WP_01… HAH_RS18… HAH_5177      PTS fruct… G            COG1299 Ph…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(peaks, hca_cogs, "G")
```

```
## # A tibble: 5 × 8
##   chr         acc    locus_tag old_locus_tag annotation COG_category Description
##   <chr>       <chr>  <chr>     <chr>         <chr>      <chr>        <chr>      
## 1 NC_015948.1 WP_01… HAH_RS06… HAH_1419      class 1 f… G            D-fructose…
## 2 NC_015948.1 WP_01… HAH_RS11… HAH_2323      phosphoen… G            Catalyzes …
## 3 NC_015948.1 WP_00… HAH_RS13… HAH_2730      type II g… G            Belongs to…
## 4 NC_015943.1 WP_01… HAH_RS16… HAH_4332      glucose 1… G            Catalyzes …
## 5 NC_015944.1 WP_01… HAH_RS17… HAH_5129      sugar por… G            COG0477 Pe…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(all_motif, hca_cogs, "E")
```

```
## # A tibble: 17 × 8
##    chr         acc   locus_tag old_locus_tag annotation COG_category Description
##    <chr>       <chr> <chr>     <chr>         <chr>      <chr>        <chr>      
##  1 NC_015948.1 WP_0… HAH_RS00… HAH_0149      alpha/bet… E            hydrolases…
##  2 NC_015948.1 WP_0… HAH_RS02… HAH_0452      oligoendo… E            oligoendop…
##  3 NC_015948.1 WP_0… HAH_RS02… HAH_0543      anthranil… E            Anthranila…
##  4 NC_015948.1 WP_0… HAH_RS04… HAH_0887      cysteine … E            Cysteine s…
##  5 NC_015948.1 WP_0… HAH_RS05… HAH_1079      aldolase   E            COG1830 Dh…
##  6 NC_015948.1 WP_0… HAH_RS08… HAH_1825      threonine… E            COG1171 Th…
##  7 NC_015948.1 WP_0… HAH_RS09… HAH_1985      amino aci… E            COG0683 AB…
##  8 NC_015948.1 WP_0… HAH_RS10… HAH_2192      amino aci… E            amino acid 
##  9 NC_015948.1 WP_0… HAH_RS12… HAH_2535      ABC trans… E            COG0747 AB…
## 10 NC_015948.1 WP_0… HAH_RS13… HAH_2729      aminopept… E            Leucyl ami…
## 11 NC_015948.1 WP_0… HAH_RS14… HAH_2900      ABC trans… E            COG0683 AB…
## 12 NC_015948.1 WP_0… HAH_RS14… HAH_3040      substrate… E            COG0683 AB…
## 13 NC_015943.1 WP_0… HAH_RS15… HAH_4082      VOC famil… E            Glyoxalase…
## 14 NC_015944.1 WP_0… HAH_RS17… HAH_5034      Glu/Leu/P… E            Belongs to…
## 15 NC_015944.1 WP_0… HAH_RS17… HAH_5080      urea ABC … E            COG0683 AB…
## 16 NC_015944.1 WP_0… HAH_RS18… HAH_5182      FAD-depen… E            COG0665 Gl…
## 17 NC_015944.1 WP_0… HAH_RS18… HAH_5197      FAD-depen… E            COG0665 Gl…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(motif_pro, hca_cogs, "E")
```

```
## # A tibble: 9 × 8
##   chr         acc    locus_tag old_locus_tag annotation COG_category Description
##   <chr>       <chr>  <chr>     <chr>         <chr>      <chr>        <chr>      
## 1 NC_015948.1 WP_01… HAH_RS02… HAH_0543      anthranil… E            Anthranila…
## 2 NC_015948.1 WP_01… HAH_RS04… HAH_0887      cysteine … E            Cysteine s…
## 3 NC_015948.1 WP_00… HAH_RS05… HAH_1079      aldolase   E            COG1830 Dh…
## 4 NC_015948.1 WP_01… HAH_RS08… HAH_1825      threonine… E            COG1171 Th…
## 5 NC_015948.1 WP_01… HAH_RS13… HAH_2729      aminopept… E            Leucyl ami…
## 6 NC_015948.1 WP_02… HAH_RS14… HAH_3040      substrate… E            COG0683 AB…
## 7 NC_015944.1 WP_01… HAH_RS17… HAH_5034      Glu/Leu/P… E            Belongs to…
## 8 NC_015944.1 WP_01… HAH_RS17… HAH_5080      urea ABC … E            COG0683 AB…
## 9 NC_015944.1 WP_01… HAH_RS18… HAH_5197      FAD-depen… E            COG0665 Gl…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(peaks, hca_cogs, "E")
```

```
## # A tibble: 4 × 8
##   chr         acc    locus_tag old_locus_tag annotation COG_category Description
##   <chr>       <chr>  <chr>     <chr>         <chr>      <chr>        <chr>      
## 1 NC_015948.1 WP_01… HAH_RS04… HAH_0887      cysteine … E            Cysteine s…
## 2 NC_015948.1 WP_01… HAH_RS13… HAH_2729      aminopept… E            Leucyl ami…
## 3 NC_015948.1 WP_02… HAH_RS14… HAH_3040      substrate… E            COG0683 AB…
## 4 NC_015944.1 WP_01… HAH_RS17… HAH_5034      Glu/Leu/P… E            Belongs to…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(all_motif, hca_cogs, "C")
```

```
## # A tibble: 17 × 8
##    chr         acc   locus_tag old_locus_tag annotation COG_category Description
##    <chr>       <chr> <chr>     <chr>         <chr>      <chr>        <chr>      
##  1 NC_015948.1 WP_0… HAH_RS01… HAH_0414      electron … C            COG2086 El…
##  2 NC_015948.1 WP_0… HAH_RS03… HAH_0647      isocitrat… C            COG0538 Is…
##  3 NC_015948.1 WP_0… HAH_RS05… HAH_1160      ABC trans… C            COG1668 AB…
##  4 NC_015948.1 WP_0… HAH_RS06… HAH_1365      FAD-depen… C            COG1018 Fl…
##  5 NC_015948.1 WP_0… HAH_RS06… HAH_1366      2-oxoacid… C            COG1014 Py…
##  6 NC_015948.1 WP_0… HAH_RS06… HAH_1373      NAD(P)/FA… C            FAD depend…
##  7 NC_015948.1 WP_0… HAH_RS06… HAH_1379      cytochrom… C            COG1290 Cy…
##  8 NC_015948.1 WP_0… HAH_RS07… HAH_1558      aldo/keto… C            oxidoreduc…
##  9 NC_015948.1 WP_0… HAH_RS08… HAH_1756      b(o/a)3-t… C            COG0843 He…
## 10 NC_015948.1 WP_0… HAH_RS08… HAH_1798      molybdopt… C            Nitrate re…
## 11 NC_015948.1 WP_0… HAH_RS12… HAH_2526      NADH-quin… C            COG0852 NA…
## 12 NC_015948.1 WP_0… HAH_RS12… HAH_2611      geranylge… C            COG0644 De…
## 13 NC_015948.1 WP_0… HAH_RS14… HAH_3045      nitrate r… C            COG0243 An…
## 14 NC_015943.1 WP_0… HAH_RS16… HAH_4313      2-oxo aci… C            COG0508 Py…
## 15 NC_015943.1 WP_0… HAH_RS17… HAH_4415      4Fe-4S di… C            COG1140 Ni…
## 16 NC_015944.1 WP_1… HAH_RS18… HAH_5259      cbb3-type… C            Belongs to…
## 17 NC_015943.1 WP_2… HAH_RS20… <NA>          halocyani… C            Copper bin…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(motif_pro, hca_cogs, "C")
```

```
## # A tibble: 3 × 8
##   chr         acc    locus_tag old_locus_tag annotation COG_category Description
##   <chr>       <chr>  <chr>     <chr>         <chr>      <chr>        <chr>      
## 1 NC_015948.1 WP_04… HAH_RS06… HAH_1365      FAD-depen… C            COG1018 Fl…
## 2 NC_015948.1 WP_01… HAH_RS06… HAH_1366      2-oxoacid… C            COG1014 Py…
## 3 NC_015948.1 WP_01… HAH_RS12… HAH_2611      geranylge… C            COG0644 De…
## # ℹ 1 more variable: Preferred_name <chr>
```

```r
cogset(peaks, hca_cogs, "C")
```

```
## # A tibble: 3 × 8
##   chr         acc    locus_tag old_locus_tag annotation COG_category Description
##   <chr>       <chr>  <chr>     <chr>         <chr>      <chr>        <chr>      
## 1 NC_015948.1 WP_01… HAH_RS06… HAH_1265      phosphoen… C            phosphoeno…
## 2 NC_015948.1 WP_04… HAH_RS06… HAH_1365      FAD-depen… C            COG1018 Fl…
## 3 NC_015948.1 WP_01… HAH_RS06… HAH_1366      2-oxoacid… C            COG1014 Py…
## # ℹ 1 more variable: Preferred_name <chr>
```

create output file:

```r
# consensus peaks
cogtest(peak_motif, hca_cogs, 1, cutoff = 1) %>%
  mutate(
    logp = -log(p_adjust),
    cat = rep("trmB_consensus", nrow(.))
  ) -> a
a
```

```
##   COG  probability    expect count     p_adjust      logp            cat
## 1   - 8.309286e-01 4.2532728     2 8.309286e-01 0.1852114 trmB_consensus
## 2   C 9.608165e-02 1.1114080     2 1.441225e-01 1.9370918 trmB_consensus
## 3   E 2.205967e-02 1.6831419     4 6.617900e-02 2.7153920 trmB_consensus
## 4   G 8.956013e-06 0.5343308     5 5.373608e-05 9.8314260 trmB_consensus
## 5   S 6.835213e-01 3.3876570     2 8.202255e-01 0.1981759 trmB_consensus
## 6   T 4.216509e-02 0.7908095     2 8.433019e-02 2.4730154 trmB_consensus
```

```r
write_csv(a, "05_functional_enrichment//fxnal_categories.csv")
```
