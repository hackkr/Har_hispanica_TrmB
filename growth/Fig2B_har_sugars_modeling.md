---
title: "Haloarcula hispanica TrmB (AKS133) growth in various carbon sources"
author: Rylee K. Hackley
output: 
  html_document:
    keep_md: true
---

```r
# load required libraries
library(growthcurver)
library(viridis)
library(scales)
library(tidyverse)
library(factoextra)
library(ggpubr)
library(rstatix)

# custom function, calculated the 95% Confidence interval
conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df = n - 1) * sd(data) / sqrt(n)
  return(error)
}
```

# Haloarcula hispanica
Load the normalized data and meta data. combine difference experiments. 

```r
expt <- "hca_sugars"

## read in files
mt1 <- read.csv("normalized_growth_data/2020_sugar1_normalized_layout.csv", stringsAsFactors = F)
dt1 <- read.csv("normalized_growth_data/2020_sugar1_normalized_data.csv")
mt2 <- read.csv("normalized_growth_data/2020_sugar2_normalized_layout.csv", stringsAsFactors = F)
dt2 <- read.csv("normalized_growth_data/2020_sugar2_normalized_data.csv")
mt3 <- read.csv("normalized_growth_data/2021_sugar3_normalized_layout.csv", stringsAsFactors = F)
dt3 <- read.csv("normalized_growth_data/2021_sugar3_normalized_data.csv")

## set well no as column names, force conseq numbering
colnames(dt1) <- c(seq(1, ncol(dt1) - 1), "time")
mt1$variable <- seq(1, ncol(dt1) - 1)

colnames(dt2) <- c(seq(1, ncol(dt2) - 1) + max(mt1$variable), "time")
mt2$variable <- seq(1, ncol(dt2) - 1) + max(mt1$variable)

colnames(dt3) <- c(seq(1, ncol(dt3) - 1) + max(mt2$variable), "time")
mt3$variable <- seq(1, ncol(dt3) - 1) + max(mt2$variable)

## combine meta and data frames across diff experiments
d <- inner_join(dt1, dt2, by = "time")
d <- inner_join(d, dt3, by = "time")
mt0 <- rbind(mt1, mt2, mt3)

# combine expt info with biorep info and remove expt column
mt0$biorep <- paste(mt0$biorep, mt0$expt, sep = "-")

## replace none with "no carb" for labeling
mt0$ID <- str_replace_all(mt0$ID, "none", "no carb")
mt0$condition <- str_replace_all(mt0$condition, "none", "no carb")
mt0
```

```
##     variable                  ID strain media      condition        biorep
## 1          1 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2020_sugar1
## 2          2      WT+25mM ribose     WT  HhCa    25mM ribose 1-2020_sugar1
## 3          3    WT+25mM fructose     WT  HhCa  25mM fructose 2-2020_sugar1
## 4          4      WT+25mM xylose     WT  HhCa    25mM xylose 1-2020_sugar1
## 5          5      WT+25mM xylose     WT  HhCa    25mM xylose 2-2020_sugar1
## 6          6  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2020_sugar1
## 7          7    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2020_sugar1
## 8          8 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2020_sugar1
## 9          9    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2020_sugar1
## 10        10      WT+25mM ribose     WT  HhCa    25mM ribose 1-2020_sugar1
## 11        11        trmB+no carb   trmB  HhCa        no carb 1-2020_sugar1
## 12        12  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2020_sugar1
## 13        13   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2020_sugar1
## 14        14 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2020_sugar1
## 15        15  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2020_sugar1
## 16        16   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2020_sugar1
## 17        17   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2020_sugar1
## 18        18   WT+25mM galactose     WT  HhCa 25mM galactose 2-2020_sugar1
## 19        19      WT+25mM xylose     WT  HhCa    25mM xylose 1-2020_sugar1
## 20        20  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2020_sugar1
## 21        21  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2020_sugar1
## 22        22   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2020_sugar1
## 23        23   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2020_sugar1
## 24        24   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2020_sugar1
## 25        25    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2020_sugar1
## 26        26      WT+25mM ribose     WT  HhCa    25mM ribose 1-2020_sugar1
## 27        27          WT+no carb     WT  HhCa        no carb 2-2020_sugar1
## 28        28      WT+25mM xylose     WT  HhCa    25mM xylose 1-2020_sugar1
## 29        29    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2020_sugar1
## 30        30          WT+no carb     WT  HhCa        no carb 1-2020_sugar1
## 31        31     WT+25mM acetate     WT  HhCa   25mM acetate 2-2020_sugar1
## 32        32   WT+25mM galactose     WT  HhCa 25mM galactose 1-2020_sugar1
## 33        33          WT+no carb     WT  HhCa        no carb 1-2020_sugar1
## 34        34  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2020_sugar1
## 35        35     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2020_sugar1
## 36        36    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2020_sugar1
## 37        37    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2020_sugar1
## 38        38     WT+25mM acetate     WT  HhCa   25mM acetate 2-2020_sugar1
## 39        39    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2020_sugar1
## 40        40        trmB+no carb   trmB  HhCa        no carb 1-2020_sugar1
## 41        41   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2020_sugar1
## 42        42    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2020_sugar1
## 43        43        trmB+no carb   trmB  HhCa        no carb 2-2020_sugar1
## 44        44    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2020_sugar1
## 45        45          WT+no carb     WT  HhCa        no carb 2-2020_sugar1
## 46        46  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2020_sugar1
## 47        47 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2020_sugar1
## 48        48   WT+25mM galactose     WT  HhCa 25mM galactose 2-2020_sugar1
## 49        49     WT+25mM glucose     WT  HhCa   25mM glucose 2-2020_sugar1
## 50        50   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2020_sugar1
## 51        51      WT+25mM xylose     WT  HhCa    25mM xylose 2-2020_sugar1
## 52        52  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2020_sugar1
## 53        53    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2020_sugar1
## 54        54   WT+25mM galactose     WT  HhCa 25mM galactose 2-2020_sugar1
## 55        55        trmB+no carb   trmB  HhCa        no carb 2-2020_sugar1
## 56        56     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2020_sugar1
## 57        57     WT+25mM acetate     WT  HhCa   25mM acetate 2-2020_sugar1
## 58        58     WT+25mM glucose     WT  HhCa   25mM glucose 1-2020_sugar1
## 59        59 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2020_sugar1
## 60        60   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2020_sugar1
## 61        61   WT+25mM galactose     WT  HhCa 25mM galactose 1-2020_sugar1
## 62        62    WT+25mM fructose     WT  HhCa  25mM fructose 1-2020_sugar1
## 63        63      WT+25mM ribose     WT  HhCa    25mM ribose 2-2020_sugar1
## 64        64     WT+25mM glucose     WT  HhCa   25mM glucose 2-2020_sugar1
## 65        65   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2020_sugar1
## 66        66    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2020_sugar1
## 67        67      WT+25mM xylose     WT  HhCa    25mM xylose 2-2020_sugar1
## 68        68  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2020_sugar1
## 69        69          WT+no carb     WT  HhCa        no carb 2-2020_sugar1
## 70        70   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2020_sugar1
## 71        71   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2020_sugar1
## 72        72     WT+25mM acetate     WT  HhCa   25mM acetate 1-2020_sugar1
## 73        73     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2020_sugar1
## 74        74  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2020_sugar1
## 75        75        trmB+no carb   trmB  HhCa        no carb 1-2020_sugar1
## 76        76    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2020_sugar1
## 77        77     WT+25mM glucose     WT  HhCa   25mM glucose 2-2020_sugar1
## 78        78  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2020_sugar1
## 79        79    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2020_sugar1
## 80        80   WT+25mM galactose     WT  HhCa 25mM galactose 1-2020_sugar1
## 81        81    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2020_sugar1
## 82        82      WT+25mM ribose     WT  HhCa    25mM ribose 2-2020_sugar1
## 83        83    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2020_sugar1
## 84        84    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2020_sugar1
## 85        85     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2020_sugar1
## 86        86    WT+25mM fructose     WT  HhCa  25mM fructose 2-2020_sugar1
## 87        87    WT+25mM fructose     WT  HhCa  25mM fructose 1-2020_sugar1
## 88        88    WT+25mM fructose     WT  HhCa  25mM fructose 1-2020_sugar1
## 89        89   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2020_sugar1
## 90        90    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2020_sugar1
## 91        91     WT+25mM acetate     WT  HhCa   25mM acetate 1-2020_sugar1
## 92        92  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2020_sugar1
## 93        93  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2020_sugar1
## 94        94  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2020_sugar1
## 95        95  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2020_sugar1
## 96        96  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2020_sugar1
## 97        97    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2020_sugar1
## 98        98     WT+25mM glucose     WT  HhCa   25mM glucose 1-2020_sugar1
## 99        99     WT+25mM acetate     WT  HhCa   25mM acetate 1-2020_sugar1
## 100      100    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2020_sugar1
## 101      101   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2020_sugar1
## 102      102   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2020_sugar1
## 103      103   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2020_sugar1
## 104      104   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2020_sugar1
## 105      105     WT+25mM glucose     WT  HhCa   25mM glucose 1-2020_sugar1
## 106      106  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2020_sugar1
## 107      107    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2020_sugar1
## 108      108      WT+25mM ribose     WT  HhCa    25mM ribose 2-2020_sugar1
## 109      109          WT+no carb     WT  HhCa        no carb 1-2020_sugar1
## 110      110        trmB+no carb   trmB  HhCa        no carb 2-2020_sugar1
## 111      111    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2020_sugar1
## 112      112 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2020_sugar1
## 113      113  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2020_sugar1
## 114      114    WT+25mM fructose     WT  HhCa  25mM fructose 2-2020_sugar1
## 115      115     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2020_sugar1
## 116      116     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2020_sugar1
## 117      117   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2020_sugar1
## 118      118    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2020_sugar1
## 119      119    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2020_sugar1
## 120      120        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar1
## 121      121        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar1
## 122      122        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar1
## 123      123        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar1
## 124      124          WT+no carb     WT  HhCa        no carb 3-2020_sugar1
## 125      125          WT+no carb     WT  HhCa        no carb 3-2020_sugar1
## 126      126          WT+no carb     WT  HhCa        no carb 3-2020_sugar1
## 127      127          WT+no carb     WT  HhCa        no carb 3-2020_sugar1
## 128      128    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2020_sugar2
## 129      129          WT+no carb     WT  HhCa        no carb 2-2020_sugar2
## 130      130  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2020_sugar2
## 131      131    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2020_sugar2
## 132      132    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2020_sugar2
## 133      133     WT+25mM glucose     WT  HhCa   25mM glucose 1-2020_sugar2
## 134      134        trmB+no carb   trmB  HhCa        no carb 1-2020_sugar2
## 135      135      WT+25mM xylose     WT  HhCa    25mM xylose 2-2020_sugar2
## 136      136      WT+25mM ribose     WT  HhCa    25mM ribose 2-2020_sugar2
## 137      137    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2020_sugar2
## 138      138  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2020_sugar2
## 139      139    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2020_sugar2
## 140      140          WT+no carb     WT  HhCa        no carb 1-2020_sugar2
## 141      141      WT+25mM xylose     WT  HhCa    25mM xylose 2-2020_sugar2
## 142      142     WT+25mM glucose     WT  HhCa   25mM glucose 1-2020_sugar2
## 143      143    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2020_sugar2
## 144      144     WT+25mM glucose     WT  HhCa   25mM glucose 2-2020_sugar2
## 145      145   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2020_sugar2
## 146      146     WT+25mM acetate     WT  HhCa   25mM acetate 1-2020_sugar2
## 147      147          WT+no carb     WT  HhCa        no carb 2-2020_sugar2
## 148      148  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2020_sugar2
## 149      149    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2020_sugar2
## 150      150          WT+no carb     WT  HhCa        no carb 1-2020_sugar2
## 151      151     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2020_sugar2
## 152      152      WT+25mM xylose     WT  HhCa    25mM xylose 1-2020_sugar2
## 153      153  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2020_sugar2
## 154      154        trmB+no carb   trmB  HhCa        no carb 2-2020_sugar2
## 155      155   WT+25mM galactose     WT  HhCa 25mM galactose 2-2020_sugar2
## 156      156    WT+25mM fructose     WT  HhCa  25mM fructose 2-2020_sugar2
## 157      157        trmB+no carb   trmB  HhCa        no carb 1-2020_sugar2
## 158      158 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2020_sugar2
## 159      159    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2020_sugar2
## 160      160     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2020_sugar2
## 161      161    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2020_sugar2
## 162      162   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2020_sugar2
## 163      163   WT+25mM galactose     WT  HhCa 25mM galactose 1-2020_sugar2
## 164      164 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2020_sugar2
## 165      165    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2020_sugar2
## 166      166    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2020_sugar2
## 167      167    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2020_sugar2
## 168      168   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2020_sugar2
## 169      169    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2020_sugar2
## 170      170  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2020_sugar2
## 171      171 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2020_sugar2
## 172      172   WT+25mM galactose     WT  HhCa 25mM galactose 2-2020_sugar2
## 173      173        trmB+no carb   trmB  HhCa        no carb 2-2020_sugar2
## 174      174  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2020_sugar2
## 175      175  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2020_sugar2
## 176      176        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar2
## 177      177          WT+no carb     WT  HhCa        no carb 3-2020_sugar2
## 178      178  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2020_sugar2
## 179      179   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2020_sugar2
## 180      180  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2020_sugar2
## 181      181     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2020_sugar2
## 182      182    WT+25mM fructose     WT  HhCa  25mM fructose 2-2020_sugar2
## 183      183   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2020_sugar2
## 184      184     WT+25mM acetate     WT  HhCa   25mM acetate 2-2020_sugar2
## 185      185     WT+25mM acetate     WT  HhCa   25mM acetate 2-2020_sugar2
## 186      186    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2020_sugar2
## 187      187  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2020_sugar2
## 188      188     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2020_sugar2
## 189      189   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2020_sugar2
## 190      190          WT+no carb     WT  HhCa        no carb 3-2020_sugar2
## 191      191     WT+25mM acetate     WT  HhCa   25mM acetate 2-2020_sugar2
## 192      192          WT+no carb     WT  HhCa        no carb 3-2020_sugar2
## 193      193    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2020_sugar2
## 194      194        trmB+no carb   trmB  HhCa        no carb 2-2020_sugar2
## 195      195   WT+25mM galactose     WT  HhCa 25mM galactose 1-2020_sugar2
## 196      196   WT+25mM galactose     WT  HhCa 25mM galactose 1-2020_sugar2
## 197      197    WT+25mM fructose     WT  HhCa  25mM fructose 2-2020_sugar2
## 198      198      WT+25mM xylose     WT  HhCa    25mM xylose 1-2020_sugar2
## 199      199  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2020_sugar2
## 200      200      WT+25mM ribose     WT  HhCa    25mM ribose 1-2020_sugar2
## 201      201   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2020_sugar2
## 202      202    WT+25mM fructose     WT  HhCa  25mM fructose 1-2020_sugar2
## 203      203          WT+no carb     WT  HhCa        no carb 2-2020_sugar2
## 204      204    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2020_sugar2
## 205      205 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2020_sugar2
## 206      206   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2020_sugar2
## 207      207 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2020_sugar2
## 208      208          WT+no carb     WT  HhCa        no carb 3-2020_sugar2
## 209      209    WT+25mM fructose     WT  HhCa  25mM fructose 1-2020_sugar2
## 210      210        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar2
## 211      211   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2020_sugar2
## 212      212     WT+25mM glucose     WT  HhCa   25mM glucose 2-2020_sugar2
## 213      213  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2020_sugar2
## 214      214      WT+25mM ribose     WT  HhCa    25mM ribose 2-2020_sugar2
## 215      215        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar2
## 216      216      WT+25mM ribose     WT  HhCa    25mM ribose 1-2020_sugar2
## 217      217   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2020_sugar2
## 218      218     WT+25mM glucose     WT  HhCa   25mM glucose 2-2020_sugar2
## 219      219  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2020_sugar2
## 220      220    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2020_sugar2
## 221      221      WT+25mM xylose     WT  HhCa    25mM xylose 1-2020_sugar2
## 222      222  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2020_sugar2
## 223      223      WT+25mM ribose     WT  HhCa    25mM ribose 2-2020_sugar2
## 224      224   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2020_sugar2
## 225      225   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2020_sugar2
## 226      226     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2020_sugar2
## 227      227    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2020_sugar2
## 228      228    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2020_sugar2
## 229      229     WT+25mM glucose     WT  HhCa   25mM glucose 1-2020_sugar2
## 230      230    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2020_sugar2
## 231      231    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2020_sugar2
## 232      232     WT+25mM acetate     WT  HhCa   25mM acetate 1-2020_sugar2
## 233      233   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2020_sugar2
## 234      234    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2020_sugar2
## 235      235   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2020_sugar2
## 236      236   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2020_sugar2
## 237      237    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2020_sugar2
## 238      238        trmB+no carb   trmB  HhCa        no carb 1-2020_sugar2
## 239      239     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2020_sugar2
## 240      240    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2020_sugar2
## 241      241      WT+25mM ribose     WT  HhCa    25mM ribose 1-2020_sugar2
## 242      242     WT+25mM acetate     WT  HhCa   25mM acetate 1-2020_sugar2
## 243      243   WT+25mM galactose     WT  HhCa 25mM galactose 2-2020_sugar2
## 244      244  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2020_sugar2
## 245      245          WT+no carb     WT  HhCa        no carb 1-2020_sugar2
## 246      246  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2020_sugar2
## 247      247   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2020_sugar2
## 248      248  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2020_sugar2
## 249      249 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2020_sugar2
## 250      250   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2020_sugar2
## 251      251      WT+25mM xylose     WT  HhCa    25mM xylose 2-2020_sugar2
## 252      252    WT+25mM fructose     WT  HhCa  25mM fructose 1-2020_sugar2
## 253      253        trmB+no carb   trmB  HhCa        no carb 3-2020_sugar2
## 254      254    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2021_sugar3
## 255      255          WT+no carb     WT  HhCa        no carb 2-2021_sugar3
## 256      256  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2021_sugar3
## 257      257    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2021_sugar3
## 258      258    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2021_sugar3
## 259      259     WT+25mM glucose     WT  HhCa   25mM glucose 1-2021_sugar3
## 260      260        trmB+no carb   trmB  HhCa        no carb 1-2021_sugar3
## 261      261      WT+25mM xylose     WT  HhCa    25mM xylose 2-2021_sugar3
## 262      262      WT+25mM ribose     WT  HhCa    25mM ribose 2-2021_sugar3
## 263      263    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2021_sugar3
## 264      264  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2021_sugar3
## 265      265    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2021_sugar3
## 266      266          WT+no carb     WT  HhCa        no carb 1-2021_sugar3
## 267      267      WT+25mM xylose     WT  HhCa    25mM xylose 2-2021_sugar3
## 268      268     WT+25mM glucose     WT  HhCa   25mM glucose 1-2021_sugar3
## 269      269    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2021_sugar3
## 270      270     WT+25mM glucose     WT  HhCa   25mM glucose 2-2021_sugar3
## 271      271   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2021_sugar3
## 272      272     WT+25mM acetate     WT  HhCa   25mM acetate 1-2021_sugar3
## 273      273          WT+no carb     WT  HhCa        no carb 2-2021_sugar3
## 274      274  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2021_sugar3
## 275      275    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2021_sugar3
## 276      276          WT+no carb     WT  HhCa        no carb 1-2021_sugar3
## 277      277      WT+25mM xylose     WT  HhCa    25mM xylose 1-2021_sugar3
## 278      278  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2021_sugar3
## 279      279        trmB+no carb   trmB  HhCa        no carb 2-2021_sugar3
## 280      280   WT+25mM galactose     WT  HhCa 25mM galactose 2-2021_sugar3
## 281      281    WT+25mM fructose     WT  HhCa  25mM fructose 2-2021_sugar3
## 282      282 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2021_sugar3
## 283      283    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2021_sugar3
## 284      284     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2021_sugar3
## 285      285    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2021_sugar3
## 286      286   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2021_sugar3
## 287      287   WT+25mM galactose     WT  HhCa 25mM galactose 1-2021_sugar3
## 288      288 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2021_sugar3
## 289      289    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2021_sugar3
## 290      290    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2021_sugar3
## 291      291    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 1-2021_sugar3
## 292      292   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2021_sugar3
## 293      293    WT+25mM pyruvate     WT  HhCa  25mM pyruvate 2-2021_sugar3
## 294      294  trmB+25mM fructose   trmB  HhCa  25mM fructose 2-2021_sugar3
## 295      295 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2021_sugar3
## 296      296   WT+25mM galactose     WT  HhCa 25mM galactose 2-2021_sugar3
## 297      297        trmB+no carb   trmB  HhCa        no carb 2-2021_sugar3
## 298      298  trmB+25mM fructose   trmB  HhCa  25mM fructose 1-2021_sugar3
## 299      299  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2021_sugar3
## 300      300        trmB+no carb   trmB  HhCa        no carb 1-2021_sugar3
## 301      301  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2021_sugar3
## 302      302          WT+no carb     WT  HhCa        no carb 1-2021_sugar3
## 303      303   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2021_sugar3
## 304      304  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2021_sugar3
## 305      305     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2021_sugar3
## 306      306    WT+25mM fructose     WT  HhCa  25mM fructose 2-2021_sugar3
## 307      307   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 1-2021_sugar3
## 308      308     WT+25mM acetate     WT  HhCa   25mM acetate 2-2021_sugar3
## 309      309     WT+25mM acetate     WT  HhCa   25mM acetate 2-2021_sugar3
## 310      310    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2021_sugar3
## 311      311  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2021_sugar3
## 312      312     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2021_sugar3
## 313      313   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2021_sugar3
## 314      314          WT+no carb     WT  HhCa        no carb 1-2021_sugar3
## 315      315     WT+25mM acetate     WT  HhCa   25mM acetate 2-2021_sugar3
## 316      316          WT+no carb     WT  HhCa        no carb 1-2021_sugar3
## 317      317    trmB+25mM ribose   trmB  HhCa    25mM ribose 1-2021_sugar3
## 318      318        trmB+no carb   trmB  HhCa        no carb 2-2021_sugar3
## 319      319   WT+25mM galactose     WT  HhCa 25mM galactose 1-2021_sugar3
## 320      320   WT+25mM galactose     WT  HhCa 25mM galactose 1-2021_sugar3
## 321      321    WT+25mM fructose     WT  HhCa  25mM fructose 2-2021_sugar3
## 322      322      WT+25mM xylose     WT  HhCa    25mM xylose 1-2021_sugar3
## 323      323  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2021_sugar3
## 324      324      WT+25mM ribose     WT  HhCa    25mM ribose 1-2021_sugar3
## 325      325   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2021_sugar3
## 326      326    WT+25mM fructose     WT  HhCa  25mM fructose 1-2021_sugar3
## 327      327          WT+no carb     WT  HhCa        no carb 2-2021_sugar3
## 328      328    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2021_sugar3
## 329      329 trmB+25mM galactose   trmB  HhCa 25mM galactose 2-2021_sugar3
## 330      330   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2021_sugar3
## 331      331 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2021_sugar3
## 332      332          WT+no carb     WT  HhCa        no carb 1-2021_sugar3
## 333      333    WT+25mM fructose     WT  HhCa  25mM fructose 1-2021_sugar3
## 334      334        trmB+no carb   trmB  HhCa        no carb 1-2021_sugar3
## 335      335   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2021_sugar3
## 336      336     WT+25mM glucose     WT  HhCa   25mM glucose 2-2021_sugar3
## 337      337  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2021_sugar3
## 338      338      WT+25mM ribose     WT  HhCa    25mM ribose 2-2021_sugar3
## 339      339        trmB+no carb   trmB  HhCa        no carb 1-2021_sugar3
## 340      340      WT+25mM ribose     WT  HhCa    25mM ribose 1-2021_sugar3
## 341      341   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2021_sugar3
## 342      342     WT+25mM glucose     WT  HhCa   25mM glucose 2-2021_sugar3
## 343      343  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 1-2021_sugar3
## 344      344    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2021_sugar3
## 345      345      WT+25mM xylose     WT  HhCa    25mM xylose 1-2021_sugar3
## 346      346  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2021_sugar3
## 347      347      WT+25mM ribose     WT  HhCa    25mM ribose 2-2021_sugar3
## 348      348   trmB+25mM sucrose   trmB  HhCa   25mM sucrose 2-2021_sugar3
## 349      349   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2021_sugar3
## 350      350     WT+25mM sucrose     WT  HhCa   25mM sucrose 1-2021_sugar3
## 351      351    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2021_sugar3
## 352      352    trmB+25mM ribose   trmB  HhCa    25mM ribose 2-2021_sugar3
## 353      353     WT+25mM glucose     WT  HhCa   25mM glucose 1-2021_sugar3
## 354      354    trmB+25mM xylose   trmB  HhCa    25mM xylose 1-2021_sugar3
## 355      355    WT+25mM glycerol     WT  HhCa  25mM glycerol 1-2021_sugar3
## 356      356     WT+25mM acetate     WT  HhCa   25mM acetate 1-2021_sugar3
## 357      357   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2021_sugar3
## 358      358    trmB+25mM xylose   trmB  HhCa    25mM xylose 2-2021_sugar3
## 359      359   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2021_sugar3
## 360      360   trmB+25mM acetate   trmB  HhCa   25mM acetate 2-2021_sugar3
## 361      361    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2021_sugar3
## 362      362   trmB+25mM glucose   trmB  HhCa   25mM glucose 1-2021_sugar3
## 363      363        trmB+no carb   trmB  HhCa        no carb 1-2021_sugar3
## 364      364     WT+25mM sucrose     WT  HhCa   25mM sucrose 2-2021_sugar3
## 365      365    WT+25mM glycerol     WT  HhCa  25mM glycerol 2-2021_sugar3
## 366      366      WT+25mM ribose     WT  HhCa    25mM ribose 1-2021_sugar3
## 367      367     WT+25mM acetate     WT  HhCa   25mM acetate 1-2021_sugar3
## 368      368   WT+25mM galactose     WT  HhCa 25mM galactose 2-2021_sugar3
## 369      369  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 1-2021_sugar3
## 370      370          WT+no carb     WT  HhCa        no carb 1-2021_sugar3
## 371      371  trmB+25mM glycerol   trmB  HhCa  25mM glycerol 2-2021_sugar3
## 372      372   trmB+25mM acetate   trmB  HhCa   25mM acetate 1-2021_sugar3
## 373      373  trmB+25mM pyruvate   trmB  HhCa  25mM pyruvate 2-2021_sugar3
## 374      374 trmB+25mM galactose   trmB  HhCa 25mM galactose 1-2021_sugar3
## 375      375   trmB+25mM glucose   trmB  HhCa   25mM glucose 2-2021_sugar3
## 376      376      WT+25mM xylose     WT  HhCa    25mM xylose 2-2021_sugar3
## 377      377    WT+25mM fructose     WT  HhCa  25mM fructose 1-2021_sugar3
## 378      378        trmB+no carb   trmB  HhCa        no carb 2-2021_sugar3
##     techrep        expt
## 1         2 2020_sugar1
## 2         1 2020_sugar1
## 3         1 2020_sugar1
## 4         2 2020_sugar1
## 5         3 2020_sugar1
## 6         2 2020_sugar1
## 7         2 2020_sugar1
## 8         3 2020_sugar1
## 9         1 2020_sugar1
## 10        3 2020_sugar1
## 11        2 2020_sugar1
## 12        2 2020_sugar1
## 13        2 2020_sugar1
## 14        1 2020_sugar1
## 15        1 2020_sugar1
## 16        2 2020_sugar1
## 17        3 2020_sugar1
## 18        3 2020_sugar1
## 19        3 2020_sugar1
## 20        3 2020_sugar1
## 21        3 2020_sugar1
## 22        2 2020_sugar1
## 23        1 2020_sugar1
## 24        1 2020_sugar1
## 25        3 2020_sugar1
## 26        2 2020_sugar1
## 27        3 2020_sugar1
## 28        1 2020_sugar1
## 29        2 2020_sugar1
## 30        1 2020_sugar1
## 31        1 2020_sugar1
## 32        3 2020_sugar1
## 33        2 2020_sugar1
## 34        1 2020_sugar1
## 35        1 2020_sugar1
## 36        1 2020_sugar1
## 37        1 2020_sugar1
## 38        2 2020_sugar1
## 39        1 2020_sugar1
## 40        3 2020_sugar1
## 41        1 2020_sugar1
## 42        3 2020_sugar1
## 43        2 2020_sugar1
## 44        3 2020_sugar1
## 45        1 2020_sugar1
## 46        3 2020_sugar1
## 47        1 2020_sugar1
## 48        2 2020_sugar1
## 49        2 2020_sugar1
## 50        3 2020_sugar1
## 51        2 2020_sugar1
## 52        3 2020_sugar1
## 53        3 2020_sugar1
## 54        1 2020_sugar1
## 55        3 2020_sugar1
## 56        2 2020_sugar1
## 57        3 2020_sugar1
## 58        2 2020_sugar1
## 59        2 2020_sugar1
## 60        2 2020_sugar1
## 61        1 2020_sugar1
## 62        3 2020_sugar1
## 63        2 2020_sugar1
## 64        1 2020_sugar1
## 65        3 2020_sugar1
## 66        2 2020_sugar1
## 67        1 2020_sugar1
## 68        2 2020_sugar1
## 69        2 2020_sugar1
## 70        3 2020_sugar1
## 71        2 2020_sugar1
## 72        2 2020_sugar1
## 73        3 2020_sugar1
## 74        1 2020_sugar1
## 75        1 2020_sugar1
## 76        2 2020_sugar1
## 77        3 2020_sugar1
## 78        3 2020_sugar1
## 79        1 2020_sugar1
## 80        2 2020_sugar1
## 81        2 2020_sugar1
## 82        3 2020_sugar1
## 83        3 2020_sugar1
## 84        2 2020_sugar1
## 85        2 2020_sugar1
## 86        2 2020_sugar1
## 87        2 2020_sugar1
## 88        1 2020_sugar1
## 89        3 2020_sugar1
## 90        2 2020_sugar1
## 91        3 2020_sugar1
## 92        2 2020_sugar1
## 93        1 2020_sugar1
## 94        2 2020_sugar1
## 95        3 2020_sugar1
## 96        1 2020_sugar1
## 97        1 2020_sugar1
## 98        3 2020_sugar1
## 99        1 2020_sugar1
## 100       1 2020_sugar1
## 101       3 2020_sugar1
## 102       2 2020_sugar1
## 103       1 2020_sugar1
## 104       1 2020_sugar1
## 105       1 2020_sugar1
## 106       1 2020_sugar1
## 107       3 2020_sugar1
## 108       1 2020_sugar1
## 109       3 2020_sugar1
## 110       1 2020_sugar1
## 111       1 2020_sugar1
## 112       3 2020_sugar1
## 113       2 2020_sugar1
## 114       3 2020_sugar1
## 115       1 2020_sugar1
## 116       3 2020_sugar1
## 117       1 2020_sugar1
## 118       3 2020_sugar1
## 119       2 2020_sugar1
## 120       1 2020_sugar1
## 121       2 2020_sugar1
## 122       3 2020_sugar1
## 123       4 2020_sugar1
## 124       1 2020_sugar1
## 125       2 2020_sugar1
## 126       3 2020_sugar1
## 127       4 2020_sugar1
## 128       3 2020_sugar2
## 129       1 2020_sugar2
## 130       2 2020_sugar2
## 131       2 2020_sugar2
## 132       1 2020_sugar2
## 133       1 2020_sugar2
## 134       3 2020_sugar2
## 135       2 2020_sugar2
## 136       2 2020_sugar2
## 137       2 2020_sugar2
## 138       3 2020_sugar2
## 139       2 2020_sugar2
## 140       2 2020_sugar2
## 141       1 2020_sugar2
## 142       2 2020_sugar2
## 143       1 2020_sugar2
## 144       1 2020_sugar2
## 145       1 2020_sugar2
## 146       1 2020_sugar2
## 147       3 2020_sugar2
## 148       3 2020_sugar2
## 149       2 2020_sugar2
## 150       3 2020_sugar2
## 151       1 2020_sugar2
## 152       1 2020_sugar2
## 153       2 2020_sugar2
## 154       3 2020_sugar2
## 155       1 2020_sugar2
## 156       2 2020_sugar2
## 157       2 2020_sugar2
## 158       3 2020_sugar2
## 159       1 2020_sugar2
## 160       1 2020_sugar2
## 161       1 2020_sugar2
## 162       1 2020_sugar2
## 163       2 2020_sugar2
## 164       1 2020_sugar2
## 165       2 2020_sugar2
## 166       2 2020_sugar2
## 167       3 2020_sugar2
## 168       3 2020_sugar2
## 169       3 2020_sugar2
## 170       1 2020_sugar2
## 171       2 2020_sugar2
## 172       3 2020_sugar2
## 173       1 2020_sugar2
## 174       1 2020_sugar2
## 175       1 2020_sugar2
## 176       1 2020_sugar2
## 177       4 2020_sugar2
## 178       2 2020_sugar2
## 179       1 2020_sugar2
## 180       1 2020_sugar2
## 181       2 2020_sugar2
## 182       1 2020_sugar2
## 183       2 2020_sugar2
## 184       2 2020_sugar2
## 185       1 2020_sugar2
## 186       2 2020_sugar2
## 187       3 2020_sugar2
## 188       3 2020_sugar2
## 189       1 2020_sugar2
## 190       2 2020_sugar2
## 191       3 2020_sugar2
## 192       3 2020_sugar2
## 193       1 2020_sugar2
## 194       2 2020_sugar2
## 195       3 2020_sugar2
## 196       1 2020_sugar2
## 197       3 2020_sugar2
## 198       2 2020_sugar2
## 199       1 2020_sugar2
## 200       1 2020_sugar2
## 201       3 2020_sugar2
## 202       3 2020_sugar2
## 203       2 2020_sugar2
## 204       3 2020_sugar2
## 205       3 2020_sugar2
## 206       2 2020_sugar2
## 207       1 2020_sugar2
## 208       1 2020_sugar2
## 209       1 2020_sugar2
## 210       4 2020_sugar2
## 211       3 2020_sugar2
## 212       3 2020_sugar2
## 213       3 2020_sugar2
## 214       1 2020_sugar2
## 215       2 2020_sugar2
## 216       2 2020_sugar2
## 217       1 2020_sugar2
## 218       2 2020_sugar2
## 219       2 2020_sugar2
## 220       2 2020_sugar2
## 221       3 2020_sugar2
## 222       3 2020_sugar2
## 223       3 2020_sugar2
## 224       2 2020_sugar2
## 225       1 2020_sugar2
## 226       3 2020_sugar2
## 227       1 2020_sugar2
## 228       1 2020_sugar2
## 229       3 2020_sugar2
## 230       3 2020_sugar2
## 231       3 2020_sugar2
## 232       3 2020_sugar2
## 233       2 2020_sugar2
## 234       3 2020_sugar2
## 235       3 2020_sugar2
## 236       3 2020_sugar2
## 237       3 2020_sugar2
## 238       1 2020_sugar2
## 239       2 2020_sugar2
## 240       1 2020_sugar2
## 241       3 2020_sugar2
## 242       2 2020_sugar2
## 243       2 2020_sugar2
## 244       2 2020_sugar2
## 245       1 2020_sugar2
## 246       2 2020_sugar2
## 247       3 2020_sugar2
## 248       1 2020_sugar2
## 249       2 2020_sugar2
## 250       2 2020_sugar2
## 251       3 2020_sugar2
## 252       2 2020_sugar2
## 253       3 2020_sugar2
## 254       3 2021_sugar3
## 255       1 2021_sugar3
## 256       2 2021_sugar3
## 257       2 2021_sugar3
## 258       1 2021_sugar3
## 259       1 2021_sugar3
## 260       3 2021_sugar3
## 261       2 2021_sugar3
## 262       2 2021_sugar3
## 263       2 2021_sugar3
## 264       3 2021_sugar3
## 265       2 2021_sugar3
## 266       2 2021_sugar3
## 267       1 2021_sugar3
## 268       2 2021_sugar3
## 269       1 2021_sugar3
## 270       1 2021_sugar3
## 271       1 2021_sugar3
## 272       1 2021_sugar3
## 273       3 2021_sugar3
## 274       3 2021_sugar3
## 275       2 2021_sugar3
## 276       3 2021_sugar3
## 277       1 2021_sugar3
## 278       2 2021_sugar3
## 279       3 2021_sugar3
## 280       1 2021_sugar3
## 281       2 2021_sugar3
## 282       3 2021_sugar3
## 283       1 2021_sugar3
## 284       1 2021_sugar3
## 285       1 2021_sugar3
## 286       1 2021_sugar3
## 287       2 2021_sugar3
## 288       1 2021_sugar3
## 289       2 2021_sugar3
## 290       2 2021_sugar3
## 291       3 2021_sugar3
## 292       3 2021_sugar3
## 293       3 2021_sugar3
## 294       1 2021_sugar3
## 295       2 2021_sugar3
## 296       3 2021_sugar3
## 297       1 2021_sugar3
## 298       1 2021_sugar3
## 299       1 2021_sugar3
## 300       4 2021_sugar3
## 301       3 2021_sugar3
## 302       4 2021_sugar3
## 303       1 2021_sugar3
## 304       1 2021_sugar3
## 305       2 2021_sugar3
## 306       1 2021_sugar3
## 307       2 2021_sugar3
## 308       2 2021_sugar3
## 309       1 2021_sugar3
## 310       2 2021_sugar3
## 311       3 2021_sugar3
## 312       3 2021_sugar3
## 313       1 2021_sugar3
## 314       5 2021_sugar3
## 315       3 2021_sugar3
## 316       6 2021_sugar3
## 317       1 2021_sugar3
## 318       2 2021_sugar3
## 319       3 2021_sugar3
## 320       1 2021_sugar3
## 321       3 2021_sugar3
## 322       2 2021_sugar3
## 323       1 2021_sugar3
## 324       1 2021_sugar3
## 325       3 2021_sugar3
## 326       3 2021_sugar3
## 327       2 2021_sugar3
## 328       3 2021_sugar3
## 329       3 2021_sugar3
## 330       2 2021_sugar3
## 331       1 2021_sugar3
## 332       7 2021_sugar3
## 333       1 2021_sugar3
## 334       5 2021_sugar3
## 335       3 2021_sugar3
## 336       3 2021_sugar3
## 337       3 2021_sugar3
## 338       1 2021_sugar3
## 339       6 2021_sugar3
## 340       2 2021_sugar3
## 341       1 2021_sugar3
## 342       2 2021_sugar3
## 343       2 2021_sugar3
## 344       2 2021_sugar3
## 345       3 2021_sugar3
## 346       3 2021_sugar3
## 347       3 2021_sugar3
## 348       2 2021_sugar3
## 349       1 2021_sugar3
## 350       3 2021_sugar3
## 351       1 2021_sugar3
## 352       1 2021_sugar3
## 353       3 2021_sugar3
## 354       3 2021_sugar3
## 355       3 2021_sugar3
## 356       3 2021_sugar3
## 357       2 2021_sugar3
## 358       3 2021_sugar3
## 359       3 2021_sugar3
## 360       3 2021_sugar3
## 361       3 2021_sugar3
## 362       2 2021_sugar3
## 363       1 2021_sugar3
## 364       2 2021_sugar3
## 365       1 2021_sugar3
## 366       3 2021_sugar3
## 367       2 2021_sugar3
## 368       2 2021_sugar3
## 369       2 2021_sugar3
## 370       1 2021_sugar3
## 371       2 2021_sugar3
## 372       3 2021_sugar3
## 373       1 2021_sugar3
## 374       2 2021_sugar3
## 375       2 2021_sugar3
## 376       3 2021_sugar3
## 377       2 2021_sugar3
## 378       8 2021_sugar3
```

```r
# bioreps with at least 3 tech reps
mt0[, c(2, 3, 4, 6, 7)] %>%
  group_by(ID, biorep) %>%
  tally() %>%
  filter(n >= 3) %>%
  group_by(ID) %>%
  tally()
```

```
## # A tibble: 20 × 2
##    ID                      n
##    <chr>               <int>
##  1 WT+25mM acetate         6
##  2 WT+25mM fructose        6
##  3 WT+25mM galactose       6
##  4 WT+25mM glucose         6
##  5 WT+25mM glycerol        6
##  6 WT+25mM pyruvate        6
##  7 WT+25mM ribose          6
##  8 WT+25mM sucrose         5
##  9 WT+25mM xylose          6
## 10 WT+no carb              8
## 11 trmB+25mM acetate       6
## 12 trmB+25mM fructose      6
## 13 trmB+25mM galactose     6
## 14 trmB+25mM glucose       5
## 15 trmB+25mM glycerol      6
## 16 trmB+25mM pyruvate      4
## 17 trmB+25mM ribose        6
## 18 trmB+25mM sucrose       6
## 19 trmB+25mM xylose        5
## 20 trmB+no carb            8
```

### Growth curve modeling

```r
# Create an output data frame to store the results in.
num_analyses <- length(names(d)) - 1
d_gc <- data.frame(
  sample = character(num_analyses),
  k = numeric(num_analyses),
  n0 = numeric(num_analyses),
  r = numeric(num_analyses),
  t_mid = numeric(num_analyses),
  t_gen = numeric(num_analyses),
  auc_l = numeric(num_analyses),
  auc_e = numeric(num_analyses),
  sigma = numeric(num_analyses),
  stringsAsFactors = FALSE
)

# Truncate or trim the input data to observations occurring in the first 80 hours.
trim_at_time <- 80.5

# Loop through all of the columns in the data frame. For each column, run Growthcurver.
pdf(paste("figures/", expt, "_growthcurver_fits.pdf", sep = ""), height = 8.5, width = 11)
par(mfcol = c(8, 12))
par(mar = c(0.25, 0.25, 0.25, 0.25))
y_lim_max <- max(d[, setdiff(names(d), "time")], na.rm = T) - min(d[, setdiff(names(d), "time")], na.rm = T)

n <- 1 # keeps track of the current row in the output data frame
for (col_name in names(d)) {

  # Don't process the time column.
  if (col_name != "time") {

    # Create a temporary data frame that contains just the time and current col
    d_loop <- d[c("time", col_name)]
    
    # Now, call Growthcurver to calculate the metrics using SummarizeGrowth
    gc_fit <- SummarizeGrowth(
      data_t = d_loop$time,
      data_n = d_loop[col_name],
      t_trim = trim_at_time,
      bg_correct = "none"
    )

    if (gc_fit$vals[[16]] == "cannot fit data") {
      d_gc$sample[n] <- col_name
      n <- n + 1
    } else {
      d_gc$sample[n] <- col_name
      d_gc[n, 2:9] <- c(
        gc_fit$vals$k,
        gc_fit$vals$n0,
        gc_fit$vals$r,
        gc_fit$vals$t_mid,
        gc_fit$vals$t_gen,
        gc_fit$vals$auc_l,
        gc_fit$vals$auc_e,
        gc_fit$vals$sigma
      )

      n <- n + 1
      # Finally, plot the raw data and the fitted curve. print some of the data points to keep the file size smaller
      n_obs <- length(gc_fit$data$t)
      idx_to_plot <- 1:20 / 20 * n_obs
      plot(gc_fit$data$t[idx_to_plot], gc_fit$data$N[idx_to_plot],
        pch = 20,
        xlim = c(0, trim_at_time),
        ylim = c(0, y_lim_max),
        cex = 0.6, xaxt = "n", yaxt = "n"
      )
      text(x = trim_at_time / 4, y = y_lim_max, labels = col_name, pos = 1)
      lines(gc_fit$data$t, predict(gc_fit$model), col = "red")
    }
  }
}
dev.off()
```

```
## png 
##   2
```

graph fitted curves:

```r
fitted_curves <- list()

n <- 1 # keeps track of the current row in the output data frame
for (col_name in names(d)) {
  if (col_name != "time") {
    d_loop <- d[c("time", col_name)]
    gc_fit <- SummarizeGrowth(
      data_t = d_loop$time,
      data_n = d_loop[col_name],
      t_trim = trim_at_time,
      bg_correct = "none"
    )
    if (gc_fit$vals[[16]] == "cannot fit data") {
      d_gc$sample[n] <- col_name
      n <- n + 1
    } else {
      fitted_curves[[n]] <- predict(gc_fit$model)
      names(fitted_curves)[n] <- col_name
      n <- n + 1
    }
  }
}

tmp <- Filter(length, fitted_curves)
fits <- as.data.frame(tmp, col.names = names(tmp))
colnames(fits) <- names(tmp)
fits <- cbind("time" = d$time[d$time < trim_at_time], fits)

# re-zero lowest value
mins <- Rfast::colMins(as.matrix(fits), value = T) # get column minimums
fits <- sweep(fits, 2, mins, "-") # subtract smallest value for each column

# Convert data from wide to long format
m_dt <- reshape2::melt(fits, id = "time")
mtdt <- merge(m_dt, mt0, by = "variable")
well <- paste(mtdt$variable, mtdt$ID, sep = " ")
mtdt <- cbind(well, mtdt)

t_mtdt2 <- mtdt[!(mtdt$ID == ""), ]

stats2 <- t_mtdt2 %>%
  group_by(ID, strain, condition, time) %>%
  summarise(
    reps = length(value),
    average = mean(value) + 1,
    CI95 = conf_int95(value)
  )
unique(stats2$strain)
```

```
## [1] "WT"   "trmB"
```

```r
# all conditions tested.
stats2 %>%
  ggplot(., aes(x = time, y = log10(average), color = strain)) +
  xlab("time (h)") +
  ylab("log10(OD)") +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(0, 0.13)) +
  scale_x_continuous(limits = c(0, 80)) +
  geom_ribbon(aes(ymin = log10(average - CI95), ymax = log10(average + CI95), fill = strain, color = NULL), alpha = 0.2) +
  scale_fill_viridis("strain", discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_color_viridis("strain", discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw() +
  facet_wrap(~condition, nrow = 2) +
  theme(
    axis.title.x = element_text(face = "bold", size = 10, angle = 0), axis.title.y = element_text(face = "bold", size = 10, angle = 90),
    axis.text.y = element_text(size = 8, angle = 0), axis.text.x = element_text(size = 8, angle = 0),
    axis.ticks.length.y.left = unit(.05, "cm"), strip.text = element_text(size = 10),
    strip.background = element_blank(), strip.placement = "outside",
    legend.position = "bottom"
  ) -> curves
plot(curves)
```

![](Fig2B_har_sugars_modeling_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
png(paste("figures/", expt, "_curves_supp.png", sep = ""), width = 6, height = 4, units = "in", res = 1200, bg = "transparent")
plot(curves)
dev.off()
```

```
## png 
##   2
```

```r
pdf(paste("figures/", expt, "_curves_supp.pdf", sep = ""), width = 6, height = 4, bg = "transparent")
plot(curves)
dev.off()
```

```
## png 
##   2
```

insets

```r
stats2$condition <- factor(stats2$condition, levels = c("25mM sucrose", "25mM glucose", "25mM glycerol", "25mM fructose", "25mM xylose", "25mM galactose", "25mM ribose", "25mM pyruvate", "25mM acetate", "no carb"))

stats2 %>%
  filter(strain == "trmB") %>%
  ggplot(., aes(x = time, y = log10(average), color = condition)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = log10(average - CI95), ymax = log10(average + CI95), fill = condition, color = NULL), alpha = 0.2) +
  xlab("time (h)") +
  ylab("log10(Optical Denisty)") +
  scale_x_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_fill_viridis("condition", discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  scale_color_viridis("condition", discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 14, angle = 0),
    axis.title.y = element_text(face = "bold", size = 14, angle = 90),
    axis.ticks = element_blank(), legend.position = "none"
  ) -> inset
plot(inset)
```

![](Fig2B_har_sugars_modeling_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
png(paste("figures/", expt, "_inset.png", sep = ""), width = 3, height = 3, units = "in", res = 1200, bg = "transparent")
plot(inset)
dev.off()
```

```
## png 
##   2
```

```r
pdf(paste("figures/", expt, "_inset.pdf", sep = ""), width = 3, height = 3, bg = "transparent")
plot(inset)
dev.off()
```

```
## png 
##   2
```

# TrmB normalized to WT by condition
take average of WT bioreps. calculate ratios of individual trmB bioreps. then average and calc SD/std err.

```r
colnames(d_gc)[1] <- "variable"
merge(d_gc, mt0) -> tmp

tmp %>%
  group_by(strain, condition, biorep, ID, media) %>%
  summarise(avg_auc = mean(auc_e)) %>% # average technical replicates
  group_by(condition, biorep) %>%
  summarise(rat_auc = avg_auc[strain == "trmB"] / avg_auc[strain == "WT"]) %>%
  group_by(condition) %>% # average ratios and calc sd
  summarise(
    avg_auc = mean(rat_auc),
    auc_sd = sd(rat_auc)
  ) -> norms.hca

# factorize for sorting
norms.hca$condition <- factor(norms.hca$condition)

# TrmB parameters normalized to WT (by condition)
norms.hca %>%
  slice(8, 4, 5, 2, 9, 3, 7, 6, 1, 10) %>%
  mutate(condition = factor(condition, levels = condition)) %>%
  ggplot(aes(x = condition, y = avg_auc)) +
  xlab("") +
  ylab(expression("AUC"[" trmB" / " parent"])) +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", alpha = 0.3) +
  geom_bar(stat = "identity", aes(fill = condition)) +
  geom_errorbar(aes(ymin = avg_auc - auc_sd, ymax = avg_auc + auc_sd),
    width = .1, size = 1, position = position_dodge(.9)
  ) +
  scale_fill_viridis("x", discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  scale_y_continuous(limits = c(0, 1.1), expand = expansion(mult = c(0, 0.02))) +
  theme_bw() +
  theme(
    axis.title.y = element_text(face = "plain", color = "#000000", size = 14, angle = 90),
    axis.text.y = element_text(face = "plain", color = "#000000", size = 12, angle = 0),
    axis.text.x = element_text(face = "plain", color = "#000000", size = 12, angle = 45, hjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  ) -> figure
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```r
plot(figure)
```

![](Fig2B_har_sugars_modeling_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
png(paste("figures/", expt, "_AUC_plot.png", sep = ""), width = 10, height = 5, units = "in", res = 1200, bg = "transparent")
plot(figure)
dev.off()
```

```
## png 
##   2
```

```r
pdf(paste("figures/", expt, "_AUC_plot.pdf", sep = ""), width = 10, height = 5, bg = "transparent")
plot(figure)
dev.off()
```

```
## png 
##   2
```

```r
# Testing for significance 
tmp %>%
  group_by(strain, condition, biorep, ID) %>%
  summarise(avg_auc = mean(auc_e)) -> tmp2

for (i in 1:length(unique(tmp2$condition))) {
  cond <- unique(tmp2$condition)[i]
  trmB <- paste("trmB", cond, sep = "+")
  wt <- paste("WT", cond, sep = "+")
  one <- filter(tmp2, ID == trmB)$avg_auc
  two <- filter(tmp2, ID == wt)$avg_auc
  test <- t.test(one, two, alternative = "less")
  print(paste("AUC difference in", cond, "is significant with p-value of ", round(test$p.value, digits = 8)))
}
```

```
## [1] "AUC difference in 25mM acetate is significant with p-value of  2.34e-06"
## [1] "AUC difference in 25mM fructose is significant with p-value of  0.00428081"
## [1] "AUC difference in 25mM galactose is significant with p-value of  2.51e-06"
## [1] "AUC difference in 25mM glucose is significant with p-value of  0.02600961"
## [1] "AUC difference in 25mM glycerol is significant with p-value of  0.01399986"
## [1] "AUC difference in 25mM pyruvate is significant with p-value of  6e-08"
## [1] "AUC difference in 25mM ribose is significant with p-value of  7.32e-06"
## [1] "AUC difference in 25mM sucrose is significant with p-value of  0.05066543"
## [1] "AUC difference in 25mM xylose is significant with p-value of  6.015e-05"
## [1] "AUC difference in no carb is significant with p-value of  4.18e-06"
```

## comparing conditions by strain

```r
# average no sugar
tmp %>%
  filter(strain == "WT" & condition == "no carb") %>%
  group_by(strain, condition) %>%
  summarise(avg_auc = mean(auc_e)) -> wt.nocarb
tmp %>%
  filter(strain == "trmB" & condition == "no carb") %>%
  group_by(strain, condition) %>%
  summarise(avg_auc = mean(auc_e)) -> trmb.nocarb

# split wt and trmB
tmp %>%
  filter(strain == "WT") %>%
  group_by(strain, condition, biorep) %>%
  summarise(
    auc = mean(auc_e),
    auc.rat = auc / wt.nocarb$avg_auc
  ) %>% # take ratio
  group_by(strain, condition) %>%
  summarise(
    avg.auc.rat = mean(auc.rat),
    sd.auc.rat = sd(auc.rat)
  ) -> reps.wt

tmp %>%
  filter(strain == "trmB") %>%
  group_by(strain, condition, biorep) %>%
  summarise(
    auc = mean(auc_e),
    auc.rat = auc / trmb.nocarb$avg_auc
  ) %>% # take ratio
  group_by(strain, condition) %>%
  summarise(
    avg.auc.rat = mean(auc.rat),
    sd.auc.rat = sd(auc.rat)
  ) -> reps.t

# combine into single DF
rbind(reps.t, reps.wt) -> test

# graph with facets
test %>%
  arrange(desc(avg.auc.rat)) %>%
  mutate(condition = factor(condition, levels = condition)) %>%
  ggplot(aes(x = condition, y = avg.auc.rat, fill = strain)) +
  ylab(expression("AUC" / "AUC"[" average no carb"])) +
  xlab("") +
  geom_hline(yintercept = 1, size = 0.5, linetype = "dashed") +
  geom_bar(stat = "identity", alpha = 0.8, position = "dodge") +
  geom_errorbar(aes(ymin = avg.auc.rat - sd.auc.rat, ymax = avg.auc.rat + sd.auc.rat),
    width = .15, size = .7, position = position_dodge(.9)
  ) +
  facet_wrap(~strain, nrow = 2) +
  scale_fill_viridis("strain", discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(limits = c(0,8), expand = expansion(mult = c(0, 0.02))) +
  theme_bw() +
  theme(
    axis.title.y = element_text(face = "bold", color = "#000000", size = 10, angle = 90),
    axis.text.y = element_text(face = "plain", color = "#000000", size = 10, angle = 0),
    axis.text.x = element_text(face = "plain", color = "#000000", size = 10, angle = 45, hjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12), strip.background = element_blank(), strip.placement = "outside",
    legend.position = "none"
  ) -> no.carb.AUC

plot(no.carb.AUC)
```

![](Fig2B_har_sugars_modeling_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
png(paste("figures/", expt, "_AUCs_supp.png", sep = ""), width = 4, height = 6, units = "in", res = 1200, bg = "transparent")
plot(no.carb.AUC)
dev.off()
```

```
## png 
##   2
```

```r
pdf(paste("figures/", expt, "_AUCs_supp.pdf", sep = ""), width = 4, height = 6, bg = "transparent")
plot(no.carb.AUC)
dev.off()
```

```
## png 
##   2
```


```r
# Testing for significance 
tmp %>%
  filter(strain == "WT") %>%
  group_by(strain, condition, biorep) %>%
  summarise(
    auc = mean(auc_e)
  ) -> tmp2

for (i in 1:9) {
  cond <- unique(tmp2$condition)[i]
  one <- filter(tmp2, condition == cond)$auc
  two <- filter(tmp2, condition == "no carb")$auc
  test <- t.test(one, two, alternative = "two.sided")
  print(paste("AUC difference in", cond, "is significant with p-value of ", round(test$p.value, digits = 10)))
}
```

```
## [1] "AUC difference in 25mM acetate is significant with p-value of  0.5406210155"
## [1] "AUC difference in 25mM fructose is significant with p-value of  0.00021249"
## [1] "AUC difference in 25mM galactose is significant with p-value of  0.0284573449"
## [1] "AUC difference in 25mM glucose is significant with p-value of  0.0007455685"
## [1] "AUC difference in 25mM glycerol is significant with p-value of  0.0010260187"
## [1] "AUC difference in 25mM pyruvate is significant with p-value of  3.842e-07"
## [1] "AUC difference in 25mM ribose is significant with p-value of  4.49534e-05"
## [1] "AUC difference in 25mM sucrose is significant with p-value of  0.0006142463"
## [1] "AUC difference in 25mM xylose is significant with p-value of  0.0075080376"
```

```r
tmp %>%
  filter(strain == "trmB") %>%
  group_by(strain, condition, biorep) %>%
  summarise(
    auc = mean(auc_e)
  ) -> tmp2

for (i in 1:9) {
  cond <- unique(tmp2$condition)[i]
  one <- filter(tmp2, condition == cond)$auc
  two <- filter(tmp2, condition == "no carb")$auc
  test <- t.test(one, two, alternative = "two.sided")
  print(paste("AUC difference in", cond, "is significant with p-value of ", round(test$p.value, digits = 10)))
}
```

```
## [1] "AUC difference in 25mM acetate is significant with p-value of  0.5148344947"
## [1] "AUC difference in 25mM fructose is significant with p-value of  8.31383e-05"
## [1] "AUC difference in 25mM galactose is significant with p-value of  0.0027970953"
## [1] "AUC difference in 25mM glucose is significant with p-value of  6.39478e-05"
## [1] "AUC difference in 25mM glycerol is significant with p-value of  4.71607e-05"
## [1] "AUC difference in 25mM pyruvate is significant with p-value of  0.2175093351"
## [1] "AUC difference in 25mM ribose is significant with p-value of  0.0315374535"
## [1] "AUC difference in 25mM sucrose is significant with p-value of  0.0001881377"
## [1] "AUC difference in 25mM xylose is significant with p-value of  0.0003298068"
```
