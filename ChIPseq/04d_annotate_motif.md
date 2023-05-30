---
title: "Annotate motifs instances identified with FIMO"
author: Rylee K. Hackley
output: 
  html_document:
    keep_md: true
---


```r
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(cowplot)
library(janitor)
```


```r
motifs <- read_tsv("04b_motif_discovery/1st_order_15peaks_19_whole_genome_FIMO.tsv") %>%
  head(-3) # XSTREME returns instances of motifs in provided sequences. these already map to peaks.
```

```
## Warning: One or more parsing issues, call `problems()` on your data frame for details,
## e.g.:
##   dat <- vroom(...)
##   problems(dat)
```

```r
gr <- makeGRangesFromDataFrame(motifs, seqnames.field = "sequence_name", keep.extra.columns = T)
peaks <- read_csv("03_peak_visualize/03_consensus_peaks.csv")[-c(1:3)] %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T)
peaks_anno <- read_csv("04a_peak_annotation/04a_consensus_genelist.csv")
peaks_anno <- makeGRangesFromDataFrame(peaks_anno, start.field = "peak_start", end.field = "peak_end", keep.extra.columns = T)

# include all peaks enriched in no glucose
peaks_DB <- rtracklayer::import("02_DiffBind/02_diffbind2022.bed")

# load genome annotation files
gff <- makeTxDbFromGFF("../00_genome_files/genomic.gff",
  format = "gff",
  dataSource = "NCBI",
  organism = "Haloarcula hispanica",
  dbxrefTag = "locus_tag"
)
```

```
## Warning in .extract_transcripts_from_GRanges(tx_IDX, gr, mcols0$type, mcols0$ID, : the transcript names ("tx_name" column in the TxDb object) imported
##   from the "Name" attribute are not unique
```

```r
gff.df <- read_csv("../00_genome_files/GCF_000223905.1_gff.key.csv")
```

Get GRanges of genes and promoters and intergenic regions

```r
genes.only <- genes(gff, columns = c("GENEID"))
genes.only$type <- rep("gene", length(genes.only))

# get promoter regions
pro250 <- promoters(genes.only, upstream = 250, downstream = 0) ## Warning: this encroaches on nearby CDS
pro250$type <- rep("promoter", length(pro250))

# non coding regions of the genome
IG <- gaps(GenomicRanges::reduce(genes.only, ignore.strand = T))
```

overlap promoter region and IG regions to get conservative promoter ranges (that respect genic regions)

```r
findOverlaps(pro250, IG) -> pros
pintersect(pro250[queryHits(pros)], IG[subjectHits(pros)]) -> trimmed.pro
trimmed.pro <- trimmed.pro[, 1:2]
names(trimmed.pro) <- NULL

# regenerate IG ranges to respect promoter regions
ig.only <- gaps(GenomicRanges::reduce(c(genes.only, trimmed.pro), ignore.strand = T))
```

How many enriched regions have motifs in them?? how many times does the motif occur in the genome?

```r
# motifs that overlap with peaks
length(subsetByOverlaps(peaks, gr))
```

```
## [1] 15
```

How many motifs are within 20bp of each other?

```r
tmp <- GenomicRanges::reduce(resize(gr, width = 25, fix = "center"), ignore.strand = T, with.revmap = T)
palindromes <- GRangesList()

for (i in 1:length(tmp)) {
  a <- length(tmp$revmap[[i]])
  if (a > 1) {
    b <- tmp$revmap[[i]][1]
    c <- tmp$revmap[[i]][2]
    rg <- width(range(gr[c(b, c)], ignore.strand = T))
    #print(paste("motif pair footprint is", rg, "bps"))
    if (rg <= 25) {
      palindromes <- append(palindromes, GRangesList(gr[c(b, c)]))
    }
  }
}
# 8 pairs of motifs within 22 bps.
as.data.frame(unlist(palindromes))
```

```
##       seqnames   start     end width strand            motif_id motif_alt_id
## 1  NC_015944.1   87858   87876    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 2  NC_015944.1   87861   87879    19      - ATTHACTMGAAWMBKWGTA       MEME-1
## 3  NC_015944.1  134417  134435    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 4  NC_015944.1  134420  134438    19      - ATTHACTMGAAWMBKWGTA       MEME-1
## 5  NC_015948.1  168916  168934    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 6  NC_015948.1  168919  168937    19      - ATTHACTMGAAWMBKWGTA       MEME-1
## 7  NC_015948.1  420780  420798    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 8  NC_015948.1  420783  420801    19      - ATTHACTMGAAWMBKWGTA       MEME-1
## 9  NC_015948.1 1358582 1358600    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 10 NC_015948.1 1358585 1358603    19      - ATTHACTMGAAWMBKWGTA       MEME-1
## 11 NC_015948.1 2237586 2237604    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 12 NC_015948.1 2237589 2237607    19      - ATTHACTMGAAWMBKWGTA       MEME-1
## 13 NC_015948.1 2492483 2492501    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 14 NC_015948.1 2492486 2492504    19      - ATTHACTMGAAWMBKWGTA       MEME-1
## 15 NC_015948.1 2701208 2701226    19      + ATTHACTMGAAWMBKWGTA       MEME-1
## 16 NC_015948.1 2701211 2701229    19      - ATTHACTMGAAWMBKWGTA       MEME-1
##       score  p.value q.value    matched_sequence
## 1  10.20410 4.43e-05  1.0000 TTGCACGCCAATTCTAGTG
## 2   9.01020 7.30e-05  1.0000 AATCACTAGAATTGGCGTG
## 3  15.44900 2.93e-06  1.0000 GTTTACTATAATTATAGTA
## 4  13.23470 1.04e-05  1.0000 AATTACTATAATTATAGTA
## 5  12.52040 1.50e-05  1.0000 GCTGACACGGAATCGAGTA
## 6  19.52040 1.46e-07  0.1620 TTTTACTCGATTCCGTGTC
## 7   9.37755 6.28e-05  1.0000 TCAGCCTCGAAACCGAGTG
## 8   9.75510 5.37e-05  1.0000 GTTCACTCGGTTTCGAGGC
## 9  23.21430 2.77e-09  0.0108 ATTCACTCGAAATGTAGTA
## 10 17.11220 9.71e-07  0.6860 ATTTACTACATTTCGAGTG
## 11 22.73470 5.26e-09  0.0136 ATTCACTCGAAATCGAGTG
## 12 12.65310 1.41e-05  1.0000 CTTCACTCGATTTCGAGTG
## 13 10.39800 4.07e-05  1.0000 TATGATTCCGAACCGAGTG
## 14  8.67347 8.35e-05  1.0000 AATCACTCGGTTCGGAATC
## 15 17.71430 6.31e-07  0.4910 AATTACTCGGTTCCGAGTC
## 16 19.04080 2.18e-07  0.2120 ATAGACTCGGAACCGAGTA
```

Annotate FIMO results: "where are the motif hits across the genome?"
identify and annotate all motifs that are in intergenic/promoter regions                

```r
GenomicRanges::findOverlaps(genes.only, gr, ignore.strand = T, minoverlap = 5) -> genes
GenomicRanges::findOverlaps(trimmed.pro, gr, ignore.strand = T, minoverlap = 5) -> promoters
GenomicRanges::findOverlaps(ig.only, gr, ignore.strand = T, minoverlap = 5) -> ig.regions

# get IRanges from hits objects and add informative metadata
genelist <- gr[subjectHits(genes)][, -3]
genelist$type <- rep("gene", length(genes))
strand(genelist) <- strand(genes.only[queryHits(genes)])
genelist$feature_start <- start(genes.only[queryHits(genes)])
genelist$feature_end <- end(genes.only[queryHits(genes)])
genelist$GENEID <- as.character(genes.only$GENEID[queryHits(genes)])

prolist <- gr[subjectHits(promoters)][, -3]
prolist$type <- rep("promoter", length(promoters))
prolist$feature_start <- start(trimmed.pro[queryHits(promoters)])
prolist$feature_end <- end(trimmed.pro[queryHits(promoters)])
strand(prolist) <- strand(trimmed.pro[queryHits(promoters)])
prolist$GENEID <- as.character(trimmed.pro$GENEID[queryHits(promoters)])

iglist <- gr[subjectHits(ig.regions)][, -3]
iglist$type <- rep("intergenic", length(ig.regions))
iglist$feature_start <- start(ig.only[queryHits(ig.regions)])
iglist$feature_end <- end(ig.only[queryHits(ig.regions)])
iglist$GENEID <- NA

# convert separate IRanges to Data frames
seqs <- seq(1, length(genes))
as.data.frame(prolist) -> one
rownames(one) <- NULL
as.data.frame(genelist, row.names(seqs)) -> two
rownames(two) <- NULL
as.data.frame(iglist, row.names(seqs)) -> three
rownames(three) <- NULL

# combine dfs (gene hits and promoter hits)
final <- rbind(one, two, three) %>% distinct(.keep_all = T)
colnames(final)[c(2, 3, 14)] <- c("motif_start", "motif_end", "locus_tag")

# merge with gff information (get NCBI annotations and locus names)
gff.df[gff.df$locus_tag %in% final$locus_tag, ] -> tmp
tmp[c(2, 3, 4, 10)] -> tmp2
(left_join(final, tmp2) %>% arrange(p.value) -> final)
```

```
## Joining with `by = join_by(locus_tag)`
```

```
##        seqnames motif_start motif_end width strand            motif_id
## 1   NC_015948.1     1315521   1315539    19      - ATTHACTMGAAWMBKWGTA
## 2   NC_015948.1     1315521   1315539    19      + ATTHACTMGAAWMBKWGTA
## 3   NC_015948.1     1358582   1358600    19      - ATTHACTMGAAWMBKWGTA
## 4   NC_015948.1     1358582   1358600    19      + ATTHACTMGAAWMBKWGTA
## 5   NC_015948.1     2237586   2237604    19      + ATTHACTMGAAWMBKWGTA
## 6   NC_015948.1      832389    832407    19      + ATTHACTMGAAWMBKWGTA
## 7   NC_015948.1     1792844   1792862    19      + ATTHACTMGAAWMBKWGTA
## 8   NC_015944.1      133790    133808    19      - ATTHACTMGAAWMBKWGTA
## 9   NC_015948.1      168919    168937    19      - ATTHACTMGAAWMBKWGTA
## 10  NC_015948.1     2701211   2701229    19      + ATTHACTMGAAWMBKWGTA
## 11  NC_015948.1     2701211   2701229    19      - ATTHACTMGAAWMBKWGTA
## 12  NC_015948.1      458400    458418    19      + ATTHACTMGAAWMBKWGTA
## 13  NC_015948.1     2701208   2701226    19      + ATTHACTMGAAWMBKWGTA
## 14  NC_015948.1     2701208   2701226    19      - ATTHACTMGAAWMBKWGTA
## 15  NC_015948.1     1358585   1358603    19      - ATTHACTMGAAWMBKWGTA
## 16  NC_015948.1     1358585   1358603    19      + ATTHACTMGAAWMBKWGTA
## 17  NC_015948.1     1211529   1211547    19      - ATTHACTMGAAWMBKWGTA
## 18  NC_015948.1     2635523   2635541    19      - ATTHACTMGAAWMBKWGTA
## 19  NC_015948.1     2635523   2635541    19      + ATTHACTMGAAWMBKWGTA
## 20  NC_015948.1     1631961   1631979    19      - ATTHACTMGAAWMBKWGTA
## 21  NC_015944.1      134417    134435    19      - ATTHACTMGAAWMBKWGTA
## 22  NC_015948.1       93670     93688    19      + ATTHACTMGAAWMBKWGTA
## 23  NC_015944.1      203356    203374    19      + ATTHACTMGAAWMBKWGTA
## 24  NC_015948.1     1367594   1367612    19      - ATTHACTMGAAWMBKWGTA
## 25  NC_015948.1     2563607   2563625    19      + ATTHACTMGAAWMBKWGTA
## 26  NC_015948.1     1235691   1235709    19      + ATTHACTMGAAWMBKWGTA
## 27  NC_015948.1      274240    274258    19      - ATTHACTMGAAWMBKWGTA
## 28  NC_015944.1       95399     95417    19      + ATTHACTMGAAWMBKWGTA
## 29  NC_015944.1       41559     41577    19      + ATTHACTMGAAWMBKWGTA
## 30  NC_015944.1       41559     41577    19      + ATTHACTMGAAWMBKWGTA
## 31  NC_015943.1       61159     61177    19      - ATTHACTMGAAWMBKWGTA
## 32  NC_015943.1       53060     53078    19      - ATTHACTMGAAWMBKWGTA
## 33  NC_015943.1       21242     21260    19      - ATTHACTMGAAWMBKWGTA
## 34  NC_015944.1      249097    249115    19      + ATTHACTMGAAWMBKWGTA
## 35  NC_015948.1     2028288   2028306    19      + ATTHACTMGAAWMBKWGTA
## 36  NC_015948.1      645912    645930    19      + ATTHACTMGAAWMBKWGTA
## 37  NC_015948.1      919561    919579    19      - ATTHACTMGAAWMBKWGTA
## 38  NC_015948.1      650367    650385    19      + ATTHACTMGAAWMBKWGTA
## 39  NC_015948.1      650367    650385    19      + ATTHACTMGAAWMBKWGTA
## 40  NC_015948.1     2011236   2011254    19      - ATTHACTMGAAWMBKWGTA
## 41  NC_015944.1      134420    134438    19      - ATTHACTMGAAWMBKWGTA
## 42  NC_015948.1     2166066   2166084    19      + ATTHACTMGAAWMBKWGTA
## 43  NC_015944.1      342697    342715    19      + ATTHACTMGAAWMBKWGTA
## 44  NC_015943.1       83719     83737    19      - ATTHACTMGAAWMBKWGTA
## 45  NC_015948.1       36857     36875    19      + ATTHACTMGAAWMBKWGTA
## 46  NC_015948.1     1615643   1615661    19      - ATTHACTMGAAWMBKWGTA
## 47  NC_015943.1      143592    143610    19      - ATTHACTMGAAWMBKWGTA
## 48  NC_015948.1     2333886   2333904    19      + ATTHACTMGAAWMBKWGTA
## 49  NC_015948.1     2381756   2381774    19      + ATTHACTMGAAWMBKWGTA
## 50  NC_015944.1       72054     72072    19      - ATTHACTMGAAWMBKWGTA
## 51  NC_015948.1      678977    678995    19      + ATTHACTMGAAWMBKWGTA
## 52  NC_015948.1     2237589   2237607    19      + ATTHACTMGAAWMBKWGTA
## 53  NC_015943.1      129326    129344    19      + ATTHACTMGAAWMBKWGTA
## 54  NC_015948.1      959945    959963    19      - ATTHACTMGAAWMBKWGTA
## 55  NC_015948.1      959945    959963    19      + ATTHACTMGAAWMBKWGTA
## 56  NC_015948.1     2731610   2731628    19      - ATTHACTMGAAWMBKWGTA
## 57  NC_015948.1      168916    168934    19      - ATTHACTMGAAWMBKWGTA
## 58  NC_015944.1      389787    389805    19      - ATTHACTMGAAWMBKWGTA
## 59  NC_015943.1      369947    369965    19      - ATTHACTMGAAWMBKWGTA
## 60  NC_015943.1      369947    369965    19      - ATTHACTMGAAWMBKWGTA
## 61  NC_015943.1      219980    219998    19      - ATTHACTMGAAWMBKWGTA
## 62  NC_015948.1      565351    565369    19      + ATTHACTMGAAWMBKWGTA
## 63  NC_015948.1      386752    386770    19      - ATTHACTMGAAWMBKWGTA
## 64  NC_015948.1     2871799   2871817    19      - ATTHACTMGAAWMBKWGTA
## 65  NC_015948.1     1108645   1108663    19      + ATTHACTMGAAWMBKWGTA
## 66  NC_015948.1     2431509   2431527    19      + ATTHACTMGAAWMBKWGTA
## 67  NC_015948.1     2700128   2700146    19      + ATTHACTMGAAWMBKWGTA
## 68  NC_015948.1      612801    612819    19      - ATTHACTMGAAWMBKWGTA
## 69  NC_015948.1     2338435   2338453    19      - ATTHACTMGAAWMBKWGTA
## 70  NC_015948.1     1770534   1770552    19      - ATTHACTMGAAWMBKWGTA
## 71  NC_015948.1     1770534   1770552    19      + ATTHACTMGAAWMBKWGTA
## 72  NC_015948.1      136839    136857    19      + ATTHACTMGAAWMBKWGTA
## 73  NC_015948.1     1478939   1478957    19      - ATTHACTMGAAWMBKWGTA
## 74  NC_015943.1       77536     77554    19      - ATTHACTMGAAWMBKWGTA
## 75  NC_015948.1      443969    443987    19      + ATTHACTMGAAWMBKWGTA
## 76  NC_015948.1      857709    857727    19      + ATTHACTMGAAWMBKWGTA
## 77  NC_015948.1     1323226   1323244    19      + ATTHACTMGAAWMBKWGTA
## 78  NC_015948.1     1570121   1570139    19      + ATTHACTMGAAWMBKWGTA
## 79  NC_015948.1     1131513   1131531    19      - ATTHACTMGAAWMBKWGTA
## 80  NC_015944.1      330326    330344    19      - ATTHACTMGAAWMBKWGTA
## 81  NC_015948.1     1196867   1196885    19      + ATTHACTMGAAWMBKWGTA
## 82  NC_015943.1       71201     71219    19      - ATTHACTMGAAWMBKWGTA
## 83  NC_015948.1     2366656   2366674    19      - ATTHACTMGAAWMBKWGTA
## 84  NC_015948.1     2366656   2366674    19      - ATTHACTMGAAWMBKWGTA
## 85  NC_015948.1     1487720   1487738    19      - ATTHACTMGAAWMBKWGTA
## 86  NC_015948.1      513562    513580    19      - ATTHACTMGAAWMBKWGTA
## 87  NC_015948.1      513562    513580    19      + ATTHACTMGAAWMBKWGTA
## 88  NC_015943.1      470064    470082    19      - ATTHACTMGAAWMBKWGTA
## 89  NC_015944.1      142296    142314    19      + ATTHACTMGAAWMBKWGTA
## 90  NC_015948.1     2543519   2543537    19      - ATTHACTMGAAWMBKWGTA
## 91  NC_015943.1      288184    288202    19      + ATTHACTMGAAWMBKWGTA
## 92  NC_015944.1      359044    359062    19      + ATTHACTMGAAWMBKWGTA
## 93  NC_015948.1     1726825   1726843    19      - ATTHACTMGAAWMBKWGTA
## 94  NC_015948.1     1726825   1726843    19      - ATTHACTMGAAWMBKWGTA
## 95  NC_015948.1      186931    186949    19      + ATTHACTMGAAWMBKWGTA
## 96  NC_015948.1     2529170   2529188    19      - ATTHACTMGAAWMBKWGTA
## 97  NC_015948.1     2529170   2529188    19      + ATTHACTMGAAWMBKWGTA
## 98  NC_015943.1      345794    345812    19      - ATTHACTMGAAWMBKWGTA
## 99  NC_015948.1     2921144   2921162    19      + ATTHACTMGAAWMBKWGTA
## 100 NC_015948.1     1280977   1280995    19      - ATTHACTMGAAWMBKWGTA
## 101 NC_015948.1     1280977   1280995    19      - ATTHACTMGAAWMBKWGTA
## 102 NC_015944.1       95098     95116    19      + ATTHACTMGAAWMBKWGTA
## 103 NC_015944.1      121472    121490    19      + ATTHACTMGAAWMBKWGTA
## 104 NC_015948.1     1489835   1489853    19      + ATTHACTMGAAWMBKWGTA
## 105 NC_015943.1       41082     41100    19      - ATTHACTMGAAWMBKWGTA
## 106 NC_015948.1     1029591   1029609    19      - ATTHACTMGAAWMBKWGTA
## 107 NC_015944.1      139092    139110    19      + ATTHACTMGAAWMBKWGTA
## 108 NC_015948.1      501532    501550    19      - ATTHACTMGAAWMBKWGTA
## 109 NC_015948.1     2186840   2186858    19      + ATTHACTMGAAWMBKWGTA
## 110 NC_015948.1     1999115   1999133    19      - ATTHACTMGAAWMBKWGTA
## 111 NC_015948.1     1877174   1877192    19      + ATTHACTMGAAWMBKWGTA
## 112 NC_015944.1      321296    321314    19      - ATTHACTMGAAWMBKWGTA
## 113 NC_015948.1      731220    731238    19      - ATTHACTMGAAWMBKWGTA
## 114 NC_015943.1      177961    177979    19      - ATTHACTMGAAWMBKWGTA
## 115 NC_015948.1     2118791   2118809    19      - ATTHACTMGAAWMBKWGTA
## 116 NC_015948.1     1020219   1020237    19      - ATTHACTMGAAWMBKWGTA
## 117 NC_015948.1     2492483   2492501    19      - ATTHACTMGAAWMBKWGTA
## 118 NC_015948.1     1020219   1020237    19      + ATTHACTMGAAWMBKWGTA
## 119 NC_015948.1     1268623   1268641    19      - ATTHACTMGAAWMBKWGTA
## 120 NC_015944.1      152398    152416    19      - ATTHACTMGAAWMBKWGTA
## 121 NC_015948.1     1534488   1534506    19      - ATTHACTMGAAWMBKWGTA
## 122 NC_015948.1     1924336   1924354    19      - ATTHACTMGAAWMBKWGTA
## 123 NC_015948.1      716018    716036    19      - ATTHACTMGAAWMBKWGTA
## 124 NC_015948.1     1486724   1486742    19      - ATTHACTMGAAWMBKWGTA
## 125 NC_015944.1       87858     87876    19      + ATTHACTMGAAWMBKWGTA
## 126 NC_015944.1      281378    281396    19      - ATTHACTMGAAWMBKWGTA
## 127 NC_015943.1      300267    300285    19      - ATTHACTMGAAWMBKWGTA
## 128 NC_015948.1     2324555   2324573    19      + ATTHACTMGAAWMBKWGTA
## 129 NC_015948.1     2751053   2751071    19      - ATTHACTMGAAWMBKWGTA
## 130 NC_015948.1     1604752   1604770    19      - ATTHACTMGAAWMBKWGTA
## 131 NC_015948.1     1510096   1510114    19      + ATTHACTMGAAWMBKWGTA
## 132 NC_015948.1      897145    897163    19      + ATTHACTMGAAWMBKWGTA
## 133 NC_015943.1      242875    242893    19      - ATTHACTMGAAWMBKWGTA
## 134 NC_015948.1     1331672   1331690    19      + ATTHACTMGAAWMBKWGTA
## 135 NC_015944.1      173025    173043    19      + ATTHACTMGAAWMBKWGTA
## 136 NC_015948.1     2791033   2791051    19      - ATTHACTMGAAWMBKWGTA
## 137 NC_015944.1      358648    358666    19      + ATTHACTMGAAWMBKWGTA
## 138 NC_015948.1     1187534   1187552    19      - ATTHACTMGAAWMBKWGTA
## 139 NC_015943.1      380467    380485    19      - ATTHACTMGAAWMBKWGTA
## 140 NC_015948.1     2713530   2713548    19      + ATTHACTMGAAWMBKWGTA
## 141 NC_015944.1      222535    222553    19      - ATTHACTMGAAWMBKWGTA
## 142 NC_015944.1      200972    200990    19      + ATTHACTMGAAWMBKWGTA
## 143 NC_015944.1      314678    314696    19      - ATTHACTMGAAWMBKWGTA
## 144 NC_015943.1      432172    432190    19      - ATTHACTMGAAWMBKWGTA
## 145 NC_015948.1      416758    416776    19      + ATTHACTMGAAWMBKWGTA
## 146 NC_015944.1      175901    175919    19      + ATTHACTMGAAWMBKWGTA
## 147 NC_015948.1     1868450   1868468    19      - ATTHACTMGAAWMBKWGTA
## 148 NC_015943.1      311832    311850    19      - ATTHACTMGAAWMBKWGTA
## 149 NC_015948.1      420783    420801    19      + ATTHACTMGAAWMBKWGTA
## 150 NC_015948.1     2532039   2532057    19      - ATTHACTMGAAWMBKWGTA
## 151 NC_015948.1     2681049   2681067    19      - ATTHACTMGAAWMBKWGTA
## 152 NC_015944.1      195237    195255    19      - ATTHACTMGAAWMBKWGTA
## 153 NC_015944.1      195237    195255    19      + ATTHACTMGAAWMBKWGTA
## 154 NC_015948.1     2460574   2460592    19      + ATTHACTMGAAWMBKWGTA
## 155 NC_015948.1     2454939   2454957    19      + ATTHACTMGAAWMBKWGTA
## 156 NC_015948.1     2278160   2278178    19      - ATTHACTMGAAWMBKWGTA
## 157 NC_015948.1     2278160   2278178    19      + ATTHACTMGAAWMBKWGTA
## 158 NC_015948.1     1794597   1794615    19      + ATTHACTMGAAWMBKWGTA
## 159 NC_015948.1      349054    349072    19      + ATTHACTMGAAWMBKWGTA
## 160 NC_015948.1     2128188   2128206    19      - ATTHACTMGAAWMBKWGTA
## 161 NC_015944.1      186779    186797    19      + ATTHACTMGAAWMBKWGTA
## 162 NC_015948.1      813021    813039    19      - ATTHACTMGAAWMBKWGTA
## 163 NC_015948.1      420780    420798    19      + ATTHACTMGAAWMBKWGTA
## 164 NC_015948.1     1125970   1125988    19      - ATTHACTMGAAWMBKWGTA
## 165 NC_015948.1     1386624   1386642    19      - ATTHACTMGAAWMBKWGTA
## 166 NC_015948.1     1057622   1057640    19      + ATTHACTMGAAWMBKWGTA
## 167 NC_015948.1     2053768   2053786    19      + ATTHACTMGAAWMBKWGTA
## 168 NC_015948.1     1925853   1925871    19      - ATTHACTMGAAWMBKWGTA
## 169 NC_015948.1     2235592   2235610    19      - ATTHACTMGAAWMBKWGTA
## 170 NC_015943.1       37215     37233    19      - ATTHACTMGAAWMBKWGTA
## 171 NC_015948.1     2182662   2182680    19      + ATTHACTMGAAWMBKWGTA
## 172 NC_015943.1      168130    168148    19      + ATTHACTMGAAWMBKWGTA
## 173 NC_015948.1     2929634   2929652    19      - ATTHACTMGAAWMBKWGTA
## 174 NC_015943.1      141612    141630    19      - ATTHACTMGAAWMBKWGTA
## 175 NC_015948.1     2344407   2344425    19      + ATTHACTMGAAWMBKWGTA
## 176 NC_015948.1     2047212   2047230    19      + ATTHACTMGAAWMBKWGTA
## 177 NC_015948.1      874985    875003    19      - ATTHACTMGAAWMBKWGTA
## 178 NC_015948.1      653501    653519    19      + ATTHACTMGAAWMBKWGTA
## 179 NC_015948.1     2153292   2153310    19      - ATTHACTMGAAWMBKWGTA
## 180 NC_015948.1     1091827   1091845    19      - ATTHACTMGAAWMBKWGTA
## 181 NC_015943.1      486890    486908    19      - ATTHACTMGAAWMBKWGTA
## 182 NC_015943.1      112048    112066    19      - ATTHACTMGAAWMBKWGTA
## 183 NC_015943.1      486890    486908    19      + ATTHACTMGAAWMBKWGTA
## 184 NC_015948.1     2072406   2072424    19      - ATTHACTMGAAWMBKWGTA
## 185 NC_015943.1       34297     34315    19      + ATTHACTMGAAWMBKWGTA
## 186 NC_015948.1     1844219   1844237    19      + ATTHACTMGAAWMBKWGTA
## 187 NC_015948.1     2107765   2107783    19      - ATTHACTMGAAWMBKWGTA
## 188 NC_015943.1      115171    115189    19      - ATTHACTMGAAWMBKWGTA
## 189 NC_015944.1       87861     87879    19      + ATTHACTMGAAWMBKWGTA
## 190 NC_015948.1      878980    878998    19      - ATTHACTMGAAWMBKWGTA
## 191 NC_015948.1     2354880   2354898    19      + ATTHACTMGAAWMBKWGTA
## 192 NC_015948.1     2876960   2876978    19      + ATTHACTMGAAWMBKWGTA
## 193 NC_015948.1     1662676   1662694    19      - ATTHACTMGAAWMBKWGTA
## 194 NC_015943.1      330415    330433    19      + ATTHACTMGAAWMBKWGTA
## 195 NC_015948.1      291956    291974    19      - ATTHACTMGAAWMBKWGTA
## 196 NC_015943.1      273634    273652    19      + ATTHACTMGAAWMBKWGTA
## 197 NC_015948.1      351329    351347    19      - ATTHACTMGAAWMBKWGTA
## 198 NC_015944.1       19035     19053    19      + ATTHACTMGAAWMBKWGTA
## 199 NC_015944.1      176354    176372    19      + ATTHACTMGAAWMBKWGTA
## 200 NC_015948.1     1177748   1177766    19      + ATTHACTMGAAWMBKWGTA
## 201 NC_015948.1     2012683   2012701    19      + ATTHACTMGAAWMBKWGTA
## 202 NC_015948.1      423608    423626    19      + ATTHACTMGAAWMBKWGTA
## 203 NC_015948.1     1701957   1701975    19      + ATTHACTMGAAWMBKWGTA
## 204 NC_015948.1     2448981   2448999    19      + ATTHACTMGAAWMBKWGTA
## 205 NC_015948.1     1244004   1244022    19      + ATTHACTMGAAWMBKWGTA
## 206 NC_015948.1     1914236   1914254    19      - ATTHACTMGAAWMBKWGTA
## 207 NC_015944.1      318318    318336    19      + ATTHACTMGAAWMBKWGTA
## 208 NC_015943.1       75739     75757    19      - ATTHACTMGAAWMBKWGTA
## 209 NC_015948.1     1582075   1582093    19      - ATTHACTMGAAWMBKWGTA
## 210 NC_015948.1     2492486   2492504    19      - ATTHACTMGAAWMBKWGTA
## 211 NC_015948.1     2168481   2168499    19      + ATTHACTMGAAWMBKWGTA
## 212 NC_015943.1      193422    193440    19      - ATTHACTMGAAWMBKWGTA
## 213 NC_015948.1     1303569   1303587    19      - ATTHACTMGAAWMBKWGTA
## 214 NC_015948.1     1980645   1980663    19      + ATTHACTMGAAWMBKWGTA
## 215 NC_015944.1      358352    358370    19      - ATTHACTMGAAWMBKWGTA
## 216 NC_015948.1     1714671   1714689    19      - ATTHACTMGAAWMBKWGTA
## 217 NC_015948.1      732642    732660    19      - ATTHACTMGAAWMBKWGTA
## 218 NC_015948.1     2148640   2148658    19      + ATTHACTMGAAWMBKWGTA
## 219 NC_015948.1     2797151   2797169    19      + ATTHACTMGAAWMBKWGTA
## 220 NC_015944.1        3212      3230    19      - ATTHACTMGAAWMBKWGTA
## 221 NC_015944.1      308379    308397    19      - ATTHACTMGAAWMBKWGTA
## 222 NC_015944.1      356860    356878    19      - ATTHACTMGAAWMBKWGTA
## 223 NC_015943.1      395959    395977    19      - ATTHACTMGAAWMBKWGTA
## 224 NC_015943.1      169286    169304    19      - ATTHACTMGAAWMBKWGTA
## 225 NC_015948.1      223769    223787    19      + ATTHACTMGAAWMBKWGTA
## 226 NC_015948.1     1573966   1573984    19      + ATTHACTMGAAWMBKWGTA
## 227 NC_015948.1     2427247   2427265    19      - ATTHACTMGAAWMBKWGTA
## 228 NC_015948.1     1741364   1741382    19      + ATTHACTMGAAWMBKWGTA
## 229 NC_015943.1      149827    149845    19      - ATTHACTMGAAWMBKWGTA
## 230 NC_015948.1     2036523   2036541    19      - ATTHACTMGAAWMBKWGTA
## 231 NC_015943.1      325242    325260    19      - ATTHACTMGAAWMBKWGTA
## 232 NC_015948.1      586831    586849    19      - ATTHACTMGAAWMBKWGTA
## 233 NC_015943.1      270411    270429    19      - ATTHACTMGAAWMBKWGTA
## 234 NC_015943.1      295925    295943    19      + ATTHACTMGAAWMBKWGTA
## 235 NC_015944.1       36156     36174    19      - ATTHACTMGAAWMBKWGTA
## 236 NC_015948.1     1605125   1605143    19      - ATTHACTMGAAWMBKWGTA
## 237 NC_015948.1     1921001   1921019    19      - ATTHACTMGAAWMBKWGTA
## 238 NC_015943.1       72895     72913    19      + ATTHACTMGAAWMBKWGTA
## 239 NC_015948.1     1381996   1382014    19      - ATTHACTMGAAWMBKWGTA
## 240 NC_015948.1     2184403   2184421    19      + ATTHACTMGAAWMBKWGTA
## 241 NC_015943.1       76904     76922    19      - ATTHACTMGAAWMBKWGTA
## 242 NC_015948.1      100274    100292    19      - ATTHACTMGAAWMBKWGTA
## 243 NC_015943.1      405373    405391    19      - ATTHACTMGAAWMBKWGTA
## 244 NC_015948.1       99983    100001    19      + ATTHACTMGAAWMBKWGTA
## 245 NC_015948.1     2222261   2222279    19      - ATTHACTMGAAWMBKWGTA
## 246 NC_015943.1       83758     83776    19      + ATTHACTMGAAWMBKWGTA
## 247 NC_015948.1     1042320   1042338    19      + ATTHACTMGAAWMBKWGTA
## 248 NC_015948.1     1759995   1760013    19      + ATTHACTMGAAWMBKWGTA
## 249 NC_015948.1     1785734   1785752    19      - ATTHACTMGAAWMBKWGTA
## 250 NC_015943.1       58497     58515    19      - ATTHACTMGAAWMBKWGTA
## 251 NC_015944.1      285361    285379    19      - ATTHACTMGAAWMBKWGTA
## 252 NC_015943.1      462932    462950    19      - ATTHACTMGAAWMBKWGTA
## 253 NC_015943.1      471606    471624    19      - ATTHACTMGAAWMBKWGTA
## 254 NC_015948.1     1327853   1327871    19      - ATTHACTMGAAWMBKWGTA
## 255 NC_015948.1      331267    331285    19      + ATTHACTMGAAWMBKWGTA
##     motif_alt_id  p.value  q.value    matched_sequence       type feature_start
## 1         MEME-1 8.92e-12 6.94e-05 ATTCACTCGAAACCGAGTA   promoter       1315395
## 2         MEME-1 8.92e-12 6.94e-05 ATTCACTCGAAACCGAGTA   promoter       1315395
## 3         MEME-1 2.77e-09 1.08e-02 ATTCACTCGAAATGTAGTA   promoter       1358460
## 4         MEME-1 2.77e-09 1.08e-02 ATTCACTCGAAATGTAGTA   promoter       1358479
## 5         MEME-1 5.26e-09 1.36e-02 ATTCACTCGAAATCGAGTG   promoter       2237553
## 6         MEME-1 1.28e-08 2.48e-02 ATTCACTCGGTTCCGAGTC   promoter        832219
## 7         MEME-1 4.26e-08 6.63e-02 TTTCACTCGATTCGGAGTC       gene       1792840
## 8         MEME-1 5.83e-08 7.56e-02 ATGAACTCGAATCTGTGTA   promoter        133780
## 9         MEME-1 1.46e-07 1.62e-01 TTTTACTCGATTCCGTGTC   promoter        168880
## 10        MEME-1 2.18e-07 2.12e-01 ATAGACTCGGAACCGAGTA   promoter       2701155
## 11        MEME-1 2.18e-07 2.12e-01 ATAGACTCGGAACCGAGTA   promoter       2701155
## 12        MEME-1 2.60e-07 2.25e-01 GTTCACTACGAACGGTGTA intergenic        458154
## 13        MEME-1 6.31e-07 4.91e-01 AATTACTCGGTTCCGAGTC   promoter       2701155
## 14        MEME-1 6.31e-07 4.91e-01 AATTACTCGGTTCCGAGTC   promoter       2701155
## 15        MEME-1 9.71e-07 6.86e-01 ATTTACTACATTTCGAGTG   promoter       1358460
## 16        MEME-1 9.71e-07 6.86e-01 ATTTACTACATTTCGAGTG   promoter       1358479
## 17        MEME-1 1.84e-06 1.00e+00 ACTTACTACAAAAGTAGTA   promoter       1211488
## 18        MEME-1 2.19e-06 1.00e+00 ATTAACTCACATAGTTGTA   promoter       2635404
## 19        MEME-1 2.19e-06 1.00e+00 ATTAACTCACATAGTTGTA   promoter       2635404
## 20        MEME-1 2.74e-06 1.00e+00 ACACACTCGAAACAGTGTA       gene       1630541
## 21        MEME-1 2.93e-06 1.00e+00 GTTTACTATAATTATAGTA   promoter        134374
## 22        MEME-1 3.62e-06 1.00e+00 ACTACCTCGATACGGTGTC       gene         93258
## 23        MEME-1 3.99e-06 1.00e+00 ACGAACTCGAACTGGAGTA       gene        201255
## 24        MEME-1 4.76e-06 1.00e+00 ATGTACTCCAGCACGAGTA       gene       1367197
## 25        MEME-1 5.25e-06 1.00e+00 GTTACTTCTATACTGAGTA intergenic       2563289
## 26        MEME-1 5.47e-06 1.00e+00 ACGCACGCGAGACCGAGTA       gene       1235432
## 27        MEME-1 5.66e-06 1.00e+00 TTGCCCTCGGTTCCGAGTA       gene        273993
## 28        MEME-1 6.93e-06 1.00e+00 GTGTATTCGAGTCGTAGTC intergenic         95250
## 29        MEME-1 7.29e-06 1.00e+00 ATTCACGAGTGACGGTGTC   promoter         41478
## 30        MEME-1 7.29e-06 1.00e+00 ATTCACGAGTGACGGTGTC   promoter         41456
## 31        MEME-1 7.49e-06 1.00e+00 ATTCAGTCAGAAACGTGTG       gene         61043
## 32        MEME-1 7.64e-06 1.00e+00 ATTAATGCTGTTTGGAGTA       gene         52829
## 33        MEME-1 7.94e-06 1.00e+00 ATTCAGTAAACACGGAGTC       gene         20769
## 34        MEME-1 8.98e-06 1.00e+00 GCGAAGTCGAAACCGTGTA       gene        248916
## 35        MEME-1 9.50e-06 1.00e+00 ACGCACTCGCTATCGAGTA       gene       2028182
## 36        MEME-1 9.50e-06 1.00e+00 AAAAATTCGGTTCGGAGTA intergenic        645894
## 37        MEME-1 9.64e-06 1.00e+00 ACTATTTCGAGACGGAGTC       gene        918846
## 38        MEME-1 9.81e-06 1.00e+00 GTATCCTCGGATATGAGTA   promoter        650321
## 39        MEME-1 9.81e-06 1.00e+00 GTATCCTCGGATATGAGTA       gene        650379
## 40        MEME-1 1.02e-05 1.00e+00 ATTACGTCCAGACCGAGTG   promoter       2011216
## 41        MEME-1 1.04e-05 1.00e+00 AATTACTATAATTATAGTA   promoter        134374
## 42        MEME-1 1.08e-05 1.00e+00 TTTCCGTCGAGAACGTGTA       gene       2165832
## 43        MEME-1 1.10e-05 1.00e+00 ATGTACTCGGCACTTTGTC       gene        342387
## 44        MEME-1 1.11e-05 1.00e+00 ATACCGTCGATACTGAGTA intergenic         83650
## 45        MEME-1 1.25e-05 1.00e+00 TTTCACGCAAACTGGAGTA       gene         35470
## 46        MEME-1 1.26e-05 1.00e+00 ACACACGATAAACCGAGTA       gene       1615461
## 47        MEME-1 1.35e-05 1.00e+00 GCTGACTCGTTACGGAGTC       gene        143197
## 48        MEME-1 1.37e-05 1.00e+00 GTGACCACGATACCGAGTA       gene       2333000
## 49        MEME-1 1.37e-05 1.00e+00 GTTCACACCGAACTTTGTC       gene       2381253
## 50        MEME-1 1.39e-05 1.00e+00 AAGAACTCGACAAGGAGTG       gene         71744
## 51        MEME-1 1.40e-05 1.00e+00 ACTCCCTCAAGACGGAGTC       gene        677626
## 52        MEME-1 1.41e-05 1.00e+00 CTTCACTCGATTTCGAGTG   promoter       2237553
## 53        MEME-1 1.42e-05 1.00e+00 GTTCTGTCGAGATGGAGTC       gene        128600
## 54        MEME-1 1.47e-05 1.00e+00 GTTATTTCGGCTCCGAGTG   promoter        959873
## 55        MEME-1 1.47e-05 1.00e+00 GTTATTTCGGCTCCGAGTG   promoter        959873
## 56        MEME-1 1.47e-05 1.00e+00 ACTACGGCGAATCCGAGTA       gene       2730814
## 57        MEME-1 1.50e-05 1.00e+00 GCTGACACGGAATCGAGTA   promoter        168880
## 58        MEME-1 1.53e-05 1.00e+00 ACTCCCTCTAACAGTAGTA intergenic        389724
## 59        MEME-1 1.56e-05 1.00e+00 ACTTTCTCGAACATTAGTA       gene        368182
## 60        MEME-1 1.56e-05 1.00e+00 ACTTTCTCGAACATTAGTA       gene        369960
## 61        MEME-1 1.58e-05 1.00e+00 GTGTATTCGACTCGTAGTG       gene        219463
## 62        MEME-1 1.65e-05 1.00e+00 GCGGATTCGAAACCGAGTC       gene        563211
## 63        MEME-1 1.70e-05 1.00e+00 ACAAACTCGACATGGAGTG       gene        386385
## 64        MEME-1 1.85e-05 1.00e+00 GTTACCTCGGACAGTTGTC intergenic       2871761
## 65        MEME-1 1.88e-05 1.00e+00 TTTGTCGCGGAACTGAGTA       gene       1108393
## 66        MEME-1 1.92e-05 1.00e+00 ATTCAGTACCTTTTGAGTA       gene       2431174
## 67        MEME-1 1.92e-05 1.00e+00 GTGTTCTCCAGAACGAGTA       gene       2699806
## 68        MEME-1 2.00e-05 1.00e+00 ATGTATTCGAGCATGAGTC       gene        612638
## 69        MEME-1 2.07e-05 1.00e+00 GTTCACACAGGTCGGTGTC       gene       2337256
## 70        MEME-1 2.14e-05 1.00e+00 GCTCTCTCGGGTAGTAGTA   promoter       1770535
## 71        MEME-1 2.14e-05 1.00e+00 GCTCTCTCGGGTAGTAGTA   promoter       1770535
## 72        MEME-1 2.16e-05 1.00e+00 ATTTCCGAGTATCGGTGTC       gene        136399
## 73        MEME-1 2.17e-05 1.00e+00 ATTTACTCGCATTCGATTA   promoter       1478796
## 74        MEME-1 2.22e-05 1.00e+00 GTTCCGTCGGGTTCGAGTC       gene         76484
## 75        MEME-1 2.23e-05 1.00e+00 ACGACTTCGAGACGGAGTA       gene        443689
## 76        MEME-1 2.28e-05 1.00e+00 ACGAACTCGACACAGAGTG       gene        857123
## 77        MEME-1 2.30e-05 1.00e+00 TCGTTCTCGAAACCGAGTC       gene       1323135
## 78        MEME-1 2.32e-05 1.00e+00 ATTAATAACAATAGTTGTC   promoter       1570087
## 79        MEME-1 2.37e-05 1.00e+00 GTAAACGCCAGTCCGAGTG       gene       1130954
## 80        MEME-1 2.39e-05 1.00e+00 GTGATTTCGAACCCGAGTC       gene        329896
## 81        MEME-1 2.45e-05 1.00e+00 ATTAACTATGAACGCTGTA   promoter       1196816
## 82        MEME-1 2.45e-05 1.00e+00 ACTTCGTCGAAATCGTGTC       gene         70939
## 83        MEME-1 2.48e-05 1.00e+00 ATTGATGCGGATTGTTGTC       gene       2366670
## 84        MEME-1 2.48e-05 1.00e+00 ATTGATGCGGATTGTTGTC intergenic       2366463
## 85        MEME-1 2.55e-05 1.00e+00 GTTAACTCGATGTCGAGTG       gene       1487653
## 86        MEME-1 2.56e-05 1.00e+00 TTGTATTCGGTTAGGAGTG   promoter        513471
## 87        MEME-1 2.56e-05 1.00e+00 TTGTATTCGGTTAGGAGTG   promoter        513471
## 88        MEME-1 2.58e-05 1.00e+00 TTTCTCTCTAGCCAGAGTA       gene        469678
## 89        MEME-1 2.63e-05 1.00e+00 AATCCCTACGATTCGAGTC       gene        140835
## 90        MEME-1 2.67e-05 1.00e+00 ATACAGTCGGCTTCGAGTG       gene       2543035
## 91        MEME-1 2.67e-05 1.00e+00 GATACGTCGGAACCGAGTC       gene        287503
## 92        MEME-1 2.68e-05 1.00e+00 GCTACTTCGGTTCCGAGTC       gene        358863
## 93        MEME-1 2.71e-05 1.00e+00 TTTGACGAGCCACCGAGTA   promoter       1726816
## 94        MEME-1 2.71e-05 1.00e+00 TTTGACGAGCCACCGAGTA   promoter       1726816
## 95        MEME-1 2.71e-05 1.00e+00 GCTTCCTCCAGAACGTGTA       gene        185043
## 96        MEME-1 2.88e-05 1.00e+00 ATTCAGACGGACTGTTGTA   promoter       2529028
## 97        MEME-1 2.88e-05 1.00e+00 ATTCAGACGGACTGTTGTA   promoter       2529028
## 98        MEME-1 3.01e-05 1.00e+00 GTTCACACGCCATGGTGTC       gene        344300
## 99        MEME-1 3.03e-05 1.00e+00 TTTAATTACAGCAGTAGTA   promoter       2921058
## 100       MEME-1 3.17e-05 1.00e+00 ACTGAGACGGATCGGTGTA       gene       1280991
## 101       MEME-1 3.17e-05 1.00e+00 ACTGAGACGGATCGGTGTA intergenic       1280959
## 102       MEME-1 3.18e-05 1.00e+00 ATGATTTCGGTTCGTAGTC       gene         94806
## 103       MEME-1 3.25e-05 1.00e+00 GCTCACGCTGAAACTAGTC       gene        121469
## 104       MEME-1 3.35e-05 1.00e+00 ATTTACAACATCAGTTGTA   promoter       1489662
## 105       MEME-1 3.35e-05 1.00e+00 TTTGAGTCCAATTCTTGTC   promoter         40975
## 106       MEME-1 3.36e-05 1.00e+00 GTTCCGTAGAACAGGAGTG       gene       1028489
## 107       MEME-1 3.43e-05 1.00e+00 TCATAGTCGAATCTGTGTA intergenic        138918
## 108       MEME-1 3.54e-05 1.00e+00 ACGAACACGAGTTCGAGTC       gene        501219
## 109       MEME-1 3.71e-05 1.00e+00 GCTGTCTCGAACAGGTGTA       gene       2186548
## 110       MEME-1 3.73e-05 1.00e+00 ATTACCGCGATCTCGTGTG       gene       1996585
## 111       MEME-1 3.74e-05 1.00e+00 TTTCAGTCCCACCCGTGTC   promoter       1877152
## 112       MEME-1 3.85e-05 1.00e+00 ATTCTTACTGGACTGAGTA intergenic        321144
## 113       MEME-1 3.89e-05 1.00e+00 GTTCCTGCTAATCTGTGTC       gene        730332
## 114       MEME-1 3.99e-05 1.00e+00 ATTGATACTGGTCGGTGTC intergenic        177621
## 115       MEME-1 4.05e-05 1.00e+00 ACGACCTCGAGTTCGAGTG       gene       2118133
## 116       MEME-1 4.07e-05 1.00e+00 ACAAACGCGAACAGTAGTA   promoter       1020167
## 117       MEME-1 4.07e-05 1.00e+00 TATGATTCCGAACCGAGTG   promoter       2492341
## 118       MEME-1 4.07e-05 1.00e+00 ACAAACGCGAACAGTAGTA   promoter       1020167
## 119       MEME-1 4.13e-05 1.00e+00 GATCACTCGTGCCTGAGTG       gene       1268625
## 120       MEME-1 4.17e-05 1.00e+00 TCTCACGCGTATCCTTGTC       gene        151583
## 121       MEME-1 4.29e-05 1.00e+00 TCGAACTCGGCCCCGAGTA       gene       1534330
## 122       MEME-1 4.33e-05 1.00e+00 TCGAAGTCGAAATGGAGTC       gene       1923612
## 123       MEME-1 4.37e-05 1.00e+00 AAGAACGCGAACAGGAGTC       gene        715933
## 124       MEME-1 4.39e-05 1.00e+00 ACTCCTTCGTTTCCGTGTC       gene       1486515
## 125       MEME-1 4.43e-05 1.00e+00 TTGCACGCCAATTCTAGTG   promoter         87852
## 126       MEME-1 4.45e-05 1.00e+00 ACACTCACGAGACCGAGTA intergenic        281339
## 127       MEME-1 4.51e-05 1.00e+00 AAATCGTCGAAAAGGAGTA       gene        299837
## 128       MEME-1 4.56e-05 1.00e+00 TCGTACTCGACCACGAGTA       gene       2323912
## 129       MEME-1 4.60e-05 1.00e+00 ATAAACGCGCCTTCGAGTC       gene       2750830
## 130       MEME-1 4.62e-05 1.00e+00 AAAAACTACGATCAGTGTA       gene       1604240
## 131       MEME-1 4.74e-05 1.00e+00 GATAACGACGAACGGTGTC       gene       1509656
## 132       MEME-1 4.80e-05 1.00e+00 GTGAACGCCGCACCGAGTC       gene        896778
## 133       MEME-1 4.81e-05 1.00e+00 GTTACTTCCCAAATGAGTC   promoter        242797
## 134       MEME-1 4.81e-05 1.00e+00 ACGACCTCTACAACGAGTA       gene       1331359
## 135       MEME-1 4.81e-05 1.00e+00 GGTCACGCGAATCGGTGTC       gene        172981
## 136       MEME-1 4.96e-05 1.00e+00 ATGAAGACGGTAACGAGTC       gene       2790324
## 137       MEME-1 5.01e-05 1.00e+00 TATACCTCGGGACTTTGTA   promoter        358613
## 138       MEME-1 5.06e-05 1.00e+00 ATTTACTCGCTACTGTGAA       gene       1187535
## 139       MEME-1 5.06e-05 1.00e+00 TTATCGTCGGAACGGAGTC       gene        379890
## 140       MEME-1 5.10e-05 1.00e+00 AATTAGTAAGATCTGAGTC       gene       2713137
## 141       MEME-1 5.12e-05 1.00e+00 TTTCAGAATAACAGGAGTA   promoter        222446
## 142       MEME-1 5.14e-05 1.00e+00 TTTTTTTACAGAATGAGTA intergenic        200343
## 143       MEME-1 5.16e-05 1.00e+00 GTACATACGGATTGTTGTA   promoter        314630
## 144       MEME-1 5.16e-05 1.00e+00 ATCAACACGAACCGGTGTA       gene        431102
## 145       MEME-1 5.23e-05 1.00e+00 ACGCCCTCCGCTCCGAGTA       gene        416547
## 146       MEME-1 5.28e-05 1.00e+00 GTAGTTTCCAATACGAGTA   promoter        175839
## 147       MEME-1 5.30e-05 1.00e+00 ATGCTCTCCACAACTAGTG   promoter       1868371
## 148       MEME-1 5.33e-05 1.00e+00 ATAAATTAGGAAAATTGTC intergenic        311212
## 149       MEME-1 5.37e-05 1.00e+00 GTTCACTCGGTTTCGAGGC       gene        420377
## 150       MEME-1 5.48e-05 1.00e+00 ACGAACACGGCACGGAGTC       gene       2531399
## 151       MEME-1 5.50e-05 1.00e+00 GCTGATTCTGGTCGGAGTG       gene       2680464
## 152       MEME-1 5.55e-05 1.00e+00 TAACACTCTCAAACGTGTA   promoter        195194
## 153       MEME-1 5.55e-05 1.00e+00 TAACACTCTCAAACGTGTA   promoter        195194
## 154       MEME-1 5.60e-05 1.00e+00 GCATCCTCGGCTCCGAGTA       gene       2459031
## 155       MEME-1 5.70e-05 1.00e+00 ATGCCCTCGAAACGCTGTA       gene       2454461
## 156       MEME-1 5.75e-05 1.00e+00 TTGTCGTCGGATCGGAGTG   promoter       2277957
## 157       MEME-1 5.75e-05 1.00e+00 TTGTCGTCGGATCGGAGTG   promoter       2278111
## 158       MEME-1 5.92e-05 1.00e+00 ATGACGTCGTGTACGAGTA       gene       1794032
## 159       MEME-1 5.98e-05 1.00e+00 GTGGTCTAGAGAAGGAGTC       gene        348424
## 160       MEME-1 6.08e-05 1.00e+00 ATTAACTCTGTACCTAGAC   promoter       2128156
## 161       MEME-1 6.08e-05 1.00e+00 GTGAACACCGGACGGTGTC intergenic        186718
## 162       MEME-1 6.19e-05 1.00e+00 TCGACTTCGAGTCGGAGTA       gene        813020
## 163       MEME-1 6.28e-05 1.00e+00 TCAGCCTCGAAACCGAGTG       gene        420377
## 164       MEME-1 6.34e-05 1.00e+00 TTGCAGACGGCACCGAGTA       gene       1124925
## 165       MEME-1 6.34e-05 1.00e+00 GCGTCCTCGATACTGAGTG       gene       1386560
## 166       MEME-1 6.34e-05 1.00e+00 ATGTCCACGGGACTGAGTC intergenic       1057536
## 167       MEME-1 6.52e-05 1.00e+00 TCTTATTCCGAAAAGTGTA   promoter       2053738
## 168       MEME-1 6.63e-05 1.00e+00 GAGTACTCAAATTGTAGTC       gene       1923612
## 169       MEME-1 6.63e-05 1.00e+00 TCAAAGTCGACTACGAGTA       gene       2233182
## 170       MEME-1 6.69e-05 1.00e+00 GATTACTCGGAACGTCGTA       gene         36558
## 171       MEME-1 6.76e-05 1.00e+00 ACGGACTCGTCTCCGAGTG       gene       2182169
## 172       MEME-1 6.79e-05 1.00e+00 GTTCTCTCAGGAAATAGTA       gene        167806
## 173       MEME-1 6.84e-05 1.00e+00 TTTCCGTCGTTTCGGTGTG       gene       2928846
## 174       MEME-1 6.84e-05 1.00e+00 TCACACTACGTACCGAGTC       gene        141313
## 175       MEME-1 6.94e-05 1.00e+00 GCGGTCTCGAAATCGAGTC       gene       2344151
## 176       MEME-1 6.94e-05 1.00e+00 GCGAACACGATAAGGAGTG       gene       2046446
## 177       MEME-1 6.97e-05 1.00e+00 TCGAACACGATACTGAGTC       gene        874870
## 178       MEME-1 7.00e-05 1.00e+00 ATTGGCTCGCGTCGGTGTA       gene        653495
## 179       MEME-1 7.00e-05 1.00e+00 ACTCTCGCAGATCCGAGTG       gene       2153113
## 180       MEME-1 7.06e-05 1.00e+00 TCGTCCTCGATTCGGTGTC       gene       1091376
## 181       MEME-1 7.06e-05 1.00e+00 ATGATCTACACTTCGAGTG       gene        486895
## 182       MEME-1 7.06e-05 1.00e+00 TCTCTTTCGCTACTGAGTA intergenic        112024
## 183       MEME-1 7.06e-05 1.00e+00 ATGATCTACACTTCGAGTG intergenic        486699
## 184       MEME-1 7.10e-05 1.00e+00 TCGGACTCGACTTCGAGTC       gene       2071907
## 185       MEME-1 7.10e-05 1.00e+00 AATTATTATCTATGTAGTA       gene         34210
## 186       MEME-1 7.13e-05 1.00e+00 TTTGATTCGATCAAGTGTC       gene       1844202
## 187       MEME-1 7.19e-05 1.00e+00 ACGAACTCGAAAAGCAGTA       gene       2106723
## 188       MEME-1 7.25e-05 1.00e+00 ACTCACAACAACTCTAGTC   promoter        115107
## 189       MEME-1 7.30e-05 1.00e+00 AATCACTAGAATTGGCGTG   promoter         87852
## 190       MEME-1 7.30e-05 1.00e+00 GCTGTCTCCGAACCGTGTG       gene        877389
## 191       MEME-1 7.30e-05 1.00e+00 AAGGCCTCGATATCGAGTC       gene       2354486
## 192       MEME-1 7.30e-05 1.00e+00 GTGTCGTCGAGTACGAGTG       gene       2876794
## 193       MEME-1 7.36e-05 1.00e+00 GCTCTCGCGACTCGGTGTC       gene       1662340
## 194       MEME-1 7.39e-05 1.00e+00 TTTGTGTCGGGACCGTGTC   promoter        330274
## 195       MEME-1 7.45e-05 1.00e+00 GATACCACGGGACGGAGTC       gene        291864
## 196       MEME-1 7.45e-05 1.00e+00 GTGTTCTCGGGCACTAGTA intergenic        272976
## 197       MEME-1 7.59e-05 1.00e+00 ATTTCCTCGAACTGGAGGA       gene        351154
## 198       MEME-1 7.62e-05 1.00e+00 AATTTCTCCAGTCATAGTC       gene         18843
## 199       MEME-1 7.62e-05 1.00e+00 ATACCCTCCGGAATGTGTC       gene        175983
## 200       MEME-1 7.68e-05 1.00e+00 ATTCTGTCCGCAATGAGTC   promoter       1177636
## 201       MEME-1 7.79e-05 1.00e+00 TCTCACTACCGAAGGAGTG       gene       2012637
## 202       MEME-1 7.86e-05 1.00e+00 AATCCCCCGAAATCGAGTA       gene        422674
## 203       MEME-1 7.86e-05 1.00e+00 ACGTTTTCGGATAGGAGTC       gene       1701945
## 204       MEME-1 7.86e-05 1.00e+00 ATTTGGTCGGTATCGAGTA       gene       2448632
## 205       MEME-1 7.89e-05 1.00e+00 TCTACCTCGAAGCGGAGTA       gene       1243499
## 206       MEME-1 8.13e-05 1.00e+00 ATTTCGACGACCCGGTGTA       gene       1913905
## 207       MEME-1 8.22e-05 1.00e+00 TCGAACGCTACACGGAGTA       gene        318146
## 208       MEME-1 8.25e-05 1.00e+00 GTGTAGTCGCAAATTAGTG intergenic         75579
## 209       MEME-1 8.31e-05 1.00e+00 GCGGCCTCGACTCCGTGTA       gene       1582038
## 210       MEME-1 8.35e-05 1.00e+00 AATCACTCGGTTCGGAATC   promoter       2492341
## 211       MEME-1 8.39e-05 1.00e+00 ACTCACTGGAGAACGTGTC       gene       2167467
## 212       MEME-1 8.39e-05 1.00e+00 GAGCAGTCGGTACAGAGTA intergenic        193417
## 213       MEME-1 8.43e-05 1.00e+00 ACGAAGTCGGGATGGAGTG       gene       1302623
## 214       MEME-1 8.43e-05 1.00e+00 ACAAACGCGGCCCGGAGTA       gene       1980221
## 215       MEME-1 8.57e-05 1.00e+00 GTTTTTGACAGTCGGTGTA   promoter        358322
## 216       MEME-1 8.60e-05 1.00e+00 AATTCCTCGAATTCGGGTA       gene       1714658
## 217       MEME-1 8.63e-05 1.00e+00 ATGGTCTCGGTCAGGAGTC       gene        732437
## 218       MEME-1 8.66e-05 1.00e+00 GCGCACGCGACTTGGAGTC       gene       2148124
## 219       MEME-1 8.69e-05 1.00e+00 ACGTCCTCGTGATGGAGTA       gene       2796643
## 220       MEME-1 8.80e-05 1.00e+00 ATTCAGGCGGAACGGATTA       gene          3135
## 221       MEME-1 8.80e-05 1.00e+00 ATAGTCTACTATCGTAGTA intergenic        308009
## 222       MEME-1 8.84e-05 1.00e+00 AGGTACTCCAGTCCGTGTA       gene        355887
## 223       MEME-1 8.88e-05 1.00e+00 TTGCATTCGTCCCCGAGTC       gene        395496
## 224       MEME-1 8.96e-05 1.00e+00 GATCAGTCAAGAAAGAGTA intergenic        169197
## 225       MEME-1 9.00e-05 1.00e+00 ACTTACGAGTGTTCGTGTC       gene        223542
## 226       MEME-1 9.03e-05 1.00e+00 GTTCACTCGGGTTGCTGTC       gene       1572793
## 227       MEME-1 9.03e-05 1.00e+00 GTTCAGTACCACTCGTGTC       gene       2427200
## 228       MEME-1 9.10e-05 1.00e+00 GTACACGCCGGACGTTGTC       gene       1740102
## 229       MEME-1 9.13e-05 1.00e+00 GTTGCTTACGGATGGAGTA   promoter        149733
## 230       MEME-1 9.16e-05 1.00e+00 ATTATTTACCGTCGTAGTG   promoter       2036306
## 231       MEME-1 9.31e-05 1.00e+00 AATTAGAATAATATGTGTA   promoter        325013
## 232       MEME-1 9.31e-05 1.00e+00 TCATCGTCGAATACGAGTA       gene        585306
## 233       MEME-1 9.31e-05 1.00e+00 ATTCGCTCCGTTTGGAGTG       gene        269270
## 234       MEME-1 9.31e-05 1.00e+00 TCTCCCACCAGACGGTGTA       gene        295600
## 235       MEME-1 9.31e-05 1.00e+00 AATACCACGCACAGGAGTA       gene         35936
## 236       MEME-1 9.31e-05 1.00e+00 ATGTTTTCGGATTATTGTA       gene       1604240
## 237       MEME-1 9.39e-05 1.00e+00 TCGAAGTCGAGATGGAGTC       gene       1920283
## 238       MEME-1 9.39e-05 1.00e+00 GTTCCCACTTGACGTTGTA       gene         72517
## 239       MEME-1 9.43e-05 1.00e+00 ACGTCCTCGCCTCCGAGTC       gene       1381152
## 240       MEME-1 9.43e-05 1.00e+00 ATAAACTCGAACAGGTGGA       gene       2182797
## 241       MEME-1 9.51e-05 1.00e+00 GTTCACTCGGGCCCGACTC       gene         76484
## 242       MEME-1 9.55e-05 1.00e+00 GTTCACTCAAGAACGTGGA       gene        100276
## 243       MEME-1 9.59e-05 1.00e+00 GTGGCGTCGGAAACGAGTC       gene        404469
## 244       MEME-1 9.62e-05 1.00e+00 TTGCGCTCGTATCGGAGTA       gene         99463
## 245       MEME-1 9.62e-05 1.00e+00 GTTTCTGCTGTTCCTTGTA       gene       2221995
## 246       MEME-1 9.62e-05 1.00e+00 ATGAATTCACGAACTAGTC intergenic         83650
## 247       MEME-1 9.66e-05 1.00e+00 TTTTCTGCCAGTACGTGTA   promoter       1042098
## 248       MEME-1 9.66e-05 1.00e+00 CCGCACTCGACACCGAGTA       gene       1757921
## 249       MEME-1 9.66e-05 1.00e+00 ATTCGCTACAACACGAGTG       gene       1785435
## 250       MEME-1 9.69e-05 1.00e+00 TTACACTCCGGATTTTGTC       gene         57493
## 251       MEME-1 9.69e-05 1.00e+00 ATGAACGAGGCCCGTAGTC       gene        284694
## 252       MEME-1 9.76e-05 1.00e+00 ACGAACAATACACCGAGTC       gene        462820
## 253       MEME-1 9.80e-05 1.00e+00 ATGCACTCTAGTTCGATTA intergenic        471502
## 254       MEME-1 9.89e-05 1.00e+00 GCTCCCGCGGCTCCTTGTA       gene       1327551
## 255       MEME-1 9.97e-05 1.00e+00 ACTCAGGCTTGACCGTGTA       gene        329515
##     feature_end   locus_tag            acc old_locus_tag
## 1       1315599 HAH_RS06645 WP_044951809.1      HAH_1365
## 2       1315599 HAH_RS06650 WP_014040223.1      HAH_1366
## 3       1358709 HAH_RS06895 WP_014040268.1      HAH_1418
## 4       1358728 HAH_RS06900 WP_014040269.1      HAH_1419
## 5       2237673 HAH_RS11270 WP_014041022.1      HAH_2323
## 6        832468 HAH_RS04315 WP_014039807.1      HAH_0887
## 7       1793430 HAH_RS08975 WP_014040638.1      HAH_1848
## 8        133893 HAH_RS17930 WP_014031120.1      HAH_5129
## 9        168953 HAH_RS19805 WP_004516449.1      HAH_0188
## 10      2701282 HAH_RS13650 WP_014041448.1      HAH_2806
## 11      2701282 HAH_RS20000 WP_008307600.1      HAH_2805
## 12       458466        <NA>           <NA>          <NA>
## 13      2701282 HAH_RS13650 WP_014041448.1      HAH_2806
## 14      2701282 HAH_RS20000 WP_008307600.1      HAH_2805
## 15      1358709 HAH_RS06895 WP_014040268.1      HAH_1418
## 16      1358728 HAH_RS06900 WP_014040269.1      HAH_1419
## 17      1211584 HAH_RS06155 WP_014040138.1      HAH_1264
## 18      2635594 HAH_RS13280 WP_014041387.1      HAH_2729
## 19      2635594 HAH_RS13285 WP_008307747.1      HAH_2730
## 20      1632523 HAH_RS08165 WP_014040499.1      HAH_1681
## 21       134531 HAH_RS17935 WP_023842995.1      HAH_5130
## 22        94133 HAH_RS00485 WP_014039124.1      HAH_0101
## 23       203789 HAH_RS18200 WP_044952813.1      HAH_5182
## 24      1367991 HAH_RS06935 WP_014040275.1      HAH_1427
## 25      2563754        <NA>           <NA>          <NA>
## 26      1237240 HAH_RS06260 WP_014040157.1      HAH_1284
## 27       275075 HAH_RS01395 WP_014039291.1      HAH_0291
## 28        95503        <NA>           <NA>          <NA>
## 29        41581 HAH_RS17470 WP_014031029.1      HAH_5034
## 30        41581 HAH_RS19660 WP_014031028.1      HAH_5033
## 31        61516 HAH_RS15415 WP_023842893.1          <NA>
## 32        53788 HAH_RS15385 WP_014030622.1      HAH_4065
## 33        23633 HAH_RS15225 WP_014030585.1      HAH_4025
## 34       249116 HAH_RS18395 WP_008312239.1      HAH_5224
## 35      2029513 HAH_RS10180 WP_014040843.1      HAH_2101
## 36       645943        <NA>           <NA>          <NA>
## 37       919616 HAH_RS04710 WP_014039884.1      HAH_0970
## 38       650378 HAH_RS03340 WP_014039640.1      HAH_0686
## 39       652307 HAH_RS03340 WP_014039640.1      HAH_0686
## 40      2011332 HAH_RS10070 WP_225306771.1      HAH_2078
## 41       134531 HAH_RS17935 WP_023842995.1      HAH_5130
## 42      2166758 HAH_RS10920 WP_044951957.1      HAH_2250
## 43       342746 HAH_RS18810 WP_014031290.1      HAH_5309
## 44        83862        <NA>           <NA>          <NA>
## 45        37515 HAH_RS00210 WP_014039074.1      HAH_0043
## 46      1616222 HAH_RS19110 WP_079891615.1          <NA>
## 47       145359 HAH_RS15890 WP_044952751.1      HAH_4167
## 48      2334523 HAH_RS11710 WP_014041106.1      HAH_2415
## 49      2382026 HAH_RS11965 WP_014041153.1      HAH_2465
## 50        72148 HAH_RS17590 WP_014031052.1      HAH_5058
## 51       679242 HAH_RS03485 WP_233425848.1      HAH_0716
## 52      2237673 HAH_RS11270 WP_014041022.1      HAH_2323
## 53       129484 HAH_RS15825 WP_014030704.1      HAH_4153
## 54       959974 HAH_RS04915 WP_004962786.1      HAH_1011
## 55       959974 HAH_RS04920 WP_014039920.1      HAH_1012
## 56      2732079 HAH_RS13815 WP_005537153.1      HAH_2841
## 57       168953 HAH_RS19805 WP_004516449.1      HAH_0188
## 58       389807        <NA>           <NA>          <NA>
## 59       369960 HAH_RS16870 WP_014030909.1      HAH_4373
## 60       371126 HAH_RS16875 WP_014030910.1      HAH_4374
## 61       220014 HAH_RS16210 WP_014030783.1      HAH_4237
## 62       566231 HAH_RS02910 WP_014039563.1      HAH_0600
## 63       387176 HAH_RS01990 WP_008311121.1      HAH_0414
## 64      2871939        <NA>           <NA>          <NA>
## 65      1108872 HAH_RS05710 WP_007187990.1      HAH_1175
## 66      2431968 HAH_RS12160 WP_023843403.1          <NA>
## 67      2701011 HAH_RS13645 WP_014041447.1      HAH_2804
## 68       613900 HAH_RS03145 WP_014039605.1      HAH_0647
## 69      2338569 HAH_RS11735 WP_014041111.1      HAH_2420
## 70      1770614 HAH_RS08865 WP_014040621.1      HAH_1825
## 71      1770614 HAH_RS08870 WP_008308731.1      HAH_1826
## 72       137241 HAH_RS00715 WP_014039167.1      HAH_0149
## 73      1479045 HAH_RS07530 WP_023843263.1      HAH_1549
## 74        77908 HAH_RS15490 WP_233425876.1      HAH_4086
## 75       444555 HAH_RS02290 WP_008311015.1      HAH_0475
## 76       857845 HAH_RS04440 WP_014039831.1      HAH_0913
## 77      1324370 HAH_RS06685 WP_014040230.1      HAH_1373
## 78      1570336 HAH_RS07930 WP_014040453.1      HAH_1633
## 79      1131628 HAH_RS05810 WP_014040076.1      HAH_1196
## 80       330996 HAH_RS18765 WP_014031281.1      HAH_5300
## 81      1196915 HAH_RS06110 WP_014040129.1      HAH_1255
## 82        71334 HAH_RS20115 WP_233425875.1          <NA>
## 83      2366903 HAH_RS11890 WP_008308192.1      HAH_2451
## 84      2366669        <NA>           <NA>          <NA>
## 85      1488642 HAH_RS07575 WP_014040386.1      HAH_1558
## 86       513579 HAH_RS02630 WP_014039518.1      HAH_0543
## 87       513579 HAH_RS02635 WP_014039519.1      HAH_0544
## 88       470145 HAH_RS17215 WP_014030979.1      HAH_4445
## 89       143240 HAH_RS17965 WP_014031127.1      HAH_5136
## 90      2543712 HAH_RS12770 WP_014041298.1      HAH_2627
## 91       288318 HAH_RS20200 WP_233425857.1      HAH_4299
## 92       362027 HAH_RS18900 WP_014031308.1      HAH_5327
## 93      1726872 HAH_RS08640 WP_014040583.1      HAH_1780
## 94      1727065 HAH_RS08645 WP_004515500.1      HAH_1781
## 95       187730 HAH_RS00985 WP_014039220.1      HAH_0207
## 96      2529226 HAH_RS12685 WP_014041281.1      HAH_2610
## 97      2529226 HAH_RS12690 WP_014041282.1      HAH_2611
## 98       347149 HAH_RS16750 WP_014030886.1      HAH_4348
## 99      2921176 HAH_RS14795 WP_023843464.1      HAH_3040
## 100     1282079 HAH_RS06485 WP_014040199.1      HAH_1331
## 101     1280990        <NA>           <NA>          <NA>
## 102       95249 HAH_RS17730 WP_014031079.1      HAH_5086
## 103      122341 HAH_RS17865 WP_014031108.1      HAH_5116
## 104     1489911 HAH_RS07585 WP_014040388.1      HAH_1560
## 105       41224 HAH_RS15315 WP_014030605.1      HAH_4046
## 106     1030498 HAH_RS05300 WP_044952190.1      HAH_1090
## 107      139107        <NA>           <NA>          <NA>
## 108      501743 HAH_RS02575 WP_014039507.1      HAH_0532
## 109     2187885 HAH_RS11020 WP_008312608.1      HAH_2271
## 110     1999683 HAH_RS10020 WP_014040815.1      HAH_2068
## 111     1877401 HAH_RS09385 WP_014040712.1      HAH_1939
## 112      321357        <NA>           <NA>          <NA>
## 113      731669 HAH_RS03765 WP_014039714.1      HAH_0772
## 114      178314        <NA>           <NA>          <NA>
## 115     2120337 HAH_RS10630 WP_014040909.1      HAH_2192
## 116     1020285 HAH_RS05250 WP_004962592.1      HAH_1079
## 117     2492506 HAH_RS12490 WP_004958397.1      HAH_2571
## 118     1020285 HAH_RS19890 WP_004962590.1      HAH_1080
## 119     1269860 HAH_RS06430 WP_014040187.1      HAH_1319
## 120      152461 HAH_RS18005 WP_044952797.1      HAH_5144
## 121     1537812 HAH_RS07790 WP_023843270.1      HAH_1601
## 122     1926380 HAH_RS09650 WP_014040753.1      HAH_1992
## 123      716673 HAH_RS03690 WP_014039700.1      HAH_0758
## 124     1487591 HAH_RS07570 WP_014040385.1      HAH_1557
## 125       88101 HAH_RS17700 WP_014031073.1      HAH_5080
## 126      281756        <NA>           <NA>          <NA>
## 127      301369 HAH_RS16575 WP_014030851.1      HAH_4313
## 128     2326344 HAH_RS11680 WP_014041100.1      HAH_2409
## 129     2751174 HAH_RS13895 WP_014041491.1      HAH_2857
## 130     1605430 HAH_RS19095 WP_079891614.1          <NA>
## 131     1510864 HAH_RS07690 WP_023843267.1      HAH_1581
## 132      898085 HAH_RS04580 WP_014039861.1      HAH_0943
## 133      243046 HAH_RS16320           <NA>      HAH_4261
## 134     1331841 HAH_RS06740 WP_008309382.1      HAH_1384
## 135      173904 HAH_RS18085 WP_225308163.1      HAH_5160
## 136     2791622 HAH_RS14100 WP_014041523.1      HAH_2900
## 137      358862 HAH_RS18900 WP_014031308.1      HAH_5327
## 138     1188572 HAH_RS06070 WP_014040121.1      HAH_1247
## 139      382298 HAH_RS16910 WP_044952599.1      HAH_4381
## 140     2714453 HAH_RS19170 WP_014041461.1      HAH_2820
## 141      222695 HAH_RS18270 WP_014031185.1      HAH_5197
## 142      201004        <NA>           <NA>          <NA>
## 143      314762 HAH_RS18685 WP_044952840.1      HAH_5284
## 144      432283 HAH_RS17075 WP_014030950.1      HAH_4415
## 145      417008 HAH_RS02130 WP_008311072.1      HAH_0443
## 146      175982 HAH_RS18095 WP_233425884.1      HAH_5162
## 147     1868620 HAH_RS09325 WP_014040703.1      HAH_1926
## 148      311900        <NA>           <NA>          <NA>
## 149      420799 HAH_RS02155 WP_014039430.1      HAH_0448
## 150     2532199 HAH_RS12700 WP_014041284.1      HAH_2613
## 151     2681207 HAH_RS13550 WP_014041429.1      HAH_2784
## 152      195437 HAH_RS18175 WP_014031166.1      HAH_5177
## 153      195437 HAH_RS18180 WP_023843005.1      HAH_5178
## 154     2460980 HAH_RS12310 WP_044952022.1      HAH_2535
## 155     2455369 HAH_RS12290 WP_014041213.1      HAH_2531
## 156     2278206 HAH_RS11460 WP_014041055.1      HAH_2362
## 157     2278360 HAH_RS11465 WP_014041057.1      HAH_2364
## 158     1795069 HAH_RS08980 WP_014040641.1      HAH_1851
## 159      349632 HAH_RS01770 WP_014039363.1      HAH_0368
## 160     2128263 HAH_RS10680 WP_230458846.1      HAH_2202
## 161      186795        <NA>           <NA>          <NA>
## 162      813424 HAH_RS04205 WP_014039788.1      HAH_0864
## 163      420799 HAH_RS02155 WP_014039430.1      HAH_0448
## 164     1126292 HAH_RS05775 WP_014040071.1      HAH_1188
## 165     1388104 HAH_RS07040 WP_014040295.1      HAH_1449
## 166     1057674        <NA>           <NA>          <NA>
## 167     2053813 HAH_RS19510 WP_144443759.1          <NA>
## 168     1926380 HAH_RS09650 WP_014040753.1      HAH_1992
## 169     2235752 HAH_RS11255 WP_014041019.1      HAH_2320
## 170       37607 HAH_RS15300 WP_014030602.1      HAH_4043
## 171     2182750 HAH_RS11000 WP_014040970.1      HAH_2267
## 172      168828 HAH_RS15980 WP_014030736.1      HAH_4186
## 173     2930951 HAH_RS14820 WP_079891627.1      HAH_3045
## 174      142110 HAH_RS15880 WP_014030716.1      HAH_4165
## 175     2346280 HAH_RS11760 WP_014041116.1      HAH_2425
## 176     2047807 HAH_RS19135 WP_014040856.1      HAH_2118
## 177      875385 HAH_RS04500 WP_014039844.1      HAH_0926
## 178      654397 HAH_RS03350 WP_014039642.1      HAH_0688
## 179     2153523 HAH_RS10835 WP_014040942.1      HAH_2233
## 180     1093250 HAH_RS05635 WP_014040049.1      HAH_1160
## 181      487719 HAH_RS17300 WP_014030996.1      HAH_4462
## 182      112077        <NA>           <NA>          <NA>
## 183      486894        <NA>           <NA>          <NA>
## 184     2072632 HAH_RS10360 WP_008312922.1      HAH_2137
## 185       35223 HAH_RS15290 WP_014030600.1      HAH_4041
## 186     1846124 HAH_RS09220 WP_014040682.1      HAH_1901
## 187     2110310 HAH_RS10580 WP_014040901.1      HAH_2182
## 188      115356 HAH_RS15730 WP_014030687.1      HAH_4135
## 189       88101 HAH_RS17700 WP_014031073.1      HAH_5080
## 190      879173 HAH_RS04510 WP_023843172.1      HAH_0928
## 191     2355511 HAH_RS11815 WP_014041127.1      HAH_2436
## 192     2877174 HAH_RS14545 WP_014041601.1      HAH_2990
## 193     1662888 HAH_RS08325 WP_014040521.1      HAH_1715
## 194      330485 HAH_RS16695 WP_014030875.1      HAH_4337
## 195      292841 HAH_RS01470 WP_014039307.1      HAH_0307
## 196      274013        <NA>           <NA>          <NA>
## 197      352209 HAH_RS01780 WP_014039365.1      HAH_0370
## 198       19133 HAH_RS17380 WP_023842982.1      HAH_5015
## 199      177341 HAH_RS18095 WP_233425884.1      HAH_5162
## 200     1177885 HAH_RS06025 WP_014040113.1      HAH_1238
## 201     2013635 HAH_RS10080 WP_014040827.1      HAH_2080
## 202      424464 HAH_RS02175 WP_014039434.1      HAH_0452
## 203     1703633 HAH_RS08520 WP_014040560.1      HAH_1756
## 204     2450269 HAH_RS12265 WP_014041208.1      HAH_2526
## 205     1244146 HAH_RS06300 WP_014040165.1      HAH_1292
## 206     1915224 HAH_RS09615 WP_014040746.1      HAH_1985
## 207      319981 HAH_RS18715 WP_014031271.1      HAH_5290
## 208       76483        <NA>           <NA>          <NA>
## 209     1583720 HAH_RS07980 WP_014040463.1      HAH_1643
## 210     2492506 HAH_RS12490 WP_004958397.1      HAH_2571
## 211     2168825 HAH_RS10930 WP_014040958.1      HAH_2252
## 212      193515        <NA>           <NA>          <NA>
## 213     1303978 HAH_RS06585 WP_014040212.1      HAH_1353
## 214     1980664 HAH_RS09935 WP_023843337.1      HAH_2049
## 215      358571 HAH_RS18895 WP_023843019.1      HAH_5326
## 216     1715131 HAH_RS08580 WP_014040571.1      HAH_1768
## 217      733147 HAH_RS03775 WP_014039715.1      HAH_0774
## 218     2149221 HAH_RS10810 WP_008312682.1      HAH_2228
## 219     2797578 HAH_RS14130 WP_014041528.1      HAH_2906
## 220        3734 HAH_RS17325 WP_008311719.1      HAH_5004
## 221      308486        <NA>           <NA>          <NA>
## 222      357677 HAH_RS18890 WP_014031306.1      HAH_5325
## 223      396788 HAH_RS16970 WP_044952679.1      HAH_4393
## 224      169942        <NA>           <NA>          <NA>
## 225      224681 HAH_RS01160 WP_014039247.1      HAH_0242
## 226     1574238 HAH_RS07940 WP_014040455.1      HAH_1635
## 227     2427349 HAH_RS19770 WP_023843402.1      HAH_2500
## 228     1742960 HAH_RS08730 WP_014040599.1      HAH_1798
## 229      149982 HAH_RS15915 WP_233425870.1          <NA>
## 230     2036555 HAH_RS10215 WP_014040848.1      HAH_2109
## 231      325262 HAH_RS16670 WP_014030870.1      HAH_4332
## 232      587231 HAH_RS02985 WP_014039577.1      HAH_0615
## 233      270481 HAH_RS16440 WP_014030825.1      HAH_4286
## 234      296214 HAH_RS16555 WP_014030847.1      HAH_4309
## 235       36475 HAH_RS17455 WP_014031024.1      HAH_5029
## 236     1605430 HAH_RS19095 WP_079891614.1          <NA>
## 237     1923039 HAH_RS09645 WP_014040752.1      HAH_1991
## 238       72966 HAH_RS15465 WP_014030637.1      HAH_4082
## 239     1382594 HAH_RS07010 WP_014040290.1      HAH_1442
## 240     2184812 HAH_RS11005 WP_014040971.1      HAH_2268
## 241       77908 HAH_RS15490 WP_233425876.1      HAH_4086
## 242      101307 HAH_RS00530 WP_014039132.1      HAH_0111
## 243      408665 HAH_RS17000 WP_044952607.1      HAH_4400
## 244      100002 HAH_RS00525 WP_014039131.1      HAH_0110
## 245     2222756 HAH_RS11205 WP_014041009.1      HAH_2310
## 246       83862        <NA>           <NA>          <NA>
## 247     1042347 HAH_RS05360 WP_014039999.1      HAH_1102
## 248     1760287 HAH_RS08805 WP_014040613.1      HAH_1813
## 249     1785956 HAH_RS08945 WP_044951885.1      HAH_1842
## 250       59886 HAH_RS15405 WP_014030626.1      HAH_4069
## 251      286400 HAH_RS18565 WP_174878676.1      HAH_5259
## 252      463671 HAH_RS17195 WP_023842944.1      HAH_4440
## 253      471834        <NA>           <NA>          <NA>
## 254     1328354 HAH_RS06715 WP_008309387.1      HAH_1379
## 255      331524 HAH_RS01685 WP_014039346.1      HAH_0351
##                                                                annotation
## 1                                            FAD-dependent oxidoreductase
## 2                         2-oxoacid:acceptor oxidoreductase subunit alpha
## 3   3-hydroxyacyl-CoA dehydrogenase NAD-binding domain-containing protein
## 4                                         class 1 fructose-bisphosphatase
## 5                                            phosphoenolpyruvate synthase
## 6                                                     cysteine synthase A
## 7                                              thioredoxin family protein
## 8                                     sugar porter family MFS transporter
## 9                                                    hypothetical protein
## 10                                                    riboflavin synthase
## 11                                                   hypothetical protein
## 12                                                                   <NA>
## 13                                                    riboflavin synthase
## 14                                                   hypothetical protein
## 15  3-hydroxyacyl-CoA dehydrogenase NAD-binding domain-containing protein
## 16                                        class 1 fructose-bisphosphatase
## 17                                               universal stress protein
## 18                                                         aminopeptidase
## 19                       type II glyceraldehyde-3-phosphate dehydrogenase
## 20                                      long-chain fatty acid--CoA ligase
## 21                                               universal stress protein
## 22                                             MBL fold metallo-hydrolase
## 23                                           FAD-dependent oxidoreductase
## 24                                      aminoglycoside phosphotransferase
## 25                                                                   <NA>
## 26                           ABC transporter ATP-binding protein/permease
## 27                                               D-xylose 1-dehydrogenase
## 28                                                                   <NA>
## 29                                          Glu/Leu/Phe/Val dehydrogenase
## 30                            rubrerythrin-like domain-containing protein
## 31                                                    PrgI family protein
## 32                                                   hypothetical protein
## 33                                UvrD-helicase domain-containing protein
## 34                                                   hypothetical protein
## 35                       nickel pincer cofactor biosynthesis protein LarC
## 36                                                                   <NA>
## 37                                            creatininase family protein
## 38                                                 threonine--tRNA ligase
## 39                                                 threonine--tRNA ligase
## 40                        phosphate ABC transporter permease subunit PstC
## 41                                               universal stress protein
## 42                                        polyprenyl diphosphate synthase
## 43                                                   hypothetical protein
## 44                                                                   <NA>
## 45                                                 DEAD/DEAH box helicase
## 46                                                   hypothetical protein
## 47                          type I-B CRISPR-associated protein Cas8b/Csh1
## 48                                     phytoene desaturase family protein
## 49                                              SDR family oxidoreductase
## 50                             helix-turn-helix domain-containing protein
## 51                                             S8 family serine peptidase
## 52                                           phosphoenolpyruvate synthase
## 53                                                       HNH endonuclease
## 54                                          CBS domain-containing protein
## 55                                     Gfo/Idh/MocA family oxidoreductase
## 56                       translation elongation factor EF-1 subunit alpha
## 57                                                   hypothetical protein
## 58                                                                   <NA>
## 59                      acetyl-CoA carboxylase biotin carboxylase subunit
## 60                                        acyl-CoA/acyl-ACP dehydrogenase
## 61                                                    diadenylate cyclase
## 62                                               PAS domain S-box protein
## 63        electron transfer flavoprotein subunit beta/FixA family protein
## 64                                                                   <NA>
## 65                                                   hypothetical protein
## 66                                                   hypothetical protein
## 67                                                phosphoglycerate kinase
## 68                                     isocitrate dehydrogenase (NADP(+))
## 69                                                        MFS transporter
## 70                                           threonine/serine dehydratase
## 71                                             class IV adenylate cyclase
## 72                                                   alpha/beta hydrolase
## 73                              ABC transporter substrate-binding protein
## 74         type IV secretion system DNA-binding domain-containing protein
## 75                                         CheF family chemotaxis protein
## 76                                                   hypothetical protein
## 77                                    NAD(P)/FAD-dependent oxidoreductase
## 78                                   glycosyltransferase family 2 protein
## 79                                         site-2 protease family protein
## 80                                               ABC transporter permease
## 81                                                        pyruvate kinase
## 82                                   halocyanin domain-containing protein
## 83                                          HTH domain-containing protein
## 84                                                                   <NA>
## 85                                                    aldo/keto reductase
## 86                                 anthranilate phosphoribosyltransferase
## 87                                  AsnC family transcriptional regulator
## 88                       pyridoxamine 5'-phosphate oxidase family protein
## 89                                      DUF2309 domain-containing protein
## 90                               phosphate signaling complex protein PhoU
## 91                                                   hypothetical protein
## 92                                          GAF domain-containing protein
## 93                                                    RNA-binding protein
## 94                                                   hypothetical protein
## 95                                                   leucine--tRNA ligase
## 96                                                   hypothetical protein
## 97                                geranylgeranyl reductase family protein
## 98                                      penicillin acylase family protein
## 99                                              substrate-binding protein
## 100                    TRC40/GET3/ArsA family transport-energizing ATPase
## 101                                                                  <NA>
## 102                                               MaoC family dehydratase
## 103                                              universal stress protein
## 104                        orotate phosphoribosyltransferase-like protein
## 105                            helix-turn-helix domain-containing protein
## 106                               TRAP transporter fused permease subunit
## 107                                                                  <NA>
## 108                                                  hypothetical protein
## 109                     hydroxymethylglutaryl-CoA synthase family protein
## 110                                              PAS domain S-box protein
## 111                                  glycosyltransferase family 4 protein
## 112                                                                  <NA>
## 113                                          sodium-dependent transporter
## 114                                                                  <NA>
## 115                                                   amino acid permease
## 116                                                              aldolase
## 117                                                DUF6432 family protein
## 118                                                  hypothetical protein
## 119                                              ABC transporter permease
## 120                                 carbohydrate ABC transporter permease
## 121                                                  hypothetical protein
## 122                                                  surface glycoprotein
## 123                                   NAD(P)/FAD-dependent oxidoreductase
## 124                                 TrmB family transcriptional regulator
## 125                        urea ABC transporter substrate-binding protein
## 126                                                                  <NA>
## 127                                   2-oxo acid dehydrogenase subunit E2
## 128                                   methyl-accepting chemotaxis protein
## 129                                                  hypothetical protein
## 130                                  glycosyltransferase family 4 protein
## 131                                                    MoxR family ATPase
## 132                                                       MFS transporter
## 133                                             IS1595 family transposase
## 134                                                  hypothetical protein
## 135                                     DUF1616 domain-containing protein
## 136                             ABC transporter substrate-binding protein
## 137                                         GAF domain-containing protein
## 138                                   tyrosine-type recombinase/integrase
## 139                                    cation-translocating P-type ATPase
## 140                                pentapeptide repeat-containing protein
## 141                                          FAD-dependent oxidoreductase
## 142                                                                  <NA>
## 143                                   carotenoid oxygenase family protein
## 144                            4Fe-4S dicluster domain-containing protein
## 145                                            30S ribosomal protein S19e
## 146                                                             sulfatase
## 147                          fumarylacetoacetate hydrolase family protein
## 148                                                                  <NA>
## 149                                       SHOCT domain-containing protein
## 150                                                     sulfurtransferase
## 151                                      DUF429 domain-containing protein
## 152                                  PTS fructose transporter subunit IIC
## 153                                    chromosome segregation protein SMC
## 154                             ABC transporter substrate-binding protein
## 155                                  serine/threonine-protein kinase RIO2
## 156                                                  hypothetical protein
## 157                                                  hypothetical protein
## 158                                                  hypothetical protein
## 159                                  glycosyltransferase family 4 protein
## 160                            helix-turn-helix domain-containing protein
## 161                                                                  <NA>
## 162                         peptide-methionine (R)-S-oxide reductase MsrB
## 163                                       SHOCT domain-containing protein
## 164                          phosphopentomutase/phosphoglucosamine mutase
## 165                                   alkaline phosphatase family protein
## 166                                                                  <NA>
## 167                                                  hypothetical protein
## 168                                                  surface glycoprotein
## 169                                                   ATP-binding protein
## 170                                                   radical SAM protein
## 171                                                DUF2150 family protein
## 172                                                  hypothetical protein
## 173                                                     nitrate reductase
## 174                              type I-B CRISPR-associated protein Cas5b
## 175                                           ATP-dependent protease LonB
## 176                                                             sulfatase
## 177                                       NUDIX domain-containing protein
## 178                                                 diacylglycerol kinase
## 179                                                   GtrA family protein
## 180                                              ABC transporter permease
## 181                                                  hypothetical protein
## 182                                                                  <NA>
## 183                                                                  <NA>
## 184                                            50S ribosomal protein L32e
## 185                                                  hypothetical protein
## 186                                         beta-CASP ribonuclease aCPSF1
## 187                                    chromosome segregation protein SMC
## 188                                                  hypothetical protein
## 189                        urea ABC transporter substrate-binding protein
## 190                                     DUF2339 domain-containing protein
## 191                                         amidohydrolase family protein
## 192                                       GNAT family N-acetyltransferase
## 193                                                       acyltransferase
## 194                                        sugar ABC transporter permease
## 195                                              ABC transporter permease
## 196                                                                  <NA>
## 197                                      DUF354 domain-containing protein
## 198                                 sulfurtransferase TusA family protein
## 199                                                             sulfatase
## 200                                 IclR family transcriptional regulator
## 201                                       phosphate uptake regulator PhoU
## 202                                                  oligoendopeptidase F
## 203                           b(o/a)3-type cytochrome-c oxidase subunit 1
## 204                                 NADH-quinone oxidoreductase subunit D
## 205                                       NUDIX domain-containing protein
## 206                  amino acid ABC transporter substrate-binding protein
## 207                                                   bacteriohemerythrin
## 208                                                                  <NA>
## 209                                                  hypothetical protein
## 210                                                DUF6432 family protein
## 211                                       DUF92 domain-containing protein
## 212                                                                  <NA>
## 213                                                  hypothetical protein
## 214                                                  hypothetical protein
## 215                                                  hypothetical protein
## 216                                      DUF555 domain-containing protein
## 217                                             SDR family oxidoreductase
## 218                                        GHKL domain-containing protein
## 219                                             serine protein kinase RIO
## 220                            TetR/AcrR family transcriptional regulator
## 221                                                                  <NA>
## 222                                                     AAA family ATPase
## 223                                              MgtC/SapB family protein
## 224                                                                  <NA>
## 225                                        LamG domain-containing protein
## 226                                       DUF58 domain-containing protein
## 227                                                  hypothetical protein
## 228                                molybdopterin-dependent oxidoreductase
## 229                           nitrate ABC transporter ATP-binding protein
## 230                                                  hypothetical protein
## 231                                               glucose 1-dehydrogenase
## 232                                                             sialidase
## 233                                exodeoxyribonuclease VII large subunit
## 234                                                  hypothetical protein
## 235                               ribbon-helix-helix protein, CopG family
## 236                                  glycosyltransferase family 4 protein
## 237                                                  surface glycoprotein
## 238                                                    VOC family protein
## 239                                        amidophosphoribosyltransferase
## 240                                         PKD domain-containing protein
## 241        type IV secretion system DNA-binding domain-containing protein
## 242                                                 RNA methyltransferase
## 243             BREX-5 system adenine-specific DNA-methyltransferase PglX
## 244                     GTP-dependent dephospho-CoA kinase family protein
## 245                                                  hypothetical protein
## 246                                                                  <NA>
## 247                                                     DNA polymerase IV
## 248                                         VWA domain-containing protein
## 249                                                  hypothetical protein
## 250                                       DUF87 domain-containing protein
## 251                              cbb3-type cytochrome c oxidase subunit I
## 252                             sulfite exporter TauE/SafE family protein
## 253                                                                  <NA>
## 254                            cytochrome bc complex cytochrome b subunit
## 255                                                  heme-binding protein
```

```r
# check that all are accounted for:
nrow(final) == nrow(one) + nrow(two) + nrow(three)
```

```
## [1] TRUE
```

```r
write_csv(final, "04d_motif_annotation/motifs_annotated.csv")
```

how about palindromic motifs?

```r
pal2 <- unlist(palindromes)

GenomicRanges::findOverlaps(genes.only, pal2, ignore.strand = T, minoverlap = 5) -> genes
GenomicRanges::findOverlaps(trimmed.pro, pal2, ignore.strand = T, minoverlap = 5) -> promoters
GenomicRanges::findOverlaps(ig.only, pal2, ignore.strand = T, minoverlap = 5) -> ig.regions

# get IRanges from hits objects and add informative metadata
genelist <- pal2[subjectHits(genes)][, -c(2:5)]
genelist$type <- rep("gene", length(genes))
strand(genelist) <- strand(genes.only[queryHits(genes)])
genelist$feature_start <- start(genes.only[queryHits(genes)])
genelist$feature_end <- end(genes.only[queryHits(genes)])
genelist$GENEID <- as.character(genes.only$GENEID[queryHits(genes)])

prolist <- pal2[subjectHits(promoters)][, -c(2:5)]
prolist$type <- rep("promoter", length(promoters))
prolist$feature_start <- start(trimmed.pro[queryHits(promoters)])
prolist$feature_end <- end(trimmed.pro[queryHits(promoters)])
strand(prolist) <- strand(trimmed.pro[queryHits(promoters)])
prolist$GENEID <- as.character(trimmed.pro$GENEID[queryHits(promoters)])

iglist <- pal2[subjectHits(ig.regions)][, -c(2:5)]
iglist$type <- rep("intergenic", length(ig.regions))
iglist$feature_start <- start(ig.only[queryHits(ig.regions)])
iglist$feature_end <- end(ig.only[queryHits(ig.regions)])

# convert separate IRanges to Data frames
seqs <- seq(1, length(genes))
as.data.frame(prolist) -> one
rownames(one) <- NULL
as.data.frame(genelist, row.names(seqs)) -> two
rownames(two) <- NULL
as.data.frame(iglist, row.names(seqs)) -> three
rownames(three) <- NULL

# combine dfs (gene hits and promoter hits)
final <- rbind(one, two, three) %>% distinct(.keep_all = T)
colnames(final)[c(2, 3, 11)] <- c("motif_start", "motif_end", "locus_tag")

# merge with gff information (get NCBI annotations and locus names)
gff.df[gff.df$locus_tag %in% final$locus_tag, ] -> tmp
tmp[c(2, 3, 4, 10)] -> tmp2
(left_join(final, tmp2) -> final)
```

```
## Joining with `by = join_by(locus_tag)`
```

```
##       seqnames motif_start motif_end width strand            motif_id
## 1  NC_015948.1     1358582   1358600    19      - ATTHACTMGAAWMBKWGTA
## 2  NC_015948.1     1358585   1358603    19      - ATTHACTMGAAWMBKWGTA
## 3  NC_015948.1     1358582   1358600    19      + ATTHACTMGAAWMBKWGTA
## 4  NC_015948.1     1358585   1358603    19      + ATTHACTMGAAWMBKWGTA
## 5  NC_015948.1     2237586   2237604    19      + ATTHACTMGAAWMBKWGTA
## 6  NC_015948.1     2237589   2237607    19      + ATTHACTMGAAWMBKWGTA
## 7  NC_015948.1     2492483   2492501    19      - ATTHACTMGAAWMBKWGTA
## 8  NC_015948.1     2492486   2492504    19      - ATTHACTMGAAWMBKWGTA
## 9  NC_015948.1     2701208   2701226    19      + ATTHACTMGAAWMBKWGTA
## 10 NC_015948.1     2701211   2701229    19      + ATTHACTMGAAWMBKWGTA
## 11 NC_015944.1       87858     87876    19      + ATTHACTMGAAWMBKWGTA
## 12 NC_015944.1       87861     87879    19      + ATTHACTMGAAWMBKWGTA
## 13 NC_015944.1      134417    134435    19      - ATTHACTMGAAWMBKWGTA
## 14 NC_015944.1      134420    134438    19      - ATTHACTMGAAWMBKWGTA
## 15 NC_015948.1      168916    168934    19      - ATTHACTMGAAWMBKWGTA
## 16 NC_015948.1      168919    168937    19      - ATTHACTMGAAWMBKWGTA
## 17 NC_015948.1     2701208   2701226    19      - ATTHACTMGAAWMBKWGTA
## 18 NC_015948.1     2701211   2701229    19      - ATTHACTMGAAWMBKWGTA
## 19 NC_015948.1      420780    420798    19      + ATTHACTMGAAWMBKWGTA
## 20 NC_015948.1      420783    420801    19      + ATTHACTMGAAWMBKWGTA
##       matched_sequence     type feature_start feature_end   locus_tag
## 1  ATTCACTCGAAATGTAGTA promoter       1358460     1358709 HAH_RS06895
## 2  ATTTACTACATTTCGAGTG promoter       1358460     1358709 HAH_RS06895
## 3  ATTCACTCGAAATGTAGTA promoter       1358479     1358728 HAH_RS06900
## 4  ATTTACTACATTTCGAGTG promoter       1358479     1358728 HAH_RS06900
## 5  ATTCACTCGAAATCGAGTG promoter       2237553     2237673 HAH_RS11270
## 6  CTTCACTCGATTTCGAGTG promoter       2237553     2237673 HAH_RS11270
## 7  TATGATTCCGAACCGAGTG promoter       2492341     2492506 HAH_RS12490
## 8  AATCACTCGGTTCGGAATC promoter       2492341     2492506 HAH_RS12490
## 9  AATTACTCGGTTCCGAGTC promoter       2701155     2701282 HAH_RS13650
## 10 ATAGACTCGGAACCGAGTA promoter       2701155     2701282 HAH_RS13650
## 11 TTGCACGCCAATTCTAGTG promoter         87852       88101 HAH_RS17700
## 12 AATCACTAGAATTGGCGTG promoter         87852       88101 HAH_RS17700
## 13 GTTTACTATAATTATAGTA promoter        134374      134531 HAH_RS17935
## 14 AATTACTATAATTATAGTA promoter        134374      134531 HAH_RS17935
## 15 GCTGACACGGAATCGAGTA promoter        168880      168953 HAH_RS19805
## 16 TTTTACTCGATTCCGTGTC promoter        168880      168953 HAH_RS19805
## 17 AATTACTCGGTTCCGAGTC promoter       2701155     2701282 HAH_RS20000
## 18 ATAGACTCGGAACCGAGTA promoter       2701155     2701282 HAH_RS20000
## 19 TCAGCCTCGAAACCGAGTG     gene        420377      420799 HAH_RS02155
## 20 GTTCACTCGGTTTCGAGGC     gene        420377      420799 HAH_RS02155
##               acc old_locus_tag
## 1  WP_014040268.1      HAH_1418
## 2  WP_014040268.1      HAH_1418
## 3  WP_014040269.1      HAH_1419
## 4  WP_014040269.1      HAH_1419
## 5  WP_014041022.1      HAH_2323
## 6  WP_014041022.1      HAH_2323
## 7  WP_004958397.1      HAH_2571
## 8  WP_004958397.1      HAH_2571
## 9  WP_014041448.1      HAH_2806
## 10 WP_014041448.1      HAH_2806
## 11 WP_014031073.1      HAH_5080
## 12 WP_014031073.1      HAH_5080
## 13 WP_023842995.1      HAH_5130
## 14 WP_023842995.1      HAH_5130
## 15 WP_004516449.1      HAH_0188
## 16 WP_004516449.1      HAH_0188
## 17 WP_008307600.1      HAH_2805
## 18 WP_008307600.1      HAH_2805
## 19 WP_014039430.1      HAH_0448
## 20 WP_014039430.1      HAH_0448
##                                                               annotation
## 1  3-hydroxyacyl-CoA dehydrogenase NAD-binding domain-containing protein
## 2  3-hydroxyacyl-CoA dehydrogenase NAD-binding domain-containing protein
## 3                                        class 1 fructose-bisphosphatase
## 4                                        class 1 fructose-bisphosphatase
## 5                                           phosphoenolpyruvate synthase
## 6                                           phosphoenolpyruvate synthase
## 7                                                 DUF6432 family protein
## 8                                                 DUF6432 family protein
## 9                                                    riboflavin synthase
## 10                                                   riboflavin synthase
## 11                        urea ABC transporter substrate-binding protein
## 12                        urea ABC transporter substrate-binding protein
## 13                                              universal stress protein
## 14                                              universal stress protein
## 15                                                  hypothetical protein
## 16                                                  hypothetical protein
## 17                                                  hypothetical protein
## 18                                                  hypothetical protein
## 19                                       SHOCT domain-containing protein
## 20                                       SHOCT domain-containing protein
```
