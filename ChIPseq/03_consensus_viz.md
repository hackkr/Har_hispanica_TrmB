---
title: "HCA trackviewer"
author: Rylee K. Hackley
output: 
  html_document:
    keep_md: true
---




```r
# BiocManager::install(c("Gviz","trackViewer", "AnnotationDbi"))
library(tidyverse)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(AnnotationDbi)
library(trackViewer)
```

Import bam and bed files:

```r
Filelist <- read_csv("01a_mosaics_peaks/01a_hca_sample_key.csv", col_names = F)

for (i in 1:nrow(Filelist)) {
  assign(paste(Filelist$X3[i], "_bam", sep = ""), importBam(file.path("00_sorted_bams", Filelist$X1[i])))
}

filenames <- list.files(path = "01a_mosaics_peaks/", pattern = ".bed")
labels <- str_replace(filenames, ".bed", "_bed")

for (i in 1:length(filenames)) {
  label <- labels[i]
  scores <- importScore(file.path("01a_mosaics_peaks/", filenames[i]), format = "BED")
  assign(label, scores)
}
```

Import peaks shared across bioreps:

```r
# differently bound regions
hca_diff <- importScore("02_DiffBind/02_diffbind2022.bed", format = "BED")
hca_diff$dat <- coverageGR(hca_diff$dat)
hca_diff$dat$score <- rep(1, length(hca_diff$dat))
strand(hca_diff$dat) <- "+"

# no glu bound regions
hca_noglu <- importScore("02_DiffBind/noglucose_peaks.bed", format = "BED")
hca_noglu$dat <- coverageGR(hca_noglu$dat)
hca_noglu$dat$score <- rep(1, length(hca_noglu$dat))
strand(hca_noglu$dat) <- "+"

# glu bound regions
hca_glu <- importScore("02_DiffBind/glucose_peaks.bed", format = "BED")
hca_glu$dat <- coverageGR(hca_glu$dat)
hca_glu$dat$score <- rep(1, length(hca_glu$dat))
strand(hca_glu$dat) <- "+"

# shared regions (with an average score of > 400)
hca_shared <- importScore("02_DiffBind/trmB_shared_peaks.bed", format = "BED")
hca_shared$dat <- hca_shared$dat[score(hca_shared$dat) > 400]
hca_shared$dat <- coverageGR(hca_shared$dat)
hca_shared$dat$score <- rep(1, length(hca_shared$dat))
strand(hca_shared$dat) <- "+"
```

Build gene model

```r
gff <- GenomicFeatures::makeTxDbFromGFF("../00_genome_files/genomic.gff", format = "gff")
gff_df <- read_csv("../00_genome_files/GCF_000223905.1_gff.key.csv")

# import motif hits
motif <- read_tsv("04b_motif_discovery/1st_order_15peaks_19_whole_genome_FIMO.tsv") %>%
  head(., -3)
motif <- makeGRangesFromDataFrame(motif, seqnames.field = "sequence_name", keep.extra.columns = T)

# import DEGs
degs <- read_csv("../RNA-seq/01_deseq2_output/group_final.csv")[, c(1, 3:5)] %>%
  left_join(., gff_df[, c(1, 3:8)])
colnames(degs)[1] <- "GENEID"
degs$log2FoldChange.pyr <- round(degs$log2FoldChange.pyr, digits = 2)
colnames(degs)[2] <- "score"
DEGs <- makeGRangesFromDataFrame(degs, keep.extra.columns = T)

# create chromosome gRanges
gr <- GRanges("NC_015948.1", IRanges(1, 2995271))
gr2 <- GRanges("NC_015943.1", IRanges(1, 488918))
gr3 <- GRanges("NC_015944.1", IRanges(1, 405816))
```

create stranded bed files

```r
# bed files of genes with strand info
genes <- GenomicFeatures::genes(gff, columns = c("GENEID"))
genes_plus <- genes[strand(genes) == "+"]
genes_plus$score <- rep(1, length(genes_plus))
genes_minus <- genes[strand(genes) == "-"]
genes_minus$score <- rep(1, length(genes_minus))

motif$score <- rep(1, length(motif))

# export beds
rtracklayer::export.bed(genes_minus, "03_peak_visualize/genes_minus.bed")
rtracklayer::export.bed(genes_plus, "03_peak_visualize/genes_plus.bed")
rtracklayer::export.bed(motif, "03_peak_visualize/motif.bed")
rtracklayer::export.bed(DEGs, "03_peak_visualize/DEGs.bed")
```

create stranded tracks

```r
# import both files to make stranded gene and motif tracks
genes_bed <- importScore(file.path("03_peak_visualize/genes_plus.bed"),
  file.path("03_peak_visualize/genes_minus.bed"),
  format = "BED"
)
strand(genes_bed$dat) <- strand(genes_bed$dat2) <- "-"

motif_bed <- importScore(file.path("03_peak_visualize/motif.bed"),
  format = "BED"
)

deg_bed <- importScore(file.path("03_peak_visualize/DEGs.bed"),
  format = "BED"
)

# add diffbind and shared peaks
all_pk <- hca_glu
all_pk$dat <- c(hca_diff$dat, hca_glu$dat)
```

graph diff bind peaks

```r
dir.create(paste("03_peak_visualize/", "DiffBind", sep = ""))
```

```
## Warning in dir.create(paste("03_peak_visualize/", "DiffBind", sep = "")):
## '03_peak_visualize\DiffBind' already exists
```

```r
for (i in 1:length(hca_diff$dat)) {
  filename <- paste("03_peak_visualize/DiffBind/", "peak_", i, ".png", sep = "")
  filename2 <- paste("03_peak_visualize/DiffBind/", "peak_", i, ".pdf", sep = "")
  begin <- start(hca_diff$dat[i])
  term <- end(hca_diff$dat[i])
  title <- paste(seqnames(hca_diff$dat[i]), ": ", start(hca_diff$dat[i]), "-", end(hca_diff$dat[i]), sep = "")

  # set figure range around peak
  if ((begin - 1000) < 1) {
    term2 <- term + 1000
    gr.tmp <- GRanges(as.vector(seqnames(hca_diff$dat[i])), IRanges(begin, term2))
  } else {
    gr.tmp <- GRanges(as.vector(seqnames(hca_diff$dat[i])), IRanges(begin, term)) + 1000
  }

  # create list of desired tracks and name
  trackList <- trackList(genes_bed, deg_bed, motif_bed, hca_diff, HCA_trmB_0.1glu_B_bam, HCA_trmB_0glu_B_bam)
  names(trackList) <- c("genes", "LFC", "motif", "consensus peaks", "glucose", "no glucose")

  optSty <- optimizeStyle(trackList, theme = "safe")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style

  # set margins and individual track heights. Tracks are plotted from bottom up.
  setTrackViewerStyleParam(viewerStyle, "margin", c(.06, .09, .01, .09))
  setTrackStyleParam(trackList[[6]], "height", .3)
  setTrackStyleParam(trackList[[5]], "height", .3)
  setTrackStyleParam(trackList[[4]], "height", .1)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[2]], "height", .1)
  setTrackStyleParam(trackList[[1]], "height", .1)

  setTrackStyleParam(trackList[[1]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[2]], "ylim", c(-7, 7))
  setTrackStyleParam(trackList[[3]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[4]], "ylim", c(-1, 1))

  # set track colors
  setTrackStyleParam(trackList[[6]], "color", "black")
  setTrackStyleParam(trackList[[5]], "color", "skyblue3")
  setTrackStyleParam(trackList[[4]], "color", "grey")
  setTrackStyleParam(trackList[[3]], "color", c("royalblue4", "royalblue4"))
  setTrackStyleParam(trackList[[2]], "color", "darkmagenta")
  setTrackStyleParam(trackList[[1]], "color", c("azure3", "azure3"))

  # save plots with title and subtitle
  png(filename = filename, width = 800, height = 1000)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()

  pdf(file = filename2, width = 8, height = 11, title = title)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()
}
```

graph no glucose peaks (7, 8, 16, 17, 20, 22, 24, 26)

```r
dir.create(paste("03_peak_visualize/", "NoGlu", sep = ""))
```

```
## Warning in dir.create(paste("03_peak_visualize/", "NoGlu", sep = "")):
## '03_peak_visualize\NoGlu' already exists
```

```r
for (i in 1:length(hca_noglu$dat)) {
  filename <- paste("03_peak_visualize/NoGlu/", "peak_", i, ".png", sep = "")
  filename2 <- paste("03_peak_visualize/NoGlu/", "peak_", i, ".pdf", sep = "")
  begin <- start(hca_noglu$dat[i])
  term <- end(hca_noglu$dat[i])
  title <- paste(seqnames(hca_noglu$dat[i]), ": ", start(hca_noglu$dat[i]), "-", end(hca_noglu$dat[i]), sep = "")

  # set figure range around peak
  if ((begin - 1000) < 1) {
    term2 <- term + 1000
    gr.tmp <- GRanges(as.vector(seqnames(hca_noglu$dat[i])), IRanges(begin, term2))
  } else {
    gr.tmp <- GRanges(as.vector(seqnames(hca_noglu$dat[i])), IRanges(begin, term)) + 1000
  }

  # create list of desired tracks and name
  trackList <- trackList(genes_bed, deg_bed, motif_bed, hca_noglu, HCA_trmB_0.1glu_B_bam, HCA_trmB_0glu_B_bam)
  names(trackList) <- c("genes", "LFC", "motif", "consensus peaks", "glucose", "no glucose")

  optSty <- optimizeStyle(trackList, theme = "safe")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style

  # set margins and individual track heights. Tracks are plotted from bottom up.
  setTrackViewerStyleParam(viewerStyle, "margin", c(.06, .09, .01, .09))
  setTrackStyleParam(trackList[[6]], "height", .3)
  setTrackStyleParam(trackList[[5]], "height", .3)
  setTrackStyleParam(trackList[[4]], "height", .1)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[2]], "height", .1)
  setTrackStyleParam(trackList[[1]], "height", .1)

  setTrackStyleParam(trackList[[1]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[2]], "ylim", c(-7, 7))
  setTrackStyleParam(trackList[[3]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[4]], "ylim", c(-1, 1))

  # set track colors
  setTrackStyleParam(trackList[[6]], "color", "black")
  setTrackStyleParam(trackList[[5]], "color", "skyblue3")
  setTrackStyleParam(trackList[[4]], "color", "grey")
  setTrackStyleParam(trackList[[3]], "color", c("royalblue4", "royalblue4"))
  setTrackStyleParam(trackList[[2]], "color", "darkmagenta")
  setTrackStyleParam(trackList[[1]], "color", c("azure3", "azure3"))

  # save plots with title and subtitle
  png(filename = filename, width = 800, height = 1000)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()

  pdf(file = filename2, width = 8, height = 11, title = title)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()
}
```

graph glucose peaks

```r
dir.create(paste("03_peak_visualize/", "Glu", sep = ""))
```

```
## Warning in dir.create(paste("03_peak_visualize/", "Glu", sep = "")):
## '03_peak_visualize\Glu' already exists
```

```r
for (i in 1:length(hca_glu$dat)) {
  filename <- paste("03_peak_visualize/Glu/", "peak_", i, ".png", sep = "")
  filename2 <- paste("03_peak_visualize/Glu/", "peak_", i, ".pdf", sep = "")
  begin <- start(hca_glu$dat[i])
  term <- end(hca_glu$dat[i])
  title <- paste(seqnames(hca_glu$dat[i]), ": ", start(hca_glu$dat[i]), "-", end(hca_glu$dat[i]), sep = "")

  # set figure range around peak
  if ((begin - 1000) < 1) {
    term2 <- term + 1000
    gr.tmp <- GRanges(as.vector(seqnames(hca_glu$dat[i])), IRanges(begin, term2))
  } else {
    gr.tmp <- GRanges(as.vector(seqnames(hca_glu$dat[i])), IRanges(begin, term)) + 1000
  }

  # create list of desired tracks and name
  trackList <- trackList(genes_bed, deg_bed, motif_bed, hca_glu, HCA_trmB_0.1glu_B_bam, HCA_trmB_0glu_B_bam)
  names(trackList) <- c("genes", "LFC", "motif", "consensus peaks", "glucose", "no glucose")

  optSty <- optimizeStyle(trackList, theme = "safe")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style

  # set margins and individual track heights. Tracks are plotted from bottom up.
  setTrackViewerStyleParam(viewerStyle, "margin", c(.06, .09, .01, .09))
  setTrackStyleParam(trackList[[6]], "height", .3)
  setTrackStyleParam(trackList[[5]], "height", .3)
  setTrackStyleParam(trackList[[4]], "height", .1)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[2]], "height", .1)
  setTrackStyleParam(trackList[[1]], "height", .1)

  setTrackStyleParam(trackList[[1]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[2]], "ylim", c(-7, 7))
  setTrackStyleParam(trackList[[3]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[4]], "ylim", c(-1, 1))

  # set track colors
  setTrackStyleParam(trackList[[6]], "color", "black")
  setTrackStyleParam(trackList[[5]], "color", "skyblue3")
  setTrackStyleParam(trackList[[4]], "color", "grey")
  setTrackStyleParam(trackList[[3]], "color", c("royalblue4", "royalblue4"))
  setTrackStyleParam(trackList[[2]], "color", "darkmagenta")
  setTrackStyleParam(trackList[[1]], "color", c("azure3", "azure3"))

  # save plots with title and subtitle
  png(filename = filename, width = 800, height = 1000)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()

  pdf(file = filename2, width = 8, height = 11, title = title)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()
}
```

graph highly scoring shared peaks (some also in DB)

```r
dir.create(paste("03_peak_visualize/", "shared", sep = ""))
```

```
## Warning in dir.create(paste("03_peak_visualize/", "shared", sep = "")):
## '03_peak_visualize\shared' already exists
```

```r
for (i in 1:length(hca_shared$dat)) {
  filename <- paste("03_peak_visualize/shared/", "peak_", i, ".png", sep = "")
  filename2 <- paste("03_peak_visualize/shared/", "peak_", i, ".pdf", sep = "")
  begin <- start(hca_shared$dat[i])
  term <- end(hca_shared$dat[i])
  title <- paste(seqnames(hca_shared$dat[i]), ": ", start(hca_shared$dat[i]), "-", end(hca_shared$dat[i]), sep = "")

  # set figure range around peak
  if ((begin - 1000) < 1) {
    term2 <- term + 1000
    gr.tmp <- GRanges(as.vector(seqnames(hca_shared$dat[i])), IRanges(begin, term2))
  } else {
    gr.tmp <- GRanges(as.vector(seqnames(hca_shared$dat[i])), IRanges(begin, term)) + 1000
  }

  # create list of desired tracks and name
  trackList <- trackList(genes_bed, deg_bed, motif_bed, hca_shared, HCA_trmB_0.1glu_B_bam, HCA_trmB_0glu_B_bam)
  names(trackList) <- c("genes", "LFC", "motif", "consensus peaks", "glucose", "no glucose")

  optSty <- optimizeStyle(trackList, theme = "safe")
  trackList <- optSty$tracks
  viewerStyle <- optSty$style

  # set margins and individual track heights. Tracks are plotted from bottom up.
  setTrackViewerStyleParam(viewerStyle, "margin", c(.06, .09, .01, .09))
  setTrackStyleParam(trackList[[6]], "height", .3)
  setTrackStyleParam(trackList[[5]], "height", .3)
  setTrackStyleParam(trackList[[4]], "height", .1)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[2]], "height", .1)
  setTrackStyleParam(trackList[[1]], "height", .1)

  setTrackStyleParam(trackList[[1]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[2]], "ylim", c(-7, 7))
  setTrackStyleParam(trackList[[3]], "ylim", c(-1, 1))
  setTrackStyleParam(trackList[[4]], "ylim", c(-1, 1))

  # set track colors
  setTrackStyleParam(trackList[[6]], "color", "black")
  setTrackStyleParam(trackList[[5]], "color", "skyblue3")
  setTrackStyleParam(trackList[[4]], "color", "grey")
  setTrackStyleParam(trackList[[3]], "color", c("royalblue4", "royalblue4"))
  setTrackStyleParam(trackList[[2]], "color", "darkmagenta")
  setTrackStyleParam(trackList[[1]], "color", c("azure3", "azure3"))

  # save plots with title and subtitle
  png(filename = filename, width = 800, height = 1000)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()

  pdf(file = filename2, width = 8, height = 11, title = title)
  viewTracks(trackList, gr = gr.tmp, viewerStyle = viewerStyle)
  grid.text(label = title, x = .5, y = .99, just = "top", gp = gpar(cex = 1))
  dev.off()
}
```
