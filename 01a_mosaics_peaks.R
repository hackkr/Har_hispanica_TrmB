### Peak calling script using MOSAiCS R package
# written by Cynthia L. Darnell, updated by Rylee K. Hackley and Amy K. Schmid

# MOSAiCS
# Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles
# https://bioconductor.org/packages/release/bioc/html/mosaics.html

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
#BiocManager::install("mosaics")

library(mosaics)
library(hexbin)
library(tidyverse)

sample_file <- read_csv("01a_mosaics_peaks/01a_hca_sample_key.csv", col_names = F)

IP_files <- sample_file$X1
WCE_files <- sample_file$X2
fragL <- sample_file$X4 #fragment length optimized using ChIPQC

for (i in 1:nrow(sample_file)) {
  constructBins(
    infile = paste("00_sorted_bams/", IP_files[i], sep = ""), #PATH/TO/BAMFILES
    fileFormat = "bam",
    outfileLoc = "01a_mosaics_peaks/bins/",
    byChr = FALSE,
    fragLen = fragL[i],
    binSize = fragL[i],
    capping = 0,
    PET = FALSE
  )

  constructBins(
    infile = paste("00_sorted_bams/", WCE_files[i], sep = ""), #PATH/TO/BAMFILES
    fileFormat = "bam",
    outfileLoc = "01a_mosaics_peaks/bins/",
    byChr = FALSE,
    fragLen = fragL[i],
    binSize = fragL[i],
    capping = 0,
    PET = FALSE
  )
}

for (i in 1:nrow(sample_file)) {
  replacements <- paste(".bam_fragL", fragL[i], "_bin", fragL[i], ".txt", sep = "")
  sample_name <- paste("01a_mosaics_peaks/bins/", sample_file[i, 1], sep = "")
  sample_name <- str_replace(string = sample_name, pattern = ".bam", replacement = replacements)
  ref_name <- paste("01a_mosaics_peaks/bins/", sample_file[i, 2], sep = "")
  ref_name <- str_replace(string = ref_name, pattern = ".bam", replacement = replacements)

  print(paste("analyzing", sample_name, "against", ref_name))

  binTest <- readBins(type = c("chip", "input"), fileName = c(sample_name, ref_name))
  fitTest <- mosaicsFit(binTest, analysisType = "IO", bgEst = "rMOM")

  peakTest <- mosaicsPeak(fitTest, signalModel = "2S", FDR = 0.01, maxgap = 220, minsize = fragL[i])
  peakTest <- extractReads(peakTest,
    chipFile = paste("00_sorted_bams/", sample_file[i, 1], sep = ""),
    chipFileFormat = "bam", chipFragLen = fragL[i],
    controlFile = paste("00_sorted_bams/", sample_file[i, 2], sep = ""),
    controlFileFormat = "bam", controlFragLen = 220
  )
  peakTest <- findSummit(peakTest)
  export(peakTest, type = "bed", filename = paste("01a_mosaics_peaks/", sample_file$X3[i], ".bed", sep = ""))
}
