## This script looks at enrichment of signals within hit enhancers

################################################################################################################################ #
## Setup ----


## Generic
# rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(readxl)
  library(tidyr)
  library(liftOver)
  library(rtracklayer)
  library(rcartocolor)

  source("../../Scripts/Functions.R")

## Load
  load("../../Data/Preprocessed/NHA Pooled.rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]

## Data information
  pos <- guides$TargetID[which(guides$TargetCat == "Promoter")]
  pos <- unique(pos)
  enh <- guides$TargetID[which(guides$TargetCat == "Enh")]
  enh <- unique(enh)
  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  sig.colours <- c("black", "firebrick1")

## Results
  # s <- read.csv("../2_DE/Enhancers - All Results Summary Filtered.csv")
  # load("../2_DE/Enhancers - All Results Summary.rda")
  res.final <- read.csv("../2_DE/Enh/Results Final.csv")
  
## Set up enhancer lists
  enh.hits <- unique(res.final$Enh.Pos[which(res.final$HitPermissive)])

  
################################################################################################################################ #
## Overlaps for validation ----


## Process input files
  ## Hnisz
    # performed by Irina in another script...
    hnisz <- read.table("/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/PROCESSED/liftOver_hg38/Astrocytes.hg38.bed")
    hnisz.annot <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/1-s2.0-S0092867413012270-mmc7/Astrocytes.csv")
    
    # Columns are: enhancer ID, chromosome, start, end, associated gene, enhancer rank, is enhancer a super-enhancer (1:yes, 0:no), H3K27ac ChIP-seq density (rpm/bp), read density in corresponding input sample (rpm/bp)
    colnames(hnisz.annot) <- c("ID", "chr", "start", "end", "gene", "rank", "super", "H3K27ac_density", "DensityInCorrespondingInputSample")
    
    hnisz <- hnisz[which(hnisz$V4 %in% hnisz.annot$ID[which(hnisz.annot$super == 1)]),] # filter to annotated superenhancers
    
    # save
    hnisz.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/PROCESSED/liftOver_hg38/Astrocytes.hg38.super.bed"
    write.bed(hnisz, dir = hnisz.in)
  
  ## ENCODE
    encode.file <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/ENCODE/ENCODE_CREs_050422_ucscDownload"
  
    # read in 
    encode <- read.delim(encode.file)
  
    # partition by ucscLabel
    table(encode$ucscLabel)
    for(j in unique(encode$ucscLabel)) {
      x <- encode[which(encode$ucscLabel == j),]
      write.bed(x[,1:4], dir = paste0(encode.file, ".", j, ".bed"))
    }
    
  ## Nott
    nott.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Nott2019/"
    nott.files <- list.files(nott.dir)[-c(1:2)]
    nott.files <- nott.files[-length(nott.files)]
    nott.files <- nott.files[-grep("sorted", nott.files)]
    # nott.files <- paste0(nott.dir, nott.files)
    
    # convert to hg38
    for (j in nott.files) {
      lift.bed(bed.in.19 = paste0(nott.dir, j), 
               bed.out.38 = paste0(nott.dir, "hg38/", j),
              isList = FALSE)
    }

## Run
  # directories
  nott.in <- paste0(nott.dir, "hg38/", nott.files)  
  encode.in <- paste0(encode.file, ".", unique(encode$ucscLabel), ".bed")
  hnisz.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/PROCESSED/liftOver_hg38/Astrocytes.hg38.super.bed"
  nha.all <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/BulkATAC_RNAseq/NHA_ATAC_S3.filtered.BAM_peaks_reformatedGJS.bed.txt"
  validation.out <- "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/BED_Annotations_ValidationSets.bed"
  validation.in <- paste0(c(nott.in, encode.in, hnisz.in), collapse = " ")
  
  # run  
  call <- paste("annotateBed",
                "-i", nha.all,
                "-files", validation.in,
                ">", validation.out)

  system(call, intern = FALSE, wait = TRUE) 
  
## Read in and wrangle
  # read
  vld <- read.table(validation.out, sep = "\t")
  vld <- vld[,-c(4:9)]
  
  # rename columns
  annot.cols <- 5:22
  colnames(vld)[-annot.cols] <- c("chr", "start", "end", "peak")
  
  x <- gsub("/mnt/Data0/PROJECTS/CROPSeq/PublicData/", "", validation.in)
  x <- gsub(".bed", "", x)
  x <- gsub("Nott2019/hg38/", "Nott_", x)
  x <- gsub("/PROCESSED/liftOver_hg38/Astrocytes.hg38.super", "", x)
  x <- gsub("/ENCODE_CREs_050422_ucscDownload.", "_", x)
  x <- gsub("optimal_peak.", "", x)
  x <- gsub("IDR_ENCODE.", "", x)
  x <- gsub("Cell2013_", "", x)
  
  x <- strsplit(x, " ")[[1]]
  
  colnames(vld)[annot.cols] <- x
  
  # binarise
  vld[annot.cols] <- vld[annot.cols] > 0
  
  # add annotation for candidate and hit peaks
  vld$peak <- paste0(vld$chr, ":", vld$start, "-", vld$end)
  vld$Candidate <- vld$peak %in% guides$TargetCoord[which(guides$TargetCat == "Enh")]
  vld$Hit <- vld$peak %in% enh.hits
  
  # save
  write.csv(vld, file = "Validation Sets - Annotations.csv")
  
## Analyse
  x <- apply(vld[,annot.cols], 2, table) %>% t() %>% as.data.frame
  colnames(x) <- c("NoOverlap", "Overlap")

  x$OverlapInCandidates <- apply(vld[which(vld$Candidate),annot.cols], 2, function(y) length(which(y)))
  x$OverlapInHits <- apply(vld[which(vld$Hit),annot.cols], 2, function(y) length(which(y)))
  
  x$FisherCanditates <- apply(vld[,annot.cols], 2, function(y) fisher.test(table(y, vld$Candidate))$estimate)
  x$FisherCanditatesP <- apply(vld[,annot.cols], 2, function(y) fisher.test(table(y, vld$Candidate))$p.value)
  
  x$FisherHits <- apply(vld[which(vld$Candidate),annot.cols], 2, function(y) {
    if (!(any(y))) return(NA) 
    fisher.test(table(y, vld$Hit[which(vld$Candidate)]))$estimate 
  })
  x$FisherHitsP <- apply(vld[which(vld$Candidate),annot.cols], 2, function(y) { 
    if (!(any(y))) return(NA) # necessary as the fisher test fails for ctcf, due to their being nothing true
    fisher.test(table(y, vld$Hit[which(vld$Candidate)]))$p.value 
  })
    
  write.csv(x, file = "Validation Sets - Statistics.csv")
  
## Plot 
  ## Fractions
    p <- data.frame(Dataset = rownames(x),
                    All.NHA.OCR = x$Overlap / nrow(vld),
                    Candidate.OCR = x$OverlapInCandidates / length(which(vld$Candidate)))
    p <- p[c(1,2,14,15,18),]
    
    p <- melt(p)
    
    pdf(file = "Validation Sets - Fractions.pdf", height = 3, width = maxw)
    ggplot(p, aes(x = Dataset, y = value, fill = variable)) +
      geom_col(position = "dodge", colour = "black") +
      theme_bw() +
      scale_fill_carto_d(palette = "Geyser") +
      theme(axis.title.x = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            panel.grid.major.x = invis, panel.border = invis, axis.line.y = element_line()) +
      scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
      labs(y = "Fraction of Overlapping Peaks")
    dev.off()
    
    # filtered version
     p <- data.frame(Dataset = rownames(x),
                    All.NHA.OCR = x$Overlap / nrow(vld),
                    Candidate.OCR = x$OverlapInCandidates / length(which(vld$Candidate)))
    p <- p[c(1,2,14,15),]
    
    p <- melt(p)
    
    p$Dataset <- factor(p$Dataset, levels = levels(as.factor(p$Dataset))[c(2,1,3,4)])
    
    
    pdf(file = "Validation Sets - Fractions (Filtered).pdf", height = 3, width = 3)
    ggplot(p, aes(x = Dataset, y = value, fill = variable, colour = variable)) +
      geom_col(position = "dodge", width = 0.7) +
      theme_bw() +
      scale_fill_carto_d(palette = "Geyser") +
      scale_colour_carto_d(palette = "Geyser") +
      coord_flip() +
      theme(axis.title.y = invis, 
            panel.grid.major.y = invis, panel.border = invis, axis.line.x = element_line(), legend.position = "none") +
      scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
      labs(y = "Fraction of NHA Peaks Overlapping With Resource")
    dev.off()
    
    # filtered version for super enhancers
     p <- data.frame(Dataset = rownames(x),
                    All.NHA.OCR = x$Overlap / nrow(vld),
                    Candidate.OCR = x$OverlapInCandidates / length(which(vld$Candidate)),
                    Hit.OCR = x$OverlapInHits / length(which(vld$Hit)))
    p <- p[c(18),]
    
    p <- melt(p)
    
    p$Dataset <- factor(p$Dataset, levels = levels(as.factor(p$Dataset))[c(2,1,3,4)])
    
    cols <- carto_pal(3, "Geyser")[c(1,3,2)]
    
    pdf(file = "Validation Sets - Fractions (Filtered, Super).pdf", height = 1.5, width = 3)
    ggplot(p, aes(x = variable, y = value, fill = variable, colour = variable)) +
      geom_col(position = "dodge", width = 0.7) +
      theme_bw() +
      scale_fill_manual(values = cols) +
      scale_colour_manual(values = cols) +
      coord_flip() +
      theme(axis.title.y = invis, 
            panel.grid.major.y = invis, panel.border = invis, axis.line.x = element_line(), legend.position = "none") +
      scale_y_continuous(expand = c(0,0), limits = c(0,0.3)) +
      labs(y = "Fraction of NHA Peaks Overlapping NHA Super-enhancers")
    dev.off()
  
  ## Odds ratios
    p <- x
    p$Dataset <- rownames(p)
    pdf(file = "Validation Sets - Odds Ratios.pdf", height = 3, width = maxw)
    ggplot(p, aes(x = Dataset, y = FisherCanditates)) +
      geom_col(colour = "black", fill = "darkorange1") +
      theme_bw() +
      theme(axis.title.x = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            panel.grid.major.x = invis, panel.border = invis, axis.line.y = element_line()) +
      scale_y_continuous(expand = c(0,0), limits = c(0,45)) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(y = "Odds Ratio")
    
    ggplot(p, aes(x = Dataset, y = log(FisherCanditates))) +
      geom_col(colour = "black", fill = "darkorange1") +
      theme_bw() +
      theme(axis.title.x = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            panel.grid.major.x = invis, panel.border = invis, axis.line.y = element_line()) +
      # scale_y_continuous(expand = c(0,0), limits = c(0,45)) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(y = "Log Odds Ratio\n(CTCF is NA, Not 0)")
    dev.off()
    
## Focused plot: ENCODE and Nott
   p <- x
    p$Dataset <- rownames(p)
    p <- p[grep("Nott", p$Dataset),]
    # p <- p[-grep("CTCF", p$Dataset),]
    pdf(file = "Validation Sets - Odds Ratios (Primary).pdf", height = 3, width = maxw)
    ggplot(p, aes(x = Dataset, y = (FisherCanditates))) +
      geom_col(colour = "black", fill = "darkorange1") +
      theme_bw() +
      theme(axis.title.x = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            panel.grid.major.x = invis, panel.border = invis, axis.line.y = element_line()) +
      # scale_y_continuous(expand = c(0,0), trans = "log2", limits = c(2^-6, 64), breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 16, 32, 64)) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(y = "Odds Ratio")
    dev.off()
    
# ## Focused plot: super enhancers
#        p <- x
#     p$Dataset <- rownames(p)
#     p <- p[grep("Hnisz", p$Dataset),]
#     # p <- p[-grep("CTCF", p$Dataset),]
#     pdf(file = "Validation Sets - Odds Ratios (Primary).pdf", height = 3, width = maxw)
#     ggplot(p, aes(x = Dataset, y = (FisherCanditates))) +
#       geom_col(colour = "black", fill = "darkorange1") +
#       theme_bw() +
#       theme(axis.title.x = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
#             panel.grid.major.x = invis, panel.border = invis, axis.line.y = element_line()) +
#       # scale_y_continuous(expand = c(0,0), trans = "log2", limits = c(2^-6, 64), breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 16, 32, 64)) +
#       geom_hline(yintercept = 0, linetype = 2) +
#       labs(y = "Odds Ratio")
#     dev.off()
    
    
################################################################################################################################ #
## Overlaps for demonstrating novelty ----
  
  

    

  
  
  

################################################################################################################################ #
## Human-gained enhancers from Vermunt et al. ----
  
## Directories
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  vermunt.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Vermun2016_HumanGainedH3K27ac/GainedBrainEnhancers.bed"
  vermunt.out <- paste0(getwd(), "/Vermunt_Annotation.bed")
  
  
## Process Supp Table 17
  x <- read_xlsx("../../../PublicData/Vermun2016_HumanGainedH3K27ac/ST17_enhancers.xlsx", sheet = 1) # human gained relative to macaque
  y <- read_xlsx("../../../PublicData/Vermun2016_HumanGainedH3K27ac/ST17_enhancers.xlsx", sheet = 2) # human gained relative to macaque and chimpanzee
  x$NotChimpanzee <- x$enhancer %in% y$enhancer 
  
  # write as BED
  write.table(x, file = vermunt.dir, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

## Run
  call <- paste("intersectBed",
                "-a", nha_hg38,
                "-b", vermunt.dir,
                "-loj", # Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B
                ">", vermunt.out)

  system(call, intern = FALSE, wait = TRUE) 
  
  
## Read in and wrangle
  vermunt <- read.table(vermunt.out, sep = "\t", header = FALSE)
  
  # column names
  colnames(vermunt) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.ID", "Vermunt.chr", "Vermunt.start", "Vermunt.end", "Vermunt.ID", "Vermunt.Chimp")
  
  # peak id
  vermunt$Peak.ID <- sub("_", ":", vermunt$Peak.ID)
  vermunt$Peak.ID <- sub("_", "-", vermunt$Peak.ID)
  vermunt$Hit <- vermunt$Peak.ID %in% enh.hits
  
## Save
  write.csv(vermunt, file = "Vermunt - Dataframe.csv")
    
################################################################################################################################ #
## Differentially-acetylated H3K27ac in ASD from  Sun et al. ----
  
  
## Directories
  sun.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Sun2016_"
  sun.out <- paste0(getwd(), "/Sun_Annotation.bed")
  nha_hg19 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19lift.bed"
  
## Process Supplementary Table 4 into beds
  sun.names <- excel_sheets("../../../PublicData/Sun2016_ST4.xlsx")
  sun.names <- sun.names[-1] # removes the notes tab
  y <- lapply(sun.names, function(y) read_xlsx("../../../PublicData/Sun2016_ST4.xlsx", sheet = y))
  sun.names <- gsub(" ", "_", sun.names)
  names(y) <- sun.names
  
  lapply(names(y), function(z) { 
    write.table(as.data.frame(y[[z]][,1:3]), file = paste0(sun.dir, z, "_hg19.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) 
  })
  
## Overlap
   call <- paste("annotateBed",
                "-i", nha_hg19,
                "-files",  paste0(sun.dir, sun.names, "_hg19.bed", collapse = " "),
                ">", sun.out)
    
   system(call, intern = FALSE, wait = TRUE) 
  
## Read in and wrangle
  sun <- read.table(sun.out, sep = "\t", header = FALSE)
  
  # column names
  annot.columns <- 5:10
  
  colnames(sun)[-annot.columns] <- c("chr.hg19", "start.hg19", "end.hg19", "ID.hg38")
  colnames(sun)[annot.columns] <- sun.names
  
  # peak id
  sun$ID.hg38 <- sub("_", ":", sun$ID.hg38)
  sun$ID.hg38 <- sub("_", "-", sun$ID.hg38)
  sun$Hit <- sun$ID.hg38 %in% enh.hits
  
## Binarise
  sun[,annot.columns] <- apply(sun[,annot.columns], 2, function(x) x > 0)
   
  table(sun$PFC_up, sun$Hit) %>% fisher.test()
  table(sun$PFC_down, sun$Hit) %>% fisher.test()
  table(sun$TC_up, sun$Hit) %>% fisher.test()
  table(sun$TC_down, sun$Hit)%>% fisher.test()
  table(sun$CB_up, sun$Hit)%>% fisher.test()
  table(sun$CB_down, sun$Hit)%>% fisher.test()
  
  
################################################################################################################################ #
## Differentially open chromatin in AD from  Bendl et al. ----
  
  
## Directories
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  bendl.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Bendl2022_AD_ATAC/Data_peaks_glia.bed"
  bendl.out <- paste0(getwd(), "/Bendl_Annotation.bed")
  
## Map bendl's peak names to our enh
  # use bedtools
  call <- paste("intersectBed",
                "-a", nha_hg38,
                "-b", bendl.dir,
                "-wb", 
                ">", bendl.out)

  system(call, intern = FALSE, wait = TRUE) 
  bendl.peaks <- read.delim(bendl.out, header = FALSE)
  
## Read in Bendl's AD-associated peaks
  # load differential openness
  bendl.hits <- read.delim("/mnt/Data0/PROJECTS/CROPSeq/PublicData/Bendl2022_AD_ATAC/Data_peaks_diff_glia.tsv")
  
  m <- match(rownames(bendl.hits), bendl.peaks$V8)
  bendl.hits <- data.frame(ScreenPeak = bendl.peaks$V4[m],
                           ScreenHit = ".",
                           ScreenGene = ".",
                           ScreenDistance = ".",
                           BendlPeak = rownames(bendl.hits),
                           bendl.hits[,1:4])
  bendl.hits <- bendl.hits[-which(is.na(bendl.hits$ScreenPeak)),]
  colnames(bendl.hits)[6:9] <- c("BraakScore", "Diagnosis", "ClinicalDementiaRating", "PlaqueLvl")
  
  # annotate with hit genes for our peaks
  bendl.hits$ScreenPeak <- sub("_", ":", bendl.hits$ScreenPeak)
  bendl.hits$ScreenPeak <- sub("_", "-", bendl.hits$ScreenPeak)
  bendl.hits$ScreenHit <- bendl.hits$ScreenPeak %in% enh.hits
  m <- match(bendl.hits$ScreenPeak, guides$TargetCoord)
  bendl.hits$ScreenPeak <- guides$TargetID[m]
  
  for (j in 1:nrow(bendl.hits)) {
    if (!(bendl.hits$ScreenHit[j])) {
      next
    } else {
      x <- res[which(res$Enh == bendl.hits$ScreenPeak[j]),]
      bendl.hits$ScreenGene[j] <- paste(x$Gene[which(x$Hit)], sep = " | ") # if there are many, but there are not...
      bendl.hits$ScreenDistance[j] <- paste(x$Gene.Distance[which(x$Hit)], sep = " | ") # if there are many, but there are not...
    }
  }
  
  write.csv(bendl.hits, file = "Bendl Results.csv")
  
  
# ################################################################################################################################ #
# ## Brain cell-type chromatin from Nott et al. ----
#   
#   
# ## Set paths
#   nha_hg19 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19lift.bed"
#   nott.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Nott2019/"
#   nott.files <- list.files(nott)[-1]
#   nott.out <- paste0(getwd(), "/Nott_Annotation.bed")
#   
# 
# ## Run
#   call <- paste("annotateBed",
#                 "-i", nha_hg19,
#                 "-files", paste0(nott.dir, nott.files, collapse = " "),
#                 ">", nott.out)
# 
#   system(call, intern = FALSE, wait = TRUE) 
#   
#   
# ## Read in and wrangle
#   nott <- read.table(nott.out, sep = "\t", header = FALSE)
#   
#   # column names
#   annot.columns <- 5:16
#   annots <- gsub("_optimal_peak", "", nott.files) 
#   annots <- gsub(".bed", "", annots)
#   annots <- gsub("_IDR_ENCODE", "", annots)
#   
#   colnames(nott)[-annot.columns] <- c("chr.hg19", "start.hg19", "end.hg19", "ID.hg38")
#   colnames(nott)[annot.columns] <- annots
#   
#   # peak id
#   nott$Peak <- sub("_", ":", nott$ID.hg38)
#   nott$Peak <- sub("_", "-", nott$Peak)
#   nott$Hit <- nott$Peak %in% enh.hits
#   
# ## Binarise
#   nott[,annot.columns] <- apply(nott[,annot.columns], 2, function(x) x > 0)
#   
# ## Categorise peaks
#   ## Function
#     cat.peak <- function(ct, dat = nott) {
#       x <- nott[,grep(ct, colnames(nott))]
#       colnames(x) <- gsub(paste0(ct, "\\."), "", colnames(x))
#       
#       x$Cat <- "N"
#       x$Cat[which(x$H3K4me3)] <- "P"
#       x$Cat[which((x$ATAC | x$H3K27) & !(x$H3K4me3))] <- "E"
#       
#       dat[,paste0(ct)] <- x$Cat
#       return(dat)
#     }
#     
#     nott <- cat.peak("LHX2")
#     nott <- cat.peak("NeuN")
#     nott <- cat.peak("Olig2")
#     nott <- cat.peak("PU1")
#   
# ## Visualise
#   p <- melt(nott[,c(18:22)], id.vars = "Hit")
#   # p$Celltype <- splitter(p$variable, "\\.", 1)
#   # p$Assay <- splitter(p$variable, "\\.", 2)
#   
#   p$value <- factor(p$value, levels = rev(c("E", "P", "N")))
#   levels(p$value) <- (c("Not in Nott et al.", "Nott et al. Promoter", "Nott et al. Enhancer"))
#   
#   pdf(file = "Nott - Barplot.pdf", height = 3.5, width = maxw)
#   ggplot(p, aes(x = Hit, fill = value)) +
#     geom_bar(position = position_fill(), colour = "black") +
#     facet_wrap(~variable, nrow = 1) +
#     theme_bw() +
#     scale_fill_manual(values = c("grey90", carto_pal(2, "Earth"))) +
#     scale_y_continuous(expand = c(0,0)) +
#     labs(y = "Fraction of Peaks", x = "Enhancer Has Significant Association")+
#     theme(legend.title = invis)
#   
#   ggplot(p, aes(x = Hit, fill = value)) +
#     geom_bar(colour = "black", position = "dodge") +
#     facet_wrap(~variable, nrow = 1) +
#     theme_bw() +
#     scale_fill_manual(values = c("grey90", carto_pal(2, "Earth"))) +
#     scale_y_continuous(expand = c(0,0)) +
#     labs(y = "Count of Peaks (Total=979)", x = "Enhancer Has Significant Association")+
#     theme(legend.title = invis)
#   dev.off()
#     
#   
# ################################################################################################################################ #
# ## ENCODE Enhancers ----
# 
# ## Paths
#   nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
#   encode.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/ENCODE_CREs_050422_ucscDownload"
#   encode.out <- paste0(getwd(), "/ENCODE_Annotation.bed")
#   
# 
# ## Run
#   call <- paste("intersectBed",
#                 "-a", nha_hg38,
#                 "-b", encode.dir,
#                 "-loj", # Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B
#                 ">", encode.out)
# 
#   system(call, intern = FALSE, wait = TRUE) 
#   
#   
# ## Read in and wrangle
#   encode <- read.table(encode.out, sep = "\t", header = FALSE)
#   
#   # Binarise
#   # encode$Overlap <- factor(encode$V5 > 0)
#   # levels(encode$Overlap) <- c("Not ")
#   
#   # rename levels of encode categories
#   encode$V17 <- factor(encode$V17, levels = c("enhD", "enhP", "K4m3", "prom", "."))
#   levels(encode$V17) <- c("Enhancer", "Enhancer", "Promoter or K4m3", "Promoter or K4m3", "No Overlap")
#   
#   # remove cases where a single NHA peak overlaps two of the same annotations
#   dup <- duplicated(encode[,c("V4", "V17")])
#   encode <- encode[which(!(dup)),]
#   
#   
#   # Add hit annotation
#   m <- match(encode$V4, nott$ID.hg38)
#   encode$Hit <- nott$Hit[m]
#   
# ## Table
#   # table(encode$Overlap, encode$Hit) # 14 of 979 are not overlapping...
# 
# ## Plot
#   # ggplot(encode, aes(x = Hit, fill = Overlap)) +
#   #   geom_bar(colour = "black") +
#   #   theme_bw() +
#   #   scale_fill_manual(values = c(carto_pal(2, "Earth"))) +
#   #   scale_y_continuous(expand = c(0,0), limits = c(0, 1000)) +
#   #   labs(y = "Count of Peaks", x = "Enhancer Has Significant Association") +
#   #   theme(panel.border = invis, axis.line.y = element_line())
#   
#   pdf(file = "ENCODE - Barplot.pdf", height = 3.5, width = maxw)
#   pA <- ggplot(encode, aes(x = Hit, fill = V17)) +
#     geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
#     theme_bw() +
#     scale_fill_manual(values = rev(c("grey90", carto_pal(2, "Earth")))) +
#     scale_y_continuous(expand = c(0,0), limits = c(0, 900)) +
#     labs(y = "Count of Peaks in Category", x = "Enhancer Has Significant Association") +
#     guides(fill = guide_legend(title = "ENCODE Annotation")) +
#     theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
#           panel.grid = invis)
#   
#   pB <- ggplot(encode, aes(x = Hit, fill = V17)) +
#     geom_bar(colour = "black", position = position_fill()) +
#     theme_bw() +
#     scale_fill_manual(values = rev(c("grey90", carto_pal(2, "Earth")))) +
#     scale_y_continuous(expand = c(0,0)) +
#     labs(y = "Fraction of Peaks in Category", x = "Enhancer Has Significant Association") +
#     theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", panel.grid = invis)
#   
#   plot_grid(pB, pA)
#   dev.off()

    
################################################################################################################################ #
## Transcription factor binding ----
  
## Binding sites from the JASPAR database (2022) were downloaded from UCSC, pre-subsetting to those overlapping our enhancers for filesize reasons
  
## Paths
  jaspar.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/JASPAR_TFBS_130422_ucscDownload"
  jaspar.out <- paste0(getwd(), "/Jaspar_Annotation.bed")
  

## Run
  call <- paste("intersectBed",
                "-a", nha_hg38,
                "-b", jaspar.dir,
                "-loj", # Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B
                ">", jaspar.out)

  system(call, intern = FALSE, wait = TRUE) 
  
## Read in
  jaspar <- read.table(jaspar.out, sep = "\t")
  colnames(jaspar) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "J.chr", "J.start", "J.end", "J.model", "J.score", "J.strand", "J.TF")

## Filter to expressed genes
  exp.thresh <- 2^-6
  means <- rowMeans(nha@assays$SCT@data)
  means <- names(means)[which(means > exp.thresh)]
  
  jaspar <- jaspar[which(jaspar$J.TF %in% means),] # 1175502 to 458188
  
## Cleaner dataframe
  renamed.hits <- sub(":", "_", enh.hits$Enh.Pos) %>% sub("-", "_", .)
  jaspar <- data.frame(Peak = jaspar$Peak.id,
                       Hit = jaspar$Peak.id %in% renamed.hits,
                       tf = jaspar$J.TF,
                       tf.start = jaspar$J.start,
                       tf.end = jaspar$J.end)

## Plot
  # simply the distribution of number of TFs per enhancer, for hit and non-hit enhancers
  p <- table(jaspar$Peak) %>% as.data.frame
  p$Hit <- p$Var1 %in% renamed.hits
    
  pdf(file = "Jaspar - Number of Peaks and TFs.pdf", height = 3, width = 3)
  ggplot(p, aes(x = ".", y = Freq, colour = Hit)) +
    geom_violin(draw_quantiles = 0.5, width = 0.8, scale = "width", fill = "white") +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), alpha = 0.1) +
    scale_colour_carto_d(palette = "Earth") +
    theme_bw() +
    labs(y = "Number of TF Binding Sites in Peak") +
    theme(panel.border = invis, axis.line.y = element_line(), axis.title.x = invis,
          axis.text.x = invis, axis.ticks.x = invis)
  
  # and the inverse: number of peaks per TFs
  p <- table(jaspar$tf) %>% as.data.frame
  
  ggplot(p, aes(x = ".", y = Freq)) +
    geom_violin(draw_quantiles = 0.5, width = 0.4, scale = "width", fill = "white") +
    geom_jitter(position = position_jitter(width = 0.1), alpha = 0.1) +
    scale_colour_carto_d(palette = "Earth") +
    theme_bw() +
    labs(y = "Number of Peak Binding Sites per TF") +
    theme(panel.border = invis, axis.line.y = element_line(), axis.title.x = invis,
          axis.text.x = invis, axis.ticks.x = invis)
  dev.off()
  
## Save
  write.csv(jaspar, file = "Jaspar - Overlaps.csv")
  
## Fisher test
  jasp <- dcast(jasp, Peak~tf) 
  
  # rename peaks
  jasp$Peak <- sub("_", ":", jasp$Peak)
  jasp$Peak <- sub("_", "-", jasp$Peak)
  
  rownames(jasp) <- jasp$Peak
  
  # hit peaks
  jasp$Hit <- jasp$Peak %in% res$Enh.Pos[which(res$Hit)]

## Overenrichment
  jasp[,-c(1,331)] <- jasp[,-c(1,331)] > 0
  jasp <- jasp[,-1]
  
  jasp.fish <- apply(jasp[,-330], 2, function(y) { 
    if (length(table(y)) == 1) return(NA) 
    f <- fisher.test(table(y, jasp$Hit))
    data.frame(OR = f$estimate, P = f$p.value)
  })
  
  jasp.fish <- do.call("rbind", jasp.fish)
  jasp.fish$Bonf <- p.adjust(jasp.fish$P, method = "bonferroni")  
  jasp.fish$FDR <- p.adjust(jasp.fish$P, method = "fdr")  
  
  write.csv(jasp.fish, file = "Jaspar - Fisher.csv")
  
## Plot
  p <- jasp.fish[which(jasp.fish$FDR < 0.1),]
  p$TF <- rownames(p)
  p$FDR <- -log10(p$FDR)
  rnk <- p$TF[order(p$FDR)]
  p$TF <- factor(p$TF, levels = rnk)
  p$Col <- p$FDR < -log10(0.05)
  
  pdf(file = "Jaspar - Fisher Barplot.pdf", height = 3, width = 5)
  ggplot(p, aes(x = TF, y = FDR, fill = Col)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values = c("darkorange1", "black")) +
    NoLegend() +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "-log10(FDR) (Fisher Test)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = invis,
          panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis)
  dev.off()  
  

################################################################################################################################ #
## Transcription factor binding V2, using ReMap ----
  
## Directories
  remap.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/TF_BindingSites/remap2022_nr_macs2_hg38_v1_0.bed"
  remap.out <- paste0(getwd(), "/ReMap.bed")
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  
## Intersect
  call <- paste("intersectBed",
                "-a", nha_hg38,
                "-b", remap.in,
                "-wb", 
                ">", remap.out)

  system(call, intern = FALSE, wait = TRUE) 
  
  
## Read in and wrangle
  remap <- read.table(remap.out, sep = "\t", header = FALSE)
  remap$V4 <- sub("_", ":", remap$V4) %>% sub("_", "-", .)
  m <- match(remap$V4, guides$TargetCoord)
  
  remap <- data.frame(Enh = guides$TargetID[m],
                      Enh.Coord = guides$TargetCoord[m],
                      TF = splitter(remap$V8, ":", 1),
                      TF.Expressed = ".",
                      TF.Coord = paste0(remap$V5, ":", remap$V6, "-", remap$V7),
                      Celltypes = splitter(remap$V8, ":", 2))
  
## Filter to those TFs expressed in NHAs
  remap$TF.Expressed <- remap$TF %in% mean
  
## Categorise peaks as: LevelA (in Ast); LevelB (2+ brain); LevelC (1 brain and 1+ other); LevelD (2+ other); LevelE (1 brain); and LevelF (1 non-brain)
  brain.ct <- c("neural-progenitor", "glioma", "tibial-nerve", "primary-glioblastoma", "glioblastoma", "NB-1643", "hippocampus",
                "cortical-interneuron", "neuron-progenitor", "neuroepithelilal-cells", "neuron", "nerve", "NPC", "dopaminergic-neuron",
                "neural", "brain-prefrontal-cortex", "SH-SY5Y", "SK-N-BE2-C", "HNPC", "SK-N-MC")
  
  split <- strsplit(remap$Celltypes, ",")
  nBrain <- sapply(split, function(x) length(which(x %in% brain.ct)))
  nNonBrain <- sapply(split, function(x) length(which(!(x %in% brain.ct))))
  
  remap$LvlA <- grepl("astrocyte", remap$Celltypes)
  remap$LvlB <- (nBrain >= 2) & !(remap$LvlA)
  remap$LvlC <- (nBrain == 1 & nNonBrain >= 1) & !(remap$LvlA)
  remap$LvlD <- (nBrain == 0 & nNonBrain >= 2) & !(remap$LvlA)
  remap$LvlE <- (nBrain == 1 & nNonBrain == 0) & !(remap$LvlA)
  remap$LvlF <- (nBrain == 0 & nNonBrain == 1) & !(remap$LvlA)
  
## Save
  write.csv(remap, file = "ReMap Processed.csv", row.names = FALSE)
  
  
## Per-enhancer statistics: total TFBS  
  remap$Category <- "."
  remap$Category[which(remap$LvlA)] <- "A: Ast"
  remap$Category[which(remap$LvlB)] <- "B: 2+ Brain"
  remap$Category[which(remap$LvlC)] <- "C: 1 Brain + 1+ NonBrain"
  remap$Category[which(remap$LvlD)] <- "D: 2+ NonBrain"
  remap$Category[which(remap$LvlE)] <- "E: 1 Brain"
  remap$Category[which(remap$LvlF)] <- "F: 1 NonBrain"
  
  tab <- table(remap$Enh, remap$Category) %>% data.frame()
  tab <- dcast(data = tab, formula = Var1~Var2)
  
  colnames(tab)[1] <- "Enh"
  tab$HitPermissive <- tab$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
  tab$HitCore <- tab$Enh %in% res.final$Enh[which(res.final$HitCore)]
  
  write.csv(tab, file = "ReMap nOverlaps Per Enhancer.csv", row.names = FALSE)
  
## Per-TF statistics: enrichment within hits
  
  
  
## What is the overlap between ReMap and JASPAR?
  # read in
  jaspar <- read.csv("Jaspar - Overlaps.csv", row.names = 1)
  remap <- read.csv("ReMap Processed.csv")
  remap.filt <- remap[which(remap$LvlA | remap$LvlB | remap$LvlC | remap$LvlD),]
  
  # intersection of tfs
  common <- intersect(jaspar$tf, remap.filt$TF)
  
  # intersect
  i <- list()
  for (j in common) {
    print(j)
    
    # jaspar data for this tf, make into bed
    x <- jaspar[which(jaspar$tf == j),]
    x <- data.frame(chr = splitter(x$Peak, "_", 1),
                    start = x$tf.start, 
                    end = x$tf.end)
    write.bed(x, "../Scratchspace/X.bed")
    
    # remap data for this tf, make into bed
    y <- remap.filt[which(remap.filt$TF == j),]
    y$TF.Coord <- sub(":", "-", y$TF.Coord)
    split <- strsplit(y$TF.Coord, "-")
    y <- data.frame(chr = sapply(split, "[", 1),
                    start = sapply(split, "[", 2),
                    end = sapply(split, "[", 3),
                    Enh = y$Enh)
    write.bed(y, "../Scratchspace/Y.bed")
    
    # bedtools
    call <- paste("intersectBed",
                "-a", "../Scratchspace/X.bed",
                "-b", "../Scratchspace/Y.bed",
                "-wb", 
                ">", "../Scratchspace/Z.bed")

  system(call, intern = FALSE, wait = TRUE) 
  
  # read in
  i[[j]] <- tryCatch(read.delim("../Scratchspace/Z.bed", header = FALSE),
                     error = function(x) return(NA))
  
  }
  
  i <- do.call("rbind", i)
  i <- i[,-4]
  colnames(i) <- c("chr", "Motif_start", "Motif_end", "Peak_start", "Peak_End", "Enh")  
  
  # further wrangling
  i$TF <- splitter(rownames(i), "\\.", 1)
  i$HitPerm <- i$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
  i$HitCore <- i$Enh %in% res.final$Enh[which(res.final$HitCore)]
  write.csv(i, file = "TF Motif and ChIP Overlap (All Enh).csv", row.names = FALSE)
  
  # filter
  i <- i[which(i$HitPerm),]
  
  # save
  write.csv(i, file = "TF Motif and ChIP Overlap.csv", row.names = FALSE)
  
## Mutate TF binding sites.
  # load in position frequency matrices, from: https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_redundant_pfms_jaspar.txt
  library(TFBSTools)
  pfm <- readJASPARMatrix("../../../PublicData/TF_BindingSites/JASPAR2022_CORE_vertebrates_redundant_pfms_jaspar.txt")
  
  ## Calculate position weight matrices from position frequency matrices
   
    make.pwm <- function(tf) {
      # collect the position frequency matrix for the tf
      x <- pfm[[tf]]  
      x <- x@profileMatrix
      
      ## Parameters
        # expected nucleotide frequencies
        enf <- pfm[[tf]]@bg # it's 0.25 for all nt across all tfs, thus arbitrary
        enf.single <- enf[1]
     
        # pseudocount constant
        pseudocount <- 0.8 # Per advice from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/  
        
      ## Get PWM by applying function across each nucleotide (column)
        y <- apply(x, 2, function(z) {
          # add pseudocount 
          z <- z + pseudocount
          
          # get frequency
          Fib <- z / sum(z) # frequency of individual bases
          CQb <- (pseudocount * length(Fib)) / sum(z) # proportion of counts that are pseudocounts
          
          # calculate
          numerator <- Fib / (sum(Fib + CQb)) # frequency of base b / sum-across-bases(frequency + (pseudocount*expected-freq)). 
          denominator <- enf.single + (enf.single * pseudocount) # expected base frequency + (pseudocount* expected base frequncy). basically, just 0.5
          
          output <- 100 * log2(numerator / denominator)
          return(output)
          
        })
        
    }
    
  pwm <- lapply(unique(i$TF), make.pwm)
  names(pwm) <- unique(i$TF)
 
  ## Get the sum of PWM for TF binding site   
    # first, get DNA sequences for motifs!
    x <- i[,1:3]
    x$ID <- paste0(x$chr, ":", x$Motif_start, "-", x$Motif_end)
    write.bed(x, "TF_Motif.bed")
    
    call <- paste("bedtools", "getfasta",
                "-fi", "/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/Genome/UCSC/hg38.fa",
                "-bed", "TF_Motif.bed", "-tab",
                ">", "TF_Motif.fa")

    
    system(call, intern = FALSE, wait = TRUE) 
  
    motif.fa <- read.delim("TF_Motif.fa", header = FALSE)
    i$Motif_seq <- motif.fa$V2
    
    pwm.length <- sapply(pwm, ncol)
    m <- match(i$TF, names(pwm.length))
    a <- cbind(i$Motif_seq, nchar(i$Motif_seq), pwm.length[m]) %>% as.data.frame()
    a$V2 <- as.numeric(a$V2); a$V3 <- as.numeric(a$V3)
    
    # function to score a sequence of nt on a given pwm
    score.nt <- function(seq, weights) {
      s <- strsplit(seq, "*")[[1]]
      
      score <- list()
      for (loop in 1:length(s)) {
        score[[loop]] <- weights[s[loop], loop]  
      }
      
      score <- do.call("c", score)
      return(score)
      
    }
    
    # for each TF
    window.scores <- window.sum <- best.window <- list()
    
    for (j in 1:nrow(i)) {
      print(j)

      # get each of its TF binding sites
      x <- i[j,]

      # get pwm
      y <- pwm[[x$TF]]

      # get the sequences corresponding to the motif
      z <- x$Motif_seq

      # slide the pwm along the window and get the best fit
      window.scores[[j]] <- window.sum[[j]] <- best.window[[j]] <- list()

      
      pwm.length <- ncol(y)
      motif.length <- nchar(z)
      
      # if the JASPAR motif is shorter than the PWM, discard (~25%)
      if (motif.length < pwm.length)  {
        window.scores[[j]] <- best.window[[j]] <- "JASPAR too short" 
        next        
      }
      
      # slide
        for (k in 1:(motif.length - pwm.length + 1)) { # start loop for each sliding window along the motif seq k

          substring <- substr(z, k, (k + pwm.length - 1)) %>% toupper()
          substring.rev <- Biostrings::complement(x = Biostrings::DNAString(substring)) %>% as.character() %>% toupper()

          # score on the pwm
          window.scores[[j]][[paste0(k, "_sense")]] <- score.nt(seq = substring, y)

          # don't forget to reverse complement!
          window.scores[[j]][[paste0(k, "_antisense")]] <- score.nt(substring.rev, y)

        } # end loop for each sliding window along the motif seq k

        # pick best
        window.sum[[j]] <- sapply(window.scores[[j]], sum)
        best.window[[j]] <- window.scores[[j]][which.max(window.sum[[j]])]

    } # end all rows of i
  
  ## Add information to i
    i$Motif_windowOptimal <- sapply(best.window, function(x) {
      if (class(x) == "character") {
         return(x)
      } else {
        return(names(x))
      }
    }) 
    
  ## Output
    save(window.scores, window.sum, best.window, file = "TF Motif and ChIP Overlap - Window Scores.rda")
    write.csv(i, file = "TF Motif and ChIP Overlap - Annotated.csv")
    
    
## Assess for TF enrichment
  x <- read.csv("TF Motif and ChIP Overlap (All Enh).csv")
  
  # table of enh has tf
  y <- table(x$Enh, x$TF) > 0
  y <- as.data.frame(y)
  y[setdiff(paste0("Enh", 1:979), rownames(y)),] <- FALSE # there are enh TFs without TFs, and these are autodropped from the df, thus return them
  y <- y[,apply(y, 2, any)] # removes 5 TFs not found in any enh?
  
  # vector of enh is hit
  hit.perm <- rownames(y) %in% res.final$Enh[which(res.final$HitPermissive)]
  hit.core <- rownames(y) %in% res.final$Enh[which(res.final$HitCore)]
  
  # fisher test
  z <- apply(y, 2, function(z) {
    test <- list()
    test$Perm <- table(z, hit.perm) %>% fisher.test()
    test$Core <- table(z, hit.core) %>% fisher.test()
    # output <- data.frame(Pval_PermCore = paste0(signif(test$Perm$p.value, 3), "/", signif(test$Core$p.value, 3)),
    #                      Odds_PermCore = paste0(signif(test$Perm$estimate, 3), "/", signif(test$Core$estimate, 3)))
    output <- data.frame(Pval_Perm = signif(test$Perm$p.value, 3),
                         Pval_Core = signif(test$Core$p.value, 3),
                         Odds_Perm = signif(test$Perm$estimate, 3),
                         Odds_Core = signif(test$Core$estimate, 3))
  })
  
  z <- do.call("rbind", z)
  z <- z[order(z$Pval_Core),]
  
  # save
  write.csv(z, "TF Motif and ChIP Overlap (All Enh) Enrichments In Hits.csv")
  
  # plot
  p <- z[which(z$FDR < 0.05),]
  p$TF <- rownames(p)
  p$FDR <- -log10(p$FDR)
  rnk <- p$TF[order(p$FDR)]
  p$TF <- factor(p$TF, levels = rnk)
  
  pdf(file = "TF Motif and ChIP Overlap (All Enh) Enrichments In Hits.pdf", height = 3, width = 5)
  ggplot(p, aes(x = TF, y = FDR)) +
    geom_col(fill = "darkorange1") +
    theme_bw() +
    # scale_fill_manual(values = c("darkorange1", "black")) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "TFBS Enrichment Within Core Hits\n(-log10(FDR), Fisher Test)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none", axis.title.x = invis,
          panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis)
  
  ggplot(p, aes(x = TF, y = Odds_Core)) +
    geom_col(fill = "darkorange1") +
    theme_bw() +
    # scale_fill_manual(values = c("darkorange1", "black")) +
    # geom_hline(yintercept = -log10(0.05), linetype = 2) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "TFBS Enrichment Within Core Hits\n(Odds Ratio)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none", axis.title.x = invis,
          panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis)
  dev.off()  
  
  
################################################################################################################################ #
## On CTCF ----
  
## In this section, I map high-confidence CTCF binding sites within NHAs, and see how they relate to enhancer-gene pairs.
  
  
## Create a bedfile where, for each enh-gene pair, it runs from tss to the nearest edge of the enhancer
  bed1 <- res.final[,c("Pair", "Enh.Pos", "Gene.TSS")]
  
  # get the chromosome
  bed1$chr <- splitter(bed1$Enh.Pos, ":", 1)
  
  # get the start, in this case the gene TSS
  bed1$Coord1 <- splitter(bed1$Gene.TSS, ":", 2) %>% as.numeric()
  
  # get the end, which is Enh's end or start (whichever is closer!)
  bed1$Coord2 <- apply(bed1, 1, function(x) {
    e.coord <- splitter(x[2], ":", 2)
    e.left <- splitter(e.coord, "-", 1) %>% as.numeric() 
    e.right <- splitter(e.coord, "-", 2) %>% as.numeric() 
    t.start <- as.numeric(x[5])
    use.left <- (abs(t.start - e.left)) < (abs(t.start - e.right))
    
    if (use.left) {
      return(e.left)
    } else {
      return(e.right)
    }
  })
  
  # the start coordinate is whichever is smaller of Coord1 and Coord2
  bed1$start <- apply(bed1[,c("Coord1", "Coord2")], 1, min)
  bed1$end <- apply(bed1[,c("Coord1", "Coord2")], 1, max) # and the opposite for the end
    
  # write
  bed1 <- bed1[,c("chr", "start", "end", "Pair")]
  write.bed(bed1, "CTCF_Bed_EnhGenePairs.bed")
  
  
## Create a bedfile with the high-confidence CTCF binding sites
  # read in file, and filter to CTCF
  
  r <- read.delim("/mnt/Data0/PROJECTS/CROPSeq/PublicData/TF_BindingSites/remap2022_nr_macs2_hg38_v1_0.bed", header = FALSE)
  r <- r[,1:4]
  r <- r[grep("CTCF", r$V4),]
  
  # annotate each peak by cell-types
  r$TF <- splitter(r$V4, ":", 1)
  r <- r[-which(r$TF == "CTCFL"),]
  r$Celltype <- splitter(r$V4, ":", 2)
  
  brain.ct <- c("neural-progenitor", "glioma", "tibial-nerve", "primary-glioblastoma", "glioblastoma", "NB-1643", "hippocampus",
                "cortical-interneuron", "neuron-progenitor", "neuroepithelilal-cells", "neuron", "nerve", "NPC", "dopaminergic-neuron",
                "neural", "brain-prefrontal-cortex", "SH-SY5Y", "SK-N-BE2-C", "HNPC", "SK-N-MC")
  
  r$Astrocyte <- grepl("astrocyte", r$Celltype) 
  r$Brain <- sapply(strsplit(r$Celltype, ","), function(x) (any(x %in% brain.ct)))
  r$Brain <- r$Brain & !(r$Astrocyte)
  r$NonBrain <- (!(r$Astrocyte) & (!(r$Brain)))
  
  # as bed
  bed2 <- r[,c(1:3,7:9)]
  write.bed(bed2, "CTCF_Bed_ChIPCTCF.bed")
  
  # overlap with JASPAR CTCF motifs for extra stringency
  call <- paste("intersectBed",
                "-a", "CTCF_Bed_ChIPCTCF.bed",
                "-b", "../../../PublicData/TF_BindingSites/MA0139.1.tsv", # access the file MA0139.1 from http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/, as that's CTCF!
                "-u", # Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B
                ">", "CTCF_Bed_ChIPCTCF_Robust.bed")
  
  system(call, intern = FALSE, wait = TRUE)
  
  robust <- read.delim("CTCF_Bed_ChIPCTCF_Robust.bed", header = FALSE) # 175008, but it was 723886 pre-filtering. So 25% survive...
  bed2$Robust <- paste0(bed2$V1, bed2$V2, bed2$V3) %in% paste0(robust$V1, robust$V2, robust$V3)

  # write
  write.bed(bed2[which(bed2$Astrocyte & bed2$Robust),], "CTCF_Bed_ChIPCTCF_Robust_Ast.bed")
  write.bed(bed2[which(bed2$Brain & bed2$Robust),], "CTCF_Bed_ChIPCTCF_Robust_Brain.bed")
  write.bed(bed2[which(bed2$NonBrain & bed2$Robust),], "CTCF_Bed_ChIPCTCF_Robust_NonBrain.bed")

  
## Intersect and read in 
  x <- list()
  for (j in c("Ast", "Brain", "NonBrain")) {
    # intersect
    files <- c(paste0("CTCF_Bed_ChIPCTCF_Robust_", j, ".bed"),
               paste0("CTCF_Bed_Intersect_", j, ".bed"))
    
    call <- paste("intersectBed",
                  "-a", "CTCF_Bed_EnhGenePairs.bed",
                  "-b", files[1],
                  "-wa", "-c", # count of b overlapping a  
                  ">", files[2])
    
    system(call, intern = FALSE, wait = TRUE) 
    
    # read in 
    x[[j]] <- read.delim(files[2], header = FALSE)
  }
  
## Augment these dataframes
  ctcf <- data.frame(res.final[,c("Pair", "Enh", "Gene")],
                     CTCF_Ast = x$Ast$V5,
                     CTCF_Brain = x$Brain$V5,
                     CTCF_Other = x$NonBrain$V5,
                     CTCF_PerKb_Ast = x$Ast$V5 / res.final$Gene.Distance * 1000,
                     CTCF_PerKb_Brain = x$Brain$V5 / res.final$Gene.Distance * 1000,
                     CTCF_PerKb_Other = x$NonBrain$V5 / res.final$Gene.Distance * 1000,
                     res.final[,c("FDR.SCEPTRE", "HitCore", "HitPermissive", "HitCategory", 
                                 "Gene.Distance", "Gene.Distance.Bin", "Gene.Nearest", "Gene.Upstream")])
  write.csv(ctcf, "CTCF - Results Summary.csv")
  
  ## The primary question: does the number of ctcf binding sites per kb per bin differ between hits and non hits?
    # plot in exploration
    p <- ctcf[,c("CTCF_PerKb_Ast", "CTCF_PerKb_Brain", "CTCF_PerKb_Other", "HitPermissive", "Gene.Distance", "Gene.Distance.Bin")]
    p$Gene.Distance.Bin <- factor(p$Gene.Distance.Bin, levels = c("2-10kb", "10-50kb", "50-100kb", "100-500kb"))
    levels(p$Gene.Distance.Bin) <- paste0(levels(p$Gene.Distance.Bin), " (Has ", table(p$HitPermissive, p$Gene.Distance.Bin)[2,], " hits)")
    p$HitPermissive <- factor(p$HitPermissive, levels = c("TRUE", "FALSE"))
    levels(p$HitPermissive) <- c("Hit", "ns")
    # p$HitPermissive <- p$HitPermissive & !(p$HitCore)
    p <- melt(p, id.vars = c("HitPermissive", "Gene.Distance", "Gene.Distance.Bin"))
    p$variable <- sub("CTCF_PerKb_", "", p$variable)
  
    pdf(file = "CTCF - Exploratory Plot.pdf", height = 6, width = 8)
    ggplot(p, aes(x = variable, y = value, fill = HitPermissive, colour = HitPermissive)) +
      geom_violin(scale = "width", draw_quantiles = c(0.5), alpha = 0.1, width = 0.75) +
      geom_jitter(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
      theme_bw() +
      facet_wrap(~Gene.Distance.Bin, scales = "free") +
      scale_fill_manual(values = c("firebrick1", "cornflowerblue")) +
      scale_colour_manual(values = c("firebrick1", "cornflowerblue")) +
      # scale_y_continuous(expand = c(0,0)) +
      guides(fill = guide_legend(title = "Enhancer-Gene Pair"), colour = guide_legend(title = "Enhancer-Gene Pair")) +
      theme(panel.grid.major.x = invis) +
      labs(x = "CTCF Datasets", y = "CTCF Sites per Kb between Enhancer-Gene Pair")
    
    ggplot(p, aes(x = variable, y = value, fill = Gene.Distance.Bin, colour = Gene.Distance.Bin)) +
      geom_violin(scale = "width", draw_quantiles = c(0.5), alpha = 0.1, width = 0.75) +
      geom_jitter(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
      theme_bw() +
      facet_wrap(~HitPermissive) +
      # scale_fill_manual(values = c("firebrick1", "cornflowerblue")) +
      # scale_colour_manual(values = c("firebrick1", "cornflowerblue")) +
      # scale_y_continuous(expand = c(0,0)) +
      guides(fill = guide_legend(title = "Enhancer-Gene Pair"), colour = guide_legend(title = "Enhancer-Gene Pair")) +
      theme(panel.grid.major.x = invis) +
      labs(x = "CTCF Datasets", y = "CTCF Sites per Kb between Enhancer-Gene Pair")
    
    ggplot(p, aes(x = variable, y = value, fill = HitPermissive, colour = HitPermissive)) +
      geom_violin(scale = "width", draw_quantiles = c(0.5), alpha = 0.1, width = 0.75) +
      geom_jitter(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
      theme_bw() +
      # facet_wrap(~Gene.Distance.Bin, scales = "free") +
      scale_fill_manual(values = c("firebrick1", "cornflowerblue")) +
      scale_colour_manual(values = c("firebrick1", "cornflowerblue")) +
      # scale_y_continuous(expand = c(0,0)) +
      guides(fill = guide_legend(title = "Enhancer-Gene Pair"), colour = guide_legend(title = "Enhancer-Gene Pair")) +
      theme(panel.grid.major.x = invis) +
      labs(x = "CTCF Datasets", y = "CTCF Sites per Kb between Enhancer-Gene Pair")
    dev.off()
 
  
## Statistics
  stats <- list()
  for (j in c(unique(ctcf$Gene.Distance.Bin)[c(1,4,2,3)], "All")) {
    # filter data
    keep.cols <- c(7:9, 12, 15)
    if (j == "All") {
      x <- ctcf[, keep.cols]
    } else {
      x <- ctcf[which(ctcf$Gene.Distance.Bin == j), keep.cols]
    }

    
    # calculate means
    means <- apply(x[,1:3], 2, function(y) {
      t.test(y ~ x$HitPermissive)$estimate
    })
    
    # calculate p
    wil <- apply(x[,1:3], 2, function(y) {
      wilcox.test(y ~ x$HitPermissive)$p.value
    })
    
    # calculate n
    n <- table(x$HitPermissive)
    
    # output
    s <- data.frame(splitter(names(wil), "_", 3), t(means), (wil)) 
    colnames(s) <- c("CTCF", "Mean.NonHit", "Mean.HitPerm", "Wilcox")
    s$Bin <- paste0(j, " (n=", paste(table (x$HitPermissive), collapse = "/"), ")")
    stats[[j]] <- s
    
  }
  
  stats <- do.call("rbind", stats)
  rownames(stats) <- 1:nrow(stats)
  stats <- stats[order(stats$CTCF),]
  write.csv(stats, file = "CTCF - Stats.csv")
  
  
## Aside: look at CTCF with promoters and enhancers
  ## Create bed file 
    # prom
    bed3 <- data.frame(V1 = (res.final$Gene.TSS), V2 = NA, V3 = NA, V4 = res.final$Gene) %>% unique()
    bed3$V2 <- as.numeric(splitter(bed3$V1, ":", 2))
    bed3$V3 <- bed3$V2 + 50
    bed3$V2 <- bed3$V2 - 1000
    
    # enh
    bed3 <- rbind(bed3, read.delim("../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed", header = FALSE))
    
    # save
    write.bed(bed3, "CTCF_Bed_PromAndEnh.bed")
    
  ## Intersect
    call <- paste("intersectBed",
                  "-a", "CTCF_Bed_PromAndEnh.bed",
                  "-b", "CTCF_Bed_ChIPCTCF_Robust_Ast.bed",
                  "-wa", "-c", # count of b overlapping a  
                  ">", "CTCF_Bed_Intersect_PromEnh.bed")
    
    system(call, intern = FALSE, wait = TRUE) 
    
    # read in 
    ctcf_PromAndEnh <- read.delim("CTCF_Bed_Intersect_PromEnh.bed", header = FALSE)
    ctcf_PromAndEnh <- ctcf_PromAndEnh[which(ctcf_PromAndEnh$V5 > 0),]
    z <- res.final
    z$CTCF_Prom <- z$Gene %in% ctcf_PromAndEnh$V4
  
################################################################################################################################ #
## On TADs ----
  
  
## Explore the TAD data
    tad <- read.delim("../../../FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38.bed", header = FALSE)
  
  # size distribution
  library(ggbeeswarm)
  p <- data.frame(Size = tad$V3- tad$V2)
  p$Size <- p$Size / 10^6
  
  pdf(file = "TAD - Size Distribution.pdf", height = 3, width = 3)
  ggplot(p, aes(x = ".", y = Size)) +
    geom_quasirandom(alpha = 0.5) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    labs(y = "Size Distribution of 1131 TADs (MB)", x = "") +
    theme(panel.border = invis, axis.line.y = element_line(), axis.text.x = invis, panel.grid.major.x = invis)
  dev.off()
  
  ## Overlaps between tads - partial and complete?
    # create bed
    x <- read.delim("../../../FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38.bed", header = FALSE)
    x$ID <- paste0("TAD_", 1:nrow(x))
    write.bed(x, "TAD_SelfIntersection_In.bed")
    
    # intersect with self
    call <- paste("intersectBed",
                  "-a", "TAD_SelfIntersection_In.bed",
                  "-b", "TAD_SelfIntersection_In.bed",
                  "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                  ">", "TAD_SelfIntersection_Out.bed")
    system(call, intern = FALSE, wait = TRUE) 
    
    self <- read.delim("TAD_SelfIntersection_Out.bed", header = FALSE) # the logic here is to intersect bed with itself
    colnames(self) <- c("A_chr", "A_start", "A_end", "A_ID", "B_chr", "B_start", "B_end", "B_ID", "Overlap")
  
    # remove duplicates
    self <- self[-which(self$A_ID == self$B_ID),] # where the a and b peaks are the same
    self <- self[-which(as.numeric(splitter(self$A_ID, "_", 2)) > as.numeric(splitter(self$B_ID, "_", 2))),] # where the a and b peak pair is the same as a b and a peak pair
    
    
    # calculate the fraction of overlap, where 1 means that a tad is entirely within another
    self$OverlapFraction <- apply(self[,-c(1,4,5,8)], 1, function (x) {
      sizeA <- x[2] - x[1]
      sizeB <- x[4] - x[3]
      frac <- x[5] / min(sizeA, sizeB)
      
    })
  
    # plot
    self$OverlapComplete <- factor(self$OverlapFraction == 1)
    levels(self$OverlapComplete) <- c("2 TAD Windows\nPartially Overlap", "Smaller TAD Wholly\nInside Larger TAD")
    pdf(file = "TAD - Self Overlap Distribution.pdf", height = 3, width = 3)
    ggplot(self, aes(x = OverlapComplete)) +
      geom_bar() +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "Count of TAD Intersections") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.title.x = invis)
    dev.off()
  
 
## Intersect EGPs with TADs
  # here, reuse the file from the CTCF section, which contains the intervening section between E and P
  call <- paste("intersectBed",
                "-a", "CTCF_Bed_EnhGenePairs.bed",
                "-b", "TAD_SelfIntersection_In.bed",
                "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                # "-wa", "-wb", "-f 1",
                ">", "TAD_Intersection.bed")
  
  system(call, intern = FALSE, wait = TRUE) 
  
## Read in
  tad <- read.delim("TAD_Intersection.bed", header = FALSE)
  colnames(tad) <- c("EGP_Chr", "EGP_Start", "EGP_End", "EGP", "TAD_Chr", "TAD_Start", "TAD_End", "TAD_ID", "Overlap_bp")
  
## How much of the EGP is contained within a tad
  tad$EGP_Distance <- (tad$EGP_End - tad$EGP_Start)
  tad$Overlap_Frac <- tad$Overlap_bp / tad$EGP_Distance
  

## Now categorise each EGP
  # first, get all EGPs
  x <- res.final[,c("Pair", "Gene.Distance.Bin", "HitPermissive", "HitCore", "HitCategory")]
  
  # define basic true/false statements to determine categories
  within <- tad$EGP[which(tad$Overlap_Frac == 1)] 
  cross <- tad$EGP[which(tad$Overlap_Frac < 1)] 
  
  # is the EGP fully contained within a single TAD?
  x$WithinTAD <- x$Pair %in% within

  # does the EGP cross a tad boundary?
  x$CrossTAD <- x$Pair %in% cross
  
  # combinations
  x$WithinAndCrossTAD <- x$WithinTAD & x$CrossTAD
  x$WithinNotCrossTAD <- x$Pair %in% setdiff(within, cross)
  x$CrossNotWithinTAD <- x$Pair %in% setdiff(cross, within)
  x$NoTAD <- rowSums(x[,6:10]) == 0 # fortunately 0, as was defined this way...

  
  table(x$Gene.Distance.Bin,x$CrossTAD)
  
## Are EGPs that have a TAD boundary between them less likely to be hits?
  # the number of cross-boundary EGPs is too low for < 50kb, so focus on those more distal
  x$Gene.Distance.Bin <- factor(x$Gene.Distance.Bin, levels = c("2-10kb", "10-50kb", "50-100kb", "100-500kb"))
  
  pdf(file = "TAD - EGPs Crossing Boundaries By Bin.pdf", height = 3, width = 5)
  ggplot(x, aes(x = Gene.Distance.Bin, fill = CrossTAD)) +
    geom_bar(position = "dodge") + 
    theme_bw() +
    labs(y = "Count of EGPs", x = "Distance Bin") +
    scale_y_continuous(expand = c(0,0))
  dev.off()  
  
## Stats
  
  
  
  p <- x[which(x$Gene.Distance.Bin == "100-500kb"),]
  table(p$HitPermissive, p$CrossTAD)
  
  
  table(x$WithinTAD, x$HitPermissive)
    
    
## Save
  write.csv(x, "TAD - EGP Annotation.csv")
  
  