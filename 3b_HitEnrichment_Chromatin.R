## This script preprocesses all counts data output from CellRanger using R

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
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist.csv", row.names = 1)
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
  load("../2_DE/Enhancers - All Results Summary.rda")
  
## Set up enhancer lists
  enh.hits <- unique(res$Enh.Pos[which(res$Hit)])

  
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
  

    