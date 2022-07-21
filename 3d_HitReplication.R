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
## Comparison to K562 screens ----
  
  
## Preprocess input files
  ## Gasperini
    gasp <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx", sheet = 2) # data from at-scale 
    gasp <- gasp[grep("candidate_enhancer", gasp$Category),] # filter to enhancers
    gasp <- gasp[,c(3,4,5,2)] # coordinates (hg19) and id
    gasp <- unique(gasp) # duplicate rows due to different guides are removed
    
    write.bed(gasp, "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/Candidates_AtScale_hg19.bed")
    
    gasp.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/Candidates_AtScale_hg38.bed"
    lift.bed(bed.in.19 = "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/Candidates_AtScale_hg19.bed", 
             bed.out.38 = gasp.in,
             isList = FALSE) # lift to hg38
    
  ## Xie
    xie <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Xie2019/mmc2.xlsx", sheet = 1) 
    x <- sub(":", "-", xie$`region pos (hg38)`)
    xie <- data.frame(chr = splitter(x, "-", 1),
                      start = splitter(x, "-", 2),
                      end = splitter(x, "-", 3), 
                      id = splitter(xie$`Position (hg19)`, "\\|", 1) %>% gsub(">", "", .))    
    xie <- unique(xie)
    xie.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Xie2019/Candidates_hg38.bed"
    write.bed(xie, xie.in)    

## Run intersection
  candidates <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  k562.out <- "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/BED_Annotations_K562.bed"
  k562.in <- paste(gasp.in, xie.in, collapse = " ")
  
  # run  
  call <- paste("intersectBed",
                "-a", candidates,
                "-b", k562.in,
                "-wa",
                "-wb", # write the original entry in B for each overlap. 
                "-names gasp xie",
                ">", k562.out)

  system(call, intern = FALSE, wait = TRUE) 
  
## Read in and wrangle
  k562 <- read.delim(k562.out, sep = "\t", header = FALSE)
  colnames(k562) <- c("chr", "start", "end", "id", "study", "k562.chr", "k562.start", "k562.end", "k562.id", "V10", "V11")
  
## Annotate with the target gene in our Astrocytes
  k562$id <- paste0(k562$chr, ":", k562$start, "-", k562$end)
  k562$HitVoineagu <- k562$id %in% enh.hits  
  
## Compare hits
  x <- res[which(res$Hit),]  
  y <- k562[which(k562$HitVoineagu),]
  
  x <- x[which(x$Enh.Pos %in% y$id),]
  
  m <- match(x$Enh.Pos, y$id)
  x$GaspEnh <- y$k562.id[m]
  
  # all are from Gasperini, so find those hits!
  gasp.hits <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx", sheet = 3) # data from at-scale 
  gasp.hits <- gasp.hits[which(gasp.hits$Target_Site %in% x$GaspEnh),] # all are unique, can use a simpler match
  m <- match(x$GaspEnh, gasp.hits$Target_Site)
  x$GaspGene <- gasp.hits$target_gene_short[m]

  # remove rows and columns for clarity
  x <- x[-which(is.na(x$GaspGene)),c(1:4, 24:25)]  
  
  # save
  write.csv(x, file = "K562 - Overlapping Peaks With Target Genes.csv")
  