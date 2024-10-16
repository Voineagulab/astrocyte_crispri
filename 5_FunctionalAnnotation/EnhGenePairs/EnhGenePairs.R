## This script looks at properties of the enhancer-gene pair, rather than either individually
## It includes looking at whether there are nearer/skipped genes than the given linkage, as well as TADs

################################################################################################################################ #
## Setup ----

## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/EnhGenePairs/")
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

  source("../../../Scripts/Functions.R")

## Load
  # expression
  # load("../../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  
  # results from earlier scripts
  res.final <- read.csv("../../2_DE/Enh/Results Final.csv")  
  annot.enh <- read.csv("../Chromatin/Final - Annotation Logical.csv")
  annot.gene <- read.csv("../Genes/Final - Annotation Logical.csv")
  

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  sig.colours <- c("black", "firebrick1")

## Setup dataframe
  ## Generic dataframe to store annotations pertinent at the pair level
    repEGP <- res.final[,c("Pair", "Enh", "Gene", "HitPermissive")]
    rownames(repEGP) <- repEGP$Pair
    colnames(repEGP)[4] <- "Hit"

  ## List of hit pairs
    hit.egp.vector <- repEGP$Pair[which(repEGP$Hit)]
    hit.egp <- repEGP[which(repEGP$Hit),]
  
## Functions
  # intersect with a “left outer join”. that is, for each feature in A report each overlap with B. if no overlaps are found, report a NULL feature for B
  # useful for binary calls
  loj <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-loj", 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
  
  # intersect with full report
  wawb <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-wa", "-wb", 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
 
  # read in a bed file
  read.bed <- function(dir) read.delim(dir, header = FALSE)
  
  # convert coordinate id to enhancer id
  convert.id <- function(x) {
    m <- match(x, guides$TargetCoord)
    y <- guides$TargetID[m]
    return(y)
  }
  
  
## Commonly used paths
  # candidate coordinates
  nha_dir_38 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed" 
  nha_dir_19 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed" 

  
################################################################################################################################ #
## On the intervening genes between enhancers and genes ----
  

## A more nuanced calculation of the nearest gene
  nrst <- res.final[,c("Pair", "Enh", "Gene", "HitPermissive"),]
  nrst[,c("nNearer", "nSkip", "nSkip_Exp", "nSkip_notExp", "Distance.Category")] <- "."
  # nrst$Distance.Category <- "."
  
  for (j in 1:nrow(nrst)) {
    
    # parameters for this run
    print(j)
    e <- nrst$Enh[j]
    g <- nrst$Gene[j]
    
    ## Check for condition 1: nearest
      if (res.final$Gene.Nearest[j]) { # note: calling res.final here, as they have the same row order
        nrst[j,c("nNearer", "nSkip", "nSkip_Exp", "nSkip_notExp")] <- 0
        nrst$Distance.Category[j] <- "Cond1_Nearest"
        next
      }
    
    ## If not nearest, which genes are nearer?
      # enhancer coordinates
      y <- guides$TargetCoord[which(guides$TargetID == e)] %>% unique()
      z <- splitter(y, ":", 2)
      e.start <- as.numeric(splitter(z, "-", 1))
      e.end <- as.numeric(splitter(z, "-", 2))
      e.mid <- round((e.start + e.end) / 2)
      
      # genes on its chromosome
      chr <- splitter(y, ":", 1)
      chr <- geneInfo[which(geneInfo$Chr == chr),]
      chr$TSS <- chr$Start; chr$TSS[which(chr$Strand == "-")] <- chr$End[which(chr$Strand == "-")]
    
      # sort
      chr$Dist <- e.mid - chr$TSS
      chr$Dist.Rank <- rank(abs(chr$Dist))
      chr$CurrentGene <- chr$Symbol == g
      chr$Side <- sign(chr$Dist)
      
      nearer <- chr[which(chr$Dist.Rank <= chr$Dist.Rank[which(chr$CurrentGene)]),]
      
    ## Are those genes on the same side of the enhancer?
      nearer$SameSide <- nearer$Side == nearer$Side[which(nearer$CurrentGene)]
    
    ## Are those genes expressed?
      nearer$Expressed <- nearer$Symbol %in% res.final$Gene
      
    ## Count
      nearer <- nearer[-which(nearer$CurrentGene),]
      
      nrst$nNearer[j] <- nrow(nearer)
      nrst$nSkip[j] <- sum(nearer$SameSide)
      nrst$nSkip_Exp[j] <- sum(nearer$SameSide & nearer$Expressed)
      nrst$nSkip_notExp[j] <- sum(nearer$SameSide & !(nearer$Expressed))
    
    ## Categorise
      # nearer <- nearer[-which(nearer$CurrentGene),] # already done above
      
      # condition 2: no skipping
      cond2 <- sum(nearer$SameSide) == 0 # that is, all of sameside is false
      
      # check condition 3: skip expressed gene
      cond3 <- any( nearer$SameSide & nearer$Expressed ) # at least one gene on the same side that is expressed
      
      # check for condition 4: skip only non-expressed genes
      if (any(nearer$SameSide)) {
        cond4 <- all(!(nearer$Expressed[which(nearer$SameSide)]))
      } else {
        cond4 <- FALSE
      }
      
      # assign
      if (cond2) nrst$Distance.Category[j] <- "Cond2_Distal_NoSkip"
      if (cond3) nrst$Distance.Category[j] <- "Cond3_Distal_SkipExp"
      if (cond4) nrst$Distance.Category[j] <- "Cond4_Distal_SkipNonExp"
      
      if(sum(cond2, cond3, cond4) > 1) nrst$Distance.Category[j] <- "Category Error" 
     
    
  }

## Save
  write.csv(nrst, file = "Intervening Gene Classification Between EGPs.csv")
  

################################################################################################################################ #
## TADs ----
  
  
## Read in
  tad <- read.delim("../../../../FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38.bed", header = FALSE)
  
    
## Some characteristics
  ## Size distribution
    p <- data.frame(Size = tad$V3 - tad$V2)
    p$Size <- p$Size / 10^6
    
    pdf(file = "TAD/Size Distribution.pdf", height = 3, width = 3)
    ggplot(p, aes(x = ".", y = Size)) +
      geom_quasirandom(alpha = 0.5) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
      labs(y = "Size Distribution of 1131 TADs (MB)", x = "") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.text.x = invis, panel.grid.major.x = invis)
    dev.off()
  
## Do TADs overlap?
  ## Setup
    x <- tad
    x$ID <- paste0("TAD_", 1:nrow(x))
    write.bed(x, "TAD/SelfIntersection_In.bed")
    
  ## Intersect with self
    call <- paste("intersectBed",
                  "-a", "TAD/SelfIntersection_In.bed",
                  "-b", "TAD/SelfIntersection_In.bed",
                  "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                  ">", "TAD/SelfIntersection_Out.bed")
    system(call, intern = FALSE, wait = TRUE) 
    
  ## Read in
    self <- read.delim("TAD/SelfIntersection_Out.bed", header = FALSE) # the logic here is to intersect bed with itself
    colnames(self) <- c("A_chr", "A_start", "A_end", "A_ID", "B_chr", "B_start", "B_end", "B_ID", "Overlap")
  
    # remove duplicates
    self <- self[-which(self$A_ID == self$B_ID),] # where the a and b peaks are the same
    self <- self[-which(as.numeric(splitter(self$A_ID, "_", 2)) > as.numeric(splitter(self$B_ID, "_", 2))),] # where the a and b peak pair is the same as a b and a peak pair
  
  ## Calculate the fraction of overlap, where 1 means that a tad is entirely within another
    self$OverlapFraction <- apply(self[,-c(1,4,5,8)], 1, function (x) {
      sizeA <- x[2] - x[1]
      sizeB <- x[4] - x[3]
      frac <- x[5] / min(sizeA, sizeB)
      
    })
  
  ## Plot
    self$OverlapComplete <- factor(self$OverlapFraction == 1)
    levels(self$OverlapComplete) <- c("2 TAD Windows\nPartially Overlap", "Smaller TAD Wholly\nInside Larger TAD")
    pdf(file = "TAD/Self Overlap Distribution.pdf", height = 3, width = 3)
    ggplot(self, aes(x = OverlapComplete)) +
      geom_bar() +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "Count of TAD Intersections") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.title.x = invis)
    dev.off()
  
 
    
    
## Intersect EGPs with TADs
  ## Create a bed file whose coordinates are the intervening section between E and P
    x <- res.final[,c("Pair", "Enh.Pos", "Gene.TSS")]
    
    # get the chromosome
    x$chr <- splitter(x$Enh.Pos, ":", 1)
    
    # get the start, in this case the gene TSS
    x$Coord1 <- splitter(x$Gene.TSS, ":", 2) %>% as.numeric()
    
    # get the end, which is Enh's end or start (whichever is closer!)
    x$Coord2 <- apply(x, 1, function(x) {
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
    x$start <- apply(x[,c("Coord1", "Coord2")], 1, min)
    x$end <- apply(x[,c("Coord1", "Coord2")], 1, max) # and the opposite for the end
    
    # write
    x <- x[,c("chr", "start", "end", "Pair")]
    write.bed(x, "TAD/Pair_Intervening_Coord.bed")
    
    
  # here, reuse the file from the CTCF section, which contains the intervening section between E and P
  call <- paste("intersectBed",
                "-a", "TAD/Pair_Intervening_Coord.bed",
                "-b", "TAD/SelfIntersection_In.bed", # despite the name, this is simply TAD coordinates with an annotation
                "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                ">", "TAD/Tad_Vs_Pair_Overlap.bed")
  
  system(call, intern = FALSE, wait = TRUE) 
  
  ## Read in
    tad_overlap <- read.delim("TAD/Tad_Vs_Pair_Overlap.bed", header = FALSE)
    colnames(tad_overlap) <- c("EGP_Chr", "EGP_Start", "EGP_End", "EGP", "TAD_Chr", "TAD_Start", "TAD_End", "TAD_ID", "Overlap_bp")
    
  ## How much of the EGP is contained within a tad?
    tad_overlap$EGP_Distance <- (tad_overlap$EGP_End - tad_overlap$EGP_Start)
    tad_overlap$Overlap_Frac <- tad_overlap$Overlap_bp / tad_overlap$EGP_Distance
    
  
  ## Now categorise each EGP
    # first, get all EGPs
    pair2tad <- res.final[,c("Pair", "Gene.Distance.Bin", "HitPermissive")]
    pair2tad$Gene.Distance.Bin <- factor(pair2tad$Gene.Distance.Bin, levels = c("2-10kb", "10-50kb", "50-100kb", "100-500kb"))
    
    # define basic true/false statements to determine categories
    within <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac == 1)] 
    cross <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac < 1)] 
    
    # is the EGP fully contained within a single TAD?
    pair2tad$WithinTAD <- pair2tad$Pair %in% within
  
    # does the EGP cross a tad boundary?
    pair2tad$CrossTAD <- pair2tad$Pair %in% cross
    
    # combinations
    pair2tad$WithinAndCrossTAD <- pair2tad$WithinTAD & pair2tad$CrossTAD
    pair2tad$WithinNotCrossTAD <- pair2tad$Pair %in% setdiff(within, cross)
    pair2tad$CrossNotWithinTAD <- pair2tad$Pair %in% setdiff(cross, within)
    pair2tad$NoTAD <- rowSums(pair2tad[,4:8]) == 0 # fortunately 0, as was defined this way...
  
  ## Save
    write.csv(pair2tad, "TAD/Pair-TAD Annotation.csv")
    
    
  
  