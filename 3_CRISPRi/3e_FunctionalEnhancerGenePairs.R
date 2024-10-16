## This script analyses SCEPTRE output for enhancer-gene pairs
## Whilst we initially defined a "core" and "permissive" set of functional hits
## Using SCEPTRE and an empirical correction, We ultimately report the latter as
## We found that they validated in an independent NanoString experiment

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Seurat)
  library(sceptre)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(tidyverse)
  library(rcartocolor)

## Load
  source("../../../Scripts/Functions.R")
  load("../../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  load("../Sceptre Input Files.rda")

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  sig.colours <- c("black", "firebrick1")

## Thresholds
  exp.thresh <- 2^-6 # this number is the minimum normalised expression in the scRNAseq data. It corresponds to a CPM of ~1.64.
  p.thresh <- 0.1 # or, more accurately, the FDR threshold
  
## Other
  used <- which(nha$AnyGuide) # defines the subset of cells on which DE was calculated
  targ.enh <- guides$TargetID[which(guides$TargetCat == "Enh")] %>% unique()



  
################################################################################################################################ #
## Loading and annotation ----
  
## Load
  load("SCEPTRE Output.rda")


## Create dataframe, based on SCEPTRE calls
  res <- de.enh
  res$Pair <- paste0(res$gRNA_id, "_", res$gene_id)
  res <- res[,c(1,2,6,4,5)]
  
  colnames(res) <- c("Gene", "Enh", "Pair", "P", "Z")
  
## Annotate with distance
  annotate.distance <- function(x) {
    # guide coord
    m <- match(x$Enh, guides$TargetID)
    x$Enh.Pos <- guides$TargetCoord[m]
    
    # gene coord
    m <- match(x$Gene, geneInfo$Symbol)
    strand <- geneInfo$Strand[m]
    strand.minus <- which(strand == "-")
    strand.plus <- which(strand == "+")
    x$Gene.Pos <- geneInfo$ID[m]
    x$Gene.TSS <- "."
    x$Gene.TSS[strand.plus] <- geneInfo$Start[m][strand.plus]
    x$Gene.TSS[strand.minus] <- geneInfo$End[m][strand.minus]
    
    # distance from enh midpoint to tss
    x$Gene.Distance <- apply(x, 1, function(y) {
      e.start <- as.numeric(splitter(splitter(y[6], "-", 1), ":", 2))
      e.end <- as.numeric(splitter(y[6], "-", 2))
      e.mid <- round((e.start + e.end) / 2)
      
      g.start <- as.numeric(y[8])
      
      dist <- abs(g.start - e.mid)
      return(dist)
    })
    
    # bin distance
    x$Gene.Distance.Bin <- cut(x$Gene.Distance, c(2000,10000,50000,100000,500000, max(x$Gene.Distance) + 1)) 
    levels(x$Gene.Distance.Bin) <- c("2-10", "10-50", "50-100", "100-500", ">500")
 
    # nearest gene: calculation for the nearest gene in the GTF for a given enh
    nearest <- list()
    for (j in targ.enh) {
      # enhancer coordinates
      y <- guides$TargetCoord[which(guides$TargetID == j)] %>% unique()
      z <- splitter(y, ":", 2)
      e.start <- as.numeric(splitter(z, "-", 1))
      e.end <- as.numeric(splitter(z, "-", 2))
      e.mid <- round((e.start + e.end) / 2)
      
      # genes on its chromosome
      chr <- splitter(y, ":", 1)
      chr <- geneInfo[which(geneInfo$Chr == chr),]
      chr$TSS <- chr$Start; chr$TSS[which(chr$Strand == "-")] <- chr$End[which(chr$Strand == "-")]
      
      # nearest
      w <- which.min(abs(e.mid - chr$TSS))
      
      nearest[[j]] <- chr$Symbol[w]
    }
    
    nearest <- do.call("c", nearest)
    x$Gene.Nearest <- "."
    for (j in 1:nrow(x)) {
      w <- which(names(nearest) == x$Enh[j])
      w <- nearest[w]
      x$Gene.Nearest[j] <- x$Gene[j] == w
    }
    
    # stream
    x$Gene.Upstream <- as.numeric(splitter(x$Enh.Pos, "-", 2)) < as.numeric(x$Gene.TSS)
    
    # finally, reformat TSS
    x$Gene.TSS <- paste0(splitter(x$Gene.Pos, ":", 1), ":", x$Gene.TSS)
    
    # return
    return(x)
    
  }

  res <- annotate.distance(res)
  
## Annotate with mean expression
  # expression normalised to library size via Seurat
  mean <- rowMeans(nha@assays$RNA@data[rownames(nha@assays$RNA@data) %in% res$Gene, used]) # rna assay, not sct
  m <- match (res$Gene, names(mean))
  res$Gene.Exp <- mean[m]
  
  # expression below threshold
  res$Gene.Exp.High <- res$Gene.Exp > exp.thresh

## Filter
  res <- res[which(res$Gene.Distance < 500000 & res$Gene.Exp.High),]
  
## FDR
  res$FDR <- p.adjust(res$P, method = "fdr")
  
## N 
  n <- rowSums(sceptre.guide.pooled)
  m <- match(res$Enh, names(n))
  res$nCells <- n[m]

## Calculate Neg Binom p
  res$P.NB <- 2 * pnorm(abs(res$Z), lower.tail = FALSE) 
  
## Fold-change
  res$logfc.vst <- res$logfc <- NaN

  norm.exp <- nha@assays$RNA@data[,used] 
  norm.exp <- norm.exp[which(rownames(norm.exp) %in% res$Gene),]
  
  norm.exp.vst <- nha@assays$VST@data[,used] 
  norm.exp.vst <- norm.exp.vst[which(rownames(norm.exp.vst) %in% res$Gene),]
  
  for (j in 1:nrow(res)) {
    print(j)

    w <- which(sceptre.guide.pooled[res$Enh[j],] == 1)

    # fold change in normalised data
    a <- norm.exp[res$Gene[j],w] # use Seurat's normalisation, natural log of ((exp / lib size) * 10000) + 1
    b <- norm.exp[res$Gene[j],-w]
    res$logfc[j] <- mean(a) - mean(b)

    # fold change in VST-normalised data
    a <- norm.exp.vst[res$Gene[j],w] # use Seurat's normalisation, natural log of ((exp / lib size) * 10000) + 1
    b <- norm.exp.vst[res$Gene[j],-w]
    res$logfc.vst[j] <- mean(a) - mean(b)

  }

  Sys.time()

# Save
  write.csv(res, file = "Results Summary.csv")
  

################################################################################################################################ #
## Empirical Correction using the negative control libraries ----
  
## Load negative control guide-gene pairs
  # we have two negative control libraries
  load("../Neg/SCEPTRE Output (NegE, Guide-level).rda") # this has 50 negative control guides. they were in the same lentiviral library as the enhancer-targeting guides
  load("../Neg/SCEPTRE Output (Neg, Guide-level).rda") # this has 250 negative control guides. they were is a separate lentiviral library
  
  neg <- list()
  neg$N250 <- de.neg.guidelevel
  neg$N50 <- de.negE
  
## Calculate negative binomial p
  neg <- lapply(neg, function(x) {
    x$NB.p <- 2 * pnorm(abs(x$z_value), lower.tail = FALSE) 
    return(x)
  })
  
## Use this to calculate an empirical P for negative binomial tests on enh-gene pairs
  total.ntc <- lapply(neg, function(x) nrow(x) + 1)
    
  # for enh-level analysis
  res$P.N50 <- NA
  res$P.N250 <- NA
  
  for (j in 1:nrow(res)) {
    print(j)
    
      # for the smaller enh pool of 50
      lower <- length(which(neg$N50$NB.p < res$P.NB[j])) # number of neg-gene pairs more significant 
      e <- (lower + 1) / total.ntc$N50 # p-value
      res$P.N50[j] <- e
      
      # for the larger neg pool of 250
      lower <- length(which(neg$N250$NB.p < res$P.NB[j])) # number of neg-gene pairs more significant 
      e <- (lower + 1) / total.ntc$N250 # p-value
      res$P.N250[j] <- e
      
  }

  res$FDR.N50 <- p.adjust(res$P.N50, method = "fdr")
  res$FDR.N250 <- p.adjust(res$P.N250, method = "fdr")

  # for guide-level analysis
  guidelvl$P.N50 <- NA
  guidelvl$P.N250 <- NA
  
  for (j in 1:nrow(guidelvl)) {
    print(j)
    
      # for the smaller enh pool of 50
      lower <- length(which(neg$N50$NB.p < guidelvl$P.NB[j])) # number of neg-gene pairs more significant 
      e <- (lower + 1) / total.ntc$N50 # p-value
      guidelvl$P.N50[j] <- e
      
      # for the larger neg pool of 250
      lower <- length(which(neg$N250$NB.p < guidelvl$P.NB[j])) # number of neg-gene pairs more significant 
      e <- (lower + 1) / total.ntc$N250 # p-value
      guidelvl$P.N250[j] <- e
      
  }

  guidelvl$FDR.N50 <- p.adjust(guidelvl$P.N50, method = "fdr")
  guidelvl$FDR.N250 <- p.adjust(guidelvl$P.N250, method = "fdr")

## Save
  write.csv(guidelvl, file = "Guide-level/Results Final (guide-level).csv")

################################################################################################################################ #
## Define core and permissive hits ----
  
## Permissive set, which is the ultimate
  permissive <- res$Pair[which(res$FDR.N50 < 0.1)] # we ultimately used the N50 library to empirically correct, rather than N250, as it show better control of p-value inflation.

## Core
  # determine which enhancers have > 2 guides significant at (bonferroni) P < 0.05 
  rep2guides <- sapply(unique(guidelvl$Pair.Enh), function(x) {
    rep <- guidelvl[which(guidelvl$Pair.Enh == x),]
    nRep <- length(which(rep$P.N50 < (0.05 / nrow(rep))))
    nRep >= 2
  })
  rep2guides <- rep2guides[which(rep2guides)] %>% names()
  
  core <- res$Pair[which((res$FDR.N50 < 0.1 & res$Pair %in% rep2guides) | (res$FDR.N50 < 0.1 & res$FDR < 0.1))]

## Add to dataframe
  res$HitPermissive <- res$Pair %in% permissive
  res$HitCore <- res$Pair %in% core
  
  res$HitCategory <- "ns"
  res$HitCategory[which((res$HitPermissive) & (res$HitCore))] <- "C"
  res$HitCategory[which((res$HitPermissive) & !(res$HitCore))] <- "P!C"
  res$HitCategory <- factor(res$HitCategory, levels = c("C", "P!C", "ns"))
  levels(res$HitCategory) <- paste0(levels(res$HitCategory), " (n=", table(res$HitCategory), ")") 
  
  res.final <- res[,c(3,2,1,6:13,17,5,18:19,4,20:22,15,23:27)]
  colnames(res.final)[c(16,20)] <- c("P.SCEPTRE", "FDR.SCEPTRE")
  res.final$Gene.Distance.Bin <- paste0(res.final$Gene.Distance.Bin, "kb")
  
  
## Save
  write.csv(res.final, file = "Results Final.csv", row.names = FALSE)
  
  
      

 
   

    