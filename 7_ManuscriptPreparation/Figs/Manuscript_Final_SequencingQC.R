## Setup
  setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_SequencingQC/")
  library(tidyverse)
  library(rcartocolor)
  library(ggsci)
  library(ggplot2)
  library(ggbeeswarm)
  library(readxl)
  source("../../../FullScale/Scripts/Functions.R")
  source("../FinalFigureFunctions.R")

## Functions to write to disk in a traceable way
  pdf_qc <- function(figNo, title, h, w) {
    pdf(file = paste0("../Final/", figNo, " - Script qc - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_qc <- function(figNo, title, toPrint) {
    sink(paste0("../Final/", figNo, " - Script qc - ", title, ".txt"))
    print(toPrint)
    sink()
  }

################################################################################################################################ #
## On Cell QC ----

## Load metadata
  load("../../../FullScale/Data/Preprocessed/Metadata.rda")
  load("../../../FullScale/Results/2_DE/Sceptre Input Files.rda")
  rm(sceptre.exp); gc()
  
  used <- which(meta$AnyGuide)
  
  p <- (meta[,c("Library", "Ribo_Pct", "Mito_Pct", "UMI_Total", "MOI", "TransductionPool")])
  p$Cell <- rownames(p)
  
  
## Write out some pertinent stats
  qc <- data.frame(nCells = nrow(p),
                   nCellsWithGuide = length(used),
                   PercentCellsWithGuide = length(used) / nrow(p),
                   MeanMOI = mean(p$MOI),
                   MeanNonZeroMOI = mean(p$MOI[used]),
                   MedianMOI = median(p$MOI),
                   MedianNonZeroMOI = median(p$MOI[used]),
                   MeanDepth = mean(p$UMI_Total),
                   MedianDepth = median(p$UMI_Total))
  qc <- t(qc) %>% as.data.frame()
  
  sink_qc(figNo = "SFig10", title = "scRNAseq qc", qc)
  

## Cells per library
  ## Original cell calls from empty drops
    load("../../../FullScale/Results/1_Processing/Library-level/EmptyDrops.rda", verbose = TRUE)
    
    keep.cells <- lapply(cell.calling, function(x) {
      rownames(x)[which(x$FDR < 0.01)]
    })
    
    keep.cells <- sapply(keep.cells, length)

  ## Final cell calls
  
    q <- table(p$Library) %>% as.data.frame.table()
    colnames(q) <- c("Batch", "Filtered")
    q$EmptyDrops <- keep.cells
    q <- melt(q)
    
    colnames(q)[2] <- "Call"
    q$Call <- factor(q$Call, levels = c("EmptyDrops", "Filtered"))
    levels(q$Call) <- c("Pre-filtering", "Post-filtering")
    
    # pdf(file = "QC - Cells per batch.pdf", height = 3, width = 2)
    pdf_qc(figNo = "SFig10A", title = "cells per batch", h = 3, w = 2)
    ggplot(q, aes(y = Batch, x = value, colour = Batch, shape = Call)) +
      geom_segment(aes(x = 0, xend = value, y = Batch, yend = Batch), linetype = 2) +
      geom_point(size = 3.5) +
      scale_colour_manual(values = rev(pals$Primary_Darker[1:7])) +
      scale_shape_manual(values = c(16, 1)) +
      theme_bw() +
      guides(colour = "none") +
      # coord_flip() +
      scale_x_continuous(expand = c(0,0), limits = c(0,9500), breaks = c(0, 2500, 5000, 7500)) +
      theme(panel.border = invis, axis.line.x = element_line(), panel.grid = invis,
            axis.title.y = invis, legend.position = "bottom", legend.title = invis) +
      labs(x = "Cell count")
    
    dev.off() 
  
  
## Library-level QC
  q <- p[,1:4]
  q$UMI_Total <- q$UMI_Total / 1000
  q <- melt(q)
  q$Library <- factor(q$Library)
  levels(q$variable) <- c("Ribosomal %", "Mitochondrial %", "Total UMIs (1000's)")
  
  # pdf(file = "QC - Cell-level expression metrics.pdf", height = 3, width = 5.5)
  pdf_qc(figNo = "SFig10B", title = "cell-level qc", h = 3, w = 5.5)
  ggplot(q, aes(x = Library, y = value, colour = Library, fill = Library)) +
    geom_violin(draw_quantiles = 0.5, alpha = 0.8, scale = "width") +
    theme_bw() +
    scale_colour_manual(values = pals$Primary_Darker) +
    scale_fill_manual(values = pals$Primary) +
    facet_wrap(~variable, scales = "free_y", strip.position = "bottom") +
    labs(y = "Value") +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    theme(panel.grid = invis, axis.text.x = invis, axis.ticks.x = invis, axis.title.x = invis,
          legend.position = "bottom", panel.border = invis, axis.line.y = element_line())
  dev.off()
  


################################################################################################################################ #
## On guide assignment ----
  
## MOI
  # there are two versions of the MOI plot: one is all pools, and the other is separately per pool
  q <- p
  q$MOIFactor <- factor(q$MOI, levels = as.character(0:40))
  q <- as.data.frame(table(q$MOIFactor))
  
  m <- mean(p$MOI)
  
  # pdf(file = "QC - MOI.pdf", height = 2.5, width = 3)
  pdf_qc(figNo = "SFig4C", title = "MOI", h = 2.5, w = 3)
  ggplot(q, aes(x = Var1, y = Freq)) +
    # geom_histogram(binwidth = 1) +
    geom_col(width = 1, fill = pals$One, colour = "white", size = 0.1) +
    theme_bw() +
    labs(x = "Number of guides per cell", y = "Number of cells") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(breaks = as.character(seq(0, 40, 2))) +
    geom_vline(xintercept = m + 1, colour = pals$Primary[8], linetype = 2) + # +1 is an offset as it is a histogram
    # coord_flip() +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = text90)
  dev.off()
  

  
## Number of cells per enhancer/guide
  # parameters
  max <- 1500
  
  ## Per enhancer
    r <- rowSums(sceptre.guide.pooled)
    r <- data.frame(Target = names(r), nCells = r)
    r <- unique(r)
    r <- r[grep("Enh", r$Target),]
    r$Bin <- cut(r$nCells, c(seq(0, 1000, 100), 2000), include.lowest = TRUE)
    levels(r$Bin) <- sub("\\(", "", levels(r$Bin)) %>% 
      sub("]", "", .) %>% 
      sub("\\[", "", .) %>% 
      sub(",", "-", .) %>% 
      gsub("e\\+03", "k", .) %>%
      sub("-2k", "+", .) 
    cats <- levels(r$Bin)
    
    # pdf(file = "QC - nCells per candidate.pdf", height = 2.5, width = 2.3)
    pdf_qc(figNo = "SFig4A", title = "nCells per candidate", h = 2.5, w = 2.3)
    ggplot(r, aes(x = Bin)) +
      geom_bar(fill = pals$One) +
      theme_bw() +
      scale_y_continuous(limits = c(0, max(table(r$Bin) + 50)), expand = c(0,0), breaks = c(0, 100, 200, 300)) +
      scale_x_discrete(limits = cats) +
      coord_flip() +
      theme(panel.border = invis, axis.line.x = element_line(), panel.grid = invis,
            axis.ticks.y = invis) +
      labs(x = "Cells with candidate targeted", y = "Number of candidates") 
    dev.off()
    
  ## Per guide
    r <- rowSums(sceptre.guide)
    r <- data.frame(Target = names(r), nCells = r)
    r <- r[grep("_", r$Target),]
    r <- unique(r)
    r <- r[grep("Enh", r$Target),]
    r$Bin <- cut(r$nCells, c(seq(0, 1000, 100), 2000), include.lowest = TRUE)
    levels(r$Bin) <- sub("\\(", "", levels(r$Bin)) %>% 
      sub("]", "", .) %>% 
      sub("\\[", "", .) %>% 
      sub(",", "-", .) %>% 
      gsub("e\\+03", "k", .) %>%
      sub("-2k", "+", .)  

    # pdf(file = "QC - nCells per guide.pdf", height = 2.5, width = 2.2)
    pdf_qc(figNo = "SFig4B", title = "nCells per guide", h = 2.5, w = 2.2)
    ggplot(r, aes(x = Bin)) +
      geom_bar(fill = pals$One) +
      theme_bw() +
      scale_y_continuous(limits = c(0, max(table(r$Bin) + 50)), expand = c(0,0)) +
      scale_x_discrete(limits = cats) +
      coord_flip() +
      theme(panel.border = invis, axis.line.x = element_line(), panel.grid = invis, axis.ticks.y = invis) +
      labs(x = "Cells assigned with sgRNA", y = "Number of sgRNAs") 
    dev.off()

################################################################################################################################ #
## Plots on Positive Control Guides ----

## Read in
  x <- read.csv("../../../FullScale/Results/2_DE/Pos/Results Matrix.csv")

## Volcano
  x$Hit <- x$FDR < 0.05
  x$Hit <- factor(x$Hit)
  levels(x$Hit) <- c("ns", "FDR < 0.05")

  x$Log2FC <- exp(x$LogFC) %>% log2()
  
  # pdf(file = "Positive Controls - Volcano.pdf", height = 2.1, width = 3)
  pdf_qc(figNo = "SFig4D", title = "positive control volcano", h = 2.1, w = 3)
  ggplot(x, aes(x = Log2FC, y = -log10(P), colour = Hit)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_colour_manual(values = pals$Hits) +
    # scale_colour_lancet() +
    scale_x_continuous(limits = c(-3,1), expand = c(0,0)) +
    geom_vline(xintercept = 0, linetype = 2) +
    # geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(expand = c(0,0.1), breaks = c(0,4,8,12)) +
    labs(y = "-Log10 Sceptre P", x = "Log2 fold-change") +
    theme(panel.border = invis, axis.line = element_line(), legend.position = "right", panel.grid = invis,
          legend.title = invis)
  dev.off()
  

    
################################################################################################################################ #
## Plots on Negative Control Guides ----
  
  
## Load data 
  # negative controls
  load("../../../FullScale/Results/2_DE/Neg/SCEPTRE Output (NegE, Guide-level).rda")
  neg <- de.negE
  neg$NB.p <- 2 * pnorm(abs(neg$z_value), lower.tail = FALSE) 
  
  # enhancer-gene pairs
  res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  
  # positive controls
  pos <- read.csv("../../../FullScale/Results/2_DE/Pos/Results Matrix.csv")
  pos$NB.p <- 2 * pnorm(abs(pos$Z), lower.tail = FALSE) 
  
## QQplots 
  ## Function, using the make_qq_plot function 
  
  revlog_trans <- function(base = exp(1)) {
    trans <- function(x) {
      -log(x, base)
    }
    inv <- function(x) {
      base^(-x)
    }
    scales::trans_new(paste("revlog-", base, sep = ""),
                      trans,
                      inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf)
    )
  }
  
 
  
## Similar qqplot to Gasperini 2019
  # on NB p-values and SCEPTRE, compare the negative control (NTC) (downsampled) and EGP libraries
  # we shall downsample to the same n of 7759 for NTC and EGP libraries
  

  
  
  # prepare negative control data
  x1 <- neg
  # x1 <- x1[sample(rownames(x1), nrow(res.final), replace = FALSE),]  # downsample
  x1 <- data.frame(Pair = paste0(x1$gRNA_id, "_", x1$gene_id),
                   Type = "NTC",
                   # SCEPTRE = x1$p_value,
                   Observed = x1$NB.p)
  x1 <- x1[order(x1$Observed),]
  x1$Expected <- (1:nrow(x1)) / nrow(x1)
  
  
  # prepare egp data
  x2 <- res.final
  x2 <- data.frame(Pair = x2$Pair,
    Type = "Inactive EGP",
    # SCEPTRE = x2$P.SCEPTRE,
    Observed = x2$P.NB)
  x2$Type[which(x2$Pair %in% res.final$Pair[which(res.final$HitPermissive)])] <- "Functional EGP"
  x2 <- x2[order(x2$Observed),]
  x2$Expected <- (1:nrow(x2)) / nrow(x2)
  
  # prepare egp data, but with empirical p
  x3 <- res.final
  x3 <- data.frame(Pair = x3$Pair,
    Type = "Inactive EGP",
    # SCEPTRE = x3$P.SCEPTRE,
    Observed = x3$P.N50)
  x3$Type[which(x3$Pair %in% res.final$Pair[which(res.final$HitPermissive)])] <- "Functional EGP"
  x3 <- x3[order(x3$Observed),]
  x3$Expected <- (1:nrow(x3)) / nrow(x3)
  
  # plot negative binomial
  x_negBinom <- rbind(x1, x2)
  # pdf(file = "Negative Control - qq V2.pdf", height = 2.5, width = 2.5)
  pdf_qc(figNo = "SFig4E", title = "qqplot", h = 2.5, w = 2.6)
  ggplot(x_negBinom, aes(x = -log10(Expected), y = -log10(Observed), colour = Type)) +
    geom_point() +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1) + 
    scale_x_continuous(expand = c(0,0), limits = c(0,5.2)) +
    # scale_y_continuous(trans = "log2") +
    labs(x = "Expected null p-value (-log10)", y = "Observed p-value (-log10)") +
    scale_colour_manual(values = c(rev(pals$Hits), "grey20")) +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line(),
          legend.position = "none")
  
  dev.off()
  