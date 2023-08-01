## This script analyses SCEPTRE analyses of the positive control guides

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
  exp.thresh <- 2^-6
  p.thresh <- 0.1
  
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
    
    # nearest gene: old calculation which is the nearest of all tested genes for a given enh
      # min.dist <- aggregate(Gene.Distance~Enh, x, min, na.rm = TRUE)
      # m <- match(x$Enh, min.dist$Enh)
      # x$Gene.Nearest <- min.dist$Gene.Distance[m]
      # x$Gene.Nearest <- x$Gene.Distance == x$Gene.Nearest
    
    # nearest gene: new calculation for the nearest gene in the GTF for a given enh
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
## Find guide-level hits ----  
  
load("Guide-level/Guide-level SCEPTRE.rda") 
guidelvl <- de.enh.guidelvl # load then rename

## Begin the plain annotation process!
  guidelvl <- data.frame(Enh = ".",
                         Guide = guidelvl$gRNA_id,
                         Gene = guidelvl$gene_id,
                         Pair.Enh = ".",
                         Pair.Guide = ".",
                         P = guidelvl$p_value,
                         Z = guidelvl$z_value)
  
  # add enhancer target
  m <- match(guidelvl$Guide, guides$GuideID)
  guidelvl$Enh <- guides$TargetID[m]
  guidelvl$Pair.Enh <- paste0(guidelvl$Enh, "_", guidelvl$Gene)
  guidelvl$Pair.Guide <- paste0(guidelvl$Enh, "_", splitter(guidelvl$Guide, "_", 2), "_", guidelvl$Gene)
  
  # filter to enhancer-gene pair tested in the pooled analysis
  # this effectively adds a consistent distance and expression filter!
  guidelvl <- guidelvl[which(guidelvl$Pair.Enh %in% res$Pair),]
  
  # add n
  n <- rowSums(sceptre.guide)
  m <- match(guidelvl$Guide, names(n))
  guidelvl$N <- n[m]
  
  # fdr correct
  guidelvl$FDR <- p.adjust(guidelvl$P, method = "fdr")
  
  # hit annotation
  guidelvl$Hit <- guidelvl$FDR < p.thresh # in this analysis, but also...
  m <- match(guidelvl$Pair.Enh, res$Pair)
  guidelvl$Hit.Pool <- res$Hit[m] # in above pooled analyses!
  
  # categorise hits
  guidelvl$Cat <- "."
  guidelvl$Cat[which(guidelvl$Hit.Pool & !(guidelvl$Hit))] <- "Enh-level Only"
  guidelvl$Cat[which(guidelvl$Hit & !(guidelvl$Hit.Pool))] <- "Guide-level Only"
  guidelvl$Cat[which(guidelvl$Hit.Pool & guidelvl$Hit)] <- "Enh- & Guide-level"
  guidelvl$Cat[which(!(guidelvl$Hit.Pool) & !(guidelvl$Hit))] <- "Non-hit"

  any <- guidelvl$Pair.Enh[which(guidelvl$Hit | guidelvl$Hit.Pool)] %>% unique()
  guidelvl$Cat[which(guidelvl$Pair.Enh %in% any & guidelvl$Cat == "Non-hit")] <- "Other Guide Only"
  

## Add nb
  guidelvl$P.NB <- 2 * pnorm(abs(guidelvl$Z), lower.tail = FALSE) 
  

## Save
  write.csv(guidelvl, file = "Guide-level/Results Summary.csv")
  
# ## Visualisation 1: replication across guides
#   p <- guidelvl[which(guidelvl$Pair.Enh %in% any),]
#   tab <- cbind(table(p$Pair.Enh, p$Cat), table(p$Pair.Enh)) %>% as.data.frame()
#   tab <- melt(tab, id.vars = "V5")
#   
#   # tab <- table(p$Pair.Enh, p$Cat)
#   # tab <- tab / rowSums(tab)
#   # tab <- cbind(tab, table(p$Pair.Enh)) %>% as.data.frame()
#   # tab <- melt(tab, id.vars = "V5")
#     
#   tab$V5 <- factor(tab$V5)
#   levels(tab$V5) <- paste0("Enh has ", levels(tab$V5), " guide(s) above threshold\n(n=", aggregate(value~V5, data = tab, FUN = sum)$value, " guides from ", table(tab$V5) / 4, "enh)")
# 
#   pdf(file = "Guide-level/Number of Hit Guides Per Enhancer.pdf", height = 7, width = 10)
#   ggplot(tab, aes(x = variable, y = value, fill = variable)) +
#     geom_col() +
#     facet_wrap(~V5, scale = "free_y") +
#     theme_bw() +
#     labs(y = "Number of Guides in Hit Category") +
#     guides(fill = guide_legend(title = "Hit Category")) +
#     theme(axis.text.x = invis, axis.title.x = invis)
#   dev.off()
#   
# ## Visualisation 2: violin plots of suppression, different files for each $Cat
#   norm.exp <- nha@assays$VST@data[,used] %>% as.data.frame()
#   norm.exp <- norm.exp[which(rownames(norm.exp) %in% res$Gene),]
#   
#   plot.guides <- function(enh, gene) {
#     # get guides for the enh
#     g <- guides$GuideID[which(guides$TargetID == enh)]
#     
#     ## Cell expression
#       e <- norm.exp[gene,] %>% as.numeric()
#       
#       x <- list()
#       for (j in g) { 
#         if (length(which(sceptre.guide[j,] ==1)) <= 30) next
#         x[[j]] <- data.frame(Guide = j, Exp = e[which(sceptre.guide[j,] == 1)]) 
#         
#         }
#       x$Control <- data.frame(Guide = "Control", Exp = e[-apply(sceptre.guide[g,], 2, any)])
#       x <- do.call("rbind", x)
#       
#     ## Annotate expression
#       # match to results dataframe
#       m <- match(paste0(x$Guide, gene), paste0(guidelvl$Guide, guidelvl$Gene))
#       
#       # n
#       x$N <- guidelvl$N[m]
#       
#       x$Category <- guidelvl$Cat[m]
#       x$Category[which(x$Guide == "Control")] <- "Non-hit"
#       x$Category <- factor(x$Category, levels = unique(guidelvl$Cat))
#       
#       # rename
#       x$Guide <- factor(x$Guide)
#       levels(x$Guide)[-1] <- paste0(enh, "_", splitter(g, "_", 2), "\nn=", table(x$Guide)[-1])
#       
#     ## Colouration
#       cols <- ggsci::pal_lancet()(4)
#       cols <- c("grey80", cols)
#       names(cols) <- unique(guidelvl$Cat)
#     
#     ## Plot
#       print(ggplot(x, aes(x = Guide, y = Exp, fill = Category, shape = Category)) +
#         geom_violin(scale = "width", colour = "black", draw_quantiles = c(0.5)) +
#         geom_jitter(width = 0.2, colour = "black", fill = "black", alpha = 0.5, show.legend = FALSE) +
#         scale_fill_manual(values = cols) +
#         scale_shape_manual(values = c(NA, 21, 21, 21, 21)) +
#         theme_bw() +
#         labs(y = paste0(gene, " Norm Exp")) +
#         theme(panel.border = invis, panel.grid.major.x = invis, axis.line = element_line(),
#               axis.title.x = invis))
# 
#     
#   }
#   
#   pdf(file = "Guide-level/Violins - Enhancer Hits.pdf", height = 3, width = 6)
#   x <- guidelvl$Pair.Enh[which(guidelvl$Hit.Pool)] %>% unique()
#   for (j in x) {
#     print(j)
#     plot.guides(enh = splitter(j, "_", 1), gene = splitter(j, "_", 2) )
#   }
#   dev.off()
#   
#   pdf(file = "Guide-level/Violins - Guide-only Hits.pdf", height = 3, width = 6)
#   x <- guidelvl$Pair.Enh[which(guidelvl$Cat == "Guide-level Only")] %>% unique()
#   for (j in x) {
#     print(j)
#     plot.guides(enh = splitter(j, "_", 1), gene = splitter(j, "_", 2) )
#   }
#   dev.off()
  
################################################################################################################################ #
## Empirical Correction ----
  
## Load negative control guide-gene pairs
  load("../Neg/SCEPTRE Output (Neg, Guide-level).rda")
  load("../Neg/SCEPTRE Output (NegE, Guide-level).rda")
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


################################################################################################################################ #
## Define core and permissive hits ----
  
## Permissive
  permissive <- res$Pair[which(res$FDR.N50 < 0.1)]

## Core
  # determine which enhancers have > 2 guides significant at (bonferroni) P < 0.05 
  rep2guides <- sapply(unique(guidelvl$Pair.Enh), function(x) {
    rep <- guidelvl[which(guidelvl$Pair.Enh == x),]
    nRep <- length(which(rep$P.N50 < (0.05 / nrow(rep))))
    nRep >= 2
  })
  rep2guides <- rep2guides[which(rep2guides)] %>% names()
  
  core <- res$Pair[which((res$FDR.N50 < 0.1 & res$Pair %in% rep2guides) | (res$FDR.N50 < 0.1 & res$FDR < 0.1))]
  # core <- res$Pair[which((res$FDR.N50 < 0.1) | (res$Pair %in% rep2guides & res$FDR < 0.1))]
  
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
  write.csv(res.final, file = "Results Final.csv", row.names = FALSE)
  
  
      
################################################################################################################################ #
## Directionality and fold-change ----

  
  
  
## Volcano
  # sig.colours2 <- c(carto_pal(2, "Earth"), carto_pal(2, "ArmyRose")[2], "grey50")
  # sig.colours2 <- c("firebrick1", carto_pal(2, "Earth")[2], "grey80")
  sig.colours2 <- c(pal_lancet()(2), "grey80")
  p <- res.final[order(res.final$HitPermissive, decreasing = FALSE),]
  
  
  pdf(file = "Volcano.pdf", height = 2.5, width = maxw)
  ggplot(p, aes(x = logfc.vst, y = -log10(P.SCEPTRE), colour = HitPermissive, fill = HitCore))  +
    geom_point(shape = 21) +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line()) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(-1,1)) +
    labs(x = "Log Fold-change (VST)", y = "-log10 Unadjusted P (SCEPTRE") +
    scale_colour_manual(values = sig.colours2[c(3,2)]) +
    scale_fill_manual(values = sig.colours2[c(3,1)])
  
  ggplot(p, aes(x = logfc.vst, y = -log10(P.N50), colour = HitPermissive, fill = HitCore))  +
    geom_point(shape = 21) +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line()) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
    scale_x_continuous(limits = c(-1,1)) +
    labs(x = "Log Fold-change (VST)", y = "-log10 Unadjusted P (N50)") +
    scale_colour_manual(values = sig.colours2[c(3,2)]) +
    scale_fill_manual(values = sig.colours2[c(3,1)])
  
  ggplot(p, aes(x = Z, y = -log10(P.SCEPTRE), colour = HitPermissive, fill = HitCore))  +
    geom_point(shape = 21) +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line()) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(-31,31)) +
    labs(x = "Z Score", y = "-log10 Unadjusted P (SCEPTRE)") +
    scale_colour_manual(values = sig.colours2[c(3,2)]) +
    scale_fill_manual(values = sig.colours2[c(3,1)])
  
  ggplot(p, aes(x = Z, y = -log10(P.N50), colour = HitPermissive, fill = HitCore))  +
    geom_point(shape = 21) +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line()) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
    scale_x_continuous(limits = c(-31,31)) +
    labs(x = "Z Score", y = "-log10 Unadjusted P (N50)") +
    scale_colour_manual(values = sig.colours2[c(3,2)]) +
    scale_fill_manual(values = sig.colours2[c(3,1)])
  
   
  
  dev.off()
    
  
## Bias in fold-change relative to all genes
  p <- table(permissivePlusCore(res.final)$HitCategory, sign(permissivePlusCore(res.final)$Z))
  # p[2,] <- table(res.final$HitPermissive, sign(res.final$Z))[2,]
  p <- p / rowSums(p)
  p <- as.data.frame(p)
  colnames(p) <- c("Hit", "Sign", "Freq")
  # levels(p$Hit) <- c("ns", "FDR < 0.1")
  levels(p$Sign) <- c("Downregulated", "Upregulated")
  # levels(p$Hit)[2] <- paste0("P (n=", length(which(res.final$HitPermissive)), ")")

  pdf(file = "Downregulation Fraction.pdf", height = 3, width = 4)
  ggplot(p, aes(x = Hit, y = Freq, fill = Sign)) +
    geom_col(colour = "black", width = 0.7) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis, axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("dodgerblue1", "firebrick1")) +
    labs(y = "Fraction of Enhancer-Gene Pairs", x = "Hit Category")
  dev.off()
  

################################################################################################################################ #
## Distance ----
 
  
## Distribution of distances
  pdf(file = "Distance Distribution.pdf", width = maxw, height = 3)
  
  p <- permissivePlusCore(res.final)
  
  # density plot, comparing hit and non-hit pairings
  m <- aggregate(p$Gene.Distance, list(p$HitCategory), median)[,2]
  m <- m / 1000
  
  ggplot(p, aes(x = Gene.Distance / 1000, colour = HitCategory, fill = HitCategory)) +
    geom_density(alpha = 0.1) +
    theme_bw() +
    scale_colour_manual(values = sig.colours2) +
    scale_fill_manual(values = sig.colours2) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 510), breaks = c(50, 100, 200, 500)) +
    geom_vline(xintercept = m, colour = sig.colours2, linetype = 2) +
    labs(x = "Distance Between Enhancer-Gene Pair (kilobases)", y = "Density") +
    theme(panel.border = invis, panel.grid = invis, legend.position = c(0.8, 0.8)) 
  
  # above, for upregulated genes
  m <- aggregate(p$Gene.Distance[which(res.final$Z > 0)], list(p$HitCategory[which(res.final$Z > 0)]), median)[,2]
  m <- m / 1000
  
  ggplot(p[which(res.final$Z > 0),], aes(x = Gene.Distance / 1000, colour = HitCategory, fill = HitCategory)) +
    geom_density(alpha = 0.1) +
    theme_bw() +
    scale_colour_manual(values = sig.colours2) +
    scale_fill_manual(values = sig.colours2) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 510), breaks = c(50, 100, 200, 500)) +
    geom_vline(xintercept = m, colour = sig.colours2, linetype = 2) +
    labs(x = "Distance Between Enhancer-Gene Pair (kilobases) (Upregulated Genes)", y = "Density") +
    theme(panel.border = invis, panel.grid = invis, legend.position = c(0.8, 0.8)) 
  dev.off()
  
  
## Annotation as nearest gene
  p <- table(permissivePlusCore(res.final)$HitCategory, permissivePlusCore(res.final)$Gene.Nearest)
  p <- p / rowSums(p)
  p <- as.data.frame(p)
  colnames(p) <- c("Hit", "NearestGene", "Freq")
  # levels(p$Hit) <- c("ns", "FDR < 0.1")
  
  pdf(file = "Nearest Gene.pdf", height = 3, width = 4)
  ggplot(p, aes(x = Hit, y = Freq, fill = NearestGene)) +
    geom_col(colour = "black", width = 0.7) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis, axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("white", "black")) +
    labs(y = "Fraction of Enhancer-Gene Pairs", x = "Hit Category")
  dev.off()  
 

  
   
################################################################################################################################ #
## Hit Interaction ----
  
## In this section, I explore situations where multiple enhancers are identified for the same gene
  
## Which?
  mult.genes <- names(which(table(hits$Gene) > 1))
  
## Function
  violin.compensation <- function(gene, enhs, neg = FALSE) {
      # get guides for the enhs
      pattern <- paste0(enhs, "_") %>% paste(collapse = "|")  
      use <- guides[grep(pattern, guides$GuideID),]
      n <- nrow(use)
      
      # name guides
      names <- splitter(use$GuideID, "_chr", 1)
      
  
      p <- data.frame(VST = nha@assays$SCT@data[gene,],
                      Neg = nha@meta.data$TransductionPool == "Neg",
                      Sum = NaN,
                      nha@meta.data[,use$GuideID])
      
      guide.cols <- 4:(n+3)
      colnames(p)[guide.cols] <- names
      
      # tabulate
      p$Sum <- rowSums(p[,guide.cols])
      p$Sum <- as.factor(p$Sum)
      levels(p$Sum) <- paste0(levels(p$Sum), "\nn=", table(p$Sum))
      
      # plot
      ggplot(p, aes(x = Sum, y = log2(exp(VST)), fill = Sum, colour = Sum)) +
        geom_violin(colour = "black", scale = "width", width = 0.8, draw_quantiles = 0.5) +
        # geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
        scale_fill_carto_d(palette = "Mint") +
        scale_colour_carto_d(palette = "Mint") +
        # scale_colour_manual(values = enh.colours) +
        theme_bw() +
        labs(y = paste0(gene, " log2-Normalised Expression"), x = "Number of Different Guides At Linked Enhancers") +
        scale_y_continuous(expand = c(0,0), limits = ylim) +
        # annotate("text", x = 1.5, y = max(p$Exp * 0.85), label = label, size = 2) +
        theme(panel.border = invis, axis.line.y = axis, panel.grid.major.x = invis,
              axis.ticks.x = invis, legend.position = "none")
  }
  
## Apply
  pdf(file = "Enhancers - Effect of Multiple Guides Across Enhancers.pdf", height = 3.5, width = maxw)
  for (j in mult.genes) {
    x <- hits[which(hits$Gene == j),]
    print(violin.compensation(gene = j, enhs = x$Enh))
  }
  dev.off()
  
  
## Extend to non-significant situations
  
  
################################################################################################################################ #
## Visualisation ----
  
## Top hits as violins
  
  ## Function
  
    vPlot <- function(row, axis = element_line(), ylim = c(NA, NA)) {
      g <- hits$Gene[row]
      e <- hits$Enh[row]
      
      # label <- paste0(hits$Pair[row], "\n",
      #                 # "FDR = ", signif(hits$FDR[row], 2), "\n",
      #                 "Fold change = ", paste0(round(exp(hits$logfc.vst[row]), 2)))
      label <- paste0("e", hits$Gene[row])
      
      p <- data.frame(Exp = nha@assays$SCT@data[g,],
                      Enh = sceptre.perturb[e,],
                      Blank = ".")
      
      p$Enh <- factor(p$Enh)
      levels(p$Enh) <- c("Enh Not Targeted", "Enh Targeted")
      
      ggplot(p, aes(x = Blank, y = log2(exp(Exp)), fill = Enh, colour = Enh)) +
        geom_violin(scale = "width", width = 0.8, draw_quantiles = 0.5, adjust = 2, fill = "white") +
        # geom_boxplot(outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
        # scale_fill_manual(values = c("")) +
        theme_bw() +
        labs(x = label, y = "log2-Normalised Expression") +
        # scale_fill_carto_d(palette = "Geyser") +
        # scale_colour_carto_d(palette = "Geyser") +
        scale_y_continuous(expand = c(0,0), limits = ylim) +
        # annotate("text", x = 1.5, y = max(p$Exp * 0.85), label = label, size = 2) +
        theme(panel.border = invis, axis.line.y = axis, panel.grid.major.x = invis, axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 8),
              axis.text.x = invis, axis.ticks.x = invis, legend.position = "none", axis.text.y = element_text(size = 6))
      
    }
  
  ## Plot for top ten hits
    hits <- hits[order(hits$logfc.vst),]
  
    # apply function
    x <- lapply((1:13)[-c(2,6,10)], vPlot, ylim = c(NA, NA))  # 2, 6, and 9 have duplicate genes
    x <- lapply(x, function(x) {
      x + theme(axis.title.y = invis)
    })
    
    # alternative to deflection
    #  x <- lapply(x, function(y) {
    #   z <- y$labels$x
    #   n <- 11 - nchar(z)
    #   z <- paste0(paste0(rep(" ", n), collapse = ""), z)
    #   y$labels$x <- z
    #           # y$labels$x <- paste0("          ", y$labels$x)
    #   # y$labels$x <- substr(y$labels$x, 1, 10)
    #   return(y)
    # })
    
    # add legend
    # x[[10]] <- x[[10]] + theme(legend.position = c(-0.9, 0.9), legend.title = invis, legend.text = element_text(size = 6), legend.key.size = unit(0.3, "cm"))
    
    # pdf(file = "Enhancers - Violin Plots V2.pdf", height = 3, width = 6)
    pdf(file = "../../../NHMRC2022_IV/Plots/GJS_Enhancer_Violins_Top10_Small.pdf", height = 2, width = 4)
    # plot_grid(plotlist = x, nrow = 1, rel_widths = c(2.8, rep(1.5, 9)))
    a <- plot_grid(plotlist = x[1:5], nrow = 1)
    b <- plot_grid(plotlist = x[6:10], nrow = 1)
    plot_grid(a,b, ncol = 1)
    dev.off()
  
## A function to visualise individual guides for a given enhancer
    # function works for any gene/guide, regardless of whether a statistical inference was performed
    violin.single.enh <- function(gene, enh, neg = FALSE, remove.multiples = TRUE, reorder = TRUE) {
      use <- guides$GuideID[grep(enh, guides$GuideID)]
      
      # reorder the above based on genomic position
      if (reorder) {
        order <- splitter(use, ":", 2) %>% splitter("-", 1) %>% as.numeric() %>% order
        use <- use[order]  
      }
      
      
      n <- length(use)
      names <- paste0("Guide", 1:n)
      enh <- gsub("_", "", enh)
      
      
      p <- data.frame(VST = nha@assays$SCT@data[gene,],
                      Neg = nha@meta.data$TransductionPool == "Neg",
                      Category = "NTC",
                      nha@meta.data[,use])
      
      guide.cols <- 4:(n+3)
      colnames(p)[guide.cols] <- names
      
      # categorise
      for (j in names) p$Category[which(p[,j])] <- j
      p$Category[which(rowSums(p[,guide.cols]) > 1)] <- "Multiple"
      
      if (neg) {
        p$Category[which(p$Neg)] <- "Pure Negative"
        
        p$Category <- factor(p$Category, levels = c("NTC", "Pure Negative", names, "Multiple"))
        # levels(p$Category) <- paste0(levels(p$Category), "\nn=", table(p$Category))
        levels(p$Category) <- paste0(levels(p$Category), " (", table(p$Category), ")")  
      } else {
        p$Category <- factor(p$Category, levels = c("NTC", names, "Multiple"))
        # levels(p$Category) <- paste0(levels(p$Category), "\nn=", table(p$Category))  
        levels(p$Category) <- paste0(levels(p$Category), " (", table(p$Category), ")")  
      }
      
      if (remove.multiples) {
        p <- p[-which(as.character(p$Category) == levels(p$Category)[length(levels(p$Category))]),]
      }
      
      
      # plot
      ggplot(p, aes(x = Category, y = log2(exp(VST)), fill = Category, colour = Category)) +
        geom_violin(colour = "black", scale = "width", width = 0.8, draw_quantiles = 0.5) +
        geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
       # geom_boxplot(fill = "white", width = 0.2, position = position_dodge(width = 0.8)) +
        # scale_fill_manual(values = c("")) +
        theme_bw() +
        labs(y = paste0("Log2-Normalised Expression"), title = paste0(gene, ":", enh)) +
        # scale_fill_carto_d(palette = "Geyser") +
        # scale_colour_carto_d(palette = "Geyser") +
        scale_y_continuous(expand = c(0,0), limits = ylim) +
        # annotate("text", x = 1.5, y = max(p$Exp * 0.85), label = label, size = 2) +
        theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis,
              axis.ticks.x = invis, legend.position = "none", axis.title.x = invis,
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    
    # singleton example
    
    
    pdf(file = "../../../NHMRC2022_IV/GJS_Guide_Violins_SERPINE1_LGALS3_small.pdf", height = 3, width = 2)
    # violin.single.enh("LGALS3", "Enh333_", reorder = FALSE)
    violin.single.enh("LGALS3", "Enh332_", reorder = TRUE)
    violin.single.enh("LGALS3", "Enh333_", reorder = TRUE)
    # violin.single.enh("LGALS3", "Enh333_", reorder = FALSE)
    violin.single.enh("LGALS3", "Enh334_", reorder = TRUE)
    
    # violin.single.enh("SERPINE1", "Enh831_", reorder = FALSE)
    violin.single.enh("SERPINE1", "Enh831_", reorder = TRUE)
    dev.off()
    
    pdf(file = "../../../NHMRC2022_IV/GJS_Guide_Violins_ID3_small.pdf", height = 3, width = 2)
    # violin.single.enh("LGALS3", "Enh333_", reorder = FALSE)
    violin.single.enh("ID3", "Enh53_", reorder = TRUE)
    violin.single.enh("ID3", "Enh54_", reorder = TRUE)
    
    dev.off()
    
    
    # en masse example
    pdf(file = "Enhancers - Guide Effect Violins.pdf", height = 3, width = 4)
    for (j in 1:nrow(hits)) {
      print(j)
      print(violin.single.enh(gene = hits$Gene[j], enh = paste0(hits$Enh[j], "_")))
    }
    dev.off()
    
    pdf(file = "Enhancers - Guide Effect Violins (With Neg).pdf", height = 3, width = 4)
    for (j in 1:nrow(hits)) {
      print(j)
      print(violin.single.enh(gene = hits$Gene[j], enh = paste0(hits$Enh[j], "_"), neg = TRUE))
    }
    dev.off()
    
    
## Plot arbitrary enhancers for a given gene, splitting by guide
    violin.any.enh <- function(gene, enhs, neg = FALSE) {
      # get guides for the enhs
      pattern <- paste0(enhs, "_") %>% paste(collapse = "|")  
      use <- guides[grep(pattern, guides$GuideID),]
      n <- nrow(use)
      
      # name guides
      names <- splitter(use$GuideID, "_chr", 1)
      
  
      p <- data.frame(VST = nha@assays$SCT@data[gene,],
                      Neg = nha@meta.data$TransductionPool == "Neg",
                      Category = "NTC",
                      nha@meta.data[,use$GuideID])
      
      guide.cols <- 4:(n+3)
      colnames(p)[guide.cols] <- names
      
      
      # categorise
      for (j in names) p$Category[which(p[,j])] <- j
      p$Category[which(rowSums(p[,guide.cols]) > 1)] <- "Multiple"
      
      if (neg) {
        p$Category[which(p$Neg)] <- "Pure Negative"
        p$Category <- factor(p$Category, levels = c("NTC", "Pure Negative", names, "Multiple"))
        
      } else {
        p$Category <- factor(p$Category, levels = c("NTC", names, "Multiple"))
      }
      
      # finally, add n
      levels(p$Category) <- paste0(levels(p$Category), "\nn=", table(p$Category))
      
      # colour by category
      p$Colour <- "_NTC"
      for (j in unique(use$TargetID)) {
        p$Colour[grep(paste0(j, "_"), p$Category)] <- j
      }
      
      enh.colours <- ggsci::pal_locuszoom()(length(unique(p$Colour)))
      names(enh.colours) <- levels(p$Colour)
      
      
      # drop cells categorised with "multiple"
      p <- p[-which(p$Category == levels(p$Category)[length(levels(p$Category))]),]
      
      # plot
      ggplot(p, aes(x = Category, y = log2(exp(VST)), fill = Colour, colour = Colour)) +
        geom_violin(colour = "black", scale = "width", width = 0.8, draw_quantiles = 0.5) +
        geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
        scale_fill_manual(values = enh.colours) +
        scale_colour_manual(values = enh.colours) +
        # scale_colour_manual(values = enh.colours) +
        theme_bw() +
        labs(y = paste0(gene, " log2-Normalised Expression")) +
        scale_y_continuous(expand = c(0,0), limits = ylim) +
        # annotate("text", x = 1.5, y = max(p$Exp * 0.85), label = label, size = 2) +
        theme(panel.border = invis, axis.line.y = axis, panel.grid.major.x = invis,
              axis.ticks.x = invis, legend.position = "none",
              axis.title.x = invis, axis.text.x = element_text(hjust = 1, vjust = 0.3, angle = 90))
    }

    
    
    pdf(file = "../../../NHMRC2022_IV/GJS_LGALS3.pdf", height = 3.5, width = 4)
    guides <- guides[rev(rownames(guides)),]
    violin.any.enh("LGALS3", which.enh)    
    guides <- guides[rev(rownames(guides)),]
    dev.off()
  
    
################################################################################################################################ #
## Transcriptome-wide DE for hit enhancers ----
    
## Load files
  load("../Sceptre Input Files.rda")
  library(sceptre)
    
  ## Including empirical p
    load("../Neg/SCEPTRE Output (NegE, Guide-level).rda")
    neg <- de.negE
    
    ## Calculate negative binomial p
      neg$NB.p <- 2 * pnorm(abs(neg$z_value), lower.tail = FALSE) 

  
    ## Function
      total.ntc <- nrow(neg) + 1 # number of tests
      calc.emp <- function(pval) {
        lower <- length(which(neg$NB.p < pval)) # number of neg-gene pairs more significant 
        empP <- (lower + 1) / total.ntc # p-value
        return(empP)
      }
      
    
## Genes to test: all expressed above threshold
  exp.thresh <- 2^-6
  exp.mean <- rowMeans(nha@assays$RNA@data[, which(nha$AnyGuide)]) 
  exp.high <- names(exp.mean)[exp.mean > exp.thresh]
  
  if (exists("nha")) rm(nha); gc() # for RAM reasons...
    
## Function
  run.sceptre.tw <- function(group, group.type = "Guide", chunkSize = 200) {
    # set tests
    pairs <- data.frame(gene_id = exp.high,
                        gRNA_id = group,
                        pair_type = "candidate")    
    
    # set groupings
    if (group.type == "Guide") guide.mx <- sceptre.guide
    if (group.type == "Enhancer") guide.mx <- sceptre.guide.pooled
    
    # setup loop
    out <- list()
    nChunks <- ceiling(nrow(pairs) / chunkSize)
    
    # run loop
    for (j in (length(out)+1):nChunks) {

      # setup
      a <- Sys.time()
      run.name <- paste0("Chunk_", j)

      # get the current chunk of a given chunkSize
      chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
      chunk <- chunk[which(chunk < nrow(pairs))] # truncate if the range is outside the full dataframe
      chunk <- pairs[chunk,]

      # run sceptre
      out[[run.name]] <- run_sceptre_high_moi(gene_matrix = sceptre.exp[chunk$gene_id,],
                                                  combined_perturbation_matrix = guide.mx, 
                                                  covariate_matrix = sceptre.covar,
                                                  gene_gRNA_group_pairs = chunk,
                                                  side = sceptre.direction,
                                                  B = 1000)
      b <- Sys.time()

      print(paste0("Chunk ", j, " of ", nChunks, ": ", b-a))
      gc()
    }
    
    out <- do.call("rbind", out)
    
    return(out)
  }
    
  neat1_sceptre_raw <- list()
  
  neat1_sceptre_raw$Enh219_g3 <- run.sceptre.tw("Enh219_g3_chr11:65419185-65419165", group.type = "Guide")
  neat1_sceptre_raw$Enh220_g1 <- run.sceptre.tw("Enh220_g1_chr11:65420132-65420112", group.type = "Guide")
  neat1_sceptre_raw$Enh219 <- run.sceptre.tw("Enh219", group.type = "Enhancer")
  neat1_sceptre_raw$Enh220 <- run.sceptre.tw("Enh220", group.type = "Enhancer")
  
  save(neat1_sceptre_raw, file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Sceptre.rda")
  

  neat1_sceptre <- lapply(neat1_sceptre_raw, function(j) {
    # rename and trim columns
    rownames(j) <- j$gene_id
    j <- j[,c("p_value", "z_value")]
    colnames(j) <- c("P_Sceptre", "Z_Sceptre")
    
    # calculate empirical p
    j$P_NB <- 2 * pnorm(abs(j$Z_Sceptre), lower.tail = FALSE) 
    j$P_Empirical <- NA
    for (k in 1:nrow(j)) {
      print(k)
      j$P_Empirical[k] <- calc.emp(j$P_NB[k])
    }
    
    # output
    return(j[,c("P_Sceptre", "P_Empirical")])
  })

  neat1_sceptre <- do.call("cbind", neat1_sceptre)
  
  # multiple testing correction
  neat1_sceptre$Enh219.FDR_Empirical <- p.adjust(neat1_sceptre$Enh219.P_Empirical, method = "fdr")
  neat1_sceptre$Enh220.FDR_Empirical <- p.adjust(neat1_sceptre$Enh220.P_Empirical, method = "fdr")
  neat1_sceptre$Enh219.FDR_Sceptre <- p.adjust(neat1_sceptre$Enh219.P_Sceptre, method = "fdr")
  neat1_sceptre$Enh220.FDR_Sceptre <- p.adjust(neat1_sceptre$Enh220.P_Sceptre, method = "fdr")
  
  neat1_sceptre$Enh219_g3.FDR_Empirical <- p.adjust(neat1_sceptre$Enh219_g3.P_Empirical, method = "fdr")
  neat1_sceptre$Enh220_g1.FDR_Empirical <- p.adjust(neat1_sceptre$Enh220_g1.P_Empirical, method = "fdr")
  neat1_sceptre$Enh219_g3.FDR_Sceptre <- p.adjust(neat1_sceptre$Enh219_g3.P_Sceptre, method = "fdr")
  neat1_sceptre$Enh220_g1.FDR_Sceptre <- p.adjust(neat1_sceptre$Enh220_g1.P_Sceptre, method = "fdr")
  
  
  # output
  write.csv(neat1_sceptre, file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Sceptre - Processed.csv")
  
## Try with Seurat's FindMarkers
  neat1_seurat <- list()
  neat1_seurat$Enh220 <- FindMarkers(nha, group.by = "Enh220", ident.1 = "TRUE",
                   logfc.threshold = 0, min.pct = 0,
                   features = exp.high)
  
  neat1_seurat$Enh219 <- FindMarkers(nha, group.by = "Enh219", ident.1 = "TRUE",
                   logfc.threshold = 0, min.pct = 0,
                   features = exp.high)
    
  save(neat1_seurat, file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Seurat.rda")
  
  neat1_seurat <- lapply(neat1_seurat, function(j) {
    j <- j[rownames(neat1_sceptre),c("p_val", "avg_log2FC")]
    colnames(j) <- c("P_Seurat", "log2fc")
    j$FDR_Seurat <- p.adjust(j$P_Seurat, method = "fdr")
    return(j)
  })
  
  neat1_seurat <- do.call("cbind", neat1_seurat)
  write.csv(neat1_seurat, file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Seurat - Processed.csv")
  
## Pooled output
    neat1_de <- neat1_sceptre[,grep("FDR", colnames(neat1_sceptre))]
    neat1_de <- cbind(neat1_de, neat1_seurat[,grep("FDR", colnames(neat1_seurat))])
    neat1_de <- neat1_de[,-grep("_g", colnames(neat1_de))]
  
    # average gene rank
    rank_fun <- function(x) {
      x <- apply(x, 2, rank)  
      return(rowMeans(x))
    }
    neat1_de$AverageRank_Enh219 <- rank_fun(neat1_de[,grep("Enh219", colnames(neat1_de))])
    neat1_de$AverageRank_Enh220 <- rank_fun(neat1_de[,grep("Enh220", colnames(neat1_de))])
    neat1_de$AverageRank_Combined <- rank_fun(neat1_de[,grep("AverageRank", colnames(neat1_de))]) 
    neat1_de$AverageRank_Combined <- rank(neat1_de$AverageRank_Combined)
    
    write.csv(neat1_de, file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - FDR.csv")
  
## A plot
  # empirical p concordance
  x <- data.frame(Enh219 = neat1_sceptre$Enh219.P_Empirical,
                  Enh220 = neat1_sceptre$Enh220.P_Empirical,
                  NEAT1 = rownames(neat1_sceptre) == "NEAT1")
  
  pdf(file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Enh219 vs Enh220 (Empirical P).pdf", height = 4, width = 4.5)
  r <- cor(-log10(x$Enh219), -log10(x$Enh220)) %>% round(2)
  ggplot(x, aes(x = -log10(Enh219), y = -log10(Enh220), colour = NEAT1, alpha = NEAT1)) +
    geom_point() +
    theme_bw() +
    theme(panel.grid.minor = invis, panel.border = invis) +
    scale_colour_lancet() +
    scale_alpha_manual(values = c(0.4, 1)) +
    labs(x = paste0("Enh219 (-log10(Empirical P))", "\nr=", r), y = "Enh220 (-log10(Empirical P))")
  dev.off()
  
  # fold-change concordance
  x <- data.frame(Enh219 = neat1_seurat$Enh219.log2fc,
                  Enh220 = neat1_seurat$Enh220.log2fc,
                  NEAT1 = rownames(neat1_sceptre) == "NEAT1")
  
  pdf(file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Enh219 vs Enh220 (Fold-change).pdf", height = 4, width = 4.5)
  r <- cor(x$Enh219, x$Enh220) %>% round(2)
  ggplot(x, aes(x = (Enh219), y = (Enh220), colour = NEAT1, alpha = NEAT1)) +
    geom_point() +
    theme_bw() +
    theme(panel.grid.minor = invis, panel.border = invis) +
    scale_colour_lancet() +
    scale_alpha_manual(values = c(0.4, 1)) +
    labs(x = paste0("Enh219 (log2FC)", "\nr=", r), y = "Enh220 (log2FC))")
  dev.off()
  
  # z concordance
  x <- data.frame(Enh219 = neat1_sceptre_raw$Enh219$z_value,
                  Enh220 = neat1_sceptre_raw$Enh220$z_value,
                  NEAT1 = neat1_sceptre_raw$Enh219$gene_id == "NEAT1")
  
  pdf(file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Enh219 vs Enh220 (Z).pdf", height = 4, width = 4.5)
  r <- cor(x$Enh219, x$Enh220) %>% round(2)
  ggplot(x, aes(x = (Enh219), y = (Enh220), colour = NEAT1, alpha = NEAT1)) +
    geom_point() +
    theme_bw() +
    theme(panel.grid.minor = invis, panel.border = invis) +
    scale_colour_lancet() +
    scale_alpha_manual(values = c(0.4, 1)) +
    labs(x = paste0("Enh219 (Z)", "\nr=", r), y = "Enh220 (Z))")
  dev.off()
  
## Looking at individual guides..
  ## Output easier csv
    x <- neat1_sceptre[,grep("_g", colnames(neat1_sceptre))]
    x <- x[,grep("FDR", colnames(x))]
  
    write.csv(x, file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - FDR (Top Guide).csv")
  
  ## Quick plots
    x <- data.frame(Enh219 = neat1_sceptre$Enh219_g3.P_Empirical,
                    Enh220 = neat1_sceptre$Enh220_g1.P_Empirical,
                    NEAT1 = rownames(neat1_sceptre) == "NEAT1")
    
    pdf(file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Enh219 vs Enh220 (Top Guide) (Empirical P).pdf", height = 4, width = 4.5)
    r <- cor(-log10(x$Enh219), -log10(x$Enh220)) %>% round(2)
    ggplot(x, aes(x = -log10(Enh219), y = -log10(Enh220), colour = NEAT1, alpha = NEAT1)) +
      geom_point() +
      theme_bw() +
      theme(panel.grid.minor = invis, panel.border = invis) +
      scale_colour_lancet() +
      scale_alpha_manual(values = c(0.4, 1)) +
      labs(x = paste0("Enh219_g3 (-log10(Empirical P))", "\nr=", r), y = "Enh220_g1 (-log10(Empirical P))")
    dev.off()
    
    # z concordance
    # load("../../Scratchspace/NEAT1/Transcriptome-wide DE - Sceptre.rda")
    
    x <- data.frame(Enh219 = neat1_sceptre_raw$Enh219_g3$z_value,
                    Enh220 = neat1_sceptre_raw$Enh220_g1$z_value,
                    NEAT1 = neat1_sceptre_raw$Enh219_g3$gene_id == "NEAT1")
    
    pdf(file = "../../Scratchspace/NEAT1/Transcriptome-wide DE - Enh219 vs Enh220 (Top Guide) (Z).pdf", height = 4, width = 4.5)
    r <- cor(x$Enh219, x$Enh220) %>% round(2)
    ggplot(x, aes(x = (Enh219), y = (Enh220), colour = NEAT1, alpha = NEAT1)) +
      geom_point() +
      theme_bw() +
      theme(panel.grid.minor = invis, panel.border = invis) +
      scale_colour_lancet() +
      scale_alpha_manual(values = c(0.4, 1)) +
      labs(x = paste0("Enh219_g3 (Z)", "\nr=", r), y = "Enh220_g1 (Z))")
    dev.off()
    