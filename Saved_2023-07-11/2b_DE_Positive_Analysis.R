## This script analyses SCEPTRE analyses of the positive control guides

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Pos/")
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

## Other
  used <- which(nha$AnyGuide) # defines the subset of cells on which DE was calculated
  
    
################################################################################################################################ #
## Format results matrix ----
  
## Load
  load("SCEPTRE Output.rda")
  
## Reformat
  colnames(de.pos) <- c("Gene", "Guides", "PairType", "P", "Z")
  
## Call hits
  de.pos$FDR <- p.adjust(de.pos$P, method = "fdr")
  de.pos$Bonf <- p.adjust(de.pos$P, method = "bonferroni")
  de.pos$Hit <- de.pos$Bonf < 0.01
  
## Add expression metrics
  # mean expression
  mean <- rowMeans(nha@assays$RNA@data[rownames(nha@assays$RNA@data) %in% de.pos$Gene, used]) # rna assay, not sct
  m <- match (de.pos$Gene, names(mean))
  de.pos$ExpMean <- mean[m]
  
  # calculate a fold-change
  de.pos$LogFC <- "."
  
  for (j in 1:nrow(de.pos)) {
    print(j)
    gene <- de.pos$Gene[j]
    guide <- de.pos$Guides[j]
    
    if (!(gene %in% rownames(nha))) next
    
    # use Seurat's normalisation, natural log of ((exp / lib size) * 10000) + 1
    k <- nha@assays$VST@data[gene, used]
    
    # fold change
    fc <- aggregate(k, list(sceptre.guide.pooled[guide,]), mean)
    fc <- fc$x[2] - fc$x[1]
    
    # return
    de.pos$LogFC[j] <- fc
  }
  
  de.pos$LogFC <- as.numeric(de.pos$LogFC)
  
## Write
  write.csv(de.pos, file = "Results Matrix.csv")
  # de.pos <- read.csv("Results Matrix.csv")
  
  
################################################################################################################################ #
## Relationship between expression and DE success ----


## Volcano
  pdf(file = "Volcano Plot.pdf", height = 2.2, width = maxw/3*2)
  ggplot(de.pos, aes(x = LogFC, y = -log10(P), colour = Hit)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_colour_manual(values = c("black", "firebrick1")) +
    scale_x_continuous(limits = c(-2,2), expand = c(0,0)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(y = "-Log10 Unadjusted P", x = "Log Fold-change") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none")
  dev.off()
  
  
## DE versus mean expression
  p <- data.frame(p = -log10(de.pos$P),
                  Hit = de.pos$Hit,
                  MeanExp = de.pos$ExpMean)
  
  pdf(file = "Mean Expression Scatterplot.pdf", height = 3, width = 3)
  ggplot(p, aes(y = p, x = MeanExp, colour = Hit)) +
    geom_point() +
    theme_bw() +
    scale_colour_manual(values = c("black", "firebrick1")) +
    scale_y_continuous(expand = c(0,0), limits = c(0,16)) +
    theme(panel.border = invis, axis.line = element_line()) +
    labs(x = "Mean Normalised Expression", y = "-Log10 Unadjusted P") +
    NoLegend() +
    scale_x_continuous(limits = c(0, 4.5), expand = c(0,0))
  dev.off()


## Violin plot
  pdf(file = "Violins.pdf", height = 2.5, width = maxw / 1.5)
  for (j in 1:nrow(de.pos)) {
    print(j)

    gene <- de.pos$Gene[j]
    guide <- de.pos$Guides[j]

    p <- data.frame(Exp = nha@assays$VST@data[gene,used],
                    Guide = sceptre.guide.pooled[guide,])

    p$Guide <- factor(p$Guide)
    levels(p$Guide) <- c("Non-targeting", "Targeting")
    tab <- (table(p$Guide))
    levels(p$Guide) <- paste0(levels(p$Guide), " (n=", tab, ")")

    pval <- signif(de.pos[j,"Bonf"], 2)
    fc <- signif(exp(de.pos[j,"LogFC"]), 2)


    print(ggplot(p, aes(y = Exp, x = Guide, fill = Guide, shape = Guide)) +
      geom_violin(draw_quantiles = 0.5, colour = "white", adjust = 2) + # note that the adjust argument specifies the extent of smoothing
      geom_jitter(width = 0.2, alpha = 0.2) +
      theme_bw() +
      scale_y_continuous() +
      scale_shape_manual(values = c(NA, 1)) +
      labs(y = "NHA Normalised Expression", x = paste0("Bonferroni=", pval, ", Relative Expression = ", fc)) +
      scale_fill_manual(values = c("black", "darkorange1")) +
      annotate("text", x = 1.5, y = max(p$Exp) * 0.95, label = gene, size = 4) +
      theme(panel.border = invis, axis.line.y = element_line(),
            legend.position = "none", axis.ticks.x = invis, panel.grid = invis))

  }

  dev.off()

## Suppression effect size
  p <- data.frame(Rel = exp(de.pos$LogFC), Ct = "NHA", Sig = de.pos$Hit)
  p$Sig <- factor(p$Sig, levels = c("TRUE", "FALSE"))
  levels(p$Sig) <- c("Hit", "ns")

  pdf(file = "Suppression Efficiency.pdf", height = 2.2, width = maxw/3)
  ggplot(p, aes(x = Sig, y = Rel, fill = Sig, colour = Sig)) +
    geom_jitter(width = 0.1, shape = 21, size = 1) +
    stat_summary(fun = median, geom = "point", colour = "black", shape = "-", size = 10) +
    scale_fill_manual(values = c("firebrick1", "black")) +
    scale_colour_manual(values = c("firebrick1", "black")) +
    theme_bw() +
    # geom_hline(yintercept = 1, linetype = 2) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
    labs(y = "Relative expression from\npositive control genes", x = "Hit (Bonferroni < 0.01)") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", panel.grid.major.x = invis,
          axis.title.x = invis)
  dev.off()
  
  
################################################################################################################################ #
## Perturbation success as a function of dCas9-KRAB expression ----
  

dcas9 <- "dCas9-KRAB-T2A-BLAST-WPRE" # name of the construct in the expression matrix
  

## Function to visualise expression across bins
  x <- list() # to store correlation coefficients
  
  pdf(file = "dCas9.pdf", height = 2.5, width = maxw)
  
  for (j in 1:nrow(de.pos)) {
    print(j)
    
    gene <- de.pos$Gene[j]
    guide <- de.pos$Guides[j]

    # setup
    p <- data.frame(Exp = nha@assays$VST@data[gene, used],
                    Guide = sceptre.guide.pooled[guide,],
                    dcas9 = nha@assays$VST@counts[dcas9,used])
    p$Bin <- cut(p$dcas9, c(-1, 0, 5, 10, 20, 150))
    levels(p$Bin) <- c("0", "1-5", "6-10", "10-20", ">20")

    w <- which(p$Guide == 0)
    lev <- c("(No Targeting Guide)", levels(p$Bin))
    p$Bin <- as.character(p$Bin)
    p$Bin[w] <- "(No Targeting Guide)"
    p$Bin <- factor(p$Bin, levels = lev)

    # stats
    cor <- round(cor(p$dcas9[w], p$Exp[w]),2)
    x[[guide]] <- cor

    print(ggplot(p, aes(y = Exp, x = Bin, fill = Bin, shape = Bin)) +
      geom_violin(draw_quantiles = 0.5, colour = "white") +
      geom_jitter(width = 0.2, alpha = 0.2) +
      theme_bw() +
      scale_y_continuous() +
      scale_shape_manual(values = c(NA, 1, 1, 1, 1, 1)) +
      labs(y = paste0("Normalised ", guide, " Expression"), x = paste0("dCas9-KRAB Expression (VST Count)\nr=", cor)) +
      scale_fill_manual(values = c("black", "darkorange1", "darkorange1", "darkorange1", "darkorange1", "darkorange1" )) +
      annotate("text", x = 3.5, y = max(p$Exp) * 0.95, label = guide, size = 4) +
        geom_vline(xintercept = 1.5) +
      theme(panel.border = invis, axis.line.y = element_line(),
            legend.position = "none", axis.ticks.x = invis, panel.grid.major.x = invis))

  }

  dev.off()
  
  # look at x: the distribution of correlation coefficients
  x <- do.call("c", x)
  x <- data.frame(Cor = x)
  lim <- max(abs(x$Cor)) * 1.1
  
  pdf(file = "dCas9 Correlations.pdf", height = 3, width = 2)
  ggplot(x, aes(y = Cor, x = ".")) +
   geom_jitter(width = 0.1, size = 1) +
    stat_summary(fun = median, geom = "point", fill = "firebrick1", colour = "firebrick1", shape = "-", size = 10) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(expand = c(0,0), limits = c(-lim, lim)) +
    labs(y = "Correlation between dCas9-KRAB and\npositive control genes", x = "Hit (Bonferroni < 0.01)") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", panel.grid.major.x = invis,
          axis.title.x = invis, axis.text.x = invis, axis.ticks.x = invis)
    dev.off()

## A summary scatterplot of expression
  # do cells with more KRAB have on average more suppression across targets?
  # x-axis: dCas9-KRAB. y-axis: suppression of positive control gene.
    
  x <- list()
  
  for (j in 1:nrow(de.pos)) {
    print(j)
    
    gene <- de.pos$Gene[j]
    guide <- de.pos$Guides[j]

    # setup
    y <- data.frame(Exp = nha@assays$VST@counts[gene, used],
                    Guide = sceptre.guide.pooled[guide,],
                    dcas9 = nha@assays$VST@counts[dcas9,used])
    
    y$Exp <- scale(y$Exp)
    
    x[[gene]] <- y
  
  }
  
  y <- do.call("rbind", x)
  y$Exp <- as.numeric(y$Exp)
  z <- y[which(y$Guide > 0), ]
  z$Density <- get_density(z$dcas9, z$Exp, n = 30)
  r <- cor(z$Exp, z$dcas9) %>% signif(2)
  p <- cor.test(z$Exp, z$dcas9)$p.value %>% signif(2)
  
  pdf(file = "dCas9 Correlations, Alternate.pdf", height = 4, width = maxw)
  ggplot(z, aes(x = dcas9, y = Exp, colour = Density)) +
    geom_point() +
    # scale_colour_viridis_c(option = "D") + 
    scale_colour_carto_c(palette = "BluYl", direction = -1) +
    theme_bw() +
    scale_x_continuous(expand = c(0,0), limits = c(0, max(z$dcas9 * 1.02))) +
    geom_smooth(method = "lm", colour = "black", se = FALSE) +
    geom_hline(yintercept = 0, linetype = 2) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(y = "Expression of 124 Postive Control Genes in Targeting Cells\n(Z-transformed Per Gene)", x = paste0("dCas9-KRAB-T2A-BLAST-WPRE\nr=", r, " p=", p)) +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis, legend.position = c(0.9, 0.75))
  dev.off()

  
################################################################################################################################ #
## Are any of the positive control genes intrinsically interesting in the brain? ----

## SFARI
  sfari.pos <- targ.pos[targ.pos %in% sfari$gene.symbol]

  
################################################################################################################################ #
## Transcriptome-wide analysis ----
  
  
## Read in
  load("Transcriptome-wide - SCEPTRE Output.rda")


## Annotate 
  norm.exp <- as.matrix(nha@assays$RNA@data[unique(de.pos.tw$gene_id),used] )
  norm.exp.vst <- as.matrix(nha@assays$VST@data[unique(de.pos.tw$gene_id),used] )

  tw <- list()
  for ( j in unique(de.pos.tw$gRNA_id)) {
    print(j)
   
   
    # get  
    x <- de.pos.tw[which(de.pos.tw$gRNA_id == j),]
    
    # fdr
    x$FDR <- p.adjust(x$p_value, method = "fdr")
    
    # bonf
    x$Bonf <- p.adjust(x$p_value, method = "bonf")
    
    # fc
    w <- which(sceptre.guide.pooled[paste0("Pos_", j),] == 1)
    
    # fold change in normalised data
    a <- rowMeans(norm.exp[x$gene_id,w]) # use Seurat's normalisation, natural log of ((exp / lib size) * 10000) + 1
    b <- rowMeans(norm.exp[x$gene_id,-w])
    x$logfc <- (a) - (b)
    
    # fold change in VST-normalised data
    a <- rowMeans(norm.exp.vst[x$gene_id,w]) # use Seurat's normalisation, natural log of ((exp / lib size) * 10000) + 1
    b <- rowMeans(norm.exp.vst[x$gene_id,-w])
    x$logfc.vst <- (a) - (b)
    
    
    # fin
    tw[[j]] <- x
  }
  
  save(tw, file = "Transcriptome-wide - Annotated.rda")
  
## Summary
   tw.sum <- lapply(names(tw), function(x) {
    # get data
    p <- tw[[x]]
    
    # target gene suppression
    t <- sceptre.pairs.pos$gene_id[which(sceptre.pairs.pos$gRNA_id == x)]
    data.frame(Target = t,
               N = sum(sceptre.guide.pooled[paste0("Pos_", x),]),
               Target.P = p$p_value[which(p$gene_id == t)],
               Target.Sig = p$p_value[which(p$gene_id == t)] < 0.05 / 125, # bonf correct for number of pos guides
               Targ.FC = p$logfc[which(p$gene_id == t)],
               Targ.FC.vst = p$logfc.vst[which(p$gene_id == t)],
               OtherHits.Bonf = length(which(p$Bonf < 0.05)),
               OtherHits.FDR = length(which(p$FDR < 0.05)))
  })
  
   tw.sum <- do.call("rbind", tw.sum)
   
   write.csv(tw.sum, file = "Transcriptome-wide - Summary.csv", row.names = FALSE)
   
## Volcano plot, highlighting the expected target
  twPlot.volc <- lapply(names(tw), function(x) {
    print(x)
    p <- tw[[x]]
    p$Highlight <- p$gene_id == sceptre.pairs.pos$gene_id[which(sceptre.pairs.pos$gRNA_id == x)]
    p$Fade <- p$FDR > 0.05 & p$gene_id != sceptre.pairs.pos$gene_id[which(sceptre.pairs.pos$gRNA_id == x)]
    
    lim <- max(abs(p$logfc.vst))
    p <- p[order(p$Highlight, decreasing = FALSE),]
    
    ggplot(p, aes(x = logfc.vst, y = -log10(p_value), colour = Highlight, alpha = Fade)) +
      geom_point(size = 1.5) +
      theme_bw() +
      labs(x = paste0("logfc (VST-normalised)\nBonf=", length(which(p$Bonf < 0.05)), ", FDR=", length(which(p$FDR<0.05))), 
           y = paste0(x, "-Suppression\n-log10(P)")) +
      theme(panel.border = invis, panel.grid = invis, axis.line.x = element_line(), legend.position = "none",
            axis.title = element_text(size = 8), axis.text = element_text(size = 6)) +
      geom_hline(yintercept = -log10(0.05 / nrow(p)), linetype = 2, colour = "black") + # bonf
      geom_hline(yintercept = -log10(max(p$p_value[which(p$FDR < 0.05)])), linetype = 2, colour = "black") + # fdr
      geom_vline(xintercept = 0, colour = "black") +
      scale_y_continuous(limits = c(0,16), expand = c(0,0)) +
      scale_x_continuous(limits = c(-lim, lim)) +
      scale_alpha_manual(values = c(1, 0.15)) +
      scale_colour_manual(values = c("black", "firebrick1"))
  })
  
  pdf(file = "Transcriptome-wide - Volcano.pdf", height = 10, width = 5)
  for (j in seq(0, 124, 12)) {
  # print((j+1):min(c(j+12,length(tw))))
    print(plot_grid(plotlist = twPlot.volc[(j+1):min(c(j+12,length(tw)))], nrow = 6, ncol = 2))
  }
  dev.off()
  
  
  
## Do these show signatures of astroyte activation?
  
  
## What about signatures of astrocyte morphology changes?
  
  
## Long genes?
  geneinfo <- read.delim("../../../Data/Whitelists/GeneInfo.txt")
  genelength <- geneinfo$End - geneinfo$Start
  names(genelength) <- geneinfo$Symbol
  
  a <- tw$TOP1
  m <- match(a$gene_id, names(genelength))  
  a$GeneLength <- genelength[m]
  
## GO
  library(gprofiler2)
  
  ## Run
  go <- lapply(tw, function(x) {
    print("1")
    if (length(which(x$FDR < 0.05)) < 5) return("Too few genes")
    g <- gost(query = x$gene_id[which(x$FDR < 0.05)], custom_bg = x$gene_id, significant = TRUE, evcodes = TRUE, organism = "hsapiens", user_threshold = 0.05, correction_method = "fdr")
    g <- as.data.frame(g$result)
    return(g)
  })
  
  save(go, file = "Transcriptome-wide - GO.rda")
    
  #   go$parents <- sapply(go$parents, function(x) paste(x, collapse = "_"))
  # 
  #   
  # ## Save
  #   g <- apply(go, 2, function(x) { # this code makes the dataframe suitable for csv format
  #     
  #     if (class(x) != "character") {
  #       return(x)
  #     } else {
  #       return(gsub(",", "_", x))
  #     }
  #     
  #   })
    
    write.csv(g, file = "GO.csv", quote = FALSE)  
  
## Analyse and plot!
  # top hits
                                                           

  