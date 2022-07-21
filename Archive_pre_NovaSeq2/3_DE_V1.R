## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_DE_V1/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(rcartocolor)

## Load
  source("../../Scripts/Functions.R")
  load("../../Data/Preprocessed/NHA Pooled.rda")
  load("../../Data/Preprocessed/SY5Y Pooled.rda")
  guides <- read.csv("../../Data/Protospacer Whitelist.csv", row.names = 1)

## Data information
  samples <- c(paste0("NHA_", 1:8),
               paste0("SY5Y_", 1:10))
  
  ct <- c("NHA", "SY5Y")
  
  ct.col <- carto_pal(n = 2, name = "Geyser")
  names(ct.col) <- ct
  
  pos <- guides$TargetID[which(guides$TargetCat == "Promoter")]
  pos <- unique(pos)
  
  de <- list()
  de$sy <- de$nha <- list()

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()


  
################################################################################################################################ #
## NHA ----
  
  # for (j in pos) {
  #   
  #   
  #   x <- FindMarkers(nha, 
  #                    group.by = j,
  #                    ident.1 = "TRUE",
  #                    logfc.threshold = 0.25,
  #                    min.pct = 0.1,
  #                    test.use = "wilcox")
  # }
  
## DE
  de$nha$pos <- lapply(pos, function(x) {
    
    counter <- paste0(x, " (", grep(x, pos)," of ", length(pos), ")")
    print(counter)
    
    if (length(which(nha@meta.data[,x])) < 5) return("< 5 cells")

    FindMarkers(nha, 
                group.by = x,
                ident.1 = "TRUE",
                logfc.threshold = 0,
                features = x,
                min.pct = 0,
                test.use = "wilcox")
  })
  
  names(de$nha$pos) <- pos


## Evaluate for Enhancers
  enh <- colnames(nha@meta.data)[grep("Enh", colnames(nha@meta.data))]

  de$nha$enh <- lapply(enh, function(x) {
    print(x)
    w <- which(ps.whitelist$TargetID[which(ps.whitelist$Celltype == "NHA")] == x)
    w <- w[1]
    w <- ps.whitelist$GuideCoord[w]
    
    nearby <- find.nearby.tss(query = w, expand.by = 10^6, reduced.output = TRUE)
    nearby <- intersect(nearby, rownames(nha))
        
    FindMarkers(nha, 
                group.by = x,
                ident.1 = "TRUE",
                logfc.threshold = 0,
                features = nearby,
                min.pct = 0,
                test.use = "wilcox")
  })
  
  names(de$nha$enh) <- enh


################################################################################################################################ #
## SY5Y ----
  

## DE
  de$sy$pos <- lapply(pos, function(x) {
        if (length(which(sy@meta.data[,x])) < 5) return("< 5 cells")

    
    counter <- paste0(x, " (", which(pos == x)," of ", length(pos), ")")
    print(counter)
    FindMarkers(sy, 
                group.by = x,
                ident.1 = "TRUE",
                logfc.threshold = 0,
                features = x,
                min.pct = 0,
                test.use = "wilcox")
  })
  
  names(de$sy$pos) <- pos



## Evaluate for Enhancers
  enh <- colnames(sy@meta.data)[grep("Enh", colnames(sy@meta.data))]

  de$sy$enh <- lapply(enh, function(x) {
    print(x)
    
    if (length(which(sy@meta.data[,x])) < 5) return("No cells")
    
    y <- ps.whitelist[which(ps.whitelist$Celltype == "SY5Y"),]
    w <- which(y$TargetID == x)
    w <- w[1]
    w <- y$GuideCoord[w]
    
    nearby <- find.nearby.tss(query = w, expand.by = 10^6, reduced.output = TRUE)
    nearby <- intersect(nearby, rownames(sy))
    
    
        
    FindMarkers(sy, 
                group.by = x,
                ident.1 = "TRUE",
                logfc.threshold = 0,
                features = nearby,
                min.pct = 0,
                test.use = "wilcox")
  })
  
  names(de$sy$enh) <- enh
  
  
## Evaluate for negative controls
  
################################################################################################################################ #
## Analyse positive controls ----
  
## Combine
  res <- lapply(de, function(x) do.call("rbind", x$pos))
  res <- do.call("cbind", res)
  res["EIF3A",] <- NA
  res <- as.data.frame(apply(res, 2, as.numeric))
  rownames(res) <- pos
  
  res$nha.logP <- -log10(res$nha.p_val)
  res$sy.logP <- -log10(res$sy.p_val)
  res$nha.sig <- res$nha.p_val < 0.05
  res$sy.sig <- res$sy.p_val < 0.05
  write.csv(res, file = "Positive Control Dataframe.csv")    
  
## Volcano
  pdf(file = "Positive Control Volcano Plots.pdf", height = 3, width = 6)
  ggplot(res, aes(x = nha.avg_log2FC, y = nha.logP, colour = nha.sig)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_colour_manual(values = c("black", "firebrick1")) +
    scale_x_continuous(limits = c(-max(abs(res$nha.avg_log2FC), na.rm = TRUE), max(abs(res$nha.avg_log2FC), na.rm = TRUE)),
                       breaks = seq(-2, 2, 0.5)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(y = "-Log10 Unadjusted P", x = "Log2 Fold-change", title = "NHA") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none")
  
  ggplot(res, aes(x = sy.avg_log2FC, y = sy.logP, colour = sy.sig)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_colour_manual(values = c("black", "firebrick1")) +
    scale_x_continuous(limits = c(-max(abs(res$sy.avg_log2FC), na.rm = TRUE), max(abs(res$sy.avg_log2FC), na.rm = TRUE)),
                       breaks = seq(-2, 2, 0.5)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(y = "-Log10 Unadjusted P", x = "Log2 Fold-change", title = "SY5Y") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none")
  dev.off()

  
  
## Violin plot
  pdf(file = "Positive Control Violins.pdf", height = 2.5, width = maxw)
  for (j in pos) {
    print(j)
    p <- data.frame(Exp = sy@assays$RNA@data[j,],
                    Guide = sy@meta.data[,j])
    
    p$Guide <- factor(p$Guide)
    tab <- (table(p$Guide))
    levels(p$Guide) <- paste0(levels(p$Guide), " (n=", tab, ")")
    
    pval <- de$sy$pos[[j]]
    if (j %in% rownames(pval)) {
      pval <- signif(pval[j,"p_val"], 2)
    } else {
      pval <- "ns"
    }
    
    pA <- ggplot(p, aes(y = Exp, x = Guide, fill = Guide, shape = Guide)) +
            geom_violin(draw_quantiles = 0.5, colour = "white") +
            geom_jitter(width = 0.2, alpha = 0.2) +
            theme_bw() +
            scale_y_continuous() +
            scale_shape_manual(values = c(NA, 1)) +
            labs(y = "SY5Y Normalised Expression", x = paste0("SY5Y cell has control guide\n(p=", pval, ")")) +
            scale_fill_manual(values = c("black", "darkorange1")) +
            annotate("text", x = 1.5, y = max(p$Exp) * 0.95, label = j, size = 4) +
            theme(panel.border = invis, axis.line.y = element_line(),
                  legend.position = "none", axis.ticks.x = invis, panel.grid.major.x = invis)
    
    p <- data.frame(Exp = nha@assays$RNA@data[j,],
                    Guide = nha@meta.data[,j])
    
    p$Guide <- factor(p$Guide)
    tab <- (table(p$Guide))
    levels(p$Guide) <- paste0(levels(p$Guide), " (n=", tab, ")")
    
    pval <- de$nha$pos[[j]]
    if (j %in% rownames(pval)) {
      pval <- signif(pval[j,"p_val"], 2)
    } else {
      pval <- "ns"
    }
    
    
    pB <- (ggplot(p, aes(y = Exp, x = Guide, fill = Guide, shape = Guide)) +
      geom_violin(draw_quantiles = 0.5, colour = "white") +
      geom_jitter(width = 0.2, alpha = 0.2) +
      theme_bw() +
      scale_y_continuous() +
      scale_shape_manual(values = c(NA, 1)) +
      labs(y = "NHA Normalised Expression", x = paste0("NHA cell has control guide\n(p=", pval, ")")) +
      scale_fill_manual(values = c("black", "darkorange1")) +
      annotate("text", x = 1.5, y = max(p$Exp) * 0.95, label = j, size = 4) +
      theme(panel.border = invis, axis.line.y = element_line(),
            legend.position = "none", axis.ticks.x = invis, panel.grid.major.x = invis))
    
    print(plot_grid(pB,pA))
  }
  
  dev.off()
  
## Comparison of repression in NHA and SY
  pdf(file = "Positive Control Relative Expression.pdf", height = 4, width = 4)
  # p <- melt(2^res[,c(2,7)])
  # levels(p$variable) <- c("NHA", "SY5Y")
  p <- melt(res[,c(1,6, 15)])
  p$value <- -log10(p$value)
  levels(p$variable) <- c("NHA", "SY5Y")
  
  ggplot(p, aes(x = CoreEssential, fill = variable, y = value)) +
    geom_violin(scale = "width", draw_quantiles = c(0.5), colour = "white") +
    # geom_boxplot() +
    facet_wrap(~variable) +
    geom_jitter(width = 0.1) +
    # scale_y_continuous(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_manual(values = ct.col) +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = 2) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(y = "Relative Expression of Positive Control Gene") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
          panel.grid.major.x = invis)
  
  p <- melt(res[,c(2,7, 16)])
  p$value <- 2^p$value
  levels(p$variable) <- c("NHA", "SY5Y")
  
  ggplot(p, aes(x = DepMap, fill = variable, y = value)) +
    geom_violin(scale = "width", draw_quantiles = c(0.5), colour = "white") +
    # geom_boxplot() +
    facet_wrap(~variable) +
    geom_jitter(width = 0.1) +
    scale_y_continuous(limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_manual(values = ct.col) +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = 2) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(y = "Relative Expression of Positive Control Gene") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
          panel.grid.major.x = invis)
  dev.off()
  
## Comparison of fold-change in NHA and SY5Y
  p <- res
  p$SignificantIn <- "."
  a <- p$nha.p_val < 0.05
  b <- p$sy.p_val < 0.05
  p$SignificantIn[which(a & b)] <- "Both"
  p$SignificantIn[which(a & !(b) )] <- "NHA"
  p$SignificantIn[which(b & !(a) )] <- "SY5Y"
  p$SignificantIn[which(!(a | b))] <- "Neither"
  
  pdf(file = "Positive Control Fold-change Scatterplot.pdf", height = 4, width = maxw)
  ggplot(p, aes(x = nha.avg_log2FC, y = sy.avg_log2FC, colour = SignificantIn)) +
    geom_point() +
    theme_bw() +
    theme(panel.border = invis, panel.grid.minor = invis, legend.position = "right") +
    labs(x = "NHA log2fc", y = "SY5Y log2fc") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(limits = c(-1.1, 0.5), expand = c(0,0)) +
    scale_x_continuous(limits = c(-2, 0.5), expand = c(0,0)) 
  dev.off()
  
  
################################################################################################################################ #
## Analyse enhancers ----
  
## Make dataframe
  res <- do.call("rbind", de$nha$enh)
  res$Enh <- splitter(rownames(res), "\\.", 1)  
  res$Gene <- splitter(rownames(res), "\\.", 2)  
  res <- res[,c(6,7,2,1,5,3,4)]

  res.filt <- res[which(res$p_val_adj < 0.05),]    
  
  write.csv(res, file = "Enh NHA DE.csv")
  write.csv(res.filt, file = "Enh NHA DE Filt.csv")
  