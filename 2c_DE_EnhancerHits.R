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
  load("SCEPTRE Output.rda")

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()

## Other
  used <- which(nha$AnyGuide) # defines the subset of cells on which DE was calculated
  

## Thresholds
  exp.thresh <- 2^-6
  p.thresh <- 0.1

  
################################################################################################################################ #
## Loading and annotation ----
  
## Load
  load("Sceptre.rda") 


## Create dataframe, based on SCEPTRE calls
  res <- de.enh
  res$Pair <- paste0(res$gRNA_id, "_", res$gene_id)
  res <- res[,c(2,1,4)]
  
  colnames(res) <- c("Enh", "Gene", "P")
  
## Annotate with distance
  annotate.distance <- function(x) {
    # guide coord
    m <- match(x$Enh, guides$TargetID)
    x$Enh.Pos <- guides$TargetCoord[m]
    
    # gene coord
    m <- match(x$Gene, geneInfo$Symbol)
    x$Gene.Pos <- geneInfo$ID[m]
   
    # distance from enh midpoint to tss
    x$Distance <- apply(x, 1, function(y) {
      e.start <- as.numeric(splitter(splitter(y[4], "-", 1), ":", 2))
      e.end <- as.numeric(splitter(y[4], "-", 2))
      e.mid <- round((e.start + e.end) / 2)
      
      g.start <- as.numeric(splitter(splitter(y[5], "-", 1), ":", 2))
      
      dist <- abs(g.start - e.mid)
      return(dist)
    })
    
    # nearest gene
    min.dist <- aggregate(Distance~Enh, x, min, na.rm = TRUE)
    m <- match(x$Enh, min.dist$Enh)
    x$Nearest <- min.dist$Distance[m]
    x$Nearest <- x$Distance == x$Nearest
    
    # bin distance
    x$Distance.Bin <- cut(x$Distance, c(2000,10000,50000,100000,500000, max(x$Distance) + 1)) 
    levels(x$Distance.Bin) <- c("2-10", "10-50", "50-100", "100-500", "500-1000")
    
    # return
    return(x)
    
  }

  res <- annotate.distance(res)
  
## Annotate with mean expression
  # expression normalised to library size via Seurat
  mean <- rowMeans(nha@assays$RNA@data[rownames(nha@assays$RNA@data) %in% res$Gene,used]) # rna assay, not sct
  m <- match (res$Gene, names(mean))
  res$Gene.Exp <- mean[m]
  
  # expression below threshold
  res$High.Expression <- res$Gene.Exp > exp.thresh
  high <- which(res$High.Expression)
  
## Fold-change
  res$logfc.vst <- res$logfc <- NaN

  for (j in 1:nrow(res)) {
    print(j)

    w <- which(sceptre.perturb[res$Enh[j],] == 1)
    
    # fold change in normalised data
    a <- nha@assays$RNA@data[res$Gene[j],w] # use Seurat's normalisation, natural log of ((exp / lib size) * 10000) + 1
    b <- nha@assays$RNA@data[res$Gene[j],-w]
    res$logfc[j] <- mean(b) - mean(a)
    
    # fold change in VST-normalised data
    a <- nha@assays$SCT@data[res$Gene[j],w] # use Seurat's normalisation, natural log of ((exp / lib size) * 10000) + 1
    b <- nha@assays$SCT@data[res$Gene[j],-w]
    res$logfc.vst[j] <- mean(b) - mean(a)

  }
  
## Add sceptre results
  res$Sceptre.P <- sceptre$enh$p_value
  res$Sceptre.Z <- sceptre$enh$z_value
  res$Sceptre.FDR <- p.adjust(res$Sceptre.P, method = "fdr")

## Annotate with negative binomial p    
    n <- nb$enh
    
    n$Enh <- splitter(rownames(n), "\\.", 1)
    n$Gene <- sub("\\.", "_", rownames(n)) %>% splitter("_", 2)
    
    n$Pair <- paste0(n$Enh, "_", n$Gene)
    
    n <- n[,c(7,6,8,2,1)]
    colnames(n) <- c("Enh", "Gene","Pair", "log2fc", "p")

  ## Annotate
    m <- match(res$Pair, n$Pair)
    res$NegBinom.P <- n$p[m]
    res$NegBinom.log2fc <- n$log2fc[m]
    res$NegBinom.FDR <- p.adjust(res$NegBinom.P, method = "fdr")

## Annotate with wilcox p-values
  ## Load  
    w <- wilcox$enh
    
    w$Enh <- splitter(rownames(w), "\\.", 1)
    w$Gene <- sub("\\.", "_", rownames(w)) %>% splitter("_", 2)
    
    w$Pair <- paste0(w$Enh, "_", w$Gene)
    
    w <- w[,c(7,6,8,2,1)]
    colnames(w) <- c("Enh", "Gene","Pair", "log2fc", "p")

  ## Annotate
    m <- match(res$Pair, w$Pair)
    res$Wilcox.P <- w$p[m]
    res$Wilcox.log2fc <- w$log2fc[m]
    res$Wilcox.FDR <- p.adjust(res$Wilcox.P, method = "fdr")



################################################################################################################################ #
## Filter tests, to better define hits ----
  
    
res.all <- res
    
  
## In this section, I explore how differential expression is influenced by mean expression
## I further use this to refine how hit calling is performed
  
  ## First, understand the distribution of expression
    pdf(file = "Enhancers - Expression Threshold Justification.pdf", height = 3, width = maxw)
    ggplot(res, aes(x = log2(Gene.ExpVST + 2^-20))) +
      geom_density(fill = "black", alpha = 0.2) +
      geom_vline(xintercept = log2(exp.thresh + 2^-20)) +
      theme_bw() +
      labs(x = "Log2 of log-normalised expression", y = "Density") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      theme(panel.border = invis, axis.line.y = element_line())
    dev.off()
    
  
  ## FDR as a function of expression  
    pdf(file = "Enhancers - Expression Threshold vs P.pdf", height = 3, width = maxw)
    yline <- max(res$Sceptre.P[which(res$Sceptre.FDR < p.thresh)])
    ggplot(res, aes(x = log2(Gene.ExpVST + 2^-20), y = -log10(Sceptre.P), colour = Sceptre.FDR < p.thresh)) +
      geom_point() +
      geom_vline(xintercept = log2(exp.thresh + 2^-20), linetype = 2, colour = "red") +
      geom_hline(yintercept = -log10(yline), linetype = 2, colour = "red") +
      scale_colour_manual(values = sig.colours) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Log2 of log-normalised expression", y = "-log10 unadjusted p") +
      theme(panel.border = invis, axis.line = element_line())
    dev.off()
    
  ## Tabulate high expression versus hit rate
    p <- table(res$Sceptre.FDR < p.thresh, res$High.Expression)
    p <- p / rowSums(p)
    p <- as.data.frame(p)
    colnames(p) <- c("Hit", "HighExpression", "Freq")
    levels(p$Hit) <- c("ns", "FDR < 0.1")
  
    pdf(file = "Enhancers - Expression Threshold vs Hit Rate.pdf", height = 3, width = 3)
    ggplot(p, aes(x = Hit, y = Freq, fill = HighExpression)) +
         geom_col(colour = "black", width = 0.7) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = c("white", "black")) +
      labs(y = "Fraction of Enhancer-Gene Pairs", x = "Statistical association")
    dev.off()
    
  ## Thus, it is evident that including lowly-expressed genes adds noise and doesn't contribute to differential expression
    
  ## Recalculate hit rate when using highly-expressed genes only
    res <- res[which(res$High.Expression),]
    res$Sceptre.FDR <- p.adjust(res$Sceptre.P, method = "fdr")
    res$NegBinom.FDR <- p.adjust(res$NegBinom.P, method = "fdr")
    res$Wilcox.FDR <- p.adjust(res$Wilcox.P, method = "fdr")
  
# ## Compare new hits and losses
#   # wrangle data
#   p <- res[which(res$Sceptre.FDR < p.thresh | res$ExpressionFilteredFDR < p.thresh),]
#   p$Hit <- p$Sceptre.FDR < p.thresh
#   p$HitV2 <- p$ExpressionFilteredFDR < p.thresh
#   p$HitV2[which(is.na(p$HitV2))] <- FALSE 
#   p$Category <- "."
#   p$Category[which(p$Hit & p$HitV2)] <- "Sig In Both"
#   p$Category[which(p$Hit & !(p$HitV2))] <- "Sig Without Thresh Only"
#   p$Category[which(!(p$Hit) & p$HitV2)] <- "Sig With Thresh Only"
#   
#   # plot
#   p <- melt(p[,c("Category", "Distance", "Z")], id.vars = "Category")
#   
#   pdf(file = "Enhancers - Expression Threshold On Pair Properties.pdf", width = maxw, height = maxw)
#   ggplot(p, aes(x = Category, y = value)) +
#     geom_violin(scale ="width") +
#     geom_jitter(width = 0.2) +
#     facet_wrap(~variable, scales = "free_y", ncol = 1, strip.position = "left") +
#     theme_bw() +
#     geom_hline(yintercept = 0) +
#     theme(axis.title = invis, panel.border = invis)
#   dev.off()
  
  
## In addition to filtering on expression, what about distance?
  tab <- table(res$Sceptre.FDR < p.thresh, res$Distance.Bin)
  
  # it is evident that half of all pairs have a distance of 500-1000MB account, yet only 5/92 hits are there!
  
  res <- res[which(res$Distance < 500000),]
      
## Now: Define hits
  ## Recalculate FDR
    res$Sceptre.FDR <- p.adjust(res$Sceptre.P, method = "fdr")
    res$NegBinom.FDR <- p.adjust(res$NegBinom.P, method = "fdr")
    res$Wilcox.FDR <- p.adjust(res$Wilcox.P, method = "fdr")
    
  ## Define hits
    res$Hit <- res$Sceptre.FDR < p.thresh
    hits <- res[which(res$Hit),]
    hits <- hits[order(hits$Sceptre.FDR),]
    

################################################################################################################################ #
## Final csv ----
  

  write.csv(res.all, file = "Enhancers - All Results Summary.csv")
  write.csv(res, file = "Enhancers - All Results Summary Filtered.csv")
  save(res, res.all, file = "Enhancers - All Results Summary.rda")
    
    
################################################################################################################################ #
## Directionality and fold-change ----

    
## Volcano
  pdf(file = "Enhancers - Volcano.pdf", height = 2.5, width = 3)
  ggplot(res, aes(x = Sceptre.Z, y = -log10(Sceptre.P), colour = Hit))  +
      geom_point() +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line()) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(-27,27)) +
    labs(x = "Z Score", y = "-log10 Unadjusted P") +
    scale_colour_manual(values = sig.colours)
  
  ggplot(res, aes(x = logfc.vst, y = -log10(Sceptre.P), colour = Hit))  +
      geom_point() +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    # scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(-0.8,0.5), expand = c(0,0)) +
    labs(x = "Log Fold-change", y = "-log10 Unadjusted P") +
    scale_colour_manual(values = sig.colours)
  dev.off()
    
  
## Bias in fold-change relative to all genes
  p <- table(res$Hit, sign(res$Sceptre.Z))
  p <- p / rowSums(p)
  p <- as.data.frame(p)
  colnames(p) <- c("Hit", "Sign", "Freq")
  levels(p$Hit) <- c("ns", "FDR < 0.1")
  levels(p$Sign) <- c("Downregulated", "Upregulated")

  pdf(file = "Enhancers - Sign of Fold-Change vs Hit Rate.pdf", height = 3, width = 3)
  ggplot(p, aes(x = Hit, y = Freq, fill = Sign)) +
       geom_col(colour = "black", width = 0.7) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("dodgerblue1", "firebrick1")) +
    labs(y = "Fraction of Enhancer-Gene Pairs", x = "Statistical association")
  dev.off()
  

################################################################################################################################ #
## Distance ----
 
  
## Distribution of distances
  pdf(file = "Enhancers - Distance Distribution.pdf", width = maxw, height = 3)
  
  # density plot, comparing hit and non-hit pairings
  m <- aggregate(res$Distance, list(res$Hit), median)[,2]
  m <- m / 1000
  
  ggplot(res, aes(x = Distance / 1000, colour = Hit, fill = Hit)) +
    geom_density(alpha = 0.1) +
    theme_bw() +
    scale_colour_manual(values = sig.colours) +
    scale_fill_manual(values = sig.colours) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 510), breaks = c(50, 100, 200, 500)) +
    geom_vline(xintercept = m, colour = sig.colours, linetype = 2) +
    labs(x = "Distance Between Enhancer-Gene Pair (kilobases)", y = "Density") +
    theme(panel.border = invis, panel.grid = invis, legend.position = c(0.8, 0.8)) 
  
  # above, for upregulated genes
  m <- aggregate(res$Distance[which(res$Sceptre.Z > 0)], list(res$Hit[which(res$Sceptre.Z > 0)]), median)[,2]
  m <- m / 1000
  
  ggplot(res[which(res$Sceptre.Z > 0),], aes(x = Distance / 1000, colour = Hit, fill = Hit)) +
    geom_density(alpha = 0.1) +
    theme_bw() +
    scale_colour_manual(values = sig.colours) +
    scale_fill_manual(values = sig.colours) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 510), breaks = c(50, 100, 200, 500)) +
    geom_vline(xintercept = m, colour = sig.colours, linetype = 2) +
    labs(x = "Distance Between Enhancer-Gene Pair (kilobases)", y = "Density") +
    theme(panel.border = invis, panel.grid = invis, legend.position = c(0.8, 0.8)) 
  
  # above,for downregulated genes
  m <- aggregate(res$Distance[which(res$Sceptre.Z < 0)], list(res$Hit[which(res$Sceptre.Z < 0)]), median)[,2]
  m <- m / 1000
  
  ggplot(res[which(res$Sceptre.Z < 0),], aes(x = Distance / 1000, colour = Hit, fill = Hit)) +
    geom_density(alpha = 0.1) +
    theme_bw() +
    scale_colour_manual(values = sig.colours) +
    scale_fill_manual(values = sig.colours) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 510), breaks = c(50, 100, 200, 500)) +
    geom_vline(xintercept = m, colour = sig.colours, linetype = 2) +
    labs(x = "Distance Between Enhancer-Gene Pair (kilobases)", y = "Density") +
    theme(panel.border = invis, panel.grid = invis, legend.position = c(0.8, 0.8)) 
  
  # histogram for hit pairings
  ggplot(res[which(res$Hit),], aes(x = Distance.Bin)) +
    geom_bar(colour = "black", fill = "firebrick1", alpha = 0.5, width = 0.7) +
    theme_bw() +
    labs(y = "Count", x = "Distance Between Enhancer-Gene Pair (kilobases)") +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis)
  
  dev.off()
  
  
## Annotation as nearest gene
  p <- table(res$Hit, res$Nearest)
  p <- p / rowSums(p)
  p <- as.data.frame(p)
  colnames(p) <- c("Hit", "NearestGene", "Freq")
  levels(p$Hit) <- c("ns", "FDR < 0.1")
  
  pdf(file = "Enhancers - Nearest Gene.pdf", height = 3, width = 3)
  ggplot(p, aes(x = Hit, y = Freq, fill = NearestGene)) +
    geom_col(colour = "black", width = 0.7) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("white", "black")) +
    labs(y = "Fraction of Enhancer-Gene Pairs", x = "Statistical association")
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
  