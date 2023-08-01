## TThis script analyses negative control guides, and the negative control population

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()

setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Neg/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Seurat)
  library(sceptre)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(rcartocolor)

## Load
  source("../../../Scripts/Functions.R")
  load("Sceptre Input Files.rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]

## Data information
  samples <- c(paste0("NHA_", c(1:5, 7:8)))
  sample.colours <- carto_pal(7, "Bold")
  names(sample.colours) <- samples
  targ.pos <- guides$TargetID[which(guides$TargetCat == "Promoter")] %>% unique()
  targ.enh <- guides$TargetID[which(guides$TargetCat == "Enh")] %>% unique()
  targ.neg <- guides$GuideID[which(guides$TargetCat == "Negative")] %>% unique()
  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  


################################################################################################################################ #
## Calculate empirical p-value ----

## Load enh data
  enh <- read.csv("../Enh/Results Summary.csv", row.names = 1)
  
  # make NB p-value from Z-score
  enh$NB.p <- 2 * pnorm(abs(enh$Z), lower.tail = FALSE) 
  
## Load pos data
  load("../Pos/SCEPTRE Output.rda")
  pos <- de.pos
  colnames(pos) <- c("Gene", "Guide", "PairType", "P", "Z")
  pos$NB.p <- 2 * pnorm(abs(pos$Z), lower.tail = FALSE) 
  
## Load neg data
  load("SCEPTRE Output (Neg, Guide-level).rda")
  load("SCEPTRE Output (NegE, Guide-level).rda")
  neg <- list()
  neg$Neg <- de.neg.guidelevel
  neg$Enh <- de.negE
  
  neg$Neg$NB.p <- 2 * pnorm(abs(neg$Neg$z_value), lower.tail = FALSE) 
  neg$Enh$NB.p <- 2 * pnorm(abs(neg$Enh$z_value), lower.tail = FALSE) 
  
## QQplot of negative control data    
  # sceptre
  pdf(file = "QQplots - SCEPTRE.pdf", height = 4, width = 4)
  make_qq_plot(neg$Neg$p_value[sample(1:nrow(neg$Neg), 100000)]) + labs(title = "Neg (100k subsample / 690k total)")
  make_qq_plot(neg$Neg$p_value[sample(1:nrow(neg$Neg), 100000)]) + labs(title = "Neg (100k subsample / 690k total)")
  make_qq_plot(neg$Neg$p_value[sample(1:nrow(neg$Neg), 100000)]) + labs(title = "Neg (100k subsample / 690k total)")
  make_qq_plot(neg$Enh$p_value[sample(1:nrow(neg$Enh), 100000)]) + labs(title = "NegE (100k subsample / 150k total)")
  make_qq_plot(neg$Enh$p_value[sample(1:nrow(neg$Enh), 100000)]) + labs(title = "NegE (100k subsample / 150k total)")
  make_qq_plot(neg$Enh$p_value[sample(1:nrow(neg$Enh), 100000)]) + labs(title = "NegE (100k subsample / 150k total)")
  make_qq_plot(enh$P) + labs(title = "Enh (7k total)")
  dev.off()  
  
  # nb
  pdf(file = "QQplots - NB.pdf", height = 4, width = 4)
  make_qq_plot(neg$Neg$NB.p[sample(1:nrow(neg$Neg), 100000)]) + labs(title = "Neg (100k subsample / 690k total)")
  make_qq_plot(neg$Neg$NB.p[sample(1:nrow(neg$Neg), 100000)]) + labs(title = "Neg (100k subsample / 690k total)")
  make_qq_plot(neg$Neg$NB.p[sample(1:nrow(neg$Neg), 100000)]) + labs(title = "Neg (100k subsample / 690k total)")
  make_qq_plot(neg$Enh$NB.p[sample(1:nrow(neg$Enh), 100000)]) + labs(title = "NegE (100k subsample / 150k total)")
  make_qq_plot(neg$Enh$NB.p[sample(1:nrow(neg$Enh), 100000)]) + labs(title = "NegE (100k subsample / 150k total)")
  make_qq_plot(neg$Enh$NB.p[sample(1:nrow(neg$Enh), 100000)]) + labs(title = "NegE (100k subsample / 150k total)")
  make_qq_plot(enh$NB.p) + labs(title = "Enh")
  dev.off()  
  
## For the enh-gene pairs, compare observed p-value to the background distribution of neg-gene pairs
  
  ## This is how Gasperini et al. (2019) calculated it:
  
    # [(the number of NTCs with a smaller P-value than that test’s raw P-value) + 1] divided by [the total number of NTCs tests + 1].
    # These empirical P-values were Benjamini-Hochberg corrected, and those < 0.1 were kept for 10% empirical FDR sets.
  
  ## Loop
    total.ntc <- lapply(neg, function(x) nrow(x) + 1)
    
    enh$EmpiricalP.Neg <- NA
    enh$EmpiricalP.Enh <- NA
  
    for (j in 1:nrow(enh)) {
      print(j)
      
      # I note an ambiguity in Gasperini et al.'s phrasing, hence I do two versions
      
      # ## Version 1: a futile endeavour, using only the neg-gene pairs for each given gene
      #   g <- bg$gene_id[j] # get gene
      #   lower <- length(which(bg.p[[g]] < s$p[j])) # number of neg-gene pairs more significant 
      #   e <- (lower + 1) / total.ntc1 # p-value
      #   s$EmpiricalP[j] <- e
      
      ## Version 2: more reasonably,using the neg-gene pairs for all genes as a background
        # for the smaller enh pool of 50
        lower <- length(which(neg$Enh$NB.p < enh$NB.p[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Enh # p-value
        enh$EmpiricalP.Enh[j] <- e
        
        # for the larger neg pool of 250
        lower <- length(which(neg$Neg$NB.p < enh$NB.p[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Neg # p-value
        enh$EmpiricalP.Neg[j] <- e
        
    }

    
    enh$EmpiricalFDR.Enh <- p.adjust(enh$EmpiricalP.Enh, method = "fdr")
    enh$EmpiricalFDR.Neg <- p.adjust(enh$EmpiricalP.Neg, method = "fdr")
    enh$NB.FDR <- p.adjust(enh$NB.p, method = "fdr")
  
    
    ## Now positive control
    pos$EmpiricalP.Neg <- NA
    pos$EmpiricalP.Enh <- NA
    
    for (j in 1:nrow(pos)) {
      print(j)
      
      # I note an ambiguity in Gasperini et al.'s phrasing, hence I do two versions

      ## Version 2: more reasonably,using the neg-gene pairs for all genes as a background
        # for the smaller enh pool of 50
        lower <- length(which(neg$Enh$NB.p < pos$NB.p[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Enh # p-value
        pos$EmpiricalP.Enh[j] <- e
        
        # for the larger neg pool of 250
        lower <- length(which(neg$Neg$NB.p < pos$NB.p[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Neg # p-value
        pos$EmpiricalP.Neg[j] <- e
        
      }
    
## Compare empirical p to raw p
    
  p <- enh
  p.thresh <- 0.1

  # p$Col <- "."
  thresh1 <- -log10(max(p$NB.p[which(p$NB.FDR < p.thresh)]))
  thresh2 <- -log10(max(p$EmpiricalP.Enh[which(p$EmpiricalFDR.Enh < p.thresh)]))
  thresh3 <- -log10(max(p$EmpiricalP.Neg[which(p$EmpiricalFDR.Neg < p.thresh)]))
  
  p$Cat <- "None"
  p$Cat[which((p$NB.FDR < p.thresh) & !(p$EmpiricalFDR.Enh) < p.thresh & !(p$EmpiricalFDR.Neg) < p.thresh)] <- "NB"
  p$Cat[which((p$NB.FDR < p.thresh) & (p$EmpiricalFDR.Enh) < p.thresh & !(p$EmpiricalFDR.Neg) < p.thresh)] <- "NB + EEnh"
  p$Cat[which((p$NB.FDR < p.thresh) & (p$EmpiricalFDR.Enh) < p.thresh & (p$EmpiricalFDR.Neg) < p.thresh)] <- "NB + EEnh + ENeg"
  
  p$Cat <- factor(p$Cat, levels = c("None", "NB", "NB + EEnh", "NB + EEnh + ENeg"))
  levels(p$Cat) <- paste0(levels(p$Cat), " (n=", table(p$Cat), ")")
  
  
   pA <- ggplot(enh[which(enh$NB.p > 1e-10),], aes(x = -log10(NB.p), y = -log10(EmpiricalP.Enh))) +
    geom_point() +
    theme_bw() + 
    labs(y = "Empirical-Enh NB -log10P", x = "Raw NB -log10P")  +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_hline(yintercept = thresh2, linetype = 2) +
    geom_vline(xintercept = thresh1, linetype = 2) +
      scale_y_continuous(limits = c(0,6)) +
    theme(legend.position = "none", panel.grid.minor = invis, panel.border = invis, axis.line.y = element_line())
  
  pB <- ggplot(enh[which(enh$NB.p < 1e-10),], aes(x = -log10(NB.p), y = -log10(EmpiricalP.Enh))) +
    geom_point() +
    theme_bw() + 
    labs(y = "Empirical-Enh NB -log10P", x = "Raw NB -log10P")  +
    geom_hline(yintercept = thresh2, linetype = 2) +
    scale_y_continuous(limits = c(0,6)) +
    theme(legend.position = "none", panel.grid.minor = invis, panel.border = invis, axis.title.y = invis,
          axis.line.y = invis, axis.text.y = invis, axis.ticks.y = invis) 
  
  pC <- ggplot(enh[which(enh$NB.p > 1e-10),], aes(x = -log10(NB.p), y = -log10(EmpiricalP.Neg))) +
    geom_point() +
    theme_bw() + 
    labs(y = "Empirical-Neg NB -log10P", x = "Raw NB -log10P")  +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_hline(yintercept = thresh3, linetype = 2) +
    geom_vline(xintercept = thresh1, linetype = 2) +
      scale_y_continuous(limits = c(0,6)) +
    theme(legend.position = "none", panel.grid.minor = invis, panel.border = invis, axis.line.y = element_line()) 
  
  pD <- ggplot(enh[which(enh$NB.p < 1e-10),], aes(x = -log10(NB.p), y = -log10(EmpiricalP.Neg))) +
    geom_point() +
    theme_bw() + 
    labs(y = "Empirical-Neg NB -log10P", x = "Raw NB -log10P")  +
    geom_hline(yintercept = thresh3, linetype = 2) +
    scale_y_continuous(limits = c(0,6)) +
    theme(legend.position = "none", panel.grid.minor = invis, panel.border = invis, axis.title.y = invis,
          axis.line.y = invis, axis.text.y = invis, axis.ticks.y = invis) 
  
  pdf(file = "Empirical P vs Raw P.pdf", height = 3, width = 10)
  plot_grid(pA, pB, pC, pD, nrow = 1)
  dev.off()
    
  pdf(file = "Empirical P vs Raw P - Sanity Checks.pdf", height = 3, width = 5)
  ggplot(p, aes(x = Cat, y = Z)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
    labs(y = "Z (from NB/SCEPTRE)") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
 
   ggplot(p, aes(x = Cat, y = Gene.Distance)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
     labs(y = "Enh-gene Distance") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
   dev.off() 
   
## Cat
  p.thresh <- 0.1
  enh$Cat <- "."
  enh$Cat[which(enh$Hit & enh$EmpiricalFDR.Enh < p.thresh & enh$EmpiricalFDR.Neg < p.thresh)] <- "3: All"
  enh$Cat[which((enh$Hit) & !(enh$EmpiricalFDR.Enh) < p.thresh & !(enh$EmpiricalFDR.Neg) < p.thresh)] <- "1: SCEPTRE"
  enh$Cat[which(!(enh$Hit) & (enh$EmpiricalFDR.Enh) < p.thresh & !(enh$EmpiricalFDR.Neg) < p.thresh)] <- "1: EEnh"
  enh$Cat[which(!(enh$Hit) & !(enh$EmpiricalFDR.Enh) < p.thresh & (enh$EmpiricalFDR.Neg) < p.thresh)] <- "1: ENeg"
  enh$Cat[which(!(enh$Hit) & (enh$EmpiricalFDR.Enh) < p.thresh & (enh$EmpiricalFDR.Neg) < p.thresh)] <- "2: EEnh+ENeg"
  enh$Cat[which((enh$Hit) & !(enh$EmpiricalFDR.Enh) < p.thresh & (enh$EmpiricalFDR.Neg) < p.thresh)] <- "2: Sceptre+ENeg"
  enh$Cat[which((enh$Hit) & (enh$EmpiricalFDR.Enh) < p.thresh & !(enh$EmpiricalFDR.Neg) < p.thresh)] <- "2: Sceptre+EEnh"
  enh$Cat[which(!(enh$Hit) & !(enh$EmpiricalFDR.Enh) < p.thresh & !(enh$EmpiricalFDR.Neg) < p.thresh)] <- "0: None"
  
  enh$Cat <- factor(enh$Cat)
  levels(enh$Cat) <- paste0(levels(enh$Cat), " (", table(enh$Cat), ")")
    
## Plots!
  cols <- pal_lancet()(6) %>% rev()
  
  ## Empirical P vs NB
  pA <- ggplot(enh, aes(x = -log10(NB.p), y = -log10(EmpiricalP.Enh), colour = Cat)) +
    geom_point() +
    theme_bw() + 
    labs(y = "Empirical-Enh NB -log10P", x = "Raw NB -log10P")  +
    scale_colour_manual(values = cols) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, colour = "black") +
     scale_y_continuous(limits = c(0, 6))
  
  pB <- ggplot(enh, aes(x = -log10(NB.p), y = -log10(EmpiricalP.Neg), colour = Cat)) +
    geom_point() +
    theme_bw() + 
    labs(y = "Empirical-Neg NB -log10P", x = "Raw NB -log10P") +
    scale_colour_manual(values = cols) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, colour = "black") +
     scale_y_continuous(limits = c(0, 6))
  
  pC <- ggplot(enh, aes(x = -log10(EmpiricalP.Enh), y = -log10(EmpiricalP.Neg), colour = Cat)) +
    geom_point() +
    theme_bw() + 
    labs(y = "Empirical-Neg NB -log10P", x = "Empirical-Enh NB -log10P") +
    scale_colour_manual(values = cols) +
    geom_abline(intercept = 0, slope = 1, colour = "black") +
     scale_y_continuous(limits = c(0, 6)) +
     scale_x_continuous(limits = c(0, 6))
  
  pD <- get_legend(pC)
  
  pdf(file = "Empirical P Scatterplots.pdf", height = 7, width = 7)
  plot_grid(pA, pB, pC + theme(legend.position = "none"), pD, ncol = 2)
  dev.off()
  

  
    
    


  
  ## Compare SCEPTRE to empirical P
  pA <- ggplot(enh, aes(x = -log10(P), y = -log10(EmpiricalP.Enh), colour = Cat)) +
    geom_point() +
    theme_bw() + 
    scale_colour_manual(values = cols) +
    labs(x = "SCEPTRE -log10P", y = "Empirical-Enh NB -log10P")  +
    geom_abline(intercept = 0, slope = 1, colour = "black") +
     scale_y_continuous(limits = c(0, 6))
    
  pB <- ggplot(enh, aes(x = -log10(P), y = -log10(EmpiricalP.Neg), colour = Cat)) +
    geom_point() +
    theme_bw() + 
    scale_colour_manual(values = cols) +
    labs(x = "SCEPTRE -log10P", y = "Empirical-Neg NB -log10P")  +
    geom_abline(intercept = 0, slope = 1, colour = "black") +
     scale_y_continuous(limits = c(0, 6))
    
  pdf(file = "Empirical P vs SCEPTRE.pdf", height = 3.5, width = 7)
  plot_grid(pA + theme(legend.position = "none"), pB, rel_widths = c(1, 1.7))
  dev.off()
  
  
  ## Properties of pairs across categories
  pdf(file = "Overlap Between Algorithms - Sanity Checks.pdf", height = 3.5, width = 5)
   ggplot(enh, aes(x = Cat, y = Z)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
    labs(y = "Z (from NB/SCEPTRE)") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
 
   ggplot(enh, aes(x = Cat, y = Gene.Distance)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
     labs(y = "Enh-gene Distance") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
   
   ggplot(enh, aes(x = Cat, y = Gene.Exp)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
     labs(y = "Gene Exp") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
   dev.off()
   
   
   

  
  
## Positive controls
  thresh <- -log10(0.05 / 125)
  pA <- ggplot(pos, aes(x = -log10(P), y = -log10(EmpiricalP.Enh))) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,17), expand = c(0,0)) +
    labs(x = "SCEPTRE -log10P", y = "Empirical Enh -log10P") +
    geom_vline(xintercept = thresh, colour = "red", linetype = 2) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.2, linetype = 2) +
    geom_hline(yintercept = thresh, colour = "red", linetype = 2)
  
  pB <- ggplot(pos, aes(x = -log10(P), y = -log10(EmpiricalP.Neg))) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,17), expand = c(0,0)) +
    labs(x = "SCEPTRE -log10P", y = "Empirical Neg -log10P") +
    geom_vline(xintercept = thresh, colour = "red", linetype = 2) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.2, linetype = 2) +
    geom_hline(yintercept = thresh, colour = "red", linetype = 2)
  
  pC <- ggplot(pos, aes(x = -log10(NB.p), y = -log10(EmpiricalP.Enh))) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,300), expand = c(0,0)) +
    labs(x = "Raw NB -log10P", y = "Empirical Enh -log10P") +
    geom_vline(xintercept = thresh, colour = "red", linetype = 2) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.2, linetype = 2) +
    geom_hline(yintercept = thresh, colour = "red", linetype = 2)
  
  pD <- ggplot(pos, aes(x = -log10(NB.p), y = -log10(EmpiricalP.Neg))) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,300), expand = c(0,0)) +
    labs(x = "Raw NB -log10P", y = "Empirical Neg -log10P") +
    geom_vline(xintercept = thresh, colour = "red", linetype = 2) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.2, linetype = 2) +
    geom_hline(yintercept = thresh, colour = "red", linetype = 2)
  
  pdf(file = "Empirical P Pos Control.pdf", height = 3, width = 8)
  plot_grid(pA, pB)
  plot_grid(pC, pD)
  dev.off()
  
  
## Simple look at Negative ontrols
  pdf(file = "Negative P - Scatterplot.pdf", height = 3, width = 5)
  xlim <- max(abs(c(neg$Neg$z_value, neg$Enh$z_value)))
  ggplot(neg$Neg, aes(x = z_value, y = -log10(p_value))) + 
    geom_point() +
    scale_x_continuous(limits = c(-xlim, xlim)) +
    labs(x = "Z", y = "SCEPTRE -log10P") +
    theme_bw() +
    geom_vline(xintercept = 0, colour = "red", linetype = 2)
  
  ggplot(neg$Enh, aes(x = z_value, y = -log10(p_value))) + 
    geom_point() +
    scale_x_continuous(limits = c(-xlim, xlim)) +
    labs(x = "Z", y = "SCEPTRE -log10P") +
    theme_bw() +
    geom_vline(xintercept = 0, colour = "red", linetype = 2)
  dev.off()
  
    pdf(file = "Negative P - Distribution.pdf", height = 3, width = 5)
    
  xlim <- max(abs(c(neg$Neg$z_value, neg$Enh$z_value)))
  ggplot(neg$Neg, aes(x = z_value)) + 
    geom_histogram(binwidth = 1) +
    scale_x_continuous(limits = c(-xlim, xlim)) +
    labs(x = "Z", y = "SCEPTRE -log10P") +
    theme_bw() +
    geom_vline(xintercept = 0, colour = "red", linetype = 2)
  
  ggplot(neg$Enh, aes(y = z_value, x = ".")) + 
    geom_point() +
    scale_y_continuous(limits = c(-xlim, xlim)) +
    labs(x = "Z", y = "SCEPTRE -log10P") +
    theme_bw() +
    geom_vline(xintercept = 0, colour = "red", linetype = 2)
  dev.off()
  
  
## Write
  write.csv(p, "Empirical P.csv")
  
  
################################################################################################################################ #
## Guide-level analyses!! ----
  
## Load
  load("../Enh/Guide-level/Guide-level SCEPTRE.rda")
  gl <- de.enh.guidelvl
  colnames(gl) <- c("Gene", "Enh", "Pair", "P", "Z")
  gl$Pair.enh <- paste0(splitter(gl$Enh, "_", 1), "_", gl$Gene)
  gl$Pair.gl <- paste0(splitter(gl$Enh, "_", 1), "_", splitter(gl$Enh, "_", 2), "_", gl$Gene)
  gl <- gl[which(gl$Pair.enh %in% enh$Pair),]
  
  
## Calculate NB p-value from Z-score
  gl$NB.p <- 2 * pnorm(abs(gl$Z), lower.tail = FALSE) 
  
## Run empirical correction
  gl$EmpiricalP.Neg <- NA
  gl$EmpiricalP.Enh <- NA
  
  for (j in 1:nrow(gl)) {
      print(j)
      
      # I note an ambiguity in Gasperini et al.'s phrasing, hence I do two versions
      
      # ## Version 1: a futile endeavour, using only the neg-gene pairs for each given gene
      #   g <- bg$gene_id[j] # get gene
      #   lower <- length(which(bg.p[[g]] < s$p[j])) # number of neg-gene pairs more significant 
      #   e <- (lower + 1) / total.ntc1 # p-value
      #   s$EmpiricalP[j] <- e
      
      ## Version 2: more reasonably,using the neg-gene pairs for all genes as a background
        # for the smaller enh pool of 50
        lower <- length(which(neg$Enh$NB.p < gl$NB.p[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Enh # p-value
        gl$EmpiricalP.Enh[j] <- e
        
        # for the larger neg pool of 250
        lower <- length(which(neg$Neg$NB.p < gl$NB.p[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Neg # p-value
        gl$EmpiricalP.Neg[j] <- e
        
    }

    
    gl$EmpiricalFDR.Enh <- p.adjust(gl$EmpiricalP.Enh, method = "fdr")
    gl$EmpiricalFDR.Neg <- p.adjust(gl$EmpiricalP.Neg, method = "fdr")
    gl$NB.FDR <- p.adjust(gl$NB.p, method = "fdr")
    gl$FDR <- p.adjust(gl$P, method = "fdr")
  
    
## Categorise
  gl$Hit <- gl$FDR < 0.1
  p.thresh <- 0.1
  gl$Cat <- "."
  gl$Cat[which(gl$Hit & gl$EmpiricalFDR.Enh < p.thresh & gl$EmpiricalFDR.Neg < p.thresh)] <- "3: All"
  gl$Cat[which((gl$Hit) & !(gl$EmpiricalFDR.Enh) < p.thresh & !(gl$EmpiricalFDR.Neg) < p.thresh)] <- "1: SCEPTRE"
  gl$Cat[which(!(gl$Hit) & (gl$EmpiricalFDR.Enh) < p.thresh & !(gl$EmpiricalFDR.Neg) < p.thresh)] <- "1: EEnh"
  gl$Cat[which(!(gl$Hit) & !(gl$EmpiricalFDR.Enh) < p.thresh & (gl$EmpiricalFDR.Neg) < p.thresh)] <- "1: ENeg"
  gl$Cat[which(!(gl$Hit) & (gl$EmpiricalFDR.Enh) < p.thresh & (gl$EmpiricalFDR.Neg) < p.thresh)] <- "2: EEnh+ENeg"
  gl$Cat[which((gl$Hit) & !(gl$EmpiricalFDR.Enh) < p.thresh & (gl$EmpiricalFDR.Neg) < p.thresh)] <- "2: Sceptre+ENeg"
  gl$Cat[which((gl$Hit) & (gl$EmpiricalFDR.Enh) < p.thresh & !(gl$EmpiricalFDR.Neg) < p.thresh)] <- "2: Sceptre+EEnh"
  gl$Cat[which(!(gl$Hit) & !(gl$EmpiricalFDR.Enh) < p.thresh & !(gl$EmpiricalFDR.Neg) < p.thresh)] <- "0: None"
  
  gl$Cat <- factor(gl$Cat)
  levels(gl$Cat) <- paste0(levels(gl$Cat), "(", table(gl$Cat), ")")
  
## Sanity
    pdf(file = "Empirical Guide-level Sanity.pdf", height = 3.5, width = 5)
   ggplot(gl, aes(x = Cat, y = Z)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
    labs(y = "Z (from NB/SCEPTRE)") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
 
   m <- match(gl$Pair.enh, res$Pair)
   gl$Gene.Distance <- res$Gene.Distance[m]
   
   ggplot(gl, aes(x = Cat, y = Gene.Distance)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
     labs(y = "Enh-gene Distance") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
   dev.off()
   
## Are these hits in the pooled analysis?
   m <- match(gl$Pair.enh, enh$Pair)
   gl$Pooled.Cat <- enh$Cat[m]
   
   p <- melt(table(gl$Cat, gl$Pooled.Cat))
   
   pdf(file = "Empirical Guide-level Overlap to Pool.pdf", height = 5, width = 5)
   ggplot(p, aes(x = Var1, y = Var2, label = value)) +
     geom_tile(fill = "white", colour = "black") +
     geom_text() +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
     labs(y = "Category of DE Replication Within Pool", x = "Hit Categories at Guide-level")
   dev.off()
   
   
################################################################################################################################ #
## Explore trends in negative  ----
   
  #  load("../../Scratchspace/temp mean exp.rda")
  #  
  #  neg <- lapply(neg, function(x) {
  #    m <- match(x$gene_id, names(mean.exp))
  #    x$gene_cpm <- mean.exp[m]
  #    x <- x[,c(5,7)]
  #    return(x)
  #  })
  #  
  #  p <- melt(neg, id.vars = "gene_cpm")
  # 
  #   p$ExpBin <- cut(p$gene_cpm, c(0,10, 25, 50, 100, 500, 1000000))
  # levels(p$ExpBin) <- c("<10cpm", "10-25cpm", "25-50cpm", "50-100cpm", "100-500cpm", ">500cpm", "NA")
  # 
  # p$L1 <- factor(p$L1)
  # levels(p$L1) <- c("N50", "N250")
  # 
  # pdf(file = "Bias In Z Across Gene Expression Bins.pdf", height = 3, width = 5)
  # ggplot(p, aes(x = L1, y = value, fill = ExpBin)) +
  #   geom_violin(draw_quantiles = 0.5) +
  #   theme_bw() +
  #   labs(x = "Negative Pool", y = "SCEPTRE Z")
  # dev.off()
  #  