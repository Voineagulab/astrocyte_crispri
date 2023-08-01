## Load
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/1_Processing/")
  # load("../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv")
  load("../../Data/Preprocessed/Guide Expression Matrix.rda")
  load("Temp Metadata.rda")
  # meta <- cbind(meta, nha@reductions$umap@cell.embeddings)
  library(tidyverse)

################################################################################################################################ #
## Test guide filtering ---- 
  
# ## Convert to fraction.
#   pos.guides <- guides$GuideID[which(guides$TargetCat == "Promoter")]
#   thresh <- 3
#   keep.cells <- colnames(guide.exp)[(colSums(guide.exp >= thresh) > 0)]
#   high.moi <- colnames(guide.exp)[(colSums(guide.exp >= thresh) >= 20)]
#   check.fraction <- apply(guide.exp, 2, function(x) {
# 
#     x[which(x < thresh)] <- 0
#     x <- x / sum(x)
#   })
#   check.fraction <- as.data.frame(check.fraction)
#   
  
## Parameter values to evaluate
  thresh <- 3
  # thresh <- 10
  par.dom <- c(0:5) / 100
  par.cpm <- c(0,10, 20, 25, 50)
  
## Expression matrices
  # cpm
  guide.cpm <- list()
  
  for (j in colnames(guide.exp)) {
    print(j)
    x <- guide.exp[,j]
    x[which(x < thresh)] <- 0
    guide.cpm[[j]] <- x / (meta[j,"UMI_Total"] / 10^6)
  }

  guide.cpm <- do.call("cbind", guide.cpm)  %>% as.data.frame()  
  rownames(guide.cpm) <- rownames(guide.exp)
  
  # dominance
  guide.fraction <- apply(guide.exp, 2, function(x) {
    
    x[which(x < thresh)] <- 0
    x <- x / sum(x)
  })
  guide.fraction <- as.data.frame(guide.fraction)
  
## Plot 1: MOI distribution as a function of CPM and dominance thresholds
  

  # note that dom is applied before cpm, in this case...
  
  # calculate the number of cells passing criteria, number of assignments, and MOI distribution for each parameter value   
  
  gRaw <- list()
  for (j in par.dom) {
    
    pd <- guide.fraction > j
    
    for (k in par.cpm){
      run.name <- paste0("Dom", j, "_CPM", k)
      print(run.name)
      
      pc <- guide.cpm > k
      
      x <- list()
      for (l in rownames(meta)) {
        x[[l]] <- length(which(pd[,l] & pc[,l]))
      }
      
      gc()
      y <- do.call("c", x)
      
      z <- data.frame(ThreshDom = j,
                      ThreshCPM = k,
                      CellHasGuide = length(which(y > 0)),
                      TotalAssignments = sum(y),
                      MeanMOI = mean(y[which(y > 0)]), 
                      MaxMOI = max(y))
      
      gRaw[[run.name]] <- list(MOI = y, Summary = z)
     
    }
  }
  
  

  gSums <- lapply(gRaw, function(q) q$Summary)  
  gSums <- do.call("rbind", gSums)

  save(gRaw, gSums, file = "Guide Expression Thresholding (UMI=10).rda")  
  
 
################################################################################################################################ #
## Plot MOI ---- 
  

## Basic MOI distribution
  p <- sapply(gRaw, function(q) q$MOI)
  p <- melt(p)
  p$DomThresh <- splitter(p$Var2, "_", 1) %>% gsub("Dom", "", .)
  p$CPMThresh <- splitter(p$Var2, "_", 2) %>% gsub("CPM", "", .)
  # p <- p[,3:5]
  p$Bin <- cut(p$value, c(0:15, 17, 19, 21, 25, 35, 50, 100, 1000), right = FALSE)
  levels(p$Bin) <- c(0:14, "15-16", "17-18", "19-20", "21-24", "25-34", "35-49", "50-99", "100+")
  bin.colour <- c("black", rep("cornflowerblue", 10), rep("goldenrod1", 7), rep("firebrick1", 5))
  
  
  pdf(file = "Final/Filtering Guide - MOI Distribution (UMI=10).pdf", height = 10, width = 14)
  ggplot(p, aes(x = Bin, fill = Bin)) +
    geom_bar() +
    facet_grid(DomThresh~CPMThresh) +
    theme_bw() +
    scale_fill_manual(values = bin.colour) +
    theme(panel.grid = invis, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") +
    labs(y = "Left: Number of Cells in Bin\nRight: Dominance Threshold", x = "Top: CPM Threshold\nBottom: MOI Bin")
  dev.off()  
  

## Comparison to expected Poisson
  # tabulate the reality
  q <- melt(table(p$Bin, p$DomThresh, p$CPMThresh))
  colnames(q) <- c("Bin", "DomThresh", "CPMThresh", "value")
  
  # calculate expected
  xpect <- (dpois(0:1000, 5) * nrow(meta)) %>% data.frame()
  xpect$Bin <- cut(as.numeric(rownames(xpect)) - 1, c(0:15, 17, 19, 21, 25, 35, 50, 100, 1000), right = FALSE)
  levels(xpect$Bin) <- levels(p$Bin)
  xpect <- aggregate(xpect$. ~ xpect$Bin, FUN = sum)
  
  # match expected to reality
  m <- match(q$Bin, xpect$`xpect$Bin`)
  q$Poisson <- xpect$`xpect$.`[m]
  
  # plot
  pdf(file = "Final/Filtering Guide - Versus Expected (UMI=10).pdf", height = 10, width = 14)
  ggplot(q, aes(x = Poisson+1, y = value+1, label = Bin )) +
    geom_text(size = 2) +
    facet_grid(DomThresh~CPMThresh) +
    labs(x = "Expected Number of Cells (+1) Per Bin at Poisson MOI = 5", y = "Observed Number of Cells (+1) Per Bin") +
    scale_y_continuous(trans = "log2") +
    scale_x_continuous(trans = "log2") +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1) +
    theme()
  
  ggplot(q, aes(x = Poisson+1, y = value+1, label = Bin )) +
    geom_text(size = 2) +
    facet_grid(DomThresh~CPMThresh) +
    labs(x = "Expected Number of Cells (+1) Per Bin at Poisson MOI = 5", y = "Observed Number of Cells (+1) Per Bin") +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1) +
    theme()
  
  dev.off()
  
  
## Repeat above but for Enh Pool only 
  # tabulate the reality
  q <- p[which(p$Var1 %in% rownames(meta)[which(meta$TransductionPool == "Primary")]),]
  q <- melt(table(q$Bin, q$DomThresh, q$CPMThresh))
  colnames(q) <- c("Bin", "DomThresh", "CPMThresh", "value")
  
  # calculate expected
  xpect <- (dpois(0:1000, 5) * length(which(meta$TransductionPool == "Primary"))) %>% data.frame()
  xpect$Bin <- cut(as.numeric(rownames(xpect)) - 1, c(0:15, 17, 19, 21, 25, 35, 50, 100, 1000), right = FALSE)
  levels(xpect$Bin) <- levels(p$Bin)
  xpect <- aggregate(xpect$. ~ xpect$Bin, FUN = sum)
  
  # match expected to reality
  m <- match(q$Bin, xpect$`xpect$Bin`)
  q$Poisson <- xpect$`xpect$.`[m]
  
  # plot
  pdf(file = "Final/Filtering Guide - Versus Expected (Enh Pool Only) (UMI=10).pdf", height = 10, width = 14)
  ggplot(q, aes(x = Poisson+1, y = value+1, label = Bin )) +
    geom_text(size = 2) +
    facet_grid(DomThresh~CPMThresh) +
    labs(x = "Expected Number of Cells (+1) Per Bin at Poisson MOI = 5", y = "Observed Number of Cells (+1) Per Bin") +
    scale_y_continuous(trans = "log2") +
    scale_x_continuous(trans = "log2") +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1) +
    theme()
  
  ggplot(q, aes(x = Poisson+1, y = value+1, label = Bin )) +
    geom_text(size = 2) +
    facet_grid(DomThresh~CPMThresh) +
    labs(x = "Expected Number of Cells (+1) Per Bin at Poisson MOI = 5", y = "Observed Number of Cells (+1) Per Bin") +
    # scale_y_continuous(trans = "log2") +
    # scale_x_continuous(trans = "log2") +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1) +
    theme()
  
  dev.off()
  
  
  
## Plot the summary information
  pdf(file = "Final/Filtering Guide - Mean vs Max MOI (UMI=10).pdf", height = 4, width = 5.5)
  ggplot(gSums, aes(x = MeanMOI, y = MaxMOI, colour = as.factor(ThreshDom), group = as.factor(ThreshDom), shape = as.factor(ThreshCPM))) +
    geom_point() +
    geom_line() +
    scale_y_continuous(trans = "log2", limits = c(0.5,1000), breaks = c(1,8,16,32,128,1024)) +
    theme_bw() +
    labs(x = "Mean MOI", y = "Max MOI") +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line())
  dev.off()  
  
  
## Plot mean versus variance of MOI...
  x <- lapply(gRaw, function(x) {
    x <- x$MOI
    zero <- which(x == 0)
    
    data.frame(Mean = mean(x),
               Var = var(x),
               MeanNoZero = mean(x[-zero]),
               VarNoZero = var(x[-zero]))
  })
  
  x <- do.call("rbind", x)
  x$DomThresh <- splitter(rownames(x), "_", 1) %>% gsub("Dom", "", .)
  x$CPMThresh <- splitter(rownames(x), "_", 2) %>% gsub("CPM", "", .)
  
  pA <- ggplot(x, aes(x = Mean, y = Var, colour = as.factor(DomThresh), group = as.factor(DomThresh), shape = as.factor(CPMThresh))) +
    geom_point() +
    geom_line() +
    # scale_y_continuous(trans = "log2", limits = c(0.5,1000), breaks = c(1,8,16,32,128,1024)) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = "Mean", y = "Var") +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line(), legend.position = "none")
  
  pB <- ggplot(x, aes(x = Mean, y = Var, colour = as.factor(DomThresh), group = as.factor(DomThresh), shape = as.factor(CPMThresh))) +
    geom_point() +
    geom_line() +
    scale_y_continuous(trans = "log2", limits = c(0.5,1000), breaks = c(1,8,16,32,128,1024)) +
    # scale_x_continuous(trans = "log2") +
    theme_bw() +
    # geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = "Mean", y = "Var") +
    guides(colour = guide_legend(title = "Dominance"), shape = guide_legend(title = "CPM")) +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line())
  
  pC <- ggplot(x, aes(x = MeanNoZero, y = VarNoZero, colour = as.factor(DomThresh), group = as.factor(DomThresh), shape = as.factor(CPMThresh))) +
    geom_point() +
    geom_line() +
    # scale_y_continuous(trans = "log2", limits = c(0.5,1000), breaks = c(1,8,16,32,128,1024)) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = "Mean (Excluding Zeroes)", y = "Var (Excluding Zeroes)") +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line(), legend.position = "none")
  
  pD <- ggplot(x, aes(x = MeanNoZero, y = VarNoZero, colour = as.factor(DomThresh), group = as.factor(DomThresh), shape = as.factor(CPMThresh))) +
    geom_point() +
    geom_line() +
    scale_y_continuous(trans = "log2", limits = c(0.5,1000), breaks = c(1,8,16,32,128,1024)) +
    # scale_x_continuous(trans = "log2") +
    theme_bw() +
    # geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = "Mean (Excluding Zeroes)", y = "Var (Excluding Zeroes)") +
    guides(colour = guide_legend(title = "Dominance"), shape = guide_legend(title = "CPM")) +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line())
  
  pdf(file = "Final/Filtering Guide - Mean vs Var.pdf", height = 3, width = 7)
  plot_grid(pA, pB, rel_widths = c(1, 1.5))
  plot_grid(pC, pD, rel_widths = c(1, 1.5))
  dev.off()

  
   
################################################################################################################################ #
## Plot suppression ---- 
   
  keep.cells <- colnames(guide.exp)[(colSums(guide.exp >= thresh) > 0)]
  pos.guides <- guides$GuideID[which(guides$TargetCat == "Promoter")]
  
  pdf(file = "Final/Filtering Guide - Suppression Scatterplots.pdf", height = 3, width = 6)
    for (j in pos.guides) {
      print(j)
      
      gene <- guides$TargetID[which(guides$GuideID == j)]
      
      if (!(gene %in% rownames(nha))) next
      if (gene == "EIF3A") next
      
      # basic dataframe with binary guide assignment
      p <- data.frame(Exp = nha@assays$VST@data[gene, keep.cells],
                      GuideUMI = t(guide.exp[j,keep.cells] >= thresh),
                      GuideCPM = t(guide.cpm[j,keep.cells]),
                      GuideFraction = t(guide.fraction[j,keep.cells]),
                      row.names = keep.cells)
      
      colnames(p) <- c("Exp", "GuideUMI", "GuideCPM", "GuideFraction")
      
      # calculate stats
      s <- summary(p$Exp[-which(p$GuideUMI)]) # stats of bg cells
      p <- p[which(p$GuideUMI),]
      
      # facilitate colouring
      p$BinFraction <- cut(p$GuideFraction, c(0, 0.02, 0.05, 1.01), right = FALSE)
      p$CPMFraction <- cut(p$GuideFraction, c(0, 25, 50, 10000), right = FALSE)
      levels(p$BinFraction) <- levels(p$CPMFraction) <- c("Low", "Medium", "High")
      cols <- c()
      
      
      # plot fraction
      r <- cor(p$Exp, p$GuideFraction) %>% round(2)
      points <- data.frame(x = c(0.01, 0.03, 0.25) * 100,
                           y = c(mean(p$Exp[which(p$GuideFraction < 0.02)]),
                                 mean(p$Exp[which(p$GuideFraction < 0.05 & p$GuideFraction > 0.02)]),
                                 mean(p$Exp[which(p$GuideFraction > 0.05)])))
      pA <- ggplot(p, aes(x = GuideFraction*100, y = Exp)) +
        geom_point(alpha = 0.5) +
        geom_point(mapping = aes(x = x, y = y), data = points, shape = 19, colour = "darkorange1", size = 5) +
        scale_x_continuous(trans = pseudo_log_trans(), breaks = c(0, 2, 5, 10, 25, 50, 100), expand = c(0,0), limits = c(0, 102)) +
        geom_smooth(se = FALSE, colour = "cornflowerblue") +
        geom_hline(yintercept = s[c(2,3,5)], linetype = c(2,1,2), colour = "red") +
        geom_vline(xintercept = c( 2,5), linetype = c(1), colour = "darkorange1") +
        labs(y = paste0(gene, " Expression\n", j), x = paste0("Dominance\nr=", r)) +
        theme_bw() +
        theme(panel.grid = invis) 
        
        
      # plot cpm
      r <- cor(p$Exp, p$GuideCPM) %>% round(2)
      points <- data.frame(x = c(20, 35, 100),
                           y = c(mean(p$Exp[which(p$GuideCPM < 25)]),
                                 mean(p$Exp[which(p$GuideCPM < 50 & p$GuideCPM > 25)]),
                                 mean(p$Exp[which(p$GuideCPM > 50)])))
      pB <- ggplot(p, aes(x = GuideCPM+1, y = Exp)) +
        geom_point(alpha = 0.5) +
        geom_point(mapping = aes(x = x, y = y), data = points, shape = 19, colour = "darkorange1", size = 5) +
        scale_x_continuous(trans = "log2", breaks = c(0, 10, 25, 50, 100, 1000)) +
        geom_smooth(se = FALSE, colour = "cornflowerblue") +
        labs(y = paste0(gene, " Expression"), x = paste0("CPM (+1 offset)\nr=", r)) +
        theme_bw() +
        theme(axis.title.y = invis, panel.grid = invis) +
        geom_vline(xintercept = c(25,50), linetype = c(1), colour = "darkorange1") +
        geom_hline(yintercept = s[c(2,3,5)], linetype = c(2,1,2), colour = "red")
      
      print(plot_grid(pA, pB), nrow = 1, rel_widths = c(1,0.8))
      
    }
  
    dev.off()
  
  # For Monday
    # Suppression
      # scatterplot gominance and cpm versus pos expression (MAYBE  UMI?)
    # MOI
      # assess fit to poisson, using scatterplots across MOI; mean = variance; and the package Irina linked
    # SCEPTRE
      # run when removing 99th and damaged cell
    

  
    
    
    
  