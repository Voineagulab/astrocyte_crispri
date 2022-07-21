## TThis script analyses negative control guides, and the negative control population

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/")
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
  source("../../Scripts/Functions.R")
  load("../../Data/Preprocessed/NHA Pooled.rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist.csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]

## Data information
  samples <- c(paste0("NHA_", c(1:5, 7:8)))
  sample.colours <- carto_pal(7, "Bold")
  names(sample.colours) <- samples
  
  pos <- guides$TargetID[which(guides$TargetCat == "Promoter")]
  pos <- unique(pos)
  enh <- guides$TargetID[which(guides$TargetCat == "Enh")]
  enh <- unique(enh)
  

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  


################################################################################################################################ #
## Analysis of the negative pool of cells ----
  
## First, use a basic wilcox test to compare Neg to Enh cells.
  # this test is run on SCT-normalised data
  # pool.de <- FindMarkers(nha, 
  #                        group.by = "TransductionPool",
  #                        ident.1 = "Neg",
  #                        ident.2 = "Enh",
  #                        logfc.threshold = log(1.1), # stringent threshold as high n makes everything significant...
  #                        min.pct = 0,
  #                        genes = unique(s$Gene[which(s$HighExpression)]),
  #                        test.use = "wilcox")
  # pool.de$Gene <- rownames(pool.de)
  # pool.de <- pool.de[which(pool.de$p_val_adj < 0.05),]
  
  write.csv(pool.de, file = "Negative Control - Effect Of Transduction Pool (FC>1.1).csv")
  
  a <- which(nha$TransductionPool == "Neg")
  b <- which(nha$TransductionPool == "Enh")
  
  dat <- as.data.frame(nha@assays$SCT@data)
  dat <- dat[unique(s$Gene[which(s$HighExpression)]),]
  counter <- 0
  y <- apply(x, 1, function(z) {
    counter <<- counter + 1
    print(counter)
    p <- wilcox.test(z[a], z[b])$p.value
    fc <- mean(z[a]) - mean(z[b])
    return(data.frame(p = p, fc = fc))
  })
  
  de.pool <- do.call("rbind", y)
  de.pool$FDR <- p.adjust(de.pool$p, method = "fdr")

  # define significance
  de.pool$Sig <- (de.pool$FDR < 0.05 & abs(de.pool$fc) > log(1.1))

  # ... finally add back gene
  de.pool$Gene <- rownames(de.pool)
  
  # is the gene a hit in the screen?
  de.pool$Hit <- (de.pool$Gene %in% hits$Gene)

## Bias in fold change
  ## Volcano plot
  p <- de.pool
  
  p$Sig <- factor(p$FDR < 0.05 & abs(p$fc) > log(1.1))
  levels(p$Sig) <- c("ns", "FDR<0.05 &\n|FC|>10%")
  
  p$Hit <- factor(p$Gene %in% hits$Gene)
  levels(p$Hit) <- c("Non-hit in Screen", "Hit in Screen")
  
  lim <- max(abs(p$fc))
  
  pdf(file = "Negative Control - Negative Pool Volcano Plot.pdf", height = 3, width = 4)
  ggplot(p, aes(x = fc, y = -log10(p), colour = Sig)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_colour_manual(values = c("black", "firebrick1")) +
    scale_x_continuous(limits = c(-lim,lim)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(y = "-Log10 Unadjusted P", x = "Log Fold-change") +
    theme(panel.border = invis, axis.line.y = element_line())
  dev.off()

  
  pdf(file = "Negative Control - Bias For Upregulation In Negative Pool.pdf", height = 3, width = maxw)
  ggplot(p, aes(x = as.factor(sign(fc)), fill = Sig)) +
    geom_bar() +
    facet_wrap(~Hit, scales = "free_y") +
    theme_bw() +
    scale_fill_manual(values = c("black", "firebrick1")) +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(),
          legend.title = invis) +
    labs(x = "Sign of Fold-change\nNegative vs. Enhancer Pools", y = "Gene Count") +
    scale_y_continuous(expand = c(0,0))
  dev.off()
  
  # this is a very worrying result, given we want to be finding downregulation!
  
## Gene ontology
  go <- gost(query = de.pool$Gene[which(de.pool$Sig)], custom_bg = de.pool$Gene,
             significant = TRUE, evcodes = TRUE, organism = "hsapiens", user_threshold = 0.05, correction_method = "fdr")
  
  go <- go$result
  
  write.csv(go[,3:11], file = "Negative Control - Pool GO.csv")
  
  # repeat, but with a more stringent fc
  go.stringent <- gost(query = de.pool$Gene[which(de.pool$Sig & abs(de.pool$fc) > log(1.2))], custom_bg = de.pool$Gene,
             significant = TRUE, evcodes = TRUE, organism = "hsapiens", user_threshold = 0.05, correction_method = "fdr")
  
  go.stringent <- go.stringent$result
  
  write.csv(go.stringent[,3:11], file = "Negative Control - Pool GO (Stronger Fold-change Subset).csv")
  
  
################################################################################################################################ #
## Effect of MOI on gene expression ----
  
## Here, I use a lm to explore how each gene's expression is influenced by the number of distinct guides it expresses
  
moi.effect <- list()
  
moi <- rowSums(nha@meta.data[,which(colnames(nha@meta.data) %in% guides$GuideID)]) # the grep finds individual guides rather than those pooled by target
log.moi <- log(moi) # recommendation per sceptre
t.dat <- t(dat)  

## Function
  run.moi.lm <- function(subset) {
    # cor
    cor <- apply(dat[,subset], 1, function(x) {
      cor(x, log.moi[subset])
    })
    
    # lm
    mod <- lm(t.dat[subset,] ~ log.moi[subset]) 
    sum <- summary(mod)
    effect <- sapply(sum, function(x) x$coefficients["log.moi[subset]", "Estimate"])
    p <- sapply(sum, function(x) x$coefficients["log.moi[subset]", "Pr(>|t|)"])
    res <- data.frame(cor = cor,
                      coef = effect,
                      p = p,
                      sig = p.adjust(p, method = "fdr") < 0.05)
    
    return(res)
  }

## Apply 
  # to all cells
  moi.effect$All <- run.moi.lm(1:ncol(nha))  
  
  # to negative control pool
  moi.effect$Neg <- run.moi.lm(which(nha$TransductionPool == "Neg"))  
  
  # to enhancer pool
  moi.effect$Enh <- run.moi.lm(which(nha$TransductionPool == "Enh"))  
  
  # bind and save
  moi.effect <- do.call("cbind", moi.effect)
  write.csv(moi.effect, file = "Negative Control - Modelling MOI vs. Log Guides.csv")
  write.csv(moi.effect[which(moi.effect$Hit),], file = "Negative Control - Modelling MOI vs. Log Guides (Hits).csv")
  
  
## Plot
  moi.effect$Hit <- rownames(moi.effect) %in% hits$Gene # then colour by hit
  
  ## Compare in Neg and Enh pools
    pdf(file = "Negative Control - MOI-Expression Correlation Scatterplot.pdf", height = 3, width = 3.5)
    Density <- get_density(moi.effect$Enh.cor, moi.effect$Neg.cor, n = 100)
    qplot(moi.effect$Enh.cor, moi.effect$Neg.cor, colour = Density) +
      scale_colour_viridis_c() +
      theme_bw() +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_vline(xintercept = 0, linetype = 2) +
      labs(x = "Expression to log(MOI) correlation in Enhancer Pool",
           y = "Expression to log(MOI) correlation in Negative Pool",) +
      theme(panel.border = invis, panel.grid = invis)
  dev.off()
  
  ## Number of significant effects
    m <- melt(moi.effect[,c(4,12,13)], id.vars = "Hit")
    levels(m$variable) <- c("MOI associated with expression\nin Enhancer Pool", "MOI associated with expression\nin Negative Pool")
    m$Hit <- factor(m$Hit)
    levels(m$Hit) <- c("Other", "Gene is Screen Hit")
    
    pdf(file = "Negative Control - MOI-Expression Significance Barplot.pdf", height = 3, width = maxw)
    ggplot(m, aes(x = Hit, fill = value)) +
      geom_bar() +
      facet_wrap(~variable, scales = "free_y") +
      theme_bw() +
      scale_fill_manual(values = c("black", "firebrick1")) +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "Number of Genes", x = "Hit In Screen") +
      theme(axis.ticks.x = invis, axis.title.x = invis)
    
    ggplot(m, aes(x = Hit, fill = value)) +
      geom_bar(position = position_fill()) +
      facet_wrap(~variable, scales = "free_y") +
      scale_fill_manual(values = c("black", "firebrick1")) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "Fraction of Genes", x = "Hit In Screen") +
      theme(axis.ticks.x = invis, axis.title.x = invis)
  dev.off()  
  
## Overlap between effect of MOI and that of the negative pool
  # that is, does MOI explain some of the observed differences between pools?
  p <- data.frame(Pool.fc = z$fc, MOI.Enh = moi.effect$Enh.cor, MOI.Neg = moi.effect$Neg.cor, row.names = rownames(moi.effect))
  p$Dens.Enh <- get_density(p$Pool.fc, p$MOI.Enh, n = 100)
  p$Dens.Neg <- get_density(p$Pool.fc, p$MOI.Neg, n = 100)
  # p <- melt(p, id.vars = "Pool.fc")
  
  pdf(file = "Negative Control - Correlation Between Gene MOI and Pool Effects.pdf", height = 4, width = 4.5)
  r <- cor(p$Pool.fc, p$MOI.Enh) %>% round(2)
  ggplot(p, aes(x = MOI.Enh, y = Pool.fc, colour = Dens.Enh)) +
    geom_point() +
    scale_colour_viridis_c(option = "D") +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = 2) +
    theme_bw() +
    theme(panel.border = invis) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    labs(x = paste0("Expression correlation with MOI\nEnhancer Pool\nr=", r), y = "Fold-change in Neg versus Enh Pool")
  
  ggplot(p[unique(hits$Gene),], aes(x = MOI.Enh, y = Pool.fc)) +
    geom_point() +
    # scale_colour_viridis_c(option = "D") +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = 2) +
    theme_bw() +
    theme(panel.border = invis) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    labs(x = paste0("Expression correlation with MOI\nEnhancer Pool\nr=", r), y = "Fold-change in Neg versus Enh Pool", title = "Hits Only")
  
  r <- cor(p$Pool.fc, p$MOI.Neg) %>% round(2)
  ggplot(p, aes(x = MOI.Neg, y = Pool.fc, colour = Dens.Neg)) +
    geom_point() +
    scale_colour_viridis_c(option = "D") +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = 2) +
    theme_bw() +
    theme(panel.border = invis) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) + 
    labs(x = paste0("Expression correlation with MOI\nNegative Pool\nr=", r), y = "Fold-change in Neg versus Enh Pool")
  
  dev.off()    
  
## Finally, log guides mapped to tSNE!
  nha$log.MOI <- log.moi

  pdf(file = "Negative Control - tSNE For log(MOI).pdf", height = 3, width = 3.5)
  FeaturePlot(nha, features = "log.MOI") + 
    theme_void() +
    scale_colour_carto_c(palette = "Geyser")
  dev.off()
  