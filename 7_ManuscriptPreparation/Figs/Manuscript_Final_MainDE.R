#This script generates figures 2A, 2C, 3A, 3B, 5C, S4G, S4H, S4I, S11 and S12
## Setup
  setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_MainDE/Final/")
  library(tidyverse)
  library(rcartocolor)
  library(ggplot2)
  invis <- element_blank()
  maxh <- 29.7 / 2.54
  maxw <- 21.0 / 2.54
  source("../../FinalFigureFunctions.R")
  source("../../../../FullScale/Scripts/Functions.R")

## Read in results dataframe
  res.final <- read.csv("../../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  
## Read in power dataframe
  pow <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable3_DE/3F_Power.csv")
  wellpowered_egps <- pow$Pair[which(pow$WellPowered015  | pow$Hit)]
  wellpowered_enh <- pow$Enhancer[which(pow$WellPowered015 | pow$Hit)] %>% unique()
  wellpowered_genes <- pow$Gene[which(pow$WellPowered015 | pow$Hit)] %>% unique()
  
## Functions to write to disk in a traceable way
  pdf_mainDE <- function(figNo, title, h, w) {
    pdf(file = paste0("../../Final/", figNo, " - Script MainDE - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_mainDE <- function(figNo, title, toPrint) {
    sink(paste0("../../Final/", figNo, " - Script MainDE - ", title, ".txt"))
    print(toPrint)
    sink()
  }
  
  

################################################################################################################################ #
## Justification of expression threshold ----
 
## Get data
  load("../../../../FullScale/Data/Preprocessed/NHA Pooled (Final).rda") # load
  
  # mean expression
  # x <- meanExp_NHA # this gets loaded in the function script
  x <- data.frame(Symbol = rownames(nha),
                  SeuratNormalised = rowMeans(nha@assays$RNA@data),
                  Pseudobulk = rowSums(nha@assays$RNA@counts),
                  Pseudobulk_CPM = NA,
                  Dropout = rowMeans(nha@assays$RNA@counts > 0))
  x$Pseudobulk_CPM <- x$Pseudobulk * (10^6 / sum(x$Pseudobulk))
  x$Expressed <- x$SeuratNormalised > meanExp_NHA_thresh
  
  # x <- x[-which(x$Mean == 0),] # removes 6 genes with no counts in any cell
  
  # filter to TSS within 500kb of target enhancers
  nby <- unique(res.final$Enh.Pos) %>%
    lapply(find.nearby.tss, expand.by = 500000) %>%
    do.call("c", .) %>%
    unique()
  
  x$Nearby <- x$Symbol %in% nby
  
  
## Plots!
  p <- x[which(x$Nearby),]
  
  # as density
  # pdf(file = "SFig11A - Expression threshold density.pdf", height = 2.5, width = 5.5)
  pdf_mainDE(figNo = "SFig11A", title = "Expression threshold density", h = 2.5, w = 5.5)
  ggplot(p, aes(x = log2(SeuratNormalised))) +
    geom_density(colour = pals$One, fill = pals$One, alpha = 0.2) +
    scale_x_continuous(breaks = c(-20, -15, -10, -8, -6, -4, -2, -1, 0, 1, 2, 3),
                        expand = c(0,0), limits = c(-20, 3)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0.02, 0.04, 0.06, 0.08)) +
    theme_bw() +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line(),
          axis.text.y = text90) +
    labs(x = "log2 Mean Seurat-normalised expression", y = "Distribution density of genes") +
    geom_vline(xintercept = log2(meanExp_NHA_thresh), linetype = 2, colour = "black")
  dev.off()
  
  
  # normalised expression versus CPM
  a <- aggregate(Pseudobulk_CPM~Expressed, data = p, FUN = summary)
  sink_mainDE(figNo = "SFig11C", title = "Expression threshold Seurat vs CPM", a)
  
  # pdf(file = "SFig11C - Expression threshold Seurat vs CPM.pdf", height = 3.5, width = 5.5)
  pdf_mainDE(figNo = "SFig11C", title = "Expression threshold Seurat vs CPM", h = 3.5, w = 5.5)
  ggplot(p, aes(x = log2(SeuratNormalised), y = (Pseudobulk_CPM))) +
    # geom_density(colour = pals$One, fill = pals$One, alpha = 0.2) +
    geom_point(colour = pals$One) +
    scale_x_continuous(breaks = c(-20, -15, -10, -8, -6, -4, -2, -1, 0, 1, 2, 3),
                        expand = c(0,0), limits = c(-20, 3)) +
    scale_y_continuous(expand = c(0,0), trans = "log10", breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000),
                       labels = c("0", "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000")) +
    # scale_colour_manual(values = pals$Hits) +
    theme_bw() +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line(),
          axis.text.y = text90, legend.position = "none") +
    labs(x = "log2 Mean Seurat-normalised expression", y = "CPM in Pseudobulk") +
    geom_vline(xintercept = log2(meanExp_NHA_thresh), linetype = 2, colour = "black")
  dev.off()
  
  # dropout rate
  a <- aggregate(Dropout~Expressed, data = p, FUN = summary)
  sink_mainDE(figNo = "SFig11B", title = "Expression threshold dropout rate", a)
  
  # pdf(file = "SFig11B - Expression threshold dropout rate.pdf", height = 3, width = 2)
  pdf_mainDE(figNo = "SFig11B", title = "Expression threshold dropout rate", h = 3, w = 2)
  ggplot(p, aes(x = Expressed, y = Dropout*100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0,100)) +
    # scale_fill_manual(values = pals$Hits) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(),
          # axis.line.x = element_line(linetype = 2),
          legend.position = "none", axis.text.y = text90) +
    geom_violin(scale = "width", alpha = 0.2, fill = pals$One) +
    labs(x = "Gene above expression\nthreshold", y = "Percent of cells expressing gene")
  dev.off()    
    
################################################################################################################################ #
## n Hits ----
  
h <- res.final[which(res.final$HitPermissive),]  

## Distribution of nGene per enhancer
  tab <- table(h$Enh) %>% table() %>% as.data.frame() # yes, table $Enh not $Gene for this
  colnames(tab) <- c("N", "Freq")
  
  yMax <- max(tab$Freq)

  # pdf(file = "nGenes Per Enhancer.pdf", height = 1.7, width = 1.6)
  pdf_mainDE(figNo = "SFig4G", title = "nGenes Per Enhancer", h = 1.7, w = 1.6)
  ggplot(tab, aes(x = N, y = Freq)) +
    geom_col(fill = pals$One, width = 0.6) +
    # theme_resizeText +
    theme_bw() +
    theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100), limits = c(0, yMax)) +
    labs(y = "Count", x = "Regulated genes\nper functional enh") +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) 
  dev.off()
  
  
  
## Distribution of nEnh per gene
  tab <- table(h$Gene) %>% table() %>% as.data.frame()
  colnames(tab) <- c("N", "Freq")
  
  # pdf(file = "nEnh Per Enhancer.pdf", height = 1.7, width = 1.5)
  pdf_mainDE(figNo = "SFig4G", title = "nEnh Per Gene", h = 1.7, w = 1.5)
  ggplot(tab, aes(x = N, y = Freq)) +
    geom_col(fill = pals$One, width = 0.7) +
    # theme_resizeText +
    theme_bw() +
    theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, yMax), breaks = c(0, 50, 100)) +
    labs(y = "Count", x = "Functional enh\nper regulated gene") +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
          axis.title.y = invis, axis.text.y = invis)
  dev.off()


  
################################################################################################################################ #
## Volcano Plot ----
  

e <- res.final[order(res.final$HitPermissive, decreasing = FALSE),]
e$log2fc.vst <- exp(e$logfc.vst) %>% log2()


# pdf(file = "EGPs - Volcano.pdf", height = 1.8, width = 3.5)
pdf_mainDE(figNo = "2A", title = "Volcano", h = 1.8, w = 3.5)
ggplot(e, aes(x = log2fc.vst, y = -log10(P.N50), colour = HitPermissive))  +
  geom_point() +
  theme_bw() +
  theme(panel.border = invis, panel.grid = invis, axis.line = element_line(), legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  scale_y_continuous(limits = c(0,6), expand = c(0,0), breaks = c(0, 2, 4, 6)) +
  scale_x_continuous(limits = c(-1.35,1.3), expand = c(0,0)) +
  # scale_color_lancet() +
  scale_colour_manual(values = pals$Hits) +
  labs(x = "Log2 fold-change", y = "-log10 empirical P")   

dev.off()

  
################################################################################################################################ #
## Bias in fold-change relative to all genes ----  


e <- res.final[order(res.final$HitPermissive, decreasing = FALSE),]
e$log2fc.vst <- exp(e$logfc.vst) %>% log2()

## Table
  p <- table(e$HitPermissive, factor(sign(e$Z), levels = c("1", "-1"))) # note that the factor call around sign makes it a test for overrepresentation of downregulation rather than upregulation
  
## Write statistics
  f <- fisher.test(p)
  sink_mainDE(figNo = "SFig4H", title = "Downregulation percentage", f)
 
  
## Plot
  p <- table(e$HitPermissive, factor(sign(e$Z), levels = c("1", "-1"))) # note that the factor call around sign makes it a test for overrepresentation of downregulation rather than upregulation
  p <- p / rowSums(p)
  p <- as.data.frame(p)
  colnames(p) <- c("Hit", "Sign", "Freq")
  levels(p$Hit) <- c("ns EGPs", "Functional\nEGPs")
  levels(p$Sign) <- c("Upregulated", "Downregulated")
  p <- p[which(p$Sign == "Downregulated"),]
  
  # pdf(file = "EGPs - Downregulation Percentage.pdf", height = 2, width = 2.2)
  pdf_mainDE(figNo = "SFig4H", title = "Downregulation Percentage", h = 2, w = 2.2)
  ggplot(p, aes(x = Hit, y = Freq*100, fill = Hit)) +
    geom_col(colour = "black", width = 0.6) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis, legend.position = "none",
          axis.title.x = invis) +
    geom_hline(yintercept = 50, linetype = 2, colour = "grey75") +
    scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
    scale_fill_manual(values = pals$Hits) +
    labs(y = "Percent of EGPs with\nnegative fold-change", x = "Hit Category")
  dev.off()
  
## Hmm, how about plotting the distribution rather than binary?
  p <- e
  p$HitPermissive <- factor(p$HitPermissive)
  levels(p$HitPermissive) <- c("ns Pairs", "Hit EGPs")
  p$Sign <- sign(p$log2fc.vst) %>% as.factor()
  levels(p$Sign) <- c("Negative\nfold-change", "Positive\nfold-change")
  
  # pdf(file = "EGPs - Downregulation Distribution.pdf", height = 2.5, width = 3)
  pdf_mainDE(figNo = "SFig4I", title = "Downregulation Distribution", h = 2.5, w = 3)
  ggplot(p, aes(x = HitPermissive, y = log2fc.vst, colour = HitPermissive)) +
    # geom_violin(draw_quantiles = 0.5, position = position_dodge(width = 0.7), scale = "width", width = 0.7) +
    facet_wrap(~Sign, ncol = 2, strip.position = "bottom", scales = "free_x") +
    geom_quasirandom(alpha = 0.7, size = 1, dodge.width = 0.7, width = 0.2) +
        stat_summary(fun = median, geom = "point", shape = "+", colour = "black", size = 5, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = pals$Hits) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = 1) + 
    geom_hline(yintercept = log2(c(1.15, 0.85)), linetype = "dotted", colour = "grey75") + 
    geom_hline(yintercept = log2(c(1.25, 0.75)), linetype = "longdash", colour = "grey75") + 
    scale_y_continuous(expand = c(0,0), limits = c(-1.5,0.6)) +
    labs(y = "Log2 Fold-change of gene") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", panel.grid = invis,
          axis.title.x = invis, axis.text.x = invis, axis.ticks.x = invis)
  
  dev.off()
  

################################################################################################################################ #
## On distance ----
 
  
## Density plot of distances
  p <- res.final
  p <- p[which(p$Pair %in% wellpowered_egps),] # well-powered
  p$HitPermissive <- factor(p$HitPermissive)
  levels(p$HitPermissive) <- c("ns", "Functional")
  
  # density plot, comparing hit and non-hit pairings
  m <- aggregate(p$Gene.Distance, list(p$HitPermissive), median)[,2]
  m <- m / 1000
  names(m) <- levels(p$HitPermissive)
  sink_mainDE(figNo = "3A", title = "Distance Density", m)
  
  # pdf(file = "EGPs - Distance Density.pdf", height = 2, width = 2.8)
  pdf_mainDE(figNo = "3A", title = "Distance Density", h = 2, w = 2.8)
  bandwidth.adjust <- c(1, 0.5, 0.25, 0.1)
  ggplot(p, aes(x = Gene.Distance / 1000, colour = HitPermissive, fill = HitPermissive)) +
    geom_density(alpha = 0.25, adjust = bandwidth.adjust[1]) +
    theme_bw() +
    scale_fill_manual(values = pals$Hits) +
    scale_colour_manual(values = pals$Hits_Darker) +
    guides(fill = guide_legend(title = "EGP"), colour = guide_legend(title = "EGP")) +
    # scale_color_lancet() +
    # scale_fill_lancet() +
    scale_y_continuous(expand = c(0,0), breaks = c(0.002, 0.004)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 510), breaks = c(0,50, 100, 250, 500)) +
    geom_vline(xintercept = m, linetype = 2, colour = pals$Hits_Darker) +
    labs(x = "Distance between EGP (kb)", y = "Density") +
    theme(panel.border = invis, panel.grid = invis, legend.position = c(0.8, 0.8),# legend.title = invis,
          axis.text.y = text90) 
  dev.off()
  
  
## Nearest genes
  ## Read in
    distCat <- read.csv(file = "../../../../FullScale/Results/3_HitEnrichment/EnhGenePairs/Intervening Gene Classification Between EGPs.csv", row.names = 1)
  
   ## Collect data
    distCat <- distCat[which(distCat$Pair %in% wellpowered_egps),] # well powered
  
    # quick test of fisher test for nearest gene
    fish <- table(distCat$HitPermissive, distCat$Distance.Category == "Cond1_Nearest") %>% fisher.test()
    sink_mainDE(figNo = "3B", title = "Nearest gene stacked barplot fisher test", fish)
    
    # now plot all categories
    p <- table(distCat$HitPermissive, distCat$Distance.Category)
    p <- p / rowSums(p)
    sink_mainDE(figNo = "3B", title = "Nearest gene stacked barplot", p)
    p <- t(p)
    p <- as.data.frame.matrix(p)
    
    colnames(p) <- c("ns", "Hit")
    p$Cat <- rownames(p)
    p <- melt(p)
    
    levels(p$variable) <- c("Inactive\nEGPs", "Functional\nEGPs")
  
  p$Cat <- factor(p$Cat)
  levels(p$Cat) <- gsub("Cond", "Type ", levels(p$Cat)) %>%
    sub("_", "\n", .) %>%
    sub("_", " \n", .) %>%
    sub("Distal", "Not nearest", .) %>%
    sub("NonExp", " off TSSs only", .) %>%
    sub("SkipExp", "Skip 1+ on TSS", .) %>%
    sub("NoSkip", "No intervening TSS", .) 
  
  legtitle <- guide_legend(title = "EGP type")
  pal_ord <- c(5,7,4,1)
  
  # pdf(file = "Nearest Gene - Stacked Barplot.pdf", height = 1.8, width = 1.8)
  pdf_mainDE(figNo = "3B", title = "Nearest gene stacked barplot", h = 1.8, w = 1.8)
  ggplot(p, aes(x = variable, fill = Cat, colour = Cat, y = value*100)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = pals$Primary[pal_ord]) +
    scale_colour_manual(values = pals$Primary_Darker[pal_ord]) +
    labs(y = "Percent of EGPs") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid = invis, axis.line.y = element_line(), panel.border = invis, 
          axis.title.x = invis, legend.position = "none") 
  dev.off()

## Compare nearer genes to distance
  # collect data
  p <- distCat[which(distCat$HitPermissive),]
  m <- match(p$Pair, res.final$Pair)
  p$Distance <- res.final$Gene.Distance[m]
  
  p$Cat <- factor(p$Distance.Category)
  levels(p$Cat) <- gsub("Cond", "Type ", levels(p$Cat)) %>%
    sub("_", "\n", .) %>%
    sub("_", " \n", .) %>%
    sub("Distal", "Not nearest", .) %>%
    sub("NonExp", " off TSSs only", .) %>%
    sub("SkipExp", "Skip 1+ on TSS", .) %>%
    sub("NoSkip", "No intervening TSS", .) 
  
  # plot number of nearer genes
  offset <- 1
  # pdf(file = "Nearest Gene - Scatterplot nNearest.pdf", height = 3.5, width = 3.5)
  pdf_mainDE(figNo = "SFig5A", title = "Nearest gene scatterplot", h = 3.5, w = 3.5)
  ggplot(p, aes(x = Distance / 1000, y = nNearer+1, colour = Cat, fill = Cat) ) +
    geom_point() +
    scale_y_continuous(trans = "log2", labels = function(x) { x - offset }, breaks = c(0, 1, 2, 4, 8, 16, 32, 64)+ offset) +
    scale_x_continuous(trans = "log10", limits = c(2, 600), 
                         expand = c(0,0), 
                         breaks = c(2, 10, 50, 200, 500)) +
    scale_colour_manual(values = pals$Primary_Darker[pal_ord]) +
    theme_bw() +
    guides(colour = guide_legend(nrow = 2)) +
    theme(panel.border = invis, axis.line = element_line(), panel.grid = invis,
          legend.position = "none", legend.title = invis) +
    labs(x = "Enhancer-to-gene distance (kb)", y = "Number of genes nearer to functional\nenhancer than target gene")
  
  ggplot(p, aes(x = Distance / 1000, y = nSkip+1, colour = Cat, fill = Cat) ) +
    geom_point() +
    scale_y_continuous(trans = "log2", labels = function(x) { x - offset }, breaks = c(0, 1, 2, 4, 8, 16, 32, 64)+ offset) +
    scale_x_continuous(trans = "log10", limits = c(2, 600), 
                         expand = c(0,0), 
                         breaks = c(2, 10, 50, 200, 500)) +
    scale_colour_manual(values = pals$Primary_Darker[pal_ord]) +
    theme_bw() +
    guides(colour = guide_legend(nrow = 2)) +
    theme(panel.border = invis, axis.line = element_line(), panel.grid = invis,
          legend.position = "none", legend.title = invis) +
    labs(x = "Enhancer-to-gene distance (kb)", y = "Number of genes skipped between\nfunctional enhancer and target gene")
  dev.off()
  
## Distribution of the number of nearer / skipped genes
  p <- melt(p[,c("nNearer", "nSkip", "Cat")])
  p$variable <- factor(p$variable)
  levels(p$variable) <- c("Nearer genes", "Skipped genes")
  
  # pdf(file = "Nearest Gene - Distribution nNearest nSkip.pdf", height = 2.8, width = 3.5)
  pdf_mainDE(figNo = "SFig5B", title = "Nearest gene distributions", h = 2.8, w = 3.5)
  ggplot(p, aes(x = variable, fill = Cat, colour = Cat, y = value+offset)) +
    geom_quasirandom(dodge.width = 0.7) +
      geom_violin(scale = "width", width = 0.7, draw_quantiles = 0.5, alpha = 0) +
    scale_fill_manual(values = pals$Primary[pal_ord]) +
    scale_colour_manual(values = pals$Primary_Darker[pal_ord]) +
    scale_y_continuous(trans = "log2", labels = function(x) { x - offset }, breaks = c(0, 1, 2, 4, 8, 16, 32, 64)+ offset) +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(), legend.position = "none",
          axis.title.x = invis) +
    labs(y = "Number of genes nearer to functional\nenhancer than target gene")
  dev.off()  
  

      
################################################################################################################################ #
  ################################################################################################################################ #
  ## Nanostring replication ----
  
  library(reshape2)
  library(Hmisc)  
  library(forcats)
  
  ##Load file
  load(file="/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Validation_RNAseq/Nanostring/Results/Final_Revised_NatNeuro/ProcessedData.rda")    
  
  
  #Remove 300ng input samples
  meta_nano_filt<-meta_nano[c(meta_nano$Input_ng=="150"),]
  
  #Remove pool samples
  meta_nano_filt<-meta_nano_filt[c(meta_nano_filt$Pool=="FALSE"),] 
  
  #match metadata and expression data files order  
  index=match(rownames(meta_nano_filt), colnames(exp_nano))
  exp_nano=exp_nano[,index]
  
  #Divide by 1000 for easy plotting axis
  exp_nano<- exp_nano/1000
  
  #Create expression data file with relevant experimental information
  nano_data<-exp_nano
  colnames(nano_data)<- meta_nano_filt$Sample 
  nano_data<-t(nano_data)
  nano_data<-as.data.frame(nano_data)
  nano_data<-cbind(Batch = meta_nano_filt$Batch,Target = meta_nano_filt$Target, Enh = meta_nano_filt$Enh, nano_data)
  nano_data[nano_data$Target == "Neg"& is.na(nano_data$Enh), "Enh"] <- "Neg"
  
  
  nano_data$Enh_ID<- sub("_.*", "",nano_data$Enh )  
  
  #Separate by batch
  nano_data_2_3<-nano_data[c(nano_data$Batch=="Run_2"|nano_data$Batch=="Run_3"),]
  nano_data_4_5<-nano_data[c(nano_data$Batch=="Run_4"|nano_data$Batch=="Run_5"),]
  
  #Plots
  pal_iv_discrete2_reorder<-c("#3a5f9a","#fc9d77","#f3ccde","#6f083d","#dd5d2b","#2c2c54","#ddf4f8")
  
  
  #plot
  
  #Batch 2+3
  for(i in ((unique(nano_data_2_3$Target[nano_data_2_3$Target != "Neg"])))){
    nano_data_2_3$Condition<-ifelse(nano_data_2_3$Target==i,TRUE, FALSE)
    nano_data_2_3$Group<- ifelse(nano_data_2_3$Condition == TRUE,nano_data_2_3$Enh_ID,"Control")
    pdf(file = paste0(i, "_no_pool.pdf"), height = 3, width = 1.4 + (length(unique(nano_data_2_3$Group))/2))
    print(ggplot(nano_data_2_3, aes_string(x="Group", y=i, fill="Group")) +
            geom_dotplot(binaxis="y", stackdir="center", stackratio=1, dotsize=1.8)+
            stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom = "linerange", size=1)+
            stat_summary(geom = "point", shape = "-", size = 8, fun = mean) +
            scale_fill_manual(values=pal_iv_discrete2_reorder[1:length(unique(nano_data_2_3$Group))])+
            theme_classic() +
            theme(legend.position="none",text = element_text(size=13),axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),plot.title = element_text(margin=margin(b = 20, unit = "pt")))+
            scale_y_continuous()+
            labs(title = i, x= invis, y = invis))
    dev.off()
  }
  
  
  #Batch 4+5
  for(i in ((unique(nano_data_4_5$Target[nano_data_4_5$Target != "Neg"])))){
    nano_data_4_5$Condition<-ifelse(nano_data_4_5$Target==i,TRUE, FALSE)
    nano_data_4_5$Group<- ifelse(nano_data_4_5$Condition == TRUE,nano_data_4_5$Enh_ID,"Control")
    pdf(file = paste0(i, "_no_pool.pdf"), height = 3, width = 1.4 + (length(unique(nano_data_4_5$Group))/2))
    print(ggplot(nano_data_4_5, aes_string(x="Group", y=i, fill="Group")) +
            geom_dotplot(binaxis="y", stackdir="center", stackratio=1, dotsize=1.8)+
            stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom = "linerange", size=1)+
            stat_summary(geom = "point", shape = "-", size = 8, fun = mean) +
            scale_fill_manual(values=pal_iv_discrete2_reorder[1:length(unique(nano_data_4_5$Group))])+
            theme_classic() +
            theme(legend.position="none",text = element_text(size=13),axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),plot.title = element_text(margin=margin(b = 20, unit = "pt")))+
            scale_y_continuous()+
            labs(title = i, x= invis, y = invis))
    dev.off()
  }
  
  
  #PCGF5
  nano_data_4_5$Condition<-ifelse(nano_data_4_5$Target=="ANKRD1",TRUE, FALSE)
  nano_data_4_5$Group<- ifelse(nano_data_4_5$Condition == TRUE,nano_data_4_5$Enh_ID,"Control")
  pdf(file = "PCGF5_no_pool.pdf", height = 3, width = 1.4 + (length(unique(nano_data_4_5$Group))/2))
  print(ggplot(nano_data_4_5, aes_string(x="Group", y="PCGF5", fill="Group")) +
          geom_dotplot(binaxis="y", stackdir="center", stackratio=1, dotsize=1.8)+
          stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom = "linerange", size=1)+
          stat_summary(geom = "point", shape = "-", size = 8, fun = mean) +
          scale_fill_manual(values=pal_iv_discrete2_reorder[1:length(unique(nano_data_4_5$Group))])+
          theme_classic() +
          theme(legend.position="none",text = element_text(size=13),axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),plot.title = element_text(margin=margin(b = 20, unit = "pt")))+
          scale_y_continuous()+
          labs(title = "PCGF5", x= invis, y = invis))
  dev.off()
  
  
  
  #Wide plot for multiple enhancers per gene
  
  #Batch 2+3
  for(i in c("LGALS3","HSPB1")){
    nano_data_2_3$Condition<-ifelse(nano_data_2_3$Target==i,TRUE, FALSE)
    nano_data_2_3$Group<- ifelse(nano_data_2_3$Condition == TRUE,nano_data_2_3$Enh_ID,"Control")
    pdf(file = paste0(i, "_no_pool.pdf"), height = 3, width = 1.4 + (5/2))
    print(ggplot(nano_data_2_3, aes_string(x="Group", y=i, fill="Group")) +
            geom_dotplot(binaxis="y", stackdir="center", stackratio=1, dotsize=1.8)+
            stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom = "linerange", size=1)+
            stat_summary(geom = "point", shape = "-", size = 8, fun = mean) +
            scale_fill_manual(values=pal_iv_discrete2_reorder[1:length(unique(nano_data_2_3$Group))])+
            theme_classic() +
            theme(legend.position="none",text = element_text(size=13),axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),plot.title = element_text(margin=margin(b = 20, unit = "pt")))+
            scale_y_continuous()+
            labs(title = i, x= invis, y = invis))
    dev.off()
  }
  
  
  #Batch 4+5 
  for(i in c("NEAT1","PTMA")){
    nano_data_4_5$Condition<-ifelse(nano_data_4_5$Target==i,TRUE, FALSE)
    nano_data_4_5$Group<- ifelse(nano_data_4_5$Condition == TRUE,nano_data_4_5$Enh_ID,"Control")
    pdf(file = paste0(i, "_no_pool.pdf"), height = 3, width = 1.4 + (5/2))
    print(ggplot(nano_data_4_5, aes_string(x="Group", y=i, fill="Group")) +
            geom_dotplot(binaxis="y", stackdir="center", stackratio=1, dotsize=1.8)+
            stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom = "linerange", size=1)+
            stat_summary(geom = "point", shape = "-", size = 8, fun = mean) +
            scale_fill_manual(values=pal_iv_discrete2_reorder[1:length(unique(nano_data_4_5$Group))])+
            theme_classic() +
            theme(legend.position="none",text = element_text(size=13),axis.text.x=element_text(size=12),
                  axis.text.y=element_text(size=12),plot.title = element_text(margin=margin(b = 20, unit = "pt")))+
            scale_y_continuous()+
            labs(title = i, x= invis, y = invis))
    dev.off()
  }
  
  
  
  #HSPB1 Plot
  pal_iv_discrete2_reorder2<-c("#3a5f9a","#2c2c54","#fc9d77","#f3ccde","#6f083d","#dd5d2b","#ddf4f8")
  
  
  nano_data_HS<-nano_data[c(nano_data$Batch=="Run_1"|nano_data$Batch=="Run_2"|nano_data$Batch=="Run_3"),]
  nano_data_HS$HS<-ifelse(nano_data_HS$Target == "HSPB1",nano_data_HS$Enh_ID,
                          ifelse(nano_data_HS$Batch=="Run_1"&nano_data_HS$Target != "HSPB1", "Control1", "Control2"))
  pdf(file = "HSPB1_Batch1-3.pdf", height = 3, width = 5.5)
  print(ggplot(nano_data_HS, aes_string(x=forcats::fct_relevel(nano_data_HS$HS,"Control1", "Enh854", "Control2", 
                                                               "Enh855", "Enh857"), y="HSPB1", fill="HS")) +
          geom_dotplot(binaxis="y", stackdir="center", stackratio=1, dotsize=1.8)+
          stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom = "linerange", size=1)+
          stat_summary(geom = "point", shape = "-", size = 8, fun = mean) +
          scale_fill_manual(values=pal_iv_discrete2_reorder2[1:length(unique(nano_data_HS$HS))])+
          theme_classic() +
          theme(legend.position="none",text = element_text(size=13),axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),plot.title = element_text(margin=margin(b = 20, unit = "pt")))+
          scale_y_continuous()+
          labs(title = "HSPB1", x= invis, y = invis))
  dev.off()
  

  
################################################################################################################################ #
## Exemplar: NEAT1 ----
  
## Setup
  # neat1's enhancers
  neat <- res.final[which(res.final$HitPermissive & res.final$Gene == "NEAT1"),]

  # gene expression
  load("../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/Pseudobulk2_ATACesquePooling.rda", verbose = TRUE)
  load("../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/scRNAseq_UMIs_Seurat.rda", verbose = TRUE)  
  
  # atac
  load("../../../FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Intersections and Coverage.rda", verbose = TRUE)

## Process
  # process expression
  enh <- herring$CoverageMean_Pooled
  enh <- enh[neat$Enh,]
  enh <- apply(enh, 1, function(x) {
    # x / max(x)
    (x - mean(x)) / sd(x)
  })
  enh <- as.data.frame(enh)
  
  # process expression
  gene <- pb2$Final  
  gene <- log2(gene + 1)
  gene <- apply(gene, 1, function(y) {
    # y / max(y)
    # scale(y)
    (y - mean(y)) / sd(y)
  })
  gene <- t(gene)

  gene <- data.frame(NEAT1 = as.numeric(gene["NEAT1",]),
                  Sample = colnames(gene),
                  row.names = colnames(gene))
  
  
  # combine
  p <- cbind(enh, gene[rownames(enh),])
  p <- melt(p)
  p$Stage <- splitter(p$Sample, "_", 2)  
  p$Ct <- splitter(p$Sample, "_", 1)  
  p$Ct <- factor(p$Ct, levels = c("Astro", "Oligo", "Micro", "Exc", "Inh"))
  p$Stage <- factor(p$Stage, levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
  levels(p$Stage) <- c("Fetal", "Neont", "Infnt", "Child", "Adols", "Adult")
  # p$Astro <- p$Ct == "Astro"
  p$variable <- factor(p$variable, levels = c("Enh219", "NEAT1", "Enh218", "Enh220"))
  
  # pdf(file = "NEAT1 - Heatmap (Z scaling).pdf", height = 2.5, width = 4.5)
  pdf_mainDE(figNo = "4F", title = "NEAT1 heatmap", h = 2.3, w = 3.2)
  ggplot(p, aes(x = Stage, y = variable, fill = value)) +
    geom_tile() +
    facet_grid(.~Ct, scales = "free", space = "free", switch = "x") +
    theme_bw() +
    scale_fill_gradient2(low = pals$grn2orng[9], mid = pals$grn2orng[5], high = pals$grn2orng[1], limits = c(-2.75, 2.75)) + # to orange
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    geom_vline(xintercept = c(0.5), linetype = 2, size = 0.5) +
    guides(fill = guide_colourbar(title = "Scaled\nvalue")) +
    theme(panel.border = invis, axis.ticks = invis, axis.line.y = element_line(linewidth = 0.6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), 
          axis.text.y = element_text(size = 7),
          legend.position = "bottom",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          strip.text = element_text(size = 8),
          axis.title = invis, panel.spacing = unit(0, "lines")) +
    labs(x = "Tissue", y = "Enhancer")   
  dev.off()
  
  
## UMAP
  # read in single-cell timecourse
  load("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/scRNAseq_UMIs_Seurat.rda", verbose = TRUE)  
  brainMat <- NormalizeData(brainMat)
  
  # read in umap
  umap <- read.csv("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/UMAP_coordinates.csv")
  # brainMat@reductions$umap <- umap
  
  # add umap to object
  umap <- as.matrix(umap[,2:3]) %>% as.data.frame()
  rownames(umap) <- colnames(brainMat)
  colnames(umap) <- paste0("UMAP_", 1:2)
  brainMat@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(brainMat))
 
  # Idents(object = brainMat) <- "major_clust"
  brainMat$Ct <- factor(brainMat$major_clust)
  levels(brainMat$Ct)[grepl("L2|L4|L5|PN_", levels(brainMat$Ct))] <- "Exc"
  levels(brainMat$Ct)[grepl("CGE|MGE|ID2|LAMP|PV|SST|VIP", levels(brainMat$Ct))] <- "Inh"
  Idents(object = brainMat) <- "Ct"
  
  p <- subset(brainMat, cells = which(brainMat$cell_type != "Poor-Quality" & brainMat$Ct != "Vas"))
  # p <- subset(p, cells = which(brainMat$Ct != "Vas"))
  
  # pdf(file = "NEAT1 - UMAP.pdf", height = 2.7, width = 3.5)
  pdf_mainDE(figNo = "4E", title = "NEAT1 UMAP", h = 2.2, w = 2.2)
  FeaturePlot(p, features = "NEAT1", raster = FALSE, label = TRUE, 
              repel = FALSE, label.size = 2.5) +
    theme_bw() +
    guides(colour = guide_colourbar(title = "NEAT1")) +
    theme(panel.grid = invis, plot.title = invis, panel.border = invis,
          axis.line = element_line(),
          legend.title = element_text(size = 10),
          legend.key.size = unit(0.4, "cm"),
          axis.text.y = text90,
          legend.position = "none") +
    scale_colour_gradientn(colours = c("#b51a00","#ee5900","#ff9d68","#feceb8","grey95"),
                         values = c(1.0,0.7,0.5,0.4,0.2,0)) +
    labs(x = "UMAP1", y = "UMAP2") 
  dev.off()
    

  
