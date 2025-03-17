#This script generates figures 2B
## Setup
  setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_SanityChecks/")
  library(tidyverse)
  library(rcartocolor)
  library(ggsci)
  library(ggplot2)
  library(ggbeeswarm)
  library(readxl)
  source("../../../FullScale/Scripts/Functions.R")
  source("../FinalFigureFunctions.R")
  
## Saving function
  write.fisher <- function(tab, name) {
    f <- fisher.test(tab)
    data.frame(OR = round(f$estimate, 2), p = f$p.value, row.names = name)
  }
  
## DE results
  res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  guides <- read.csv("../../../FullScale/Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  
  candidate.annot <- read.csv("../../../FullScale/Results/3_HitEnrichment/Chromatin/Final - Annotation Logical.csv", row.names = 1)
  candidate.enrich <- read.csv("../../../FullScale/Results/3_HitEnrichment/Chromatin/Final - Enrichments.csv", row.names = 1)
  
## Read in power dataframe
  pow <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable3_DE/3F_Power.csv")
  wellpowered_egps <- pow$Pair[which(pow$WellPowered015  | pow$Hit)]
  wellpowered_enh <- pow$Enhancer[which(pow$WellPowered015 | pow$Hit)] %>% unique()
  wellpowered_genes <- pow$Gene[which(pow$WellPowered015 | pow$Hit)] %>% unique()

  
## Plotting parameters
  invis <- element_blank()
  maxh <- 29.7 / 2.54
  maxw <- 21.0 / 2.54
  
  theme_gjs <- theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) 
  NoExpand <- scale_y_continuous(expand = c(0,0))

## Functions to write to disk in a traceable way
  pdf_sanity <- function(figNo, title, h, w) {
    pdf(file = paste0("../Final/", figNo, " - Script SanityCheck - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_sanity <- function(figNo, title, toPrint) {
    sink(paste0("../Final/", figNo, " - Script SanityCheck - ", title, ".txt"))
    print(toPrint)
    sink()
  }  
  
  
################################################################################################################################ #
## Plots on Hit Enhancers ----
  

## Plot 0: ATAC pileup
  p <- read.csv("../../../FullLibrary_Selection/Results/Peaks_Annotated.csv")
  # p <- p[which(p$id %in% res.final$Enh.Pos),]
  p <- p[which(p$id %in% res.final$Enh.Pos[which(res.final$Enh %in% wellpowered_enh)]),]
  
  # annotate as hit
  p$Hit <- p$id %in% res.final$Enh.Pos[which(res.final$HitPermissive)]
  
  # reformat
  p$Hit <- factor(p$Hit)
  levels(p$Hit) <- c("Inactive\ncandidates", "Functional\nenhancers")
  
  # statistic
  sink_sanity(figNo = "3G", title = "ATAC pileup", wilcox.test(p$cpm ~ p$Hit))
  
  # plot
  # pdf(file = "Inhouse ATAC Pileup.pdf", height = 2, width = 1.6)
  pdf_sanity(figNo = "3G", title = "ATAC pileup", h = 2, w = 1.6)
  ggplot(p, aes(y = cpm, x = Hit, fill = Hit, colour = Hit)) +
    geom_violin(scale = "width", alpha = 0.75) +
    geom_boxplot(width = 0.2, outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = pals$Hits) +
    scale_colour_manual(values = pals$Hits_Darker) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), limits = c(0, 35)) +
    theme(axis.title.x = invis, panel.grid = invis, panel.border = invis,
          axis.text.y = text90, axis.text.x = invis, axis.ticks.x = invis,
          axis.line.y = element_line(), legend.position = "none") +
    labs(y = "ATAC pileup (CPM)")
  dev.off()




################################################################################################################################ #
## A combined plot ----
        
## Of superenhancer, TAD, and K562 enhancers
  yMax <- 25    

## Collect data
  # tads
  tad <- read.csv("../../../FullScale/Results/3_HitEnrichment/EnhGenePairs/TAD/Pair-TAD Annotation.csv", row.names = 1)
  tad <- tad[which(tad$Pair %in% wellpowered_egps),]
  x <- table(tad$HitPermissive, tad$CrossTAD)
  fish <- fisher.test(x)
  x <- x / rowSums(x) * 100
  x <- as.data.frame(x)
  x <- x[which(as.logical(x$Var2)),]
  # levels(x$Var1) <- c("Inactive\nEGP links", "Functional\nEGP links")
  x$Resource <- "TAD Boundary"
  
  x_ttl <- paste0("p=", signif(fish$p.value, 1), ", ", 
                "OR=", signif(fish$estimate, 2))
  sink_sanity(figNo = "2C", title = "TAD", x_ttl)
  
  # superenhancers
  y <- candidate.annot[which(candidate.annot$Tested),]
  y <- y[which(rownames(y) %in% wellpowered_enh),]
  y <- table(y$Hit, y$Superenhancer)
  
  y_ttl <- fisher.test(y)
  
  y <- y / rowSums(y) * 100
  y <- as.data.frame(y)
  y <- y[which(as.logical(y$Var2)),]
  # levels(y$Var1) <- c("Inactive\ncandidates", "Functional\nenhancers")
  y$Resource <- "Superenhancer"
  
  # y_ttl <- candidate.enrich["Superenhancer",]
  # y_ttl <- paste0("p=", signif(y_ttl$p, 2), ", ", 
  #               "OR=", signif(y_ttl$OR, 2))
  sink_sanity(figNo = "2C", title = "Superenhancer", y_ttl)
  
  # k562 enhancers
  z <- candidate.annot[which(candidate.annot$Tested),]
  z <- z[which(rownames(z) %in% wellpowered_enh),]
  z <- table(z$Hit, z$ValidatedEnh_K562_Yao2022)
  z_ttl <- fisher.test(z)
  z <- z / rowSums(z) * 100
  z <- as.data.frame(z)
  z <- z[which(as.logical(z$Var2)),]
  # levels(z$Var1) <- c("Inactive\ncandidates", "Functional\nenhancers")
    
  z$Resource <- "K562 enhancer"
  
  # z_ttl <- candidate.enrich["ValidatedEnh_K562_Yao2022",]
  # z_ttl <- paste0("p=", signif(z_ttl$p, 2), ", ", 
  #               "OR=", signif(z_ttl$OR, 2))
  sink_sanity(figNo = "2B", title = "K562", z_ttl)
  
## Plot
  p <- rbind(x,y,z)
  p$Var1 <- factor(p$Var1)
  levels(p$Var1) <- c("Inactive\ncandidates", "Functional\nenhancers")
  p$Resource <- gsub(" ", "\n", p$Resource)
  
  # pdf(file = "Combined sanity check barplot.pdf", height = 2, width = 3)
  pdf_sanity(figNo = "2B", title = "K562, Superenhancer, TAD", h = 2, w = 3)
  ggplot(p, aes(x = Resource, y = Freq, colour = Var1, fill = Var1)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_colour_manual(values = pals$Hits_Darker) +
      scale_fill_manual(values = pals$Hits) +
      # facet_wrap(~Resource, scales = "free_x") +
      theme_bw() +
      theme_gjs +
      scale_y_continuous(limits = c(0, yMax), expand = c(0,0), breaks = c(0, 10, 20)) +
      theme(axis.title.x = invis, legend.position = "none", axis.ticks.x = invis,
            panel.grid = invis, plot.title = element_text(size = 6)) +
      labs(y = "Percentage")
  dev.off()
  
