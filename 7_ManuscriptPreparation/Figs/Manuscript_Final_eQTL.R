#This script generates figures 6A, 6B, 6C and S8A
## Setup
  setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_eQTL/")
  library(tidyverse)
  library(rcartocolor)
  library(ggsci)
  library(ggplot2)
  library(ggbeeswarm)
  library(readxl)
  source("../../../FullScale/Scripts/Functions.R")
  source("../FinalFigureFunctions.R")
  
## Load eQTL results
  snp.annot <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/Final - SNP Annotation.csv")
  interesting <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/Final - SNPs of Interest.csv")
  load("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL - All Overlaps.rda", verbose = TRUE)
  # rep_perPair <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/")
  
## DE results
  res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  guides <- read.csv("../../../FullScale/Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
 
## Plotting parameters
  invis <- element_blank()
  maxh <- 29.7 / 2.54
  maxw <- 21.0 / 2.54
  
  theme_gjs <- theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) 
  NoExpand <- scale_y_continuous(expand = c(0,0))

## Functions to write to disk in a traceable way
  pdf_eqtl <- function(figNo, title, h, w) {
    pdf(file = paste0("../Final/", figNo, " - Script eqtl - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_eqtl <- function(figNo, title, toPrint) {
    sink(paste0("../Final/", figNo, " - Script eqtl - ", title, ".txt"))
    print(toPrint)
    sink()
  }
  
################################################################################################################################ #
## Count your SNPs ----
  

## Number of SNPs in hit enhancers 
  
  x <- snp.annot
  x <- x[which(x$Hit),]

  sink_eqtl(figNo = "6A", title = "SNPs per enhancer", table(x$Category))
  
  # table(x$Category)
  
  # LD      Nearby Overlapping 
  # 13931        1454         301 
  
  # distributions
  p <- snp.annot
  p <- p[which(p$Hit),]
  
  tab <- table(p$Enh, p$Category)
  tab <- as.data.frame.matrix(tab)
  tab$Total <- rowSums(tab)
  tab$Enh <- rownames(tab)
  tab$Enh <- factor(tab$Enh, levels = tab$Enh[order(tab$Total)])
  tab <- melt(tab, id.vars = c("Enh", "Total"))
  tab$variable <- factor(tab$variable, levels = rev(levels(tab$variable)))
  levels(tab$variable) <- gsub("Overlapping", "Overlaps", levels(tab$variable))
  
  offset <- 1
  
  # as violin
  # pdf(file = "Count of SNPs per enhancer.pdf", height = 2.3, width = 2.1)
  pdf_eqtl(figNo = "6A", title = "SNPs per enhancer", h = 2.3, w = 2.1)
  ggplot(tab, aes(x = variable, y = value+offset, fill = variable, colour = variable)) +
    geom_violin(draw_quantiles = 0.5, scale = "width", width = 0.7) +
    labs(x = "SNP position", y = "SNPs per functional enhancer") +
    theme_bw() +
    scale_fill_manual(values = pals$Primary[c(4,5,7)]) +
    scale_colour_manual(values = pals$Primary_Darker[c(4,5,7)]) +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(), axis.title.x = invis,
          legend.position = "none", axis.text.y = text90) +
    scale_y_continuous(expand = c(0,0), limits = c(1, 1025),
                       trans = "log2", labels = function(x) { comma(x - offset) },
                       breaks = c(0, 1, 2, 8, 32, 128, 512) + offset)
    dev.off()
 

################################################################################################################################ #
## Are eQTLs supported by experimental data  ---- 
  
## Analyses focusing on the pooled data
  # read in
  sigeQTL_pair_pooled <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs pooled.csv") # from significance call
  fmeQTL_pair_pooled <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs pooled, Fine-mapped.csv") # from significance call
  asteQTL_pair_pooled <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs pooled, Ast-specific.csv") # from significance call
  
  # plot
  p <- data.frame(Sig = sigeQTL_pair_pooled$eQTL_Category,
                  FM = fmeQTL_pair_pooled$eQTL_Category,
                  Ast = asteQTL_pair_pooled$eQTL_Category,
                  Pair = sigeQTL_pair_pooled$Pair)
  p <- melt(p, id.vars = "Pair")

  p$variable <- factor(p$variable, levels = c("Sig", "FM", "Ast"))
  levels(p$variable) <- c("Brain\nSignificant", "Brain\nFine-mapped", "Astrocyte\nSignificant")
  
  p$value <- gsub(" non-expressed", "", p$value) %>%
    splitter(., " \\(", 1)
  p$value <- factor(p$value, levels = unique(p$value)[c(2,1,3)])
  
  pal <- c("grey95", pals$Primary[c(2,6)])
  
  p <- table(p$value, p$variable) %>%
    proportions(2) %>% 
    as.data.frame()
  
  levels(p$Var1) <- gsub(" ", "\n", levels(p$Var1)) %>%
    gsub("gene", "gene\ntargeted", .)
  
  p$Label2 <- p$Label <- paste0(round(p$Freq*100), "%")
  toDodge <- which(p$Label == "4%") 
  p$Label2[-toDodge] <- NA
  p$Label[toDodge] <- NA

  # pdf(file = "Consistency rate (pooled) (V2).pdf", height = 3, width = 3)
  pdf_eqtl(figNo = "6B", title = "eqtl vs crispri (pooled)", h = 3, w = 3)
  ggplot(p, aes(x = ".", fill = Var1, y = Freq * nrow(sigeQTL_pair_pooled))) +
    geom_col(position = "stack", width = 0.7, colour = "black", linewidth = 0.5) +
    theme_bw() +
    geom_text(mapping = aes(label = Label), position = "stack", vjust = 2, size = 3) +
    geom_text(mapping = aes(label = Label2), position = "stack", vjust = -1, size = 3) +
    facet_wrap(~Var2, strip.position = "bottom", scales = "free_x") +
    guides(fill = guide_legend(title = "eQTL target gene", title.position = "top", title.hjust = 0.5)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 158), limits = c(0,165)) +
    labs(y = "CRISPRi EGPs", x = "Pooled eQTL datasets") +
    # scale_x_discrete(expand = c(0.02, 0.02)) +
    scale_fill_manual(values = pal) +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
          axis.text.y = text90,
          # legend.title = invis,
          axis.text.x = invis,
          # axis.title.x = invis,
          # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid = invis, axis.ticks.x = invis)
    dev.off()
    
    
## Looking within each dataset
  eQTL_pair_across <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs across studies.csv", row.names = 1) # from significance call
  
  p <- eQTL_pair_across[,-c(2:3)]
  p <- melt(p, id.vars = "Pair")
  colnames(p) <- c("Pair", "Study", "Cat")
  # p$Study <- factor(p$Study, levels = unique(p$Study)[c(1,2,3,4,6,5)])
  p$Astro <- grepl("Ast", p$Study) %>% factor()
  p$Finemapped <- grepl("FM", p$Study) %>% factor()
  levels(p$Astro) <- c(" Brain eQTL", "Astro eQTL")
  levels(p$Finemapped) <- c("Significant", "Finemapped")
  levels(p$Study) <- gsub("Ast|_FM", "", levels(p$Study))
  
  p$Cat <- splitter(p$Cat, " \\(", 1)
  p$Cat <- factor(p$Cat, levels = unique(p$Cat)[1:3])
  pal <- c("grey95", pals$Primary[c(2,6)])
  p$Group <- paste0(p$Finemapped, "\n", p$Astro)
  p$Group <- gsub(" eQTL", "", p$Group)
  p$Group <- factor(p$Group)
  levels(p$Group) <- c("Brain\nFine-mapped", "Brain\nSignificant", "Astrocyte\nSignificant")
  
  # pdf(file = "Consistency rate (across).pdf", height = 3, width = 4)
  pdf_eqtl(figNo = "SFig8", title = "eqtl vs crispri (per study)", h = 3, w = 4)
  ggplot(p, aes(x = Study, fill = Cat)) +
    geom_bar(position = "stack", width = 0.7, linewidth = 0.5) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 158), limits = c(0,165)) +
    facet_grid(.~Group, scales = "free", space = "free", switch = "x") +
    labs(y = "CRISPRi EGPs", x = "eQTL Resource") +
    scale_fill_manual(values = pal) +
    scale_alpha_manual(values = c(0.85, 0.85, 1)) +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
          axis.text.y = text90,
          legend.title = invis, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = invis, panel.grid = invis)
  dev.off()


################################################################################################################################ #
## How replicable are the eQTL target genes across studies? How does this relate to having experimental evidence? ----  
  
## Read in significant eQTL calls for four studies
  sigeQTL_pair_pooled <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs pooled.csv") # from significance call
  x <- sigeQTL_pair_pooled
  x <- x[-which(x$eQTL_Category == "No eQTL"),]

## Are they replicable?
  x$Rep <- grepl("(replicably)", x$eQTL_Category, fixed = TRUE) # this works for when "different" gene identified
  x$Rep[which(x$nRep >= 2)] <- TRUE # and works for the "same" gene identified
  x$Same <- x$eQTL_Category == "Same gene"
  
## Check the replication rate versus same gene, at the eQTL level
  rep <- list()
 
  # when an "other" gene is detected for the enhancer, is it replicated?
  rep$Expressed <- list()
  
  rep$Expressed <- strsplit(x$OtherGenes_Expressed, "|", fixed = TRUE) # creates a list 
  rep$Expressed <- do.call("c", rep$Expressed) # concatenates the list into a vector
  rep$Expressed <- data.frame(Gene = rep$Expressed, Same = FALSE) # converts the list to a dataframe, which states that the same gene is "FALSE" (by definion we selected these)
  rep$Expressed$Rep <- splitter(rep$Expressed$Gene, "\\(", 2) # get number in brackets after the gene, which denotes how many studies it's significant in
  rep$Expressed <- rep$Expressed[-which(is.na(rep$Expressed$Gene)),] # remove rows with NA (when there is no expressed other genes for the EGP)
  rep$Expressed$Rep <- !(grepl("1", rep$Expressed$Rep)) # it is replicated if a "1" is not found. That is 2/3/4 is the number of replications
  
  # as above, but for non-expressed genes
  rep$NotExpressed <- strsplit(x$OtherGenes_NotExpressed, "|", fixed = TRUE)
  rep$NotExpressed <- do.call("c", rep$NotExpressed)
  rep$NotExpressed <- data.frame(Gene = rep$NotExpressed, Same = FALSE)
  rep$NotExpressed$Rep <- splitter(rep$NotExpressed$Gene, "\\(", 2)
  rep$NotExpressed <- rep$NotExpressed[-which(is.na(rep$NotExpressed$Gene)),]
  rep$NotExpressed$Rep <- !(grepl("1", rep$NotExpressed$Rep))
  
  # what about when the same gene is detected in eQTL as the CRISPRi EGP?
  rep$Same <- list()
  rep$Same <- x[which(x$Same),]
  rep$Same <- data.frame(Gene = rep$Same$Gene,
                  Same = TRUE,
                  Rep = rep$Same$nRep >= 2)
  
  # bind
  p <- rbind(rep$Expressed,
             rep$NotExpressed,
             rep$Same)
  
  f <- table(p$Same, p$Rep) %>% fisher.test()
  sink_eqtl(figNo = "6C", title = "eqtl reproducibility", f)
  
  
  # plot
  tab <- table(p$Same, p$Rep) %>%
    proportions(1) %>%
    as.data.frame.table()
  
  colnames(tab) <- c("Same", "eQTL", "Fraction")
  tab$Same <- factor(tab$Same, levels = c("TRUE", "FALSE"))
  levels(tab$Same) <- c("eQTL targets\nsame gene\nas CRISPRi", "eQTL targets\ndifferent gene\nto CRISPRi")
  # levels(tab$eQTL) <- c("1", "2+")
  levels(tab$eQTL) <- c("False", "True")
  

  # pdf(file = "eQTL Cross-study Reproducibility.pdf", height = 2.2, width = 2.4)
  pdf_eqtl(figNo = "6C", title = "eqtl reproducibility", h = 2.2, w = 2.4)
  ggplot(tab, aes(y = Same, fill = eQTL, colour = eQTL, x = Fraction * 100)) +
    geom_col(width = 0.5) +
    theme_bw() +
    scale_x_continuous(expand = c(0,0), breaks = c(0, 50, 100), limits = c(-2,105)) +
    # guides(fill = guide_legend(nrow = 2, title = "Number of studies\nwhere eQTL is called"), colour = "none") +
    guides(fill = guide_legend(nrow = 2, title = "eQTL reproduced\nacross studies"), colour = "none") +
    labs(x = "Percent", y = "eQTL vs. enhancer target gene") +
    scale_fill_manual(values = c("grey95", pals$One)) +
    scale_colour_manual(values = c("grey80", "black")) +
    theme(panel.border = invis, axis.line.x = element_line(), legend.position = "none", axis.title.y = invis,
          panel.grid = invis, axis.ticks.y = invis)
  
  dev.off()
  
  
