## Setup
  setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_TTseq/")
  library(tidyverse)
  library(rcartocolor)
  library(ggsci)
  library(ggplot2)
  library(ggbeeswarm)
  library(reshape2)
  library(readxl)
  library(VennDiagram)
  source("../../../FullScale/Scripts/Functions.R")
  source("../FinalFigureFunctions.R")
  
## Saving function
  write.fisher <- function(tab, name) {
    f <- fisher.test(tab)
    data.frame(OR = round(f$estimate, 2), p = f$p.value, row.names = name)
  }
  
  write.fisher.full <- function(tab, name) {
    f <- fisher.test(tab)
    data.frame(OR = round(f$estimate, 2), 
               p = f$p.value, 
               lower = f$conf.int[1],
               upper = f$conf.int[2],
               row.names = name)
  }
  
  convert.id <- function(x) {
    m <- match(x, guides$TargetCoord)
    y <- guides$TargetID[m]
    return(y)
  }
  
## DE results
  res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  guides <- read.csv("../../../FullScale/Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
 
  invis <- element_blank()
  maxh <- 29.7 / 2.54
  maxw <- 21.0 / 2.54
  
  theme_gjs <- theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) 
  NoExpand <- scale_y_continuous(expand = c(0,0))

## TTseq data
  tt <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/TTseq/Results Table.csv", row.names = 1)
  tt <- tt[which(rownames(tt) %in% res.final$Enh),]
  transcribed <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/TTseq/Transcriptional classification.csv")
  
  # reorder and filter the latter
  m <- match(rownames(tt), transcribed$Enh)
  transcribed <- transcribed[m,] # just tested enhancers (957/979)
  
## Functions to write to disk in a traceable way
  pdf_tt <- function(figNo, title, h, w) {
    pdf(file = paste0("../Final/", figNo, " - Script tt - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_tt <- function(figNo, title, toPrint) {
    sink(paste0("../Final/", figNo, " - Script tt - ", title, ".txt"))
    print(toPrint)
    sink()
  }
  
## Read in power dataframe
  pow <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable3_DE/3F_Power.csv")
  wellpowered_egps <- pow$Pair[which(pow$WellPowered015  | pow$Hit)]
  wellpowered_enh <- pow$Enhancer[which(pow$WellPowered015 | pow$Hit)] %>% unique()
  wellpowered_genes <- pow$Gene[which(pow$WellPowered015 | pow$Hit)] %>% unique()

################################################################################################################################ #
## Comparison of TTseq to RNAseq ----
    
    
## Expression distribution
  p <- data.frame(Hit = tt$Hit,
                  TTseq = (tt$TTseq_Total),
                  RNAseq = (tt$RNAseq_Total))

  # p <- melt(p)
  p$Hit <- factor(p$Hit)
  levels(p$Hit) <- c("Inactive candidate", "Functional Enhancer")
  offset <- 1
  p$Density <- get_density(p$TTseq, p$RNAseq)
  
  # pdf(file = "RNAseq vs TTseq scatterplot.pdf", height = 3, width = 3)
  pdf_tt(figNo = "SFig6B", title = "Scatterplot RNAseq vs TTseq", h = 3, w = 3)
  ggplot(p, aes(x = RNAseq+offset, y = TTseq+offset )) +
    geom_point(colour = pals$One ) +
    geom_rug(col = rgb(0.5, 0, 0, alpha = 0.2)) +
    theme_bw() +
    scale_y_continuous(trans = "log2", labels = function(x) { x - offset }, breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128)+ offset) +
    scale_x_continuous(trans = "log2", labels = function(x) { x - offset }, breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128)+ offset) +
    geom_hline(yintercept = 3+offset, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = 3+offset, linetype = 2, alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
    theme(panel.grid = invis, #panel.border = invis,
          axis.line.y = element_line(), legend.position = "none") +
    labs(x = "RNAseq read count in candidates", y = "TTseq read count in candidates")
  dev.off()  
  

################################################################################################################################ #
## Read distribution in hits and non-hits ----  
  
  
  p <- data.frame(Enh = rownames(tt),
                  Hit = tt$Hit,
                  TTseq = (tt$TTseq_Total),
                  RNAseq = (tt$RNAseq_Total))
  m <- match(p$Enh, transcribed$Enh)
  p$Cat <- transcribed$Category_3[m]
  p <- p[-which(p$Cat == "Not transcribed"),]
  
  p <- p[which(p$Enh %in% wellpowered_enh),] # 424 to 350


p <- melt(p)
p$Hit <- factor(p$Hit)
levels(p$Hit) <- c("Inactive\ncandidates", "Functional\nenhancers")

offset <- 0.1
  
# pdf(file = "Counts in hits versus non-hits.pdf", height = 2, width = 2)
pdf_tt(figNo = "3E", title = "Reads in hits versus non-hits", h = 2, w = 2)
ggplot(p[which(p$variable == "TTseq"),], aes(x = Hit, colour = Hit, fill = Hit, y = value+offset)) +
  geom_violin(scale = "width", position = position_dodge(width = 0.7), width = 0.7, alpha = 0.7) +
  geom_boxplot(colour = "black", width = 0.2, outlier.shape = NA, show.legend = FALSE) +
  theme_bw() +
  scale_fill_manual(values = pals$Hits) +
  scale_colour_manual(values = pals$Hits_Darker) +
  scale_y_continuous(trans = "log2", 
                     expand = c(0,0),
                     limits = c(1, NA),
                     labels = function(x) { x - offset }, 
                     breaks = c(0, 1, 4, 16, 64, 256)+ offset) +
  theme(axis.title.x = invis, panel.grid = invis, panel.border = invis,
        axis.line.y = element_line(), legend.position = "none",
        axis.text.y = text90, axis.text.x = invis, axis.ticks.x = invis) +
  labs(y = "TTseq read count at\ntranscribed enhancers")

dev.off()
  

## Statistics
q <- p[which(p$variable == "TTseq"),]
sink_tt(figNo = "4E", title = "Reads in hits versus non-hits", wilcox.test(q$value ~ q$Hit))

 # p = 0.0446


################################################################################################################################ #
## New definition of expressed ----  
  
  
## Here, a peak is expressing eRNA if > X reads on at least one strand, and the total reads are greater in TT than RNA
  

## Plot
  final.thresh <- 3
  
  ## Functions
    plot.cats <- function(thresh, complex = FALSE) {
      
      transcribed_filt <- transcribed[which(transcribed$Enh %in% wellpowered_enh),]
      
      # stats
      fish_trns <- table(transcribed_filt$Hit, transcribed_filt[,paste0("Category_", thresh)] != "Not transcribed") %>% fisher.test()
      lab_trns <- paste0("Transcribed: p=", 
                         signif(fish_trns$p.value, 2),
                         ", OR=", signif(fish_trns$estimate, 2))
      
      fish_uni <- table(transcribed_filt$Hit, transcribed_filt[,paste0("Category_", thresh)] == "Unidirectional") %>% fisher.test()
      lab_uni <- paste0("Unidirectional: p=", 
                         signif(fish_uni$p.value, 2),
                         ", OR=", signif(fish_uni$estimate, 2))
      
      fish_bi <- table(transcribed_filt$Hit, transcribed_filt[,paste0("Category_", thresh)] == "Bidirectional") %>% fisher.test()
      lab_bi <- paste0("Bidirectional: p=", 
                         signif(fish_bi$p.value, 2),
                         ", OR=", signif(fish_bi$estimate, 2))
    
      # title
      ttl <- paste0("Threshold: ", thresh,
                    "\n", lab_trns,
                    "\n", lab_uni,
                    "\n", lab_bi)
      
      # plot values
      tab <- table(transcribed_filt$Hit, transcribed_filt[,paste0("Category_", thresh)])
      tab <- tab / rowSums(tab)
      tab <- t(tab)
      tab <- as.data.frame.matrix(tab)
      tab$Transcribed <- rownames(tab)
      
      p <- melt(tab)
      colnames(p)[2] <- "Hit"
      p$Hit <- factor(p$Hit, levels = c("FALSE", "TRUE"))  
      levels(p$Hit) <- (c("Inactive\ncandidates", "Functional\nenhancers")  )
      p$Transcribed <- factor(p$Transcribed, levels = c("Not transcribed", "Unidirectional", "Bidirectional"))
  
      if (complex) {
        pal <- c("grey90", pals$Hits[1], pals$Hits_Darker[1], # non-hit colours
               "grey91", pals$Hits[2], pals$Hits_Darker[2]) # hit colours
        
         ggplot(p, aes(x = Hit, fill = interaction(Transcribed, Hit), y = value*100, alpha = Transcribed)) +
          geom_col(width = 0.7) +
          theme_bw() +
           scale_alpha_manual(values = c(1, 0.8, 1)) +
          theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
                axis.title.x = invis, plot.title = element_text(size = 6),
                axis.text.y = text90) +
          scale_fill_manual(values = pal) +
          scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100)) +
          labs(y = "Percent of peaks", title = ttl)
        
        
      } else {
        # pal <- c("grey90", pals$Primary[1], pals$Primary[8])
        pal <- c("grey90", pals$Primary[2], pals$Primary_Darker[2])
        ggplot(p, aes(x = Hit, fill = Transcribed, y = value*100, linetype = Transcribed, alpha = Transcribed)) +
          geom_col(colour = "black", width = 0.7) +
          theme_bw() +
          scale_alpha_manual(values = c(1, 0.9, 1)) +
          theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
                axis.title.x = invis, plot.title = element_text(size = 6),
                axis.text.y = text90) +
          scale_fill_manual(values = pal) +
          scale_linetype_manual(values = c("blank", "dashed", "solid")) +
          scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100)) +
          labs(y = "Percent of peaks", title = ttl)
      }
      
      
      
    }
    

  ## Apply function
    # pdf(file = "Transcriptional Categories (V2.2).pdf", height = 2.2, width = 2)  
    pdf_tt(figNo = "3D", title = "eRNA", h = 2.2, w = 4)
    plot.cats(final.thresh) 
    dev.off()
 
## Test whether functional EGPs are more likely to be bidirectional (given transcription is assumed)
    p <- transcribed[which(transcribed$Enh %in% wellpowered_enh),]
    p <- p[which(p$Category_3 != "Not transcribed"),]
    fish <- table(p$Hit, p$Bidirectional_3)  %>% fisher.test()
    sink_tt(figNo = "3D", title = "eRNA bidirectional if transcribed", fish)
    
## Plot across thresholds
  p <- transcribed[which(transcribed$Enh %in% wellpowered_enh),]
  p <- apply(p[, grep("Category", colnames(p)),], 2, function(x) {
    y <- table(p$Hit, x)
    y <- y / rowSums(y)
    y <- as.data.frame.matrix(y)
  })
  
  p <- do.call("rbind", p)
  p$Threshold <- splitter(rownames(p), "\\.", 1) %>%
    splitter("_", 2) 
  p$Hit <- splitter(rownames(p), "\\.", 2) %>% factor()
  levels(p$Hit) <- c("Inactive candidates", "Functional Enhancers")  

  p <- melt(p, id.vars = c("Hit", "Threshold"))
  p$variable <- factor(p$variable, levels = c("Not transcribed", "Unidirectional", "Bidirectional"))
 
  
  p <- p[-which(p$variable == "Not transcribed"),]
  pal_complex <- c(pals$Hits[1], pals$Hits_Darker[1], # non-hit colours
                   pals$Hits[2], pals$Hits_Darker[2]) # hit colours
  p$Threshold <- factor(p$Threshold, levels = as.character(c(2, 3, 5, 10)))
  labels <- c("Inactive candidates / Unidirectional",
              "Inactive candidates / Bidirectional",
              "Functional enhancers / Unidirectional",
              "Functional enhancers / Bidirectional")
  p$Annot <- round(p$value * 100) %>% paste0("%") 
  
  # pdf(file = "Transcriptional categories across thresholds (TTseq).pdf", height = 2.5, width = 3.5)
  pdf_tt(figNo = "SFig5C", title = "eRNA across thresholds", h = 2.5, w = 3.5)
  ggplot(p, aes(x = Hit, fill = interaction(variable, Hit), y = value*100, label = Annot)) +
    geom_col(alpha = 0.8, width = 0.9, position = "stack") +
    theme_bw() +
    facet_wrap(~Threshold, scales = "free_x", nrow = 1, strip.position = "bottom") +
    labs(y = "Percent of peaks transcribed", x = "TTseq minimum read threshold") +
    scale_fill_manual(values = pal_complex, labels = labels) +
    scale_y_continuous(limits = c(0,105), expand = c(0,0)) +
    geom_text(position = "stack", vjust = -0.5, size = 2.5) +
    guides(fill = guide_legend(ncol = 2)) +
    theme(panel.border = invis, axis.line.y = element_line(), 
          legend.position = "none",
          panel.grid = invis, axis.text.x = invis, axis.ticks.x = invis,
          strip.background = element_rect(fill = "white", linetype = 0),
          axis.text.y = text90,
          legend.title = invis)
  dev.off()

    
################################################################################################################################ #
## Correlate TTseq to other variables ----
    
 
# ## Scatterplot transcription versus distance for EGPs
#   ## Setup
#     p <- res.final
#     distCat <- read.csv("../../../FullScale/Results/3_HitEnrichment/EnhGenePairs/Intervening Gene Classification Between EGPs.csv", row.names = 1)
#     
#     m <- match(p$Enh, rownames(tt))
#     p$TT <- tt$TTseq_Total[m]
#     
#     m <- match(p$Enh, transcribed$Enh)
#     p$Transcribed <- transcribed$Category_3[m]  != "Not transcribed"
#     p$Transcribed <- factor(p$Transcribed, levels = c("FALSE", "TRUE"))
#     levels(p$Transcribed) <- c("eRNA -", "eRNA +")
#     
#     m <- match(p$Pair, distCat$Pair)
#     p$DistCat <- distCat$Distance.Category[m]
#     p$DistCat <- factor(p$DistCat)
#   
#     p$HitPermissive <- factor(p$HitPermissive, levels = c("FALSE", "TRUE"))
#     levels(p$HitPermissive) <- c("Inactive EGP", "Functional EGP")
#     
#   ## Plot
#     offset <- 1 # on TTseq counts
# 
#     # pdf(file = "Distance versus TT scatterplot (V2).pdf", height = 1.8, width = 2)
#     pdf_tt(figNo = "3F", title = "Distance versus eRNA", h = 1.8, w = 2)
#     pal <- pals$Primary[c(5,7,4,1)]
#     
#     levels(p$DistCat) <- c("Type A", "Type B", "Type C", "Type D")
#     ggplot(p[which(p$HitPermissive == "Functional EGP"),], aes(y = Gene.Distance / 1000, x = Transcribed, colour = DistCat)) +
#       geom_quasirandom(dodge.width = 0.25) +
#       geom_violin(colour = "black", scale = "width", width = 0.7, draw_quantiles = 0.5, alpha = 0) +
#       scale_colour_manual(values = pal) +
#       theme_bw() +
#       theme(panel.grid = invis, legend.position = "none", #  axis.text.y = text90,
#             panel.border = invis, axis.line.y = element_line(), axis.title.x = invis) +
#       labs(y = "Distance (kb)") +
#       scale_y_continuous(trans = "log10", limits = c(2, 510), 
#                          expand = c(0,0), 
#                          breaks = c(2, 5, 10, 20, 50, 100, 200, 500)) 
#     dev.off()
#     
#   ## Statistics
#     # within hits, what is the relationship between distance and transcription
#     h <- which(p$HitPermissive == "Functional EGP")
#    
#     w <- wilcox.test((p$Gene.Distance[h]) ~ (p$Transcribed[h] )) # p = 0.02
#     sink_tt(figNo = "3F", title = "Distance versus eRNA", w)
#   
################################################################################################################################ #
## FANTOM5 ----
  

## Read in
  f5 <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/FANTOM5/Peak Annotation.csv")
  f5 <- f5[which(f5$Enh %in% res.final$Enh),]
  
  f5 <- f5[which(f5$Enh %in% wellpowered_enh),]
  
## Plot 1: Binary calls at enhancers as independent validation
  p <- apply(f5[,-c(1:2, 5)], 2, function(x) {
    tab <- table(f5$Hit, x > 0) # note this coerces the logical, thus is valid
    tab <- tab / rowSums(tab)
    
    tab <- as.data.frame(tab)
    colnames(tab) <- c("Hit", "FANTOM5", "Freq")
    
    return(tab)
    
  })
  
  p <- do.call("rbind", p)
  p$Sample <- splitter(rownames(p), "\\.", 1)
  p$Sample <- gsub("Used", "", p$Sample) %>% 
    gsub("FANTOM5", "Any Sample", .) %>%
    gsub("Ast", "Astrocytes", .)
  p <- p[-which(p$FANTOM5 == "FALSE"),]
  levels(p$Hit) <- c("Inactive candidates", "Functional enhancers")

  
  # pdf(file = "FANTOM5 - Overlap Barplot.pdf", height = 2, width = 3)
  pdf_tt(figNo = "SFig5D", title = "FANTOM5", h = 2, w = 3)
  ggplot(p[which(p$Sample == "Astrocytes"),], aes(x = Sample, y = Freq*100, fill = Hit, colour = Hit)) +
    geom_col(position = "dodge", width = 0.75) +
    theme_bw() +
    scale_fill_manual(values = pals$Hits) +
    scale_colour_manual(values = pals$Hits_Darker) +
    scale_y_continuous(limits = c(0,30), expand = c(0,0)) +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
          legend.position = "right", legend.title = invis, axis.title.x = invis, axis.text.x = invis,
          axis.ticks.x = invis, plot.title = element_text(size = 6)) +
    labs(y = "Percent with bidirectional\nCAGE peak in FANTOM5", x = "FANTOM5 Sample Set")
  dev.off()
  
  # stat
  f <- table(f5$UsedAst > 0, f5$Hit) %>% fisher.test()
  sink_tt(figNo = "SFig5D", title = "FANTOM5", f)
  
   
   
    
## Plot 2: Consistency between FANTOM5 and TTseq calls, across TTseq thresholds
  ## Read in TTseq data (again)
    # tt <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/TTseq/Results Table.csv", row.names = 1)
    # tt <- tt[which(rownames(tt) %in% res.final$Enh),]
    # x <- data.frame(f5, Pos = tt$TTseq_Pos, Neg = tt$TTseq_Neg)
    # x$UsedAst <- x$UsedAst > 0
    # x$UsedGBM <- x$UsedGBM > 0
    
    

## Pooled categorisation plot
  # get data for TTseq category
  p <- data.frame(Enh = transcribed$Enh,
                  TTseq = transcribed$Category_3,
                  Hit = transcribed$Hit)
  
  p <- p[which(p$Enh %in% wellpowered_enh),]
  
  # add in FANTOM5 calls
  m <- match(p$Enh, f5$Enh)
  p$FANTOM5 <- f5$UsedAst > 0
  p$FANTOM5 <- factor(p$FANTOM5)
  levels(p$FANTOM5) <- c("CAGE-", "CAGE+")
  
  # aside: fisher test comparing fantom5 to ttseq
  fish <- table(p$TTseq != "Not transcribed", p$FANTOM5) %>% fisher.test()
  sink_tt(figNo = "SFig6E", title = "TTseq vs FANTOM5", fish)
  
  # stratify FANTOM5 calls by transcription
  tab <- table(p$TTseq, p$FANTOM5)
  tab <- proportions(tab, 2)
  tab <- as.data.frame(tab)
  tab$Var1 <- factor(tab$Var1, levels = levels(tab$Var1)[c(2,3,1)])
  colnames(tab) <- c("TTseq", "CAGE", "Freq")
  pal <- c("grey90", pals$Primary[2], pals$Primary_Darker[2])
  
  # pdf(file = "FANTOM5 - Proportions of TTseq Calls.pdf", height = 2.5, width = 4)
  # ggplot(tab, aes(x = CAGE, y = Freq*100, fill = TTseq, linetype = TTseq)) +
  #   geom_col(colour = "black") +
  #   theme_bw() +
  #   scale_linetype_manual(values = c("blank", "dashed", "solid")) +
  #   theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(), axis.text.y = text90) +
  #   scale_fill_manual(values = pal) +
  #   labs(x = "Candidate overlaps FANTOM5 CAGE peak", y = "Percent") +
  #   scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100), limits = c(0,105))
  # dev.off()  
  
  # a combined annotation for CAGE/TTseq
  p$Cat <- paste0(p$TTseq, "_", p$FANTOM5)
  p$Cat <- gsub("Uni|Bi", "", p$Cat) %>%
    gsub("Not transcribed", "TTseq-", .) %>%
    gsub("directional", "TTseq+", .) %>%
    gsub("_", " / ", .)
  p$Cat <- factor(p$Cat)
  levels(p$Cat) <- c("Neither", "CAGE+", "TTseq+", "CAGE+ / TTseq+")
  
  # pal <- c("grey95", pals$Primary[c(2,5)], pals$Primary_Darker[2])
  
  tab <- table(p$Cat, p$Hit)
  tab <- proportions(tab, 2)
  tab <- as.data.frame(tab)
  # tab$Var1 <- factor(tab$Var1, levels = levels(tab$Var1)[c(2,3,1)])
  # colnames(tab) <- c("TTseq", "CAGE", "Freq")
  # pal <- c("grey90", pals$Primary[2], pals$Primary_Darker[2])
  tab$Var2 <- factor(tab$Var2)
  levels(tab$Var2) <- c("Inactive\ncandidates", "Functional\nenhancers")
  
  pal <- c("grey95", pals$Primary[c(4,5,7)])
  
  # pdf("FANTOM5 - Combined category.pdf", height = 3, width = 3.8)
  pdf_tt(figNo = "SFig5E", title = "FANTOM5 vs TTseq", h = 3, w = 3.8)
  ggplot(tab, aes(x = Var2, fill = Var1, y = Freq*100)) +
    geom_col(colour = "black", width = 0.7) +
    theme_bw() +
    # scale_linetype_manual(values = c("blank", "dashed", "solid")) +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(), axis.text.y = text90,
          axis.title.x = invis, legend.title = invis) +
    scale_fill_manual(values = pal) +
    labs(y = "Percent of candidates") +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100), limits = c(0, 105))
  dev.off()
  
  