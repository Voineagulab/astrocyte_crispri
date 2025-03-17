#This script generates figures 4A and 5A
## Setup
  setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_RelevanceToHumanBiology/")
  library(tidyverse)
  library(rcartocolor)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  
  invis <- element_blank()
  maxh <- 29.7 / 2.54
  maxw <- 21.0 / 2.54
  options(stringsAsFactors = FALSE)
  source("../../../FullScale/Scripts/Functions.R")
  source("../FinalFigureFunctions.R")
  
## Read in power dataframe
  pow <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable3_DE/3F_Power.csv")
  wellpowered_egps <- pow$Pair[which(pow$WellPowered015)]
  wellpowered_enh <- pow$Enhancer[which(pow$WellPowered015)] %>% unique()
  wellpowered_genes <- pow$Gene[which(pow$WellPowered015)] %>% unique()

## Read in results dataframe
  res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  res.final <- res.final[which(res.final$Pair %in% wellpowered_egps),]
  gene.hit <- unique(res.final$Gene[which(res.final$HitPermissive)])
  enh.hit <- unique(res.final$Enh[which(res.final$HitPermissive)])
  gene.bg <- unique(res.final$Gene)
  
## Plotting
  cols.perm <- pal_lancet()(2)
  
  invis <- element_blank()
  maxh <- 29.7 / 2.54
  maxw <- 21.0 / 2.54
  
  theme_gjs <- theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) 
  NoExpand <- scale_y_continuous(expand = c(0,0))
  
## Functions to write to disk in a traceable way
  pdf_biol <- function(figNo, title, h, w) {
    pdf(file = paste0("../Final/", figNo, " - Script HumBiol - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_biol <- function(figNo, title, toPrint) {
    sink(paste0("../Final/", figNo, " - Script HumBiol - ", title, ".txt"))
    print(toPrint)
    sink()
  }
  
  # sink_qc(figNo = "SFig10", title = "scRNAseq qc", qc)

  

################################################################################################################################ #
## A combined plot for functional annotations ----
    
  
## Setup GO
  # p_go <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/GO - Clusterprofiler.csv", row.names = 1)
  p_go <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Wellpowered/GO - Clusterprofiler.csv", row.names = 1)
  
  # filter
  keepGO <- c("tissue development", "response to external stimulus",
              "chemotaxis", "regulation of cell migration", "apoptotic process",
              "response to wounding", "response to growth factor",  
              "regulation of response to stress", "cell surface receptor signaling pathway",
              "response to cytokine")
  
  p_go <- p_go[which(p_go$Description %in% keepGO),c("ONTOLOGY", "Description", "GeneRatio", "BgRatio", "pvalue", "Count")]
  
  # calculate odds ratio (note: this retains the original P)
  p_go$OR <- apply(p_go, 1, function(y) {
    # columns
    h <- grep("GeneRatio", colnames(p_go))
    bg <- grep("BgRatio", colnames(p_go))
    
    # numbers for fisher test
    mx <- matrix(nrow = 2, ncol = 2)
    
    a <- splitter(y[h], "/", 1) %>% as.numeric()
    b <- splitter(y[h], "/", 2) %>% as.numeric()
    c <- splitter(y[bg], "/", 1) %>% as.numeric()
    d <- splitter(y[bg], "/", 2) %>% as.numeric()
    
    # convert to matrix
    mx[1,1] <- d - c # not hit, and not in GO category
    mx[1,2] <- b - a # hit, and not in GO category
    mx[2,1] <- c # not hit, and in GO category
    mx[2,2] <- a # hit, and in GO category
    print(mx)
    # test
    fish <- fisher.test(mx)
    or <- fish$estimate
    return(or)
    
    
  })
  
  p_go <- data.frame(Type = "GO",
                     Name = p_go$Description,
                     P = p_go$pvalue,
                     Count = p_go$Count,
                     OR = p_go$OR)
  
  

## Read in other categories
  # use combined data for ageing, maturation, and markers
  load("../../../FullScale/Results/3_HitEnrichment/Genes/Wellpowered/Final.rda", verbose = TRUE)
  
  p_ast <- enrichments.combined
  p_ast <- p_ast[c(1,4,5,6,9),]
  p_ast <- data.frame(Type = "DE",
                      Name = p_ast$Resource,
                      P = p_ast$p,
                      Count = p_ast$Total_Hit_TRUE,
                      OR = p_ast$OR)
  p_ast$Name <- gsub("Ast|_HumanOnly", "", p_ast$Name) %>%
    gsub("Marker", "Astrocyte marker", .)
  
  
  
## Combine
  p <- rbind(p_go, p_ast)
  p$Name <- str_to_sentence(p$Name)
 
  p$Name <- gsub("Regulation of response ", "Regulation of response\n", p$Name) %>%
    gsub("Cell surface receptor ", "Cell surface receptor\n", .) %>%
    gsub("Regulation of cell ", "Regulation of\ncell ", .) %>%
    gsub("Response to ", "Rsp. to ", .) %>%
    # gsub("stimulus ", "Rsp. to", .) %>%
    gsub("external", "ext.", .) %>%
    gsub("stimulus", "stim.", .) %>%
    gsub("Regulation", "Reg.", .) %>%
    gsub("Reg. of response\nto stress", "Reg. of rsp. to stress", .)
  
  p$Type <- factor(p$Type, levels = c("GO", "DE"))
  levels(p$Type) <- c("Gene ontology", "Functional sets")
  p$mlog10p <- -log10(p$P)
  # p$pBin <- cut(p$mlog10p, c(0, -log10(0.05), 3, 5, 100))
  p$pBin <- cut(p$P, c(1, 0.05, 1e-3, 1e-5, 0))
  levels(p$pBin) <- c("1e-5", "1e-3", "0.05", "ns")
  
  
## Plot
  p$Name <- factor(p$Name, levels = p$Name[order(p$OR)])
  
  # pdf(file = "Genes/Latest annotation (V3 - OR Bins).pdf", height = 6, width = 2.6)
  pdf_biol(figNo = "4A", title = "Functional annotation", h = 4.8, w = 2.6)
  
  ggplot(p, aes(y = Name, x = OR, size = Count, fill = pBin)) +
    geom_segment(mapping = aes(xend = OR, x = 1, yend = Name), size = 0.2, colour = pals$One) +
    geom_point(colour = "black", shape = 21) +
    scale_fill_manual(values = pals$grn2orng[c(1,3,5)]) +
    # facet_wrap(~Type, ncol = 1, scales = "free_y", strip.position = "left") +
    facet_grid(Type~., scales = "free_y", space = "free_y", switch = "y") +
    scale_size_continuous(trans = "log2") +
    guides(fill = guide_legend(title = "P-value\nthreshold", nrow = 3), size = guide_legend(ncol = 1, title.position = "top")) +
    scale_x_continuous(limits = c(0.1, 11), expand = c(0,0), trans = "log2", 
                       breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8), labels = c("0.125", "", "0.5", "1", "2", "4", "8")) +
    labs(x = "Odds ratio") +
    theme_bw() +
    theme(panel.border = invis, axis.line.x = element_line(),
          panel.grid = invis, legend.position = "bottom", # legend.position = "bottom",
          axis.text.y = element_text(size = 7, lineheight = 0.75),
          strip.text = element_text(size = 8), 
          legend.box = "vertical", axis.title.y = invis) +
    geom_vline(xintercept = 1, linetype = 2, alpha = 0.5, colour = "grey80")
  dev.off()
  

################################################################################################################################ #
## A combined plot for disease ----
  
## The aim of this plot is to summarise disgenet and disease associations
  # for disgenet: the same categories
  # for disease: de
  
## Process disgenet data
  disg <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Wellpowered/Disgenet - Enrichment.csv", row.names = 1)

  # filter
  disg <- disg[which(disg$FDR < 0.05),]
  
  # adjust the "total" columns to be total - hits
  disg$Total_count <- disg$Total_count - disg$Hit_count
  disg$Total_ratio <- disg$Total_count / length(unique(res.final$Gene))
  
  # clean
  disg$Disease_Type <- str_to_sentence(disg$Disease_Type)
  disg$Disease <- str_to_sentence(disg$Disease)
  
  # load variable disg above
  # keep <- c(4,5,21,64,115,250,256,281,296,408)
  #p1 <- disg[keep,]
  keep <-  c("Malignant glioma","Glioma","Neuroblastoma","Astrocytoma","Alzheimer's disease","Multiple sclerosis","Amyotrophic lateral sclerosis 1", "Depressive disorder","Cerebral palsy","Parkinson disease") 
  p1 <- disg[which(disg$Disease %in% keep),]
  
  # recategorise
  p1$Cat <- "Disgenet\nNeoplasms"
  p1$Cat[-which(p1$Disease_SemanticType == "Neoplastic Process")] <- "Disgenet\nCNS Diseases"
  
  # filter columns
  p1 <- p1[,c("Disease", "p", "OR", "Hit_count", "Cat")]
  
  
## Process DE data
  # read in
  p2 <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Wellpowered/Final - Enrichments.csv", row.names = 1)
  
  # filter
  p2 <- p2[which(p2$Signed == "Unsigned" & grepl("Disease", p2$Resource)),]
  
  # clean
  # p2$Resource <- gsub("AstDisease_", "", p2$Resource) %>%
  #   gsub("_", "\n", .) %>%
  #   gsub("20", " 20", .) %>%
  #   gsub("FocalCorticalDysplasia", "FCD", .) %>%
  #   gsub("20 20", "2020", .)
  
  p2$Resource <- gsub("AstDisease_", "", p2$Resource) %>%
    gsub("_", " (", .) %>%
    gsub("$", ")", .) %>%
    gsub("20", " 20", .) %>%
    gsub("FocalCorticalDysplasia", "FCD", .) %>%
    gsub("20 20", "2020", .)
  
  p2$Cat <- "Astrocyte\nDE resources"
  
  # filter columns
  p2 <- p2[,c("Resource", "p", "OR", "Total_Hit_TRUE", "Cat")]
  colnames(p2)[c(1,4)] <- c("Disease", "Hit_count")
  
## Get ROSMAP data
  p3 <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.Ast.FisherTests.csv")
  p3 <- p3[53,]
  p3 <- data.frame(Disease = "AD (ROSMAP 2023)",
                   p = p3$loc.fisher.p,
                   OR = p3$loc.fisher.or,
                   Hit_count = p3$loc.Nintersect,
                   Cat = "Astrocyte\nDE resources")
  
## Combine
  p <- rbind(p1, p2, p3)  
  
  p$OR_bin <- cut(p$OR, c(0,1,2,3,5,10), include.lowest = TRUE)
  levels(p$OR_bin) <- c("<1", "1-2", "2-3", "3-5", "5+")
  
  colnames(p)[4] <- "Count"
  p <- p[-grep("Habib", p$Disease),]
  p <- p[-which(p$Count == 0),]
  
  # p$Disease <- gsub("lateral", "lateral\n", p$Disease)
  p$Disease <- factor(p$Disease, levels = p$Disease[order(-p$p)])
  p$Cat <- factor(p$Cat, levels = unique(p$Cat)[c(2,1,3)])
  
  levels(p$Disease) <- gsub(" 1", "", levels(p$Disease))

## Plot
  common_size <- scale_size_continuous(trans = "log2", breaks = c(4, 8, 16, 32))
  common_x <- scale_x_continuous(limits = c(0, 10.5), expand = c(0,0), breaks = c(0, 1, 3, 5, 7, 9)) 
  common_lab <- labs(x = "-log10(P)") 
  common_theme <- theme_bw(base_size = 7)
  common_vline <- geom_vline(xintercept = -log10(0.01), linetype = 2, alpha = 1, colour = "grey80")
  pal <- c(rev(pals$grn2orng[c(1,2,3,5)])) 
  
  # pdf(file = "Genes/Disease annotation (V3, Wide).pdf", height = 2, width = 3)
  pdf_biol(figNo = "5A", title = "Disease annotation (wellpowered)", h = 2, w = 3)
  p1 <- ggplot(p[-grep("Disgenet", p$Cat),], aes(y = Disease, x = -log10(p), size = Count, fill = OR_bin)) +
    common_size + common_x + common_lab + common_theme + common_vline +
    geom_segment(mapping = aes(xend = -log10(p), x = 0, yend = Disease), size = 0.2, colour = pals$One) +
    geom_point(colour = "black", shape = 21) +
    scale_fill_manual(values = pal) +
    facet_grid(Cat~., scales = "free", space = "free", switch = "y", drop = TRUE) +
    guides(fill = guide_legend(title = "Odds ratio", nrow = 2)) +
    guides(size = guide_legend(nrow = 2)) +
    theme(panel.border = invis, axis.line.x = element_line(),
          panel.grid = invis,
          legend.position = "right",
          # legend.position = c(0.75, 0.15),
          legend.box = "horizontal", axis.title.y = invis) 
  
  p2 <- ggplot(p[grep("Disgenet", p$Cat),], aes(y = Disease, x = -log10(p), size = Count, fill = OR_bin)) +
    common_size + common_x + common_lab + common_theme + common_vline +
    geom_segment(mapping = aes(xend = -log10(p), x = 0, yend = Disease), size = 0.2, colour = pals$One) +
    geom_point(colour = "black", shape = 21) +
    scale_fill_manual(values = pal) +
    facet_grid(Cat~., scales = "free", space = "free", switch = "y", drop = TRUE) +
    
    guides(fill = guide_legend(title = "Odds ratio", nrow = 2)) +
    guides(size = guide_legend(nrow = 2)) +
    theme(panel.border = invis, axis.line.x = element_line(),
          panel.grid = invis, 
          legend.position = "none",
          # legend.position = c(0.75, 0.15), 
          legend.box = "horizontal", axis.title.y = invis) 
  
  # output plots
  p3 <- get_legend(p1)
  p1 + theme(legend.position = "none")
  p2
  grid.newpage()
  grid.draw(p3)
  # print(p3)


    
  dev.off()
       
# 
# ################################################################################################################################ #
# ## Gene regulation within the astrocytes ----
#     
# ## Dotplot of signed enrichment within individual DE resources
#   e <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Final - Enrichments.csv", row.names = 1)
#   
#   e$p_bin <- cut(e$p, c(0, 1e-5, 1e-3, 1e-2, 1))
#   levels(e$p_bin) <- c("< 1e-5", "< 1e-3", "< 0.01", "ns")
# 
#   e$or_bin <- cut(e$OR, c(0, 1, 2, 3, 5, 100))
#   levels(e$or_bin) <- c("<1", "1-2", "2-3", "3-5", "5+")
# 
#   p_colours <- or_colours <- carto_pal(7, "Geyser")[(c(7,6,4,1))]
#   names(p_colours) <- levels(e$p_bin)
#   or_colours <- c("grey50", rev(or_colours))
#   names(or_colours) <- levels(e$or_bin)
#   
#   e$Class <- splitter(e$Resource, "_", 1)
#   
#   e <- e[grepl("Ageing|Maturation|Markers", e$Class),]
#   e$Resource <- gsub(paste(e$Class, collapse = "|"), "", e$Resource) %>%
#     gsub("^_", "", .) %>%
#     gsub("20", " 20", .) %>%
#     gsub("_", "\n", .)
#   g <- grep("Maturation|Ageing", e$Class)
#   e$Resource[g] <- paste0(e$Class[g], "\n", e$Resource[g]) %>% gsub("Ast", "", .)
#   e$Resource <- gsub("Maturation\n|Ageing\n", "", e$Resource)
#   e$Resource <- gsub("AllAges", "Across development", e$Resource)
#   e$Class <- gsub("Ast", "", e$Class)
# 
#   # reorder
#   e$Class <- factor(e$Class, levels = c("Ageing", "Maturation", "Markers"))
#   
#   res.order <- aggregate(p~Resource, data = e, FUN = min) 
#   res.order <- res.order$Resource[order(res.order$p)]
#   e$Resource <- factor(e$Resource, levels = res.order)
#   
#   # pdf(file = "Genes/Signed Dotplot Enrichment (v2).pdf", height = 4, width = 7)
#   pdf_biol(figNo = "SFig7A", title = "Ageing, maturation, markers", h = 4, w = 7)
#   ggplot(e, aes(x = Resource, y = -log10(p), yend = -log10(p), xend = Resource, fill = or_bin, shape = Signed, colour = or_bin)) +
#     geom_segment(aes(y = 0), colour = "black", size = 0.2, linetype = 2, alpha = 0.5) +
#     geom_point(size = 2) + 
#     geom_point(colour = "black", size = 2.5, show.legend = FALSE) + 
#     
#     scale_shape_manual(values = c(25,21,24)) +
#     facet_grid(. ~ Class, scales = "free_x", space = "free", switch = "x") +
#     # facet_wrap(~Class, nrow = 1, scales = "free_x", strip.position = "left") +
#     scale_fill_manual(values = or_colours, limits = names(or_colours)) +
#     scale_colour_manual(values = or_colours, limits = names(or_colours)) +
#     theme_bw() +
#     theme_gjs +
#     guides(fill = guide_legend(title = "Odds Ratio", nrow = 2), 
#            colour = guide_legend(title = "Odds Ratio", nrow = 2),
#            shape = guide_legend(ncol = 1)) +
#     scale_y_continuous(limits = c(0, 7.5), expand = c(0,0)) +
#     geom_hline(yintercept = 2, linetype = 2) +
#     theme(axis.title.x = invis, legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.box = "vertical", panel.grid = invis) +
#     labs(y = "-log10(P)")
#     
#   dev.off()
#   
#   
# ## On astrocyte activation
#   f <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Final - Enrichments.csv", row.names = 1)
#   f$or_bin <- cut(f$OR, c(0, 1, 1.5, 2, 3, 100))
#   levels(f$or_bin) <- c("<1", "1-1.5", "1.5-2", "2-3", "3+")
# 
#   or_colours <- rev(c(carto_pal(7, "Geyser")[(c(7,6,3,1))], "grey50"))
#   names(or_colours) <- levels(f$or_bin)
#   
#   f <- f[grepl("Activation", f$Resource),]
#   f <- f[grepl("hiPSC", f$Resource),]
#   f$Class <- NA
#   f$Class[grep("Rep", f$Resource)] <- "Replication N"
#   f$Class[grep("IRAS", f$Resource)] <- "Response Subtype"
#   f$Class[which(is.na(f$Class))] <- "hiPSC Protocol"
#   # f <- f[-grep("Velm|Focal", f$Resource),]
#   f$Resource <- gsub("AstActivation", "", f$Resource) %>%
#     gsub("^_", "", .) %>%
#     gsub("TIC_hiPSC_", "", .) %>%
#     gsub("Leng", "iAstrocyte", .) %>%
#     gsub("_", "\n", .)
#   f$Resource[grep("Rep", f$Resource)] <- gsub("Rep", "", f$Resource[grep("Rep", f$Resource)]) %>%
#     paste0(., " hiPSC line(s)")
#   
#   # f$Resource <- gsub("IRAS1", "IL-1/IL-6-responsive\nSubtype Markers", f$Resource)
#   # f$Resource <- gsub("IRAS2", "TNF/IFN-responsive\nSubtype Markers", f$Resource)
#   
#   f <- f[-which(f$Class == "Response Subtype"),]
#   
#   
#   # plot 
#   # pdf(file = "Genes/Signed Dotplot Enrichment (Activation).pdf", height = 2.8, width = 7.5)
#   pdf_biol(figNo = "SFig7B", title = "Astrocyte activation", h = 2.8, w = 7.5)
#   ggplot(f, aes(x = Resource, y = -log10(p), yend = -log10(p), xend = Resource, colour = or_bin, fill = or_bin, shape = Signed)) +
#     geom_segment(aes(y = 0), colour = "black", size = 0.2, linetype = 2, alpha = 0.5) +
#     # geom_point(shape = 21, colour = "black") +
#     geom_point(size = 2.5) +
#     scale_shape_manual(values = c(25,21,24)) +
#     scale_size_continuous(limits = c(0,85), breaks = c(0, 10, 25, 40, 60, 85)) +
#     facet_grid(. ~ Class, scales = "free_x", space = "free", switch = "x") +
#     # facet_wrap(~Class, nrow = 1, scales = "free_x", strip.position = "left") +
#     scale_fill_manual(values = or_colours, limits = names(or_colours)) +
#     scale_colour_manual(values = or_colours, limits = names(or_colours)) +
#     theme_bw() +
#     theme_gjs +
#     guides(fill = guide_legend(title = "Odds Ratio", nrow = 2), size = guide_legend(title = "% Hits in Set", ncol = 2),
#            shape = guide_legend(ncol = 1), colour = guide_legend(title = "Odds Ratio", nrow = 2)) +
#     scale_y_continuous(limits = c(0, 9), expand = c(0,0)) +
#     # scale_y_continuous(limits = c(0.5, NA), expand = c(0,0), trans = "log2", breaks = c(0.5, 1, 2, 4, 8, 16)) +
#     geom_hline(yintercept = 2, linetype = 2) +
#     theme(axis.title.x = invis, legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.box = "vertical", panel.grid.minor.y = invis, panel.grid = invis) +
#     labs(y = "-log10(P)")
#   dev.off()
#   
#     
#   
