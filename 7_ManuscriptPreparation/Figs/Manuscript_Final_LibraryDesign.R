#This script generates figures 1D, 4F, S2A, S2B, S3A and S3B.
## Setup
  setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_LibraryDesign/")
  source("../../../FullScale/Scripts/Functions.R")
  source("../FinalFigureFunctions.R")
  
  options(stringsAsFactors = FALSE)
  library(tidyverse)
  library(rcartocolor)
  library(ggplot2)
  library(ggdendro)
  library(egg)
  invis <- element_blank()
  maxh <- 29.7 / 2.54
  maxw <- 21.0 / 2.54

## Plotting function
  ggDendro_gjs <- function(clust) {
    x <- as.dendrogram(clust) %>% 
      dendro_data(type = "rectangle") %>%
      segment()
    
    ggplot(x) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      theme_void() +
      scale_x_continuous(expand = c(0,0), limits = c(0.5, max(clust$order) + 0.5)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, max(x$y) * 1.05)) 
  }
  

## Functions to write to disk in a traceable way
  pdf_LibDesign <- function(figNo, title, h, w) {
    pdf(file = paste0("../Final/", figNo, " - Script LibraryDesign - ", title, ".pdf"), height = h, width = w)
  }
  
  sink_LibDesign <- function(figNo, title, toPrint) {
    sink(paste0("../Final/", figNo, " - Script LibraryDesign - ", title, ".txt"))
    print(toPrint)
    sink()
  }
  
################################################################################################################################ #
## Characteristics of our chosen candidate enhancers ----
  
  guides <- read.csv("../../../FullScale/Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[,-which(colnames(guides) == "Celltype")] # extraneous
  

## Enhancer length distribution, split by number of guides
  e <- guides[which(guides$TargetCat == "Enh"),]
  e <- e[-which(duplicated(e$TargetID)),]
  p <- strsplit(e$TargetCoord, ":") %>% sapply("[", 2)
  q <- strsplit(p, "-") %>% sapply(., function(x) as.numeric(x[[2]]) - as.numeric(x[[1]]))
  
  q <- data.frame(q)
  q$Bin <- cut(q$q, seq(0, 1200, 200), dig.lab = 10)
  levels(q$Bin) <- gsub("\\(", "", levels(q$Bin)) %>% gsub("]", "", .)  %>% gsub(",", "-", .)
  rownames(q) <- e$TargetID
  tab <- table(guides$TargetID)
  m <- match(rownames(q), names(tab))
  q$N <- as.numeric(tab[m])
  q$N[which(q$N >= 5)] <- "5+"
  
  pal <- pals$Primary[c(3,4,7,8)]
 
  # pdf(file = "Enhancer Width Distribution, With nGuides (Wider Bins).pdf", height = 2.2, width = 3.1)
  pdf_LibDesign(figNo = "ExtFig3A", title = "Candidate length histogram", h = 2.2, w = 3.1)
  ggplot(q, aes(x = Bin, fill = N, colour = N)) +
    geom_bar(position = "stack") +
    theme_bw() +
    # theme_resizeText +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 200, 400)) +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal) +
    
    guides(fill = guide_legend(title = "sgRNAs", ncol = 2), colour = "none") +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis, axis.ticks.x = invis,
          legend.position = c(0.8, 0.75), axis.text.x = element_text(angle = 30, hjust = 0.9),
          axis.text.y = text90, legend.key.size = unit(0.5, "cm"),
          legend.title = element_text(size = 10)) +
    labs(y = "Number of candidates", x = "Candidate enhancer length (bp)")
  dev.off()

  
  
## Plot 3: number of final number of tested genes per enhancer
  res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  p <- factor(res.final$Enh, levels = unique(guides$TargetID[which(guides$TargetCat == "Enh")]))
  p <- table(p) %>% as.data.frame()
  p$Factor <- factor(p$Freq, levels = c(0:max(p$Freq)))
  p$Cut <- cut(p$Freq, c(-1:9, 11, 13, 15, 17, 19, 24, 29, 39, 49), include.lowest = TRUE)
  levels(p$Cut) <- c(0:9, "10-11", "12-13", "14-15", "16-17", "18-19", "20-24", "25-29", "30-39", "40-49")
  
  # pdf(file = "Genes (Tested) Per Enhancer.pdf", height = 2.1, width = 3.8)
  pdf_LibDesign(figNo = "ExtFig3B", title = "Genes tested per candidate", h = 2.1, w = 3.8)
  ggplot(p, aes(x = Cut)) +
    geom_bar(fill = pals$One, width = 0.8) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis, axis.ticks.x = invis,
          axis.text.y = text90,
          axis.text.x = element_text(angle = 90, hjust = 0.9, vjust = 0.5)) +
    # scale_fill_manual(values = pal_iv_single) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Number of tested genes per candidate enhancer", y = "Number of candidates")
  dev.off()
        
      
################################################################################################################################ #
## Upset plot of peak annotation ----
  
  
## Load data
  atac <- read.csv("../../../FullLibrary_Selection/Results/Peaks_Annotated.csv") # annotation for each of the MACS2 peak calls
  
## Plot
  p <- atac[, c("GENCODE32", "TAD", "PE", "high.cpm", "nRep")]
  p$GENCODE32 <- !(p$GENCODE32) # converts to intergenic for TRUE/FALSE rather than genic
  p$nRep <- p$nRep >= 2 # binarise for >= 2 
  colnames(p) <- c("Intergenic", "TAD", "PsychENCODE", "CPM", "Replicated 2+")
  
  p$List <- apply(p, 1, function(x) {
    if (any(x)) {
      names(x)[which(x)]  
    } else {
      x <- "None"
    }
    
  })
  
  # plot
  library(ggupset)
  # pdf(file = "Candidate upset.pdf", height = 3.5, width = 7)
  pdf_LibDesign(figNo = "ExtFig2A", title = "Candidate annotation upset", h = 3.5, w = 7)
  ggplot(p, aes(x = List)) +
    geom_bar(fill = pals$One, colour = pals$One) +
    geom_text(stat = "Count", aes(label = after_stat(count)), vjust = 0.5, hjust = -0.2, size = 2.5,
              angle = 90) +
    scale_x_upset(order_by = "degree", reverse = TRUE) +
    theme_bw() +
    scale_y_continuous(position = "left", expand = expansion(mult = c(0, 0.2)),
                       labels = function(x) x / 1000) +
    labs(y = "\n\n\nNumber of peaks (1000's)", x = "Set") +
    theme(panel.grid = invis, panel.border = invis, axis.line = element_line(),
          legend.position = c(0.7, 0.7)) 
  dev.off()
  
## Barplot of totals
  p <- p[,1:5]
  p$Coord <- atac$id
  p <- melt(p, id.vars = "Coord")
  p <- p[which(p$value),]

  # pdf(file = "Candidate total.pdf", height = 2.5, width = 4)
  pdf_LibDesign(figNo = "ExtFig2B", title = "Candidate annotation barplot", h = 2.5, w = 4)
  ggplot(p, aes(x = variable)) +
    geom_bar(fill = pals$One, colour = pals$One, width = 0.7) +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis,
          axis.title.x = invis, axis.line.y = element_line()) +
    labs(y = "Number of peaks (1000's)") +
    scale_y_continuous(labels = function(x) x / 1000, expand = c(0,0))
  dev.off()
    
    
  
  
################################################################################################################################ #
## Validation of the enhancer set ----
  
  
## Using ENCODE data

  ## Get data
    encode <- read.csv("../../../FullLibrary_Selection/Results/Candidate_Validation/Annotations.csv", row.names = 1)

  ## Run stats
    g <- grep("ENCODE", colnames(encode))
    
    x <- apply(encode[,g], 2, table) %>% t() %>% as.data.frame
    colnames(x) <- c("NoOverlap", "Overlap")

    x$OverlapInCandidates <- apply(encode[which(encode$Candidate),g], 2, function(y) length(which(y)))
    x$FisherCanditates <- apply(encode[,g], 2, function(y) fisher.test(table(y, encode$Candidate))$estimate)
    x$FisherCanditatesP <- apply(encode[,g], 2, function(y) fisher.test(table(y, encode$Candidate))$p.value)

    sink_LibDesign(figNo = "1B", title = "Encode annotation of candidates", x)
    
    # write.csv(x, file = "ENCODE - Stats.csv")
    
  ## Plot
    p <- data.frame(Dataset = rownames(x),
                    All.NHA.OCR = x$Overlap / nrow(encode),
                    Candidate.OCR = x$OverlapInCandidates / length(which(encode$Candidate)))
  
    p <- melt(p)
    # p$Dataset <- factor(p$Dataset, levels = levels(as.factor(p$Dataset))[c(2,1,3,4)])
    p$Dataset <- sub("ENCODE_", "", p$Dataset)
    
    p <- p[-which(p$Dataset %in% c("CTCF", "K4m3")),]
    p$Dataset <- factor(p$Dataset)
    levels(p$Dataset) <- c("Distal\nEnhancer", "Proximal\nEnhancer", "Promoter")
    
    levels(p$variable) <- c("All peaks\n(260,193)", "Candidate\nenhancers\n(979)")
    
    # pdf(file = "ENCODE - Barplot V4.1.pdf", height = 2.5, width = 4)
    pdf_LibDesign(figNo = "1B", title = "Encode annotation of candidates", h = 2.5, w = 4)
    pal <- rev(pals$Primary[1:3])
    ggplot(p, aes(x = variable, y = value*100, fill = Dataset, colour = Dataset)) +
      geom_col(position = "dodge", width = 0.65) +
      theme_bw() +
      scale_fill_manual(values = pal) +
      scale_colour_manual(values = pal) +
  
      coord_flip() +
      guides(fill = guide_legend(nrow = 1), colour = FALSE) +
      theme(axis.title.y = invis, 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = invis, panel.border = invis, axis.line.x = element_line(), 
            legend.position = "bottom", legend.title = invis) +
      scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
      labs(y = "Percent intersecting ENCODE annotation")
    dev.off()
    
  
  
## Plot 3: Nott. 
  ## Get data
    nott <- read.csv("../../../FullLibrary_Selection/Results/Candidate_Validation/Annotations.csv", row.names = 1)
    
    # filter to intergenic
    intergenic <- read.csv("../../../FullLibrary_Selection/Results/Peaks_Annotated.csv")
    intergenic <- intergenic[which(!(intergenic$GENCODE32)),]
    nott <- nott[which(nott$peak %in% intergenic$id),]
    
    # filter to nott only, and remove H3K4me3
    nott <- nott[,-grep("H3K4me3|Hnisz|ENCODE", colnames(nott))]
    
     g <- grep("Nott", colnames(nott))
    
    x <- apply(nott[,g], 2, table) %>% t() %>% as.data.frame
    colnames(x) <- c("NoOverlap", "Overlap")

    x$OverlapInCandidates <- apply(nott[which(nott$Candidate),g], 2, function(y) length(which(y)))
    x$FisherCanditates <- apply(nott[,g], 2, function(y) fisher.test(table(y, nott$Candidate))$estimate)
    x$FisherCanditatesP <- apply(nott[,g], 2, function(y) fisher.test(table(y, nott$Candidate))$p.value)
    
    # save
    sink_LibDesign(figNo = "SFig2C", title = "Nott annotation of candidates", x[grep("LHX2", rownames(x)),])
    
  ## Plot
    p <- data.frame(Dataset = rownames(x),
                    All.NHA.OCR = x$Overlap / nrow(nott),
                    Candidate.OCR = x$OverlapInCandidates / length(which(nott$Candidate)))
  
    p <- melt(p)
    p$Dataset <- sub("Nott_", "", p$Dataset)
    
    p <- p[-grep("Olig2|PU1|NeuN", p$Dataset),] # filters out other celltypes
    p$Dataset <- factor(p$Dataset)
    levels(p$Dataset) <- c("Astrocyte ATAC\n(Nott 2019)", "Astrocyte H3K27ac\n(Nott 2019)")
    
    levels(p$variable) <- c("Intergenic\nPeaks\n(64,519)", "Candidate\nEnhancers\n(979)")
  
    # pdf(file = "Nott - Barplot.pdf", height = 2.5, width = 3.5)
    pdf_LibDesign(figNo = "SFig2C", title = "Nott annotation of candidates", h = 2.5, w = 3.5)
    pal <- pals$Primary[1:2]
    ggplot(p, aes(x = variable, y = value, fill = Dataset, colour = Dataset)) +
      geom_col(position = "dodge", width = 0.7) +
      theme_bw() +
      scale_fill_manual(values = pal) +
      scale_colour_manual(values = pal) +
      coord_flip() +
      guides(fill = guide_legend(nrow = 1)) +
      theme(axis.title.y = invis, 
            panel.grid = invis, panel.border = invis, axis.line.x = element_line(), legend.position = c("bottom"), legend.title = invis) +
      scale_y_continuous(expand = c(0,0), limits = c(0,1.05)) +
      labs(y = "Fraction intersecting holdout annotation")
    dev.off()
    
################################################################################################################################ #
## Enhancer trends in the ENCODE V3 DHS dataset ----
    
## Variables will refer to this as "Altius" due to the institution in which it was generated


## Read in
  load("../../../FullScale/Results/3_HitEnrichment/Chromatin/Altius - Pooled Expression Matrices.rda", verbose = TRUE)
  p <- altius_pool$Mean
  
  altius_stats <- read.csv("../../../FullScale/Results/3_HitEnrichment/Chromatin/Altius - CV vs Onrate.csv", row.names = 1)
  
  
  m <- melt(p)
  ggplot(m, aes(y = variable, x = value)) +
    geom_boxplot() 
  
  m <- apply(p, 2, median) %>% sort()
  
## Wrangle
  p$Enh <- splitter(rownames(p), "_", 1)
  p$Classification <- altius_stats$Classification
  p <- melt(p)
  p$value[which(p$value < 0.5)] <- NA # remove values below expression threshold
  colnames(p)[4] <- "DHS"
  # p$value <- log2(p$value)
  
  # reorder enhancers: by classification, and then on-rate within classification
  e <- data.frame(Enh = splitter(rownames(altius_stats), "_", 1),
                  Classification = altius_stats$Classification,
                  On = altius_stats$OnRate)
  e$Classification <- factor(e$Classification, levels = unique(e$Classification)[c(3,2,4,1)])
  e <- e[order(e$On),]
  e <- e[order(e$Classification),]
  p$Enh <- factor(p$Enh, levels = e$Enh)
  
  # reorder tissues: by clustering
  x <- altius_pool$Mean
  x <- x[-which(splitter(rownames(x), "_", 1) %in% e$Enh[which(e$Classification == "Not detected")]),] # just detected ones
  hc <- hclust(dist(t(x)), method = "complete")
  ord <- as.dendrogram(hc)
  colnames(p)[3] <- "Tissue"
  p$Tissue <- factor(p$Tissue, levels = colnames(x)[hc$order])
  
  # vertical lines to delineate the brain
  vline <- grep("Brain", levels(p$Tissue))
  vline[1] <- vline[1] - 0.5
  vline[2] <- vline[2] + 0.5
  
  # filter out the undetected peaks
  # p <- p[-which(p$Classification == "Not detected"),]
  
  # heatmap
  ggHeatmap <- ggplot(p, aes_string(x = "Tissue", y = "Enh", fill = "DHS")) +
    geom_tile() +
    theme_bw() +
    # scale_fill_carto_c(palette = "Geyser", na.value = "grey95") +
    # scale_fill_scico(palette = "lapaz", trans = "log2", na.value = "grey95", direction = -1) +
    scale_fill_gradient2(low = "white", mid = "#ffe2d6", high = "#b51a00", trans = "log2", na.value = "white") +
    geom_vline(xintercept = vline) +
    theme(panel.border = invis, panel.grid = invis, axis.text.y = invis, axis.ticks.y = invis,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "bottom",
          axis.title.x = invis) +
    # guides(fill = guide_legend(title = "log2 DHS")) +
    labs(x = "Tissue", y = "Candidate Enhancers") 
  
  # empty plot
  ggEmpty <- ggplot(p, aes(x = Tissue, y = Enh, fill = DHS)) +
    geom_blank() +
    theme_void()
    
  # column dendrogram
  library(ggdendro)
  ggColDendro <- as.dendrogram(hc) %>% 
    dendro_data(type = "rectangle") %>%
    segment %>%
    ggplot(.) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_void() +
    scale_x_continuous(expand = c(0,0), limits = c(0.5, 48.5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 240)) 
  
  # row marginal distribution
  tab <- table(p$Enh, !(is.na(p$DHS)))
  tab <- as.data.frame.matrix(tab)
  tab <- tab[which(rownames(tab) %in% p$Enh),]
  colnames(tab) <- c("Off", "On")
  tab$Enh <- rownames(tab)
  m <- match(tab$Enh, splitter(rownames(altius_stats), "_", 1))
  tab$Classification <- factor(altius_stats$Classification[m], levels = levels(e$Classification))
  # levels(tab$Classification) <- paste0(levels(tab$Classification), "\n(", table(tab$Classification), ")")
  tab$Enh <- factor(tab$Enh, levels = levels(p$Enh))
  
  rowCols <- c("grey90", pals$Primary[c(2,6,7)])
  showColours(rowCols)
  ggRowCounts <- ggplot(tab, aes(x = Enh, y = On, colour = Classification, fill = Classification)) +
    geom_col() +
    theme_void() +
    coord_flip() +
    scale_colour_manual(values = rowCols) +
    scale_fill_manual(values = rowCols) +
    theme(legend.position = c(0.7, 0.1), legend.title = invis, legend.text = element_text(size = 8),
          legend.direction = "vertical", legend.key.size = unit(3, "mm"))
  
  # legends
  leg1 <- get_legend(ggHeatmap)
  leg2 <- get_legend(ggRowCounts)
  
  # output
  # pdf(file = "Altius - Final Heatmap.pdf", height = 5, width = 7)
  pdf_LibDesign(figNo = "SFig3A", title = "Altius annotation of candidates heatmap", h = 5, w = 7)
  ggarrange(ggColDendro, ggEmpty, ggHeatmap + theme(legend.position = "none"), ggRowCounts + theme(legend.position = "none"),
            heights = c(0.2, 1), widths = c(1, 0.2))
  plot_grid(leg1, leg2)
  dev.off()
  
## Tabulate classifications
  altius_stats$Classification <- factor(altius_stats$Classification, levels = levels(e$Classification))
  
  # pdf(file = "Altius - Classification barplot.pdf", height = 2.5, width = 3.5)
  pdf_LibDesign(figNo = "SFig3C", title = "Altius annotation of candidates barplot", h = 2.5, w = 3.5)
  ggplot(altius_stats, aes(x = Classification, fill = Classification, colour = Classification)) +
    geom_bar() +
    scale_fill_manual(values = c("grey90", pals$Primary[c(2,6,7)])) +
    scale_colour_manual(values = c("grey50", pals$Primary_Darker[c(2,6,7)])) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
          axis.title.x = invis, axis.line.x = invis, axis.text.x = invis, axis.ticks.x = invis) +
    labs(y = "Count of candidates")
  dev.off()
  
## Scatterplot cv versus on-rate
  # pdf(file = "Altius - Onrate vs CV.pdf", height = 3, width = 3)
  pdf_LibDesign(figNo = "SFig3B", title = "Altius annotation of candidates scatterplot", h = 3, w = 3)
  ggplot(altius_stats[-which(altius_stats$OnRate == 0),], aes(x = OnRate, y = CV, fill = Classification, colour = Classification)) +
    # geom_point(shape = 21, colour = "black") +
    geom_point(size = 1) +
    theme_bw() +
    scale_y_continuous(trans = "log2") +
    scale_fill_manual(values = pals$Primary[c(2,6,7)]) +
    scale_colour_manual(values = pals$Primary[c(2,6,7)]) +
    labs(x = "Tissue on-rate (Mean DHS > 0.5)", y = "Coefficient of variation") +
    geom_hline(yintercept = 1, colour = "red") +
    geom_vline(xintercept = 0.5, colour = "red") +
    theme(panel.border = invis, axis.line = element_line(), panel.grid = invis,
          legend.position = "none")
  dev.off()
  
################################################################################################################################ #
## Enhancer trends in the Herring dataset ----
  
## Load
  load("../../../FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Intersections and Coverage.rda")
  hMods <- read.csv("../../../FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Cell-type Specificity Models.csv")

## Heatmap in ggplot2
  ## Wrangle for ggplot2
    p <- herring$CoverageMean_Pooled
   
    p <- apply(p, 1, function(x) {
      if (max(x) > 0 ){
        # x / max(x)
        scale(x)
      } else {
        x
      }

    })
    
    rownames(p) <- colnames(herring$CoverageMean_Pooled)
    p <- t(p)
    p <- melt(p)
    colnames(p) <- c("Enh", "Tissue", "Coverage")

  ## Clustering
    # log transform to stabilise data pre clustering
    x <- herring$CoverageMean_Pooled
    min_exp <- 0.1 # at this cut-off, 80% of all peaks but 11% of non-peaks are "expressed"
    x <- log2(x + min_exp)
    
    # reorder columns
    hc_col <- hclust(dist(t(x)), method = "complete")
    # p$Tissue <- factor(p$Tissue, levels = colnames(x)[hc_col$order]) # instead of running this line, the below will rotate the far left branch - a purely visual change without affecting clustering
    
    dd_col <- as.dendrogram(hc_col)
    weights_col <- ifelse(colnames(x) == "Astro_Neonatal", yes = 50, no = 950) # + ifelse(colnames(x) == "Oligo_Adolescence", yes = 100, no = 0)
    dd_col <- reorder(dd_col, wts = weights_col, agglo.FUN = mean) %>% dendro_data(type = "rectangle")
    p$Tissue <- factor(p$Tissue, levels = dd_col$labels$label)
   
    levels(p$Tissue) <- gsub("Infancy", "Infnt", levels(p$Tissue)) %>% # to rename stages
      gsub("Adolescence", "Adlsc", .) %>%
      gsub("Neonatal", "Neont", .) %>%
      gsub("Childhood", "Child", .) %>%
      gsub("Adolescence", "Adlsc", .) %>%
      gsub("_", " ", .) 

    # reorder rows
    e <- data.frame(Enh = hMods$X,
                    Classification_Ct = "Other",
                    # Classification_Dev = "Other",
                    Mean = rowMeans(herring$CoverageMean_Pooled),
                    Height = 1)
    e$Classification_Ct[which(hMods$Glia.Enriched)] <- "Glial"
    e$Classification_Ct[which(hMods$Astro.Enriched)] <- "Astrocytic" # due to ordering of lines, (astro & glial) outputs as astro
    e$Classification_Ct[which(!(hMods$Astro.Enriched | hMods$Glia.Enriched))] <- "Shared" 
    e$Classification_Ct <- factor(e$Classification_Ct)
    
    e <- e[rev(order(e$Classification_Ct, -e$Mean)),]
    p$Enh <- factor(p$Enh, levels = e$Enh)

  # heatmap
  text.col <- levels(p$Tissue)
  text.col[grep("Oligo|Micro", text.col)] <- "#3a5f9a"
  text.col[grep("Astro", text.col)] <- "#6f083d"
  text.col[grep("^Exc|^Inh", text.col)] <- "#e4784e"
  
  ggHeatmap <- ggplot(p, aes(x = Tissue, y = Enh, fill = Coverage)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient2(low = "white", mid = "grey97", high = "#b51a00") + # to orange
    # scale_fill_gradient2(low = pals$grn2orng[9], mid = pals$grn2orng[5], high = pals$grn2orng[1], midpoint = 0) +
    theme(panel.border = invis, panel.grid = invis, axis.text.y = invis, axis.ticks = invis,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "left",
          axis.title.x = invis) +
    theme(axis.text.x = element_text(colour = text.col)) +
    labs(x = "Tissue", y = "Candidate enhancers") 
  
  # empty plot
  ggEmpty <- ggplot(p, aes(x = Tissue, y = Enh, fill = Coverage)) +
    geom_blank() +
    theme_void()
  
  # rowbars, as astrocyte specific
  e$Enh <- factor(e$Enh, levels = e$Enh)
  # rowLabCol <- c("#6f083d", "#3a5f9a", "#e4784e")
  rowLabCol <- pals$Primary[c(7, 6, 1)]
  
  ggRowLabs <- ggplot(e, aes(x = Enh, y = Height, fill = Classification_Ct, colour = Classification_Ct)) +
    geom_col() +
    scale_fill_manual(values = rowLabCol) +
    scale_colour_manual(values = rowLabCol) +
    theme_void() +
    theme(legend.position = "none") +
    theme(legend.title = invis) +
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1)) 
  
    # dendrogram for column
  # ggColDendro <- ggDendro_gjs(hc_col)
  
  # rotated dendrogram
  ggColDendro <- segment(dd_col) %>%
     ggplot(.) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      theme_void() +
      scale_x_continuous(expand = c(0,0), limits = c(0.5, 27 + 0.5)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 60))
  
  # output
  # pdf(file = "Herring - Final Heatmap (ggplot) (v11).pdf", height = 4.5, width = 4.5)
  pdf_LibDesign(figNo = "1C", title = "Herring heatmap", h = 4.5, w = 4.5)
  ggarrange(ggColDendro, ggEmpty, ggHeatmap  + theme(legend.position = "none"), ggRowLabs,
            heights = c(0.2, 1), widths = c(1, 0.1))
  
 
  leg0 <- get_legend(ggHeatmap + theme(legend.position = "bottom") + guides(fill = guide_colorbar(title = "Coverage", title.position = "top", title.hjust = 0.5)))
  leg1 <- get_legend(ggRowLabs + theme(legend.position = "bottom"))
  plot_grid(leg0, leg1, ncol = 1)
  
  dev.off()

 
################################################################################################################################ #
## Validation of our astrocyte scRNAseq using transcriptional correlation ----
    
## Load NHA
  load("../../../FullScale/Data/Preprocessed/NHA Pooled (Final).rda")
  keepgenes <- rowMeans(nha@assays$RNA@data) > 2^-6
  # nha_filt <- subset(nha, features = keepgenes)

## Load pseudobulks from Herring 2022
  load("../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/Pseudobulk_byGJS.rda", verbose = TRUE)
    
  # now pool samples within each stage/cell-type combination 
  pb$Meta$Category <- paste0(pb$Meta$Celltype, "_", pb$Meta$Stage)
  tc <- list()
  for (j in unique(pb$Meta$Category)) { tc[[j]] <- rowSums(pb$Exp[,which(pb$Meta$Category == j)])  }
  tc <- do.call("cbind", tc)
  tc <- apply(tc, 2, function(x) { x / (sum(x) / 10^6) })
  tc <- as.data.frame(tc)

## Run singleR
  # common genes
  common <- intersect(rownames(nha), rownames(tc))
  tc <- tc[common,]  
  nha_filt <- subset(nha, features = common)
  
  # function
  singleR_annots <- SingleR(test = nha_filt@assays$RNA@data,
                            ref = log2(tc+0.5),
                            labels = colnames(tc),
                            genes = "de")

  save(singleR_annots, file = "SingleR.rda")
  
## Barplot
  p <- singleR_annots$pruned.labels
  p[is.na(p)] <- "NA_NA"
  p <- factor(p, levels = c(colnames(tc), "NA_NA"))
  # p <- factor(p, levels = c(colnames(tc), "NA_NA", paste0(unique(splitter(colnames(tc), "_", 1)), "_NA"), paste0("NA_", unique(splitter(colnames(tc), "_", 2)))) )

  # barplot
  q <- as.data.frame(table(p))
  q$Stage <- splitter(q$p, "_", 2)
  q$Stage <- factor(q$Stage, levels = unique(q$Stage)[c(5,1,3,2,4,6,7)])
  levels(q$Stage)[7] <- "Unassigned"
  q$Celltype <- splitter(q$p, "_", 1)
  q$Celltype <- factor(q$Celltype, levels = unique(q$Celltype)[c(3,1,2,4,5,6,7)])
  levels(q$Celltype) <- gsub("Astro", "Astrocyte", levels(q$Celltype)) %>%
    gsub("Exc", "Excitatory\nNeurons", .) %>%
    gsub("Inh", "Inhibitory\nNeurons", .) %>%
    gsub("Micro", "Microglia", .) %>%
    gsub("Oligo", "Oligodendrocytes", .) 
  q$Label <- q$Freq
  q$Label[which(q$Label == 0)] <- NA
  
  pal <- c(pals$Primary[c(6:1)], "grey50")
  pal2 <- c(pals$Primary_Darker[c(6:1)], "grey50")
 
  # pdf(file = "Celltype validation - SingleR Barplot.pdf", height = 3, width = 7.5)
  pdf_LibDesign(figNo = "SFig1A", title = "SingleR", h = 3, w = 7.5)
  offset <- 1
  ggplot(q, aes(x = Celltype, fill = Stage, y = Freq+offset, colour = Stage)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = Label), position = position_dodge(width = 1), vjust = -0.25,
              hjust = 0.4, size = 2.5, show.legend = FALSE) +
    theme_bw() +
    facet_grid(.~Celltype, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal2) +
    guides(fill = guide_legend(ncol = 2)) +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
          legend.position = c(0.7, 0.7), strip.background = invis,
          axis.text.x = invis, axis.ticks.x = invis) +
    scale_y_continuous(expand = c(0,0), limits = c(1, 100001),
                       trans = "log10", labels = function(x) { comma(x - offset) },
                       breaks = c(0, 10, 1000, 100000) + offset) +
    scale_x_discrete(expand = c(0.15,0)) +
    labs(y = "Number of cells", x = "Celltype annotation from brain maturation timecourse")
  dev.off()  
  
  
################################################################################################################################ #
## Validation of our astrocyte ATAC-seq ----
  
## Here, the intention is to show that our NHA ATAC-seq more closely aligns to astrocyte/brain than other samples
  
## In Herring et al.
  # load
  load("../../../FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Intersections and Coverage.rda", verbose = TRUE)
  p <- herring$CoverageMean_Pooled
  # p <- p > 0.1
  p <- colSums(p) / nrow(p)
  p <- as.data.frame(p)
  p$Celltype <- splitter(rownames(p), "_", 1)
  p$Stage <- splitter(rownames(p), "_", 2)
  p$Stage <- factor(p$Stage, levels = unique(p$Stage)[c(4,6,5,3,1,2)])

  
  ctOrd <- aggregate(p~Celltype, FUN = mean, data = p)
  p$Celltype <- factor(p$Celltype, levels = ctOrd$Celltype[order(-ctOrd$p)])

  # pdf(file = "Celltype validation (ATAC) - Herring open pct.pdf", height = 2.2, width = 3)
  pdf_LibDesign(figNo = "SFig1B", title = "NHA ATAC vs Herring snATACseq", h = 2.2, w = 3)
  ggplot(p, aes(x = Celltype, y = p*100, colour = Stage)) +
    geom_point(position = position_dodge(width = 0.5)) +
    stat_summary(fun = mean, geom = "point", shape = "+", colour = "black", size = 5) +
    scale_colour_manual(values = pals$Primary_Darker[6:1]) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 105)) +
    theme_bw() +
    labs(y = "Percentage of candidates\nin open ATAC") +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
          axis.title.x = invis, axis.text.y = text90, legend.position = "none")
  dev.off()
  
  
  

################################################################################################################################ #
## dCas9-KRAB expression ----
  
## Load
  if (!(exists("nha"))) {
    load("../../../FullScale/Data/Preprocessed/NHA Pooled (Final).rda", verbose = TRUE)    
  }
  
## Get dCas9-KRAB expression
  g <- "dCas9-KRAB-T2A-BLAST-WPRE"
  
  krabStats <- data.frame(Counts = nha@assays$RNA@counts[g,], 
                          Norm = nha@assays$RNA@data[g,])
  krabStats$CPM <- krabStats$Counts / (nha$UMI_Total / 10^6)
  
  
  means <- rowMeans(nha@assays$RNA@data)
  expStats <- data.frame(Gene = names(means), Mean = means, Rank = rank(means))

## Plot rank
  # expStats[,c("ShowShape", "ShowLabel")] <- "Yes"
  expStats$ShowShape <- expStats$Gene == g
  
  expStats$ShowLabel <- expStats$Gene
  expStats$ShowLabel[which(expStats$Gene != g)] <- NA
  pct <- (expStats[g,"Rank"] / max(expStats$Rank)) %>% 
    round(3)
  expStats$ShowLabel <- gsub(g, paste0("dCas9-KRAB\n(", 100*pct, "th percentile)"), expStats$ShowLabel)
  
  expStats$Expressed <- expStats$Mean < 2^-6
  expStats$Expressed[which(expStats$Gene == g)] <- g
  expStats$Expressed <- factor(expStats$Expressed)
  
  levels(expStats$Expressed) <- c(paste0("dCas9-KRAB"), 
                                  "Expressed", "Not expressed")
  expStats <- expStats[rev(order(expStats$Expressed)),]
  
  vlines <- expStats[g, "Mean"]
  
  expStats <- expStats[-which(expStats$Mean == 0),]
  
  # pdf(file = "dCas9-KRAB expression rank.pdf", height = 2.2, width = 2.4)
  pdf_LibDesign(figNo = "SFig4F", title = "dCas9KRAB", h = 2.2, w = 2.4)
  ggplot(expStats, aes(y = Rank, x = log2(Mean), label = ShowLabel, colour = Expressed, size = Expressed)) +
    # geom_hline(yintercept = 0) +
    # geom_col(colour = "black") +
    geom_point(show.legend = TRUE) +
    scale_colour_manual(values = c(pals$Primary[8], "black", "grey90")) +
    scale_size_manual(values = c(5, 1, 1)) +
    # geom_text(nudge_y = -3000, nudge_x = 1.5, size = 3, show.legend = FALSE) +
    geom_text(nudge_y = -1000, nudge_x = -9, size = 3, show.legend = FALSE) +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.line = element_line(),
          axis.text.y = text90, legend.position = "none", legend.title = invis,
          legend.key.size = unit(0.5, "cm")) +
    scale_x_continuous(limits = c(-20,3), expand = c(0,0)) +
    scale_shape_manual(values = c(NA, 21)) +
    labs(y = "Expression rank", x = "log2 mean norm expression") +
    scale_y_continuous(limits = c(0, max(expStats$Rank+1000)), expand = c(0,0),
                       breaks = c(10000, 20000))
  dev.off()

  #####
  
  
    

  
