## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----


## Generic
  rm(list = ls()); gc()
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/1_Processing/")
  options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(tricycle)
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(tidyverse)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(rcartocolor)
  library(sceptre)
  library(glmGamPoi)
  
  
  
## Load
  source("../../Scripts/Functions.R")

## Data information
  samples <- c(paste0("NHA_", c(1:5, 7:8)))
  sample.colours <- carto_pal(7, "Bold")
  names(sample.colours) <- samples
            
## Plotting parameters
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  

################################################################################################################################ #
## Sequencing summary statistics ----
  

## Read in Cellranger sequencing summaries
  seqsum <- lapply(samples, function(j) {
    
    x <- read.csv(paste0("../../Data/Sequencing/FinalCount/", j, "/outs/metrics_summary.csv"))  
    x <- as.data.frame(t(x))
    x$Metric.Name <- rownames(x)
    
    # convert to numeric
    g <- grep("%", x$V1)
    x$V1 <- gsub(",", "", x$V1)
    x$V1 <- gsub("%", "", x$V1)
    x$V1 <- as.numeric(x$V1)
    x$V1[g] <- x$V1[g] / 100
    
    # rename metrics
    x$Metric.Name <- gsub("\\.", " ", x$Metric.Name)
    x$Metric.Name <- gsub("Fraction", "", x$Metric.Name)
    x$Metric.Name <- gsub("Confidently to", "Confidently to\n ", x$Metric.Name)
    x$Metric.Name <- gsub("Antisense to", "Antisense\nto", x$Metric.Name)
    x$Metric.Name[g] <- paste0(x$Metric.Name[g], " (%)")
    
    # output
    x$Sample <- j
    return(x)
    })
  
  seqsum <- do.call("rbind", seqsum)
 
## Plot
 pdf(file = "Library-level/Pre-filter/CellRanger Report.pdf", height = maxh, width = maxw)
 ggplot(seqsum[grep("%", seqsum$Metric.Name),], aes(x = Sample, y = V1*100, fill = Sample)) +
   geom_col(colour = "black") +
   scale_fill_manual(values = sample.colours) +
   facet_wrap(~Metric.Name, ncol = 3) +
   scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
   theme_bw() +
   guides(fill = guide_legend(nrow = 3)) +
   theme(axis.text.x = invis, axis.ticks.x = invis, axis.title = invis, legend.position = c(0.7, 0.1))
 dev.off()
 
 
################################################################################################################################ #
## Read in CellRanger count output into Seurat ----
 
 obj <- lapply(samples, function(j) {
   as.Seurat(file = paste0("../../Data/Sequencing/FinalCount/", j,"/outs/filtered_feature_bc_matrix/"),
             sample.id = j, 
             h5 = FALSE,
             min.cells = 3, 
             min.features = 200)
 })
 
 names(obj) <- samples
  
 
 
################################################################################################################################ #
## QC of counts matrices ----
 
## Basic qc using Seurat
  obj <- lapply(obj, qc.features, cc = TRUE)
 

## Plot
  # extract data for plotting
  p <- lapply(obj, function(x) x@meta.data)
  p <- do.call("rbind", p)
  
  # format
  p <- melt(p, id.vars = "orig.ident")
  p$value <- as.numeric(p$value)
  levels(p$variable) <- c("UMIs Per Nucleus", "Genes Per Nucleus", 
                          "Mito %", "Small Ribo %", "Large Ribo %", "MALAT1 %", 
                          "S Score", "G2M Score", "Phase")
  
 
  # ggplot
  qc.plot <- function(vars, trans = "identity", lim = c(0,NA)) {
    x <- p[which(p$variable %in% vars),]
    
    ggplot(x, aes(x = orig.ident, y = value, fill = orig.ident, colour = orig.ident)) +
      geom_violin(colour = "black", scale = "width") +
      # geom_boxplot(fill = "white", outlier.shape = NA, width = 0.20) +
      scale_fill_manual(values = sample.colours) +
      scale_y_continuous(limits = lim, expand = c(0,0), labels = scales::comma, trans = trans) +
      scale_colour_manual(values = sample.colours) +
      facet_wrap(~variable, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
            axis.title = element_blank(), panel.grid.major.x = element_blank(),
            legend.position = "bottom")
  }
  
  pdf("Library-level/Pre-filter/Expression Per Nucleus.pdf", height = 6, width = maxw)
  qc.plot(c("UMIs Per Nucleus", "Genes Per Nucleus"), trans = "log10", lim = c(100, NA)) 
  dev.off()
  
## Explore relationship between depth and QC
    # extract data for plotting
  p <- lapply(obj, function(x) x@meta.data)
  p <- do.call("rbind", p)
  
   # annotate cells as low or high depth
  p$Depth <- cut(p$nCount_RNA, c(0, 3000, 5000, 10000, 1000000))
  levels(p$Depth) <- c("<3000", "3000-5000", "5000-10000", ">10000")
  
  # format
  p <- melt(p, id.vars = c("orig.ident", "Depth"))
  p$value <- as.numeric(p$value)
  levels(p$variable) <- c("UMIs Per Nucleus", "Genes Per Nucleus", 
                          "Mito %", "Small Ribo %", "Large Ribo %", "MALAT1 %", 
                          "S Score", "G2M Score", "Phase")
  
  
  pdf("Library-level/Pre-filter/Other Transcripts Per Nucleus (Bin By Total UMIs).pdf", height = maxh, width = maxw)
  x <- p[which(p$variable %in% c( "Mito %", "Small Ribo %", "Large Ribo %", "MALAT1 %")),]
  ggplot(x, aes(x = orig.ident, y = value, fill = Depth)) +
    geom_violin(colour = "black", scale = "width") +
    scale_y_continuous(limits = c(0,NA), expand = c(0,0), labels = scales::comma, trans = "identity") +
    facet_wrap(~variable, scales = "free_y", ncol = 1) +
    theme_bw() +
    geom_hline(yintercept = 10, linetype = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.title = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom")
  
  dev.off()
  
  # now plot numbers in each bin
  x <- as.data.frame(table(p$Depth, p$orig.ident) / 6)
  pdf(file = "Library-level/Pre-filter/Number of cells per expression bin.pdf", height = 4, width = maxw)
  ggplot(x, aes(x = Var2, y = Freq, fill = Var1)) +
    geom_col(colour = "black", position = "dodge") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.title.x = invis, panel.grid.major = invis, panel.border = invis, axis.line.y = element_line())
  dev.off()  
  
  
################################################################################################################################ #
## Normalise and filter each library separately  ----
  
  
## Quick and dirty threshold: 10000 UMIs per cell
  umi.thresh <- 10000
  
## Normalise
  obj <- lapply(obj, norm.bySeurat, min.nUMI = umi.thresh) # normalise and threshold 
  
## PCA
  obj <- lapply(obj, RunPCA) # PCA
  
## tSNE/UMAP
  # heuristic for number of dimensions to use in tSNE
  p <- lapply(obj, ElbowPlot, ndims = 50)
  p <- lapply(p, function(x) x + geom_hline(yintercept = c(1, 1.5, 2, 2.5)))
  
  pdf(file = "Library-level/Pre-filter/Elbow Plots for PCA Dimensions.pdf", height = maxh, width = maxw*2)
  plot_grid(plotlist = p, ncol = 3, labels = samples, align = "hv") # based on this, set to 25 as it is a common plateau point
  dev.off()
  
## Run basic UMAP 
  obj <- lapply(obj, clust.bySeurat, dim = 25, resolution = 1, method = "UMAP")
  
  # visualise
  p <- lapply(obj, function(x) {
    data.frame(x@meta.data, x@reductions$umap@cell.embeddings)
  })
  
  p <- do.call("rbind", p)
  
  umap.plot <- function(meta, facet = TRUE) {
    plot <- ggplot(p, aes_string(x = "UMAP_1", y = "UMAP_2", colour = meta)) +
      geom_point(size = 1) +
      theme_void() +
      scale_colour_carto_c(palette = "Geyser") +
      
      theme(legend.position = c(0.7, 0.15)) 
    
    if (facet) {
      plot + facet_wrap(~orig.ident) 
    } else {
      plot
    }
  }
  
  pdf(file = "Library-level/UMAP Plots.pdf", height = maxh-1.5, width = maxw)
  umap.plot("nCount_RNA")
  umap.plot("percent.mito")
  umap.plot("percent.rps")
  umap.plot("percent.rpl")
  umap.plot("percent.malat1")
  umap.plot("S.Score")
  umap.plot("G2M.Score")
  umap.plot("Phase") + scale_colour_carto_d(palette = "Geyser") 
  dev.off()
  
  
################################################################################################################################ #
## Save ----
    
save(obj, file = "../../Data/Preprocessed/NHA Separate Libraries.rda")
 

################################################################################################################################ #
## Create pooled object ----
    

## Merge
    nha <- merge(x = obj$NHA_1, y = obj[2:7], add.cell.ids = samples[1:7])
  
## Normalise
  # restrict to cells with a guide
  # nha <- subset(nha, cells = which(nha$AnyGuide)) # update: will be applied in a later script

  # normalise expression
  nha <- norm.bySeurat(nha, min.nUMI = umi.thresh)
  
## Cluster
  # pca
  nha <- RunPCA(nha)
  pdf(file = "Pooled/Elbow Plot.pdf", width = maxw, height = 3)
  ElbowPlot(nha, ndims = 50) + geom_hline(yintercept = c(1, 1.5, 2, 2.5))
  dev.off()
  
  # UMAP
  nha <- clust.bySeurat(nha, dim = 25, resolution = 2, method = "UMAP") # formerly 20, now 25
  
  
################################################################################################################################ #
## Final QC on pooled cells ----
 
  
## Remove cluster annotation
  nha@meta.data <- nha@meta.data[,-grep("RNA_snn", colnames(nha@meta.data))]
  
## Recalculate metadata that change depending on relative properties
  # recalculate cell-cycle score using Seurat
  nha <- CellCycleScoring(nha, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)

  # calculate cell-cycle score using Tricycle
  x <- normalizeCounts(as.matrix(nha@assays$RNA@counts)) # per recommendations, use this function
  nha$Tricycle <- estimate_cycle_position(x = x,
                                          species = "human",
                                          gname.type = "SYMBOL")
  nha$Tricycle.Stage <- cut(nha$Tricycle, pi * (c(0, 0.5, 1, 1.75, 2)))
  levels(nha$Tricycle.Stage) <- c("G1/0", "S", "G2M", "G1/0")
    
    # per Tricycle documentation:
      # 0.5pi to be the start of S stage,   
      # pi to be the start of G2M
      # 1.5pi to be the middle of M stage, 
      # 1.75pi-0.25pi to be G1/G0 stage             
    

  # calculate cell-cycle score using the Tricycle implementation of the Schwabb approach  
  nha$Schwabe.Stage <- estimate_Schwabe_stage(x = x,
                         species = "human",
                         gname.type = "SYMBOL")
            
## Redo columnames
  colnames(nha@meta.data) <- c("Library", "UMI_Total", "Gene_Total", 
                               "Mito_Pct", "RiboS_Pct", "RiboL_Pct", "MALAT1_Pct", 
                               "Cycle_Seurat_S", "Cycle_Seurat_G2M", "Cycle_Seurat_Phase", 
                               "Cluster",
                               "Cycle_Tricycle", "Cycle_Tricycle_Phase", "Cycle_Schwabe_Phase")
  
  
## Plot UMAP
  p <- data.frame(nha@meta.data, nha@reductions$umap@cell.embeddings)
  
  pdf(file = "Pooled/UMAP - Technical.pdf", height = maxh, width = maxw)
  umap.plot("UMI_Total", facet = FALSE)
  umap.plot("Mito_Pct", facet = FALSE)
  umap.plot("RiboS_Pct", facet = FALSE)
  umap.plot("RiboL_Pct", facet = FALSE)
  umap.plot("MALAT1_Pct", facet = FALSE)
  dev.off()
  
  pdf(file = "Pooled/UMAP - Cell Cycle (Seurat).pdf", height = maxh, width = maxw)
  umap.plot("Cycle_Seurat_S", facet = FALSE)
  umap.plot("Cycle_Seurat_G2M", facet = FALSE)
  umap.plot("Cycle_Seurat_Phase", facet = FALSE) + scale_colour_lancet()
  umap.plot("Cycle_Seurat_Phase", facet = FALSE) + scale_colour_lancet() + facet_wrap(~Cycle_Seurat_Phase) + NoLegend()
  dev.off()
  
  pdf(file = "Pooled/UMAP - Cell Cycle (Tricycle).pdf", height = maxh, width = maxw)
  umap.plot("Cycle_Tricycle", facet = FALSE)
  umap.plot("Cycle_Tricycle_Phase", facet = FALSE) + scale_colour_lancet()
  umap.plot("Cycle_Tricycle_Phase", facet = FALSE) + scale_colour_lancet() + facet_wrap(~Cycle_Tricycle_Phase) + NoLegend()
  dev.off()
  
  pdf(file = "Pooled/UMAP - Library.pdf", height = maxw, width = maxw)
  DimPlot(nha, group.by = "Library", pt.size = 0.2, shuffle = TRUE) + scale_colour_manual(values = sample.colours) + theme_void()
  dev.off()
  


################################################################################################################################ #
## Save ---- 
    
## Save
  save(nha, file = "../../Data/Preprocessed/NHA Pooled (No Guide Annotation).rda")
  
  ## Friday
    # run guide-to-cell assignment
    # check MAGMA for Sam
    # other plots for Irina? Check Slack!


  
      
   
    