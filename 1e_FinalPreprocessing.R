setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/1_Processing/")

# Make hda5
# Make plots
# Make guide expression matrix (then check ambient assignments and ambiguous as a function of thresholding)
# dCas9-KRAB
# gene summary matrix (mean exp, type, etc)

################################################################################################################################ #
## Cell-level QC Plots ----


basic.col <- carto_pal(2, "Magenta")[1]

## UMIs per cell
  # version 1
  x <- data.frame(UMI = nha$nCount_RNA, Dummy = ".")
  pdf(file = "Final/UMIs Per Cell.pdf", height = 2.5, width = maxw/3)
  m <- median(x$UMI)
  ggplot(x, aes(y = UMI, x = Dummy)) +
    geom_violin(fill = basic.col, draw_quantiles = 0.5) +
    theme_bw() +
    # coord_flip() +
    theme(panel.border = invis, axis.line.y = element_line(), axis.text.x = invis, axis.ticks.x = invis,
          axis.title.x = invis, panel.grid.major.x = invis) +
    # scale_y_continuous(expand = c(0,0), limits = c(1000, 250000), trans = "log10", labels = scales::comma) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 240000)) +
    # geom_hline(yintercept = m, linetype = 2) +
    labs(y = "UMIs Per Cell")
  dev.off()    
  
  # version 2: includes comparison to Gasperini's data
  gasp <- readRDS("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/GSE120861_at_scale_screen.cds.rds")
  x <- data.frame(UMI = nha$nCount_RNA, Study = "Current")
  y <- data.frame(UMI = colSums(gasp@assayData$exprs), Study = "Gasperini")
  x <- rbind(x, y)
  x$UMI <- x$UMI / 1000
  pdf(file = "Final/UMIs Per Cell (vs. Gasperini).pdf", height = 2, width = maxw/3)
  ggplot(x, aes(x = Study, y = UMI)) +
    geom_violin(scale = "width", draw_quantiles = c(0.5), fill = basic.col) +
    theme_bw() +
    # coord_flip() +
    theme(panel.border = invis, axis.line = element_line(), axis.title.x = invis,
          panel.grid.major.x = invis, panel.grid.minor.y = invis) +
    # geom_hline(yintercept = 915, linetype = 2) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 240)) +
    # scale_y_continuous(expand = c(0,0), limits = c(100, 250000),trans = "log10") +
    labs(y = "UMIs Per Cell (1000s)")
  dev.off()    
  
  
## UMAP plots
  # functions
  filter.plot <- function(var, min, max, title) {
    x <- p
    x$Var <- (min < x[,var]) & (x[,var] < max)
    x <- x[order(x$Var, decreasing = FALSE),]
    title <- paste0(title, " (n=", length(which(x$Var)), ")")
    ggplot(x, aes_string(x = "UMAP_1", y = "UMAP_2", colour = "Var", alpha = "Var")) +
      geom_point(size = 1) +
      theme_void() +
      NoLegend() +
      scale_colour_manual(values = c("grey80", "darkorange1")) +
      scale_alpha_manual(values = c(0.25, 1)) +
      labs(title = title)
  }
  
  umap.plot <- function(meta, facet = TRUE) {
    plot <- ggplot(p, aes_string(x = "UMAP_1", y = "UMAP_2", colour = meta)) +
      geom_point(size = 1) +
      theme_void() +
      scale_colour_carto_c(palette = "Geyser") +
      theme(legend.position = leg.pos) 
    
    if (facet) {
      plot + facet_wrap(~orig.ident) 
    } else {
      plot
    }
  }
  
  # plots
  p <- data.frame(nha@meta.data, nha@reductions$umap@cell.embeddings)
  leg.pos <- c(0.8, 0.15)
  
  pdf(file = "Final/UMAP - Technical.pdf", height = maxh, width = maxw)
  umap.plot("UMI_Total", facet = FALSE)
  umap.plot("Mito_Pct", facet = FALSE)
  umap.plot("RiboS_Pct", facet = FALSE)
  umap.plot("RiboL_Pct", facet = FALSE)
  umap.plot("Ribo_Pct", facet = FALSE)
  umap.plot("MALAT1_Pct", facet = FALSE)
  dev.off()
  
  pdf(file = "Final/UMAP - Cell Cycle (Seurat).pdf", height = maxh, width = maxw)
  umap.plot("Cycle_Seurat_S", facet = FALSE)
  umap.plot("Cycle_Seurat_G2M", facet = FALSE)
  umap.plot("Cycle_Seurat_Phase", facet = FALSE) + scale_colour_lancet()
  umap.plot("Cycle_Seurat_Phase", facet = FALSE) + scale_colour_lancet() + facet_wrap(~Cycle_Seurat_Phase) + NoLegend()
  dev.off()
  
  pdf(file = "Final/UMAP - Cell Cycle (Tricycle).pdf", height = maxh, width = maxw)
  umap.plot("Cycle_Tricycle", facet = FALSE)
  umap.plot("Cycle_Tricycle_Phase", facet = FALSE) + scale_colour_lancet()
  umap.plot("Cycle_Tricycle_Phase", facet = FALSE) + scale_colour_lancet() + facet_wrap(~Cycle_Tricycle_Phase) + NoLegend()
  dev.off()
  
  pdf(file = "Final/UMAP - Library.pdf", height = maxw, width = maxw)
  DimPlot(nha, group.by = "Library", pt.size = 0.2, shuffle = TRUE) + scale_colour_manual(values = sample.colours) + theme_void()
  dev.off()
  
  pdf(file = "Final/UMAP - Check On Variables To Threshold.pdf", height = 12, width = 8)
  plot_grid(filter.plot("Mito_Pct", 10, 100, "Mito > 10"),
            filter.plot("Mito_Pct", 15, 100, "Mito > 15"),
            filter.plot("Ribo_Pct", 10, 100, "Ribo > 10"),
            filter.plot("Ribo_Pct", 15, 100, "Ribo > 15"),
            filter.plot("MALAT1_Pct", 10, 100, "MALAT1 > 10"),
            filter.plot("MALAT1_Pct", 15, 100, "MALAT1 > 15"), ncol = 2)
  plot_grid(filter.plot("UMI_Total", 0, 10000, "nUMI < 10000"),
            filter.plot("UMI_Total", 0, 15000, "nUMI < 15000"),
            filter.plot("UMI_Total", 0, 20000, "nUMI < 20000"),
            filter.plot("UMI_Total", 0, 25000, "nUMI < 25000"),
            filter.plot("UMI_Total", 150000, 1000000, "nUMI > 150000"),
            filter.plot("UMI_Total", 200000, 1000000, "nUMI > 200000"), ncol = 2)
  
  dev.off()
  


################################################################################################################################ #
## Guide-level QC Plots ----

  
## Dataframes
  load("../2_DE/Sceptre Input Files.rda")
  
  # at the guide level
  dat.guide <- as.data.frame(sceptre.guide)
  # dat.guide <- dat.guide[grep("Enh", rownames(dat.guide)),]
  
  # pooled at enhancer level
  dat.pool <- as.data.frame(sceptre.guide.pooled)
  # dat.enh <- dat.enh[grep("Enh", rownames(dat.enh)),]

    
## Cells per guide
  # basic function
  cat.guide <- function(x) {
    x[grep("Enh", x)] <- "Enhancer"
    x[grep("_A|_B", x)] <- "Positive"
    x[grep("Neg", x)] <- "Negative"
    
    return(x)
  }
  
  ## Enhancer transduction pool
    # collect n
    x <- dat.guide[,colnames(nha)[which(nha$TransductionPool == "Enh")]] # columns: "used" is the first filter for 
    y <- rowSums(x)
    y <- data.frame(Guide = names(y), nCells = y)
    
    # categorise guides
    y$Category <- cat.guide(y$Guide)
    e <- read.csv("../../../TransferFromRNA/FullLibrary_Selection/Results/Final_List/SR_final_library_design/NHA Enh Library_with Neg Cont.csv")
    e <- e$GuideID
    y$InPool <- (y$Guide %in% e) | grepl("Neg_E", y$Guide)
    
    x <- x[-which(rownames(x) %in% e),]
    y <- apply(x, 2, function(z) any(z == 1))
    
    # plot for guides in pool
    pdf(file = "Final/Number of Cells Per Guide (Enhancer Pool, Expected Guides) (Hist).pdf", height = 2.5, width = maxw)
    ggplot(y[which(y$InPool),], aes(x = nCells)) +
      geom_histogram(colour = "black", fill = basic.col) +
      facet_wrap(~Category, scales = "free_y") +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Expressing Guide", y = "Count of Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line(), legend.position = "none")
    dev.off()
    
    pdf(file = "Final/Number of Cells Per Guide (Enhancer Pool, Expected Guides) (Density).pdf", height = 2.2, width = maxw/2)
    ggplot(y[which(y$InPool),], aes(x = nCells, fill = Category, colour = Category)) +
      geom_density(alpha = 0.2) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Per Guide", y = "Fraction Of Cells") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line(), legend.position = "bottom")
    dev.off()
     
    
     # plot for guides out of pool
     pdf(file = "Final/Number of Cells Per Guide (Enhancer Pool, Unexpected Guides).pdf", height = 2.5, width = maxw)
    ggplot(y[-which(y$InPool),], aes(x = nCells, fill = Category)) +
      geom_histogram(colour = "black", fill = basic.col) +
      facet_wrap(~Category, scales = "free_y") +
      # geom_density() +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(limits = c(0,max(y$nCells))) +
      labs(x = "Cells Expressing Guide", y = "Count of Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
    
    ggplot(y[-which(y$InPool),], aes(x = nCells, fill = Category, colour = Category)) +
      # geom_histogram(colour = "black") +
       geom_density(alpha = 0.2) +
      # facet_wrap(~Category, scales = "free_y") +
      # geom_density() +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(limits = c(0,max(y$nCells))) +
      labs(x = "Cells Expressing Guide", y = "Density") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
    dev.off()
    
  ## Positive transduction pool
    pdf(file = "Final/Number of Cells Per Guide (Positive Pool).pdf", height = 2.5, width = 2.5)
    
    # collect n
    x <- dat.guide[,colnames(nha)[which(nha$TransductionPool == "Pos")]] # columns: "used" is the first filter for 
    y <- rowSums(x)
    y <- data.frame(Guide = names(y), nCells = y)
    y <- y[which(y$nCells > 0),]
    
    # plot for guides in pool
    ggplot(y, aes(x = nCells)) +
      geom_histogram(colour = "black", fill = basic.col) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Expressing Guide", y = "Count of Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line(), legend.position = "none")
    
    # repeat above for pooling guides
    x <- dat.pool[,colnames(nha)[which(nha$TransductionPool == "Pos")]] # columns: "used" is the first filter for 
    y <- rowSums(x)
    y <- data.frame(Guide = names(y), nCells = y)
    y <- y[which(y$nCells > 0),]
    
    ggplot(y, aes(x = nCells)) +
      geom_histogram(colour = "black", fill = basic.col) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Expressing Guide Pair", y = "Count of Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line(), legend.position = "none")
    dev.off()
    
  ## Negative transduction pool
    
    
    # collect n
    x <- dat.guide[,colnames(nha)[which(nha$TransductionPool == "Neg")]] # columns: "used" is the first filter for 
    y <- rowSums(x)
    y <- data.frame(Guide = names(y), nCells = y)
    y <- y[-which(y$nCells == 0 | grepl("_E_", y$Guide)),]
    
    pdf(file = "Final/Number of Cells Per Guide (Negative Pool).pdf", height = 2.5, width = 2.5)
    ggplot(y, aes(x = nCells)) +
      geom_histogram(colour = "black", fill = basic.col) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Expressing Guide", y = "Count of Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line(), legend.position = "none")
    dev.off()
    
    
## Cells per guide: just by category, then coloured by pool
  pool.colours <- pal_lancet()(3)
  names(pool.colours) <- c("Primary", "Negative", "Positive")
  
  ## Enhancer guides
   x <- dat.guide[,colnames(nha)[which(nha$TransductionPool == "Enh")]] 
  y <- rowSums(dat.guide)
  y <- data.frame(Guide = names(y), nCells = y)
    
    # categorise guides
    y$Category <- cat.guide(y$Guide)
    y <- y[which(y$Category == "Enhancer"),]
    # e <- read.csv("../../../TransferFromRNA/FullLibrary_Selection/Results/Final_List/SR_final_library_design/NHA Enh Library_with Neg Cont.csv")
    # e <- e$GuideID
    # y$InPool <- (y$Guide %in% e) | grepl("Neg_E", y$Guide)  
    
    # plot 
     pdf(file = "Final/Number of Cells Per Guide (By Guide Class).pdf", height = 2.2, width = maxw/3)
    ggplot(y, aes(x = nCells)) +
      geom_histogram(colour = "black", fill = pool.colours["Primary"]) +
      # facet_wrap(~Category, scales = "free_y") +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Expressing Guide", y = "Count of Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line(), legend.position = "none")
    dev.off()
    
  ## Positive
     x <- dat.guide[,colnames(nha)[which(nha$TransductionPool == "Enh")]]  %>% rowSums()
     y <- dat.guide[,colnames(nha)[which(nha$TransductionPool == "Pos")]] %>% rowSums()
     z <- data.frame(Guide = names(y), Primary = x, Positive = y)
     z$Category <- cat.guide(z$Guide)
     z <- z[which(z$Category == "Positive"),]
     z2 <- z; z2$Primary[21:250] <- NA
     z <- melt(z[,-1])
     z2 <- melt(z2[,-1])
     
    #
     pdf(file = "Final/Number of Cells Per Guide (By Guide Class).pdf", height = 2.2, width = maxw)
    ggplot(z, aes(x = value, fill = variable)) +
      geom_histogram(colour = "black", position = "stack") +
      facet_wrap(~variable, scales = "free_y") +
      scale_fill_manual(values = pool.colours[c(1,3)]) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Pers Positive Guide", y = "Count of Positive Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
    
    ggplot(z, aes(x = value, fill = variable)) +
      geom_histogram(colour = "black", position = "stack") +
      # facet_wrap(~variable, scales = "free_y") +
      scale_fill_manual(values = pool.colours[c(1,3)]) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Pers Positive Guide", y = "Count of Positive Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
    
    
    ggplot(z2, aes(x = value, fill = variable)) +
      geom_histogram(colour = "black", position = "stack") +
      # facet_wrap(~variable, scales = "free_y") +
      scale_fill_manual(values = pool.colours[c(1,3)]) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Cells Pers Positive Guide", y = "Count of Positive Guides") +
      theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
    dev.off()
  
## Cells per enhancer
  x <- dat.pool[grep("Enh", rownames(dat.pool)),]
  y <- rowSums(x)
  y <- data.frame(Enh = names(y), nCells = y)
  
  m <- median(y$nCells)
  
  pdf(file = "Final/Number of Cells Per Enhancer.pdf", height = 2, width = maxw*2/3)
  ggplot(y, aes(x = nCells)) +
    geom_histogram(fill = basic.col, colour = "black", bins = 30) +
    theme_bw() +
    geom_vline(xintercept = m, linetype = 2) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Cells Per Enhancer Target", y = "Count of Enhancers") +
    theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
  dev.off()
    
    
## Guides per cell
  x <- dat.guide[,colnames(nha)[which(nha$TransductionPool == "Enh")]]
  y <- colSums(x)
  y <- data.frame(Cell = names(y), nGuides = y)
  
  pdf(file = "Final/Number of Guides Per Cell (Unfiltered).pdf", height = 2, width = maxw/2)
  ggplot(y, aes(x = nGuides)) +
    theme_bw() +
    geom_histogram(fill = basic.col, colour = "black", bins = 50) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "Number of Cells", x = "Guides Per Cell") +
    theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
  dev.off()
  
  y$Extreme <- y$nGuides > 50
  pdf(file = "Final/Number of Guides Per Cell (Filtered).pdf", height = 2, width = maxw/2)
  ggplot(y[-which(y$Extreme),], aes(x = nGuides)) +
    geom_histogram(fill = basic.col, colour = "black", bins = 50) +
    theme_bw() +
    # facet_wrap(~Extreme, scales = "free") +
    # geom_vline(xintercept = m, linetype = 2) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "Number of Cells", x = "Guides Per Cell") +
    theme(panel.grid = invis, panel.border = invis, axis.line.y = element_line())
  dev.off()
  
  table(y$Extreme) 
  
    # FALSE  TRUE 
    # 35240  1056 
  
  pdf(file = "Final/Number of Guides Per Cell (MOI UMAP).pdf", height = 4, width = 4)
  filter.plot("MOI", nrow(guides) / 100, 1000, "MOI > 31")
  filter.plot("MOI", 50, 1000, "MOI > 50")
  filter.plot("MOI", 60, 1000, "MOI > 60")
  filter.plot("MOI", 70, 1000, "MOI > 70")
  filter.plot("MOI", 80, 1000, "MOI > 80")
  filter.plot("MOI", 90, 1000, "MOI > 90")
  filter.plot("MOI", 100, 1000, "MOI > 100")
  filter.plot("UMAP_1", 5, 1000, "UMAP1 > 5")
  dev.off()
  
## NGuides per enhancer length
  table(p$UMAP_1 < 5, p$MOI > 50) %>% fisher.test()
  table(p$UMAP_1 > 5, p$MOI > 25) %>% fisher.test()
  table(p$UMAP_1 > 5, p$MOI > 50) %>% fisher.test()
  table(p$UMAP_1 > 5, p$MOI > 100) %>% fisher.test()
  
################################################################################################################################ #
## dCas9-KRAB ----

dcas9 <- "dCas9-KRAB-T2A-BLAST-WPRE" # name of the construct in the expression matrix

## Plot DE as a function of dCas9-KRAB bin
  # what is the distribution of expression?
  pdf(file = "Final/dCas9 Distribution.pdf", height = 3, width = maxw)
  p <- data.frame(Raw_Count = nha@assays$RNA@counts[dcas9,],
                  VST_Count = nha@assays$VST@counts[dcas9,],
                  LibSize = nha$nCount_RNA)
  
  p$Raw_Density <- get_density(p$Raw_Count, p$LibSize, n = 30)
  p$VST_Density <- get_density(p$VST_Count, p$LibSize, n = 30)
  
  # p <- melt(p, id.vars = "LibSize")
  cor <- round(cor(p$Raw_Count, p$LibSize), 2)
  ggplot(p, aes(x = LibSize, y = Raw_Count, colour = Raw_Density)) +
    geom_point() +
    scale_colour_viridis_c() + 
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = paste0("Library Size Per Cell\nr=", cor), y = "dCas9-KRAB-T2A-BLAST-WPRE Expression\n(Raw Count)") +
    theme(panel.border = invis, axis.line = element_line())
  
  cor <- round(cor(p$VST_Count, p$LibSize), 2)
  ggplot(p, aes(x = LibSize, y = VST_Count, colour = VST_Density)) +
    geom_point() +
    scale_colour_viridis_c() + 
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = paste0("Library Size Per Cell\nr=", cor), y = "dCas9-KRAB-T2A-BLAST-WPRE Expression\n(VST Count)") +
    theme(panel.border = invis, axis.line = element_line())
  
  # use VST to define bins, as it is not driven by library size
  p$Bin <- cut(p$VST_Count, c(-1, 0, 5, 10, 20, 150))
  levels(p$Bin) <- c("0", "1-5", "6-10", "10-20", ">20")
  
  ggplot(p, aes(x = Bin)) +
    geom_bar(fill = "black") + 
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "dCas9-KRAB-T2A-BLAST-WPRE Expression\n(VST Count)", y = "Number of Cells in Bin") +
    # labs(x = paste0("Library Size Per Cell\nr=", cor), y = "dCas9-KRAB-T2A-BLAST-WPRE Expression\n(VST Count)") +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis)
  dev.off()
  
  
  

  
  
################################################################################################################################ #
## Normalisation of expression for covariates, focusing particularly on depth ----
  
## Use SCTransform to normalise these data
  nha <- SCTransform(nha, 
                     vars.to.regress = c("Library", "UMI_Total", "Mito_Pct", "Cycle_Seurat_S", "LogMOI"), 
                     vst.flavor = "v2", 
                     method = "glmGamPoi",
                     new.assay.name = "VST",
                     return.only.var.genes = FALSE)

    
## Reset default assay
  DefaultAssay(nha) <- "RNA"
  

  
################################################################################################################################ #
## Save ----
  
  
## Metadata
  # all
  x <- as.data.frame(nha@meta.data)
  write.csv(x, file = "Final/Metadata.csv")
  save(x, file = "Final/Metadata.rda")
  
  
## As Seurat object in rda
  save(nha, file = "../../Data/Preprocessed/NHA Pooled (Final).rda")
  

## As hda5
  library(SeuratDisk)
  
  # note: this uses information from the output of 2_DE, thus is run out of order
  load("../2_DE/Archive_V1/Enhancers - All Results Summary.rda")
  
  # output when pooling guides to targets
  x <- nha
  y <- x@meta.data
  hit <- res$Enh[which(res$Hit)]
  y <- y[,c(colnames(y)[1:25], unique(hit))]
  x@meta.data <- y
  
  SaveH5Seurat(x, filename = "../../Data/Preprocessed/H5AD/NHA_TargetLevel.h5Seurat")
  Convert("../../Data/Preprocessed/H5AD/NHA_TargetLevel.h5Seurat", dest = "h5ad")
  
   # output for guides in hit targets
  x <- nha
  y <- x@meta.data
  hit <- guide.list$GuideID[which(guide.list$TargetID %in% hit)]
  y <- y[,c(colnames(y)[1:25], unique(hit))]
  x@meta.data <- y
  
  SaveH5Seurat(x, filename = "../../Data/Preprocessed/H5AD/NHA_GuideLevel.h5Seurat")
  Convert("../../Data/Preprocessed/H5AD/NHA_GuideLevel.h5Seurat", dest = "h5ad")
  
