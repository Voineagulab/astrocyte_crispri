setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/1_Processing/")
load("../../Data/Preprocessed/NHA Pooled (Final).rda")
guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv")

# Make hda5
# Make plots
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
  leg.pos <- c(0.2, 0.15)
  
  pdf(file = "Final/UMAP - Technical.pdf", height = maxh, width = maxw)
  umap.plot("UMI_Total", facet = FALSE)
  umap.plot("Mito_Pct", facet = FALSE)
  umap.plot("Ribo_Pct", facet = FALSE)
  umap.plot("MALAT1_Pct", facet = FALSE)
  umap.plot("NuclearFraction", facet = FALSE) 
  umap.plot("DropletQC", facet = FALSE) + scale_colour_lancet()
  dev.off()
  
  pdf(file = "Final/UMAP - Cell Cycle (Continuous).pdf", height = 12, width = 10)
  plot_grid(umap.plot("Cycle_Seurat_S", facet = FALSE),
            umap.plot("Cycle_Seurat_G2M", facet = FALSE),
            umap.plot("Cycle_Tricycle_PC1", facet = FALSE),
            umap.plot("Cycle_Tricycle_PC2", facet = FALSE),
            umap.plot("Cycle_Tricycle", facet = FALSE),
            nrow = 3)
  dev.off()
  
  pdf(file = "Final/UMAP - Cell Cycle (Phases).pdf", height = maxh, width = maxw)
  umap.plot("Cycle_Seurat_Phase", facet = FALSE) + scale_colour_lancet()
  umap.plot("Cycle_Seurat_Phase", facet = FALSE) + scale_colour_lancet() + facet_wrap(~Cycle_Seurat_Phase) + NoLegend()
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
  
## Focus on DropletQC
  p <- nha@meta.data
  p2 <- ggplot(p, aes(x = NuclearFraction, colour = DropletQC)) +
    theme_bw() + 
    theme() +
    scale_x_continuous(limits = c(0,1)) +
    labs(x = "Nuclear Fraction") 
  
  pdf(file = "Final/DropletQC - Nuclear Fraction.pdf", height = 3, width = 5)
  p2 + geom_point(mapping = aes(y = log10(UMI_Total))) 
  p2 + geom_point(mapping = aes(y = Mito_Pct))
  p2 + geom_point(mapping = aes(y = Ribo_Pct))
  p2 + geom_point(mapping = aes(y = MALAT1_Pct))
  p2 + geom_point(mapping = aes(y = log10(MOI+1))) 
  dev.off()
  
  p3 <- ggplot(p, aes(x = DropletQC, colour = DropletQC)) +
    theme_bw() + 
    theme(axis.title.x = invis, legend.position = "none")

  pdf(file = "Final/DropletQC - Violins.pdf", height = 3, width = 4)
  p3 + geom_violin(mapping = aes(y = log10(UMI_Total)), draw_quantiles = 0.5, scale = "width") 
  p3 + geom_violin(mapping = aes(y = Mito_Pct), draw_quantiles = 0.5, scale = "width") 
  p3 + geom_violin(mapping = aes(y = Ribo_Pct), draw_quantiles = 0.5, scale = "width") 
  p3 + geom_violin(mapping = aes(y = MALAT1_Pct), draw_quantiles = 0.5, scale = "width") 
  p3 + geom_violin(mapping = aes(y = log10(MOI+1)), draw_quantiles = 0.5, scale = "width")
  dev.off()

## Assess for doublets using simulation
  # get umis of cells, and simulated doublets
  get.info <- function(w, label) {
    data.frame(UMI = p$UMI_Total[w],
               MOI = p$MOI[w],
               Label = label)
  }

  a <- get.info(which(p$MOI == 1), "MOI=1")
  b <- get.info(which(p$MOI %in% c(1:10)), "MOI<10")
  bb <- get.info(which(p$MOI > 20), "MOI>20")
  c <- a[sample(rownames(a), 300, FALSE),1:2] + a[sample(rownames(a), 300, FALSE),1:2]
  c$Label <- "Doublet1"
  
  d <- b[sample(rownames(b), 300, FALSE),1:2] + b[sample(rownames(b), 300, FALSE),1:2]
  d$Label <- "Doublet10"
  
  e <- rbind(a, b, bb, c, d)
  e$Bin <- cut(e$MOI, c(1:10, 15, 20, 1000), right = FALSE)
  levels(e$Bin) <- c(1:10, "11-15", "16-20", ">20")
  e$Label <- factor(e$Label, levels = levels(as.factor(e$Label))[c(4,3,5,1,2)])
  
  ggplot(e, aes(x = Label, y = UMI, fill = Label)) +
    geom_violin(draw_quantiles = c(0.5, 0.75, 0.95)) +
    coord_flip() +
    scale_fill_manual(values = c("white", "white", "white", "grey70", "grey70"))
  
  ggplot(e, aes(x = Bin, y = UMI)) +
    geom_violin(scale = "width", draw_quantiles = 0.5) +
    facet_wrap(~Label, nrow = 1) 
  
  table(p$MOI > 20, p$UMI_Total > 100000)
  table(p$MOI > 20, p$UMI_Total > 125000)
  table(p$MOI > 20, p$UMI_Total > 150000)
  
  
  # calculate 99th% of UMIs under null, per library
  l <- data.frame(Lib = p$Library, UMI = p$UMI_Total, MOI = p$MOI)
  l$Bin <- cut(l$MOI, c(0,1,10,20,1000), right = FALSE)
  levels(l$Bin) <- c(0, "1-10", "11-20", "21+")
  m <- median(l$UMI)
  # l <- l[which(l$MOI == 1),]
  
  pdf(file = "Final/Doublets - Distribution of UMI Across Libs and MOI.pdf", height = 4, width = 10)
  ggplot(l, aes(x = Bin, y = UMI, fill = Lib)) +
    geom_violin(draw_quantiles = c(0.5, 0.75, 0.95), scale = "width", alpha = 0.4) +
    theme() +
    geom_hline(yintercept = c(m,125000), linetype = 2, colour = "red") +
    scale_fill_lancet()
  
  ggplot(l, aes(x = Lib, y = UMI, fill = Bin)) +
    geom_violin(draw_quantiles = c(0.5, 0.75, 0.95), scale = "width", alpha = 0.4) +
    theme() +
    geom_hline(yintercept = c(m,125000), linetype = 2, colour = "red") +
    scale_fill_aaas()
  
  ggplot(l, aes(x = Lib, y = UMI)) +
    geom_violin(draw_quantiles = c(0.5, 0.75, 0.95), scale = "width") +
    theme() +
    geom_hline(yintercept = c(m,125000), linetype = 2, colour = "red") +
    facet_wrap(~Bin) 
  dev.off()
  
  ## Generate doublets per library
    # the sim
    doub.sim <- list() 
    
    for (j in unique(p$Library)) {
      print(j)
      
      a <- get.info(which(p$MOI == 1 & p$Library == j), "MOI=1")
      b <- get.info(which(p$MOI %in% c(1:10) & p$Library == j), "MOI<10")
      bb <- get.info(which(p$MOI > 20 & p$Library == j), "MOI>20")
      c <- a[sample(rownames(a), 300, FALSE),1:2] + a[sample(rownames(a), 300, FALSE),1:2]
      c$Label <- "Doublet1"
      d <- b[sample(rownames(b), 300, FALSE),1:2] + b[sample(rownames(b), 300, FALSE),1:2]
      d$Label <- "Doublet10"
    
      e <- rbind(a, b, bb, c, d)
      e$Bin <- cut(e$MOI, c(1:10, 15, 20, 1000), right = FALSE)
      levels(e$Bin) <- c(1:10, "11-15", "16-20", ">20")
      e$Label <- factor(e$Label, levels = levels(as.factor(e$Label))[c(4,3,5,1,2)])
      
      e$Lib <- j
      
      doub.sim[[j]] <- e
    }
    
    doub.sim <- do.call("rbind", doub.sim)
  
  ## Plot 
    pdf(file = "Final/Doublets - Doublet Simulations.pdf", height = 8, width = 9)
     ggplot(doub.sim, aes(x = Label, y = UMI / 1000, fill = Label)) +
      geom_violin(draw_quantiles = c(0.5, 0.75, 0.95), scale = "width") +
      facet_wrap(~Lib) +
       labs(y = "UMIs (thousands)") +
      scale_fill_manual(values = c(carto_pal(7, "Earth")[5:7], "grey70", "grey50")) +
       theme(axis.text.x = invis, axis.title.x = invis, axis.ticks.x = invis) +
       geom_hline(yintercept = c(median(p$UMI_Total / 1000), 125), linetype = 2, colour = "red")
    dev.off()
    
  ## Stats
    # the strategy: check every cell with MOI 1. 
    
    # 99th percentile of MOI 1 
  
 ## Hypothesis:
   # NHA_8 has a lot of ambient RNA.
   # hard to test using traditional tools, as no cell-types!
   # instead, use guides, predicting: 1) more "impossible" guide assignments, 2) greater balance between guides
   
  
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
## Guide expression! ----
  
## Load
  load("../../Data/Preprocessed/Guide Expression Matrix.rda")
  load("../../Data/Preprocessed/Metadata.rda")

  
## Quick validation...
  validating <- list()
  for (j in colnames(guide.exp)) {
    print(j)
    a <- guide.exp[,j] >= 3 %>% as.numeric()
    b <- nha@meta.data[j,rownames(guide.exp)] %>% as.numeric()
    
    c <- data.frame(N = length(which(b == 1)), Mx_Not_Meta = length(which((a & ! (b)))), Meta_Not_Mx = length(which((b & !(a)))))
    validating[[j]] <- c
  }
  
  save(validating, file = "Final/Guide Expression Matrix - Temp.rda")
  
  x <- do.call("rbind", validating)
  
## Plot the number of guides per cell as a function of threshold
  fun1 <- function(thresh, x = guide.exp) {
    # x <- guide.exp > 3
    colSums(x >= thresh)
  }
  
  test.thresholds <- c(1, 2, 3, 4, 5, 10, 20, 25, 50)
  p <- sapply(test.thresholds, fun1) %>% as.data.frame()
  colnames(p) <- paste0("MinUMI_", test.thresholds)
  
  p <- melt(p)
  p$Bin <- cut(p$value, c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 25, 50, 100, 4000), right = FALSE)
  levels(p$Bin) <- c(as.character(0:10), "11-14", "15-24", "25-49", "50-99", ">100")
  levels(p$variable) <- gsub("MinUMI_", "", levels(p$variable))
  
  # ggplot(p, aes(x = variable, y = value, colour = variable)) +
  #   geom_violin() +
  #   scale_y_continuous(limits = c(0, 50))
  # ggplot(p, aes(x = value, colour = variable)) +
  #   geom_density() +
  #   facet_wrap(~variable, scales = "free") +
  #   scale
  #   
  
  pdf(file = "Final/Guide Expression Matrix - MOI Across UMI Thresholds (Bin).pdf", height = 9.5, width = 8)
  q <- table(p$Bin, p$variable) %>% as.data.frame()
  ggplot(q, aes(x = Var2, fill = Var1, y = Freq)) +
    geom_col(colour = "black") +
    # scale_fill_viridis_d()
    scale_fill_manual(values = c("grey90", carto_pal(11, "Prism")[1:10], carto_pal(6, "Pastel")[1:5])) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    guides(fill = guide_legend(title = "MOI", ncol = 3)) +
    theme(panel.border = invis, axis.line.y = element_line()) +
    labs(y = "Number of Cells", x = "Minimum Guide UMI")
  dev.off()
    
  
  # and the hypothetical under a Poisson distribution?
  mean.moi <- 1:10
  test.moi.prob <- 0:30
  
  p <- sapply(mean.moi, function(y) {
    dpois(x = test.moi.prob, lambda = y)
  })
  
  p <- as.data.frame(p)
  colnames(p) <- paste0("TransductionMOI", mean.moi)
  p$CellHasNGuides <- test.moi.prob
  rownames(p) <- paste0("CellHas", test.moi.prob, "Guides")
  
  
  p <- melt(p, id.vars = "CellHasNGuides")
  p$value <- p$value * ncol(guide.exp)
  levels(p$variable) <- gsub("TransductionMOI", "", levels(p$variable))
  
  library(ggsci)
  pdf(file = "Final/Guide Expression Matrix - Expected (Poisson) Number of Guides Per Cell Across Transduction MOIs.pdf", height = 4, width = 8)
  ggplot(p, aes(x = CellHasNGuides, y = value, colour = variable)) +
    geom_line() +
    geom_point() +
    scale_colour_d3() +
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(xintercept = 5, linetype = 2) +
    guides(colour = guide_legend(title = "Transduction\nMOI", ncol = 2)) +
    labs(y = "Expected Number of Cells With nGuide", x = "nGuides Per Given Cell") +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.8, 0.6))
  dev.off()
  
  
  
  p$Bin <- cut(p$CellHasNGuides, c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 25, 50, 100, 4000), right = FALSE)
  levels(p$Bin) <- c(as.character(0:10), "11-14", "15-24", "25-49", "50-99", ">100")
  
  pdf(file = "Final/Guide Expression Matrix - Expected (Poisson) Number of Guides Per Cell Across Transduction MOIs (Bin).pdf", height = 9.5, width = 8)
  ggplot(p, aes(x = as.factor(variable), fill = Bin, y = value)) +
    geom_col(colour = "black") +
    # scale_fill_viridis_d()
    scale_fill_manual(values = c("grey90", carto_pal(11, "Prism")[1:10], carto_pal(6, "Pastel")[1:5])) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    guides(fill = guide_legend(title = "MOI", ncol = 3)) +
    theme(panel.border = invis, axis.line.y = element_line()) +
    labs(y = "Expected Number of Cells Under Poisson", x = "Poisson Mean MOI Parameter")
  dev.off()
  
  ## Scatterplot expected vs observed cells above a given nGuide as a function of UMI and MOI
  
    fun2 <- function(UMI.thresh, ExpectedMOI, labs = FALSE) {
      # collect the number of guides per cell at the given UMI threshold (y)
      a <- fun1(UMI.thresh)
      
      # calculate the proportion of cells with at least test.moi.prob guides
      b <- sapply(test.moi.prob, function(z) length(which(a >= z)) / length(a))
      
      # calculate the expected proportion of cells with at least test.moi.prob guides, under the poisson assumption
      c <- ppois(q = test.moi.prob, lambda = ExpectedMOI, lower = FALSE) # ppois, with lower = FALSE, gives the probability of a the poisson distribution being < q
      
      # plot
      d <- data.frame(Observed = b, Expected = c, nGuides = test.moi.prob, UMI.thresh = UMI.thresh, PoissonMOI = ExpectedMOI)
      scale.factor <- ncol(guide.exp) / 1000
      d$Observed <- d$Observed * scale.factor
      d$Expected <- d$Expected * scale.factor
      return(d)
      
      # e <- ggplot(d, aes(x = Expected, y = Observed, label = nGuides)) + 
      #   # geom_point(shape = as.character(test.moi.prob)) +
      #   geom_text() +
      #   geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.8, colour = "red") +
      #   scale_x_continuous(limits = c(0,scale.factor)) +
      #   scale_y_continuous(limits = c(0,scale.factor)) +
      #   theme_bw() +
      #   theme(panel.border = invis, panel.grid.minor = invis) +
      #   labs(x = paste0("Expected nCells (1000s) With nGuides > Label:\nAssuming MOI=", ExpectedMOI),
      #        y = paste0("Observed nCells (1000s) With nGuides > Label:\nThresholding At ", UMI.thresh, " UMIs"))
      # 
      # if (labs) {
      #   xlab
      # } else {
      #   e + theme()
      # }
        
      
  }
  
  res <- list()
    for (j in test.thresholds) {
      print(j)
      for (k in mean.moi) {
        res[[paste0("Thresh", j, "_MOI", k)]] <- fun2(UMI.thresh = j, ExpectedMOI = k)
      }
    }
    
  res <- do.call("rbind", res)
  # res <- res[which(res$PoissonMOI < 11),]
  res$UMI.thresh <- factor(res$UMI.thresh)
  res$PoissonMOI <- factor(res$PoissonMOI)
  # levels(res$UMI.threshlevels) <- paste0(levels(res$UMI.thresh), " UMI Threshold")
  # levels(res$PoissonMOI) <- paste0("Poisson MOI ", levels(res$PoissonMOI))
  levels(res$UMI.thresh) <- paste0("Min UMI:\n", levels(res$UMI.thresh))
  levels(res$PoissonMOI) <- paste0("Poisson MOI:\n", levels(res$PoissonMOI))
  
  plot <- ggplot(res, aes(x = Expected, y = Observed, label = nGuides)) + 
    # geom_point(shape = as.character(test.moi.prob)) +
    geom_text(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.8, colour = "red") +
    facet_grid(PoissonMOI~UMI.thresh) +
    scale_x_continuous(limits = c(0,scale.factor)) +
    scale_y_continuous(limits = c(0,scale.factor)) +
    theme_bw() +
    labs(x = "Thousands of Expected nCells With nGuides > Label", 
         y = "Thousands of Observed nCells With nGuides > Label") +
    theme(panel.grid = invis) 
        
  pdf(file = "Final/Guide Expression Matrix - Expected Versus Observed MOI.pdf", height = 12, width = 10)
  plot
  plot + scale_y_continuous(trans = pseudo_log_trans(sigma = 5, base = 10)) + scale_x_continuous(trans = pseudo_log_trans(sigma = 5, base = 10))
  dev.off()
  
  
## Plot the fraction of guide reads occupied by each guide
  # get information
  guide.frac <- apply(guide.exp, 2, function(x) {
    if (max(x) < 3) return(0)
    x <- x[which(x >= 3)]
    x <- x / sum(x)
    
    data.frame(AboveThresh = length(x), 
               Percent05 = length(which(x > 0.05)),
               Percent02 = length(which(x > 0.02)),
               Percent01 = length(which(x > 0.01)))
    })
  
  guide.frac <- do.call("rbind", guide.frac)
  
  # bin MOI
  guide.frac$Bin <- cut(guide.frac$AboveThresh, c(0,1,2,3,4,5,6,8,10,13,16,21,26,31,41,51,1000), right = FALSE)
  levels(guide.frac$Bin) <- c("0", "1", "2", "3", "4", "5", "6-7", "8-9", "10-12", "13-15", "16-20", "21-25", "26-30", "31-40", "41-50", "51+")
  
  # bin change in MOI
  guide.frac$Loss <- cut(guide.frac$AboveThresh - guide.frac$Percent05, c(0,1,2,3,4,5,6,8,10,13,16,21,26,31,41,51,1000), right = FALSE)
  levels(guide.frac$Loss) <- c("0", "1", "2", "3", "4", "5", "6-7", "8-9", "10-12", "13-15", "16-20", "21-25", "26-30", "31-40", "41-50", "51+")
  
  # ggplot(guide.frac, aes(x = AboveThresh, y = AboveThresh-Percent05)) +
  #   # geom_violin() +
  #   geom_col() +
  #   facet_wrap(~Bin, ncol = 4, scales = "free")

  pdf(file = "Final/Guide Expression Matrix - 5pct Dominance.pdf", height = 12, width = 14)
    ggplot(guide.frac, aes(x = Loss)) +
    # geom_violin() +
    geom_bar() +
    facet_wrap(~Bin, ncol = 4, scales = "free") +
      labs(x = paste0("Number of Guides Lost Per Cell (<5% Dominance)\n\n(Each facet title notes the original cell MOI)"),
           y = "Number of Cells") +
      theme_bw()
  dev.off()
  
## Check the effect of these guides on positive controls!
  use <- guides$GuideID[which(guides$TargetCat == "Promoter")]
  thresh <- 3
  keep.cells <- colnames(guide.exp)[(colSums(guide.exp >= thresh) > 0)]
  # high.moi <- colnames(guide.exp)[(colSums(guide.exp >= thresh) >= 20)]
  check.fraction <- apply(guide.exp, 2, function(x) {

    x[which(x < thresh)] <- 0
    x <- x / sum(x)
  })
  check.fraction <- as.data.frame(check.fraction)
  
  pdf(file = "Final/Guide Expression Matrix - High MOI On PosCon Suppression.pdf", height = 3, width = 4)
  for (j in use) {
    print(j)
    
    gene <- guides$TargetID[which(guides$GuideID == j)]
    
    if (!(gene %in% rownames(nha))) next
    if (gene == "EIF3A") next

    # basic dataframe with binary guide assignment
    p <- data.frame(Exp = nha@assays$VST@data[gene, keep.cells],
                    GuideUMI = t(guide.exp[j,keep.cells] >= thresh),
                    GuideFraction = t(check.fraction[j,keep.cells]),
                    GuideDominant = t(check.fraction[j,keep.cells] >= 0.02))
    colnames(p) <- c("Exp", "GuideUMI", "GuideFraction", "GuideDominant")
    # add cell moi
    p$MOI <- "<20"
    p[high.moi,"MOI"] <- ">20"
    # p$MOI[which(p$GuideFractionHigh & p$MOI == ">20")] <- ">20 (but dominant)"
    # p$MOI[which(p$GuideFractionHigh & p$MOI == "<20")] <- "<20 (but minor)"
    p$MOI <- factor(p$MOI)
    levels(p$MOI) <- paste0(levels(p$MOI), paste0(" (n=", table(p$MOI, p$GuideUMI)[,2], ")"))
    
    downsample <- sample(c(1:nrow(p))[-which(p$GuideUMI)], nrow(p)-1000, FALSE)
    p <- p[-downsample,]
    
    print(ggplot(p, aes(y = Exp, x = GuideUMI, colour = MOI)) +
      geom_violin(draw_quantiles = 0.5, fill = "white", width = 0.7, scale = "width") +
      geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.4) +
      theme_bw() +
      scale_colour_carto_d(palette = "Earth") +
      scale_y_continuous() +
      labs(y = paste0("Normalised ", gene, " Expression"), x = paste0("Has ", j)) +
      theme(panel.border = invis, axis.line.y = element_line(),
            axis.ticks.x = invis, panel.grid.major.x = invis))
  }
  dev.off()
  
## Check the effect dominant guides on positive controls!
  use <- guides$GuideID[which(guides$TargetCat == "Promoter")]
  thresh <- 3
  keep.cells <- colnames(guide.exp)[(colSums(guide.exp >= thresh) > 0)]
  high.moi <- colnames(guide.exp)[(colSums(guide.exp >= thresh) >= 20)]
  check.fraction <- apply(guide.exp, 2, function(x) {

    x[which(x < thresh)] <- 0
    x <- x / sum(x)
  })
  check.fraction <- as.data.frame(check.fraction)
  
  
  # 210fth1, 333lgals3,id353
  
  pdf(file = "Final/Guide Expression Matrix - Guide Dominance (02) On Enh Suppression.pdf", height = 3, width = 5)
  # for (j in guides$GuideID[grep("Enh333", guides$GuideID)]) {
  for (j in use) {
    print(j)
    
    gene <- guides$TargetID[which(guides$GuideID == j)]
    # gene <- "LGALS3"
    
    if (!(gene %in% rownames(nha))) next
    if (gene == "EIF3A") next

    # basic dataframe with binary guide assignment
    p <- data.frame(Exp = nha@assays$VST@data[gene, keep.cells],
                    GuideUMI = t(guide.exp[j,keep.cells] >= thresh),
                    GuideFraction = t(check.fraction[j,keep.cells]),
                    GuideDominant = t(check.fraction[j,keep.cells] >= 0.02))
    colnames(p) <- c("Exp", "GuideUMI", "GuideFraction", "GuideDominant")
    
    
    # add cell category column
    m <- rownames(p) %in% high.moi
    p$Category <- "."
    p$Category[which(!(p$GuideUMI))] <- "No Guide"
    p$Category[which((p$GuideDominant & !(m)))] <- "Dominant, Low MOI"
    p$Category[which(!(p$GuideDominant) & p$GuideUMI  & !(m))] <- "Minor, Low MOI"
    p$Category[which((p$GuideDominant & (m)))] <- "Dominant, High MOI"
    p$Category[which(!(p$GuideDominant) & p$GuideUMI  & (m))] <- "Minor, High MOI"
    p$Category <- factor(p$Category, levels = c("No Guide", "Dominant, Low MOI",
                                                "Minor, Low MOI", "Dominant, High MOI", "Minor, High MOI"))
    levels(p$Category) <- paste0(levels(p$Category), " (n=", table(p$Category), ")")
    
    cols <- c("black", "firebrick1", "firebrick3", "dodgerblue1", "dodgerblue3")
    fills <- c("white", "grey80", "white", "grey80", "white")
    names(cols) <- names(fills) <- levels(p$Category)
    
    downsample <- sample(c(1:nrow(p))[-which(p$GuideUMI)], nrow(p)-1000, FALSE)
    p <- p[-downsample,]
    
    print(ggplot(p, aes(y = Exp, x = Category, fill = Category, colour = Category)) +
      geom_violin(draw_quantiles = 0.5, width = 0.8, scale = "width") +
      # geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.4) +
        geom_jitter(width = 0.15, alpha = 0.2) +
      theme_bw() +
      # scale_colour_carto_d(palette = "Earth") +
      scale_colour_manual(values = cols) +
        scale_fill_manual(values = fills) +
      scale_y_continuous() +
      labs(y = paste0("Normalised ", gene, " Expression"), x = j) +
      theme(panel.border = invis, axis.line.y = element_line(),
            axis.ticks.x = invis, panel.grid.major.x = invis, axis.text.x = invis ))
    
    ggplot(p, aes(x = GuideFraction*100, y = Exp)) +
      geom_point() +
      scale_x_continuous(trans = pseudo_log_trans())
  }
  dev.off()
  
## How many guide assignments would be lost by removing high MOI cells?
  x <- guide.exp >= 3
  remove <- which(colSums(x) >= 20)
  enh.guides <- guides$GuideID[which(guides$TargetCat == "Enh")]
  sum(x[enh.guides,]) # 438448
  sum(x[enh.guides,-remove]) # 197366
  
  length(which(rowSums(x[enh.guides,]) < 30)) # 36
  length(which(rowSums(x[enh.guides,-remove]) < 30)) # 70
  
## Plot the relationship between guide UMI expression and MOI
  # create dataframe
  x <- apply(guide.exp, 2, function(x) sum(x[which(x >= 3)]))
  
  p <- data.frame(GuideUMIs = apply(guide.exp, 2, function(x) sum(x[which(x >= 3)])), # sum total of UMIs for guides, where said guide has 3+ UMIs
                  MaxGuideUMIs = apply(guide.exp, 2, function(x) max(x[which(x >= 3)])), # value of most highly expressed guide
                  MeanGuideUMIs = apply(guide.exp, 2, function(x) mean(x[which(x >= 3)])), # value of most highly expressed guide
                  GeneUMIs = nha$UMI_Total, # total UMIs mapping to genes
                  MOI = nha$MOI, # moi
                  Pool = nha$TransductionPool) # inferred transduction pool
  p <- p[-which(p$GuideUMIs == 0),]
  
  # plot guide UMIs vs. gene UMIs
  q <- melt(p[,-5], id.vars = c("GeneUMIs", "Pool"))
  q$Pool <- factor(q$Pool)
  levels(q$Pool) <- paste0(levels(q$Pool), " Pool")
  
  pdf(file = "Final/Guide Expression Matrix - Gene vs Guide UMIs 1.pdf", height = 7, width = 10)
  ggplot(q, aes(x = (GeneUMIs), y = (value))) +
    geom_point() +
    scale_x_continuous(trans = "log10", labels = scales::comma) +
    scale_y_continuous(trans = "log10") +
    facet_grid(variable~Pool, scales = "free_y", switch = "y") +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(x = "Total (Gene) UMIs", y = "Guide UMI Statistic")
  dev.off()
  
  q <- p
  q$GuideUMIs <- log10(q$GuideUMIs)
  q$GeneUMIs <- log10(q$GeneUMIs)
  q$Density <- get_density(q$GuideUMIs, q$GeneUMIs, n = 100)
  r <- cor(q$GuideUMIs, q$GeneUMIs) %>% round(2)
  
  pdf(file = "Final/Guide Expression Matrix - Gene vs Guide UMIs 2.pdf", height = 3, width = 4)
  ggplot(q[which(q$Pool == "Primary"),], aes(x = GeneUMIs, y = GuideUMIs, colour = Density)) +
    geom_point() +
    # scale_x_continuous(trans = "log10", labels = scales::comma) +
    scale_colour_viridis_c() +
    geom_smooth(method = "lm") +
    # scale_y_continuous(trans = "log10") +
    theme_bw() +
    labs(x = paste0("Log10 Total (Gene) UMIs\nr=", r), y = "Log10 Total Guide UMIs")
  dev.off()

  # plot guide UMIs vs MOI  
  q <- melt(p[,c(1,3,5,6)], id.vars = c("MOI", "Pool"))
  levels(q$variable) <- c("Total Guide Expression Per Cell", "Mean Guide Expression Per Cell")
  
  pdf(file = "Final/Guide Expression Matrix - Guide UMIs vs MOI.pdf", height = 8, width = 5)
  ggplot(q, aes(x = MOI, y = value)) +
    geom_point() +
    facet_grid(Pool~variable, scales = "free", switch = "y") +
    scale_x_continuous(trans = "log10", labels = scales::comma) +
    scale_y_continuous(trans = "log10", labels = scales::comma) +
    theme_bw() +
    theme() +
    labs(x = "MOI", y = "Guide Expression (UMIs)")
  dev.off()  
  
## Distribution of guide UMIs
  y <- melt(check.fraction)
  y <- y[-which(y$value == 0 | is.nan(y$value)),]
  
  moi <- colSums(guide.exp >= thresh)
  m <- match(y$variable, names(moi))
  y$MOI <- moi[m]
  y$Bin <- cut(y$MOI, c(0,5,10,15,20,50,1000), right = FALSE)
  levels(y$Bin) <- gsub("\\[", "", levels(y$Bin)) %>% gsub("\\)", "", .) %>% gsub(",", "-", .)
  levels(y$Bin)[6] <- "50+"
  
  pdf(file = "Final/Guide Expression Matrix - Distribution of Guide Percent.pdf", height = 6, width = 9)
  ggplot(y, aes(x = value*100)) +
    geom_histogram(bins = 100) +
    facet_wrap(~Bin) +
    scale_x_continuous(trans = "log2", breaks = c(0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1)*100) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = c(2, 5), colour = "red", linetype = 2) +
    labs(x = "Guide UMI %", y = "Number of Guide-Cell Assignments")
  dev.off()
  
## Distribution of guide UMIs, at the CPM level!
  z <- melt(guide.exp)
  z <- z[which(z$value >= 3),]
  m <- match(z$variable, colnames(nha))
  z$LibSize <- nha$UMI_Total[m]
  z$CPM <- z$value * ((10^6) / z$LibSize)
  z$Log2CPM <- log2(z$CPM)
  # z$Bin <- 
  z$Lib <- substr(z$variable, 1, 5)
  z$MOI <- nha$MOI[m]
  z$MOIBin <- cut(z$MOI, c(0,5,10,15,20,50,1000), right = FALSE)
  levels(z$MOIBin) <- gsub("\\[", "", levels(z$MOIBin)) %>% gsub("\\)", "", .) %>% gsub(",", "-", .)
  levels(z$MOIBin)[6] <- "50+"
  z$Fraction <- y$value
  z$DropletQC <- meta$DropletQC[m]
  
  pdf(file = "Final/Guide Expression Matrix - CPM.pdf", height = 6, width = 8)
  
  # plot distribution of CPM
  ggplot(z, aes(x = MOIBin, y = Log2CPM)) +
    geom_violin(draw_quantiles = 0.5) +
    # facet_wrap(~MOIBin) +
    # scale_y_continuous(trans = "log2") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "MOI Bin", y = "log2(CPM)")
  dev.off()
  # plot CPM vs fraction
  
  # colour by density
  # draw thresholds, with n
  
  min(z$Fraction)
  
  library(ggExtra)
  plot.fun1 <- function(MOIBin) {
    # filter
    zz <- z[which(z$MOIBin == MOIBin),]
    
    # stats
    zz$Density <- get_density(log10(zz$Fraction), zz$Log2CPM)
    r <- cor(log10(zz$Fraction), zz$Log2CPM) %>% round(2)
    
    # plot
    a <- ggplot(zz, aes(x = Fraction, y = Log2CPM, colour = Density)) +
      geom_point() +
      scale_x_continuous(trans = "log10", limits = c(min(z$Fraction), 1)) +
      scale_y_continuous(limits = c(min(z$Log2CPM), max(z$Log2CPM))) +
      scale_colour_viridis_c() +
      labs(y = "Guide log2(CPM)", x = paste0("Guide Fraction\nr=", r), title = paste0("MOI: ", MOIBin)) +
      geom_vline(xintercept = c(0.02, 0.05), linetype = 2, colour = "red")  +
      geom_hline(yintercept = c(log2(25)), linetype = 2, colour = "red")  +
      theme_bw() +
      theme(legend.position = "none", panel.border = invis, plot.title = element_text(hjust = 0.5))
    
    ggMarginal(a, type = "density")
  }

  pdf(file = "Final/Guide Expression Matrix - CPM vs. Fraction.pdf", height = 8, width = 10)
  plot_grid(plot.fun1(levels(z$MOIBin)[1]),
            plot.fun1(levels(z$MOIBin)[2]),
            plot.fun1(levels(z$MOIBin)[3]),
            plot.fun1(levels(z$MOIBin)[4]),
            plot.fun1(levels(z$MOIBin)[5]),
            plot.fun1(levels(z$MOIBin)[6]),
            nrow = 2)
  dev.off()
  
## Compare each guide species for its expression patterns
  guide.cpm <- list()
  
  for (j in colnames(guide.exp)) {
    print(j)
    x <- guide.exp[,j]
    x[which(x < 3)] <- 0
    guide.cpm[[j]] <- x / (nha@meta.data[j,"UMI_Total"] / 10^6)
  }

  guide.cpm <- do.call("cbind", guide.cpm) %>% t() %>% as.data.frame()  
  colnames(guide.cpm) <- rownames(guide.exp)
  
  
  # stats
  x <- apply(guide.cpm, 2, function(y) {
    summary(y[which(y != 0)])
  })
  
  x <- t(x) %>% as.data.frame()
  x$N <- apply(guide.cpm, 2, function(x) length(which(x > 0)))
  x$Category <- guides$TargetCat
  
  # plot
  # p <- melt(guide.cpm)  
  # p <- p[-which(p$value == 0),]
  # 
  # p$variable <- splitter(p$variable, "_chr", 1)
  
  pdf(file = "Final/Guide Expression Matrix - Cross-Guide Variability.pdf", height = 3, width = 4)
  ggplot(x, aes(x = Category, y = Median)) +
    geom_violin(draw_quantiles = 0.5) +
    scale_y_continuous(limits = c(0,NA)) +
    labs(y = "Median Expression of Guide", x = "Guide Category")
  dev.off()  
  
  