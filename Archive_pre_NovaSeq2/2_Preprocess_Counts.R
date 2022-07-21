## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----


## Generic
  rm(list = ls()); gc()
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Data/Preprocessed/")
  options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(rcartocolor)

## Load
  # load("getProtospacers/barcodes2ps.rda")
  source("../../Scripts/Functions.R")


## Data information
  samples <- c(paste0("NHA_", 1:8),
               paste0("SY5Y_", 1:10))
  
  ct <- c("NHA", "SY5Y")
  
  ct.col <- carto_pal(n = 2, name = "Geyser")
  names(ct.col) <- ct
  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  
  

################################################################################################################################ #
## Sequencing summary statistics ----
  

## Read in Cellranger sequencing summaries
  seqsum <- list()
  for (j in samples) {
    x <- read.csv(paste0("../GIMR_GWCCG_200802_IRIVOI_10X/211215_A00152_0509_BHWJ5MDSX2_summary/", j,"/summary/metrics_summary.csv"))  
    x <- x[c(1,13,18,19,20,21,22:27,29:33),c(1,5,6)]
    g <- grep("%", x$Metric.Value); x$Metric.Name[g] <- paste0(x$Metric.Name[g], " (%)")
    x$Sample <- j
    x$Celltype <- splitter(j, "_", 1)
    
    colnames(x)[1] <- "Metric.Category"
    seqsum[[j]] <- x
  }
  
 seqsum <- do.call("rbind", seqsum)
 
## Convert character to numeric
 seqsum$Metric.Value <- gsub(",", "", seqsum$Metric.Value)
 
 seqsum$Metric.Value <- gsub("%", "", seqsum$Metric.Value)
 seqsum$Metric.Value <- as.numeric(seqsum$Metric.Value)
 
## Refactor
 seqsum$Metric.Name <- gsub("Confidently mapped", "Confidently mapped\n", seqsum$Metric.Name)
 seqsum$Metric.Name <- gsub("Number of reads assigned", "Number of reads assigned\n", seqsum$Metric.Name)
 seqsum$Metric.Name <- factor(seqsum$Metric.Name)
 seqsum$Metric.Name <- factor(seqsum$Metric.Name, levels = levels(seqsum$Metric.Name)[c(12, 1, 15,  # reads and cell total
                                                                                        13, 14, 8, # reads assigned, sequencing saturation, fric
                                                                                        2:7, # confident map
                                                                                        9:11, # median per cell
                                                                                        16:17 # valid barcodes and umis
 )]
 )
 
 
## Plot
 pdf(file = "QC/Sequencing Statistics.pdf", height = maxh, width = 7)
 ggplot(seqsum, aes(x = Sample, y = Metric.Value, fill = Celltype)) +
   geom_col(colour = "black") +
   scale_fill_manual(values = ct.col) +
   facet_wrap(~Metric.Name, scale = "free_y", ncol = 3) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
         axis.title = element_blank(), legend.position = c(0.8, 0.05))
 dev.off()
 
 
################################################################################################################################ #
## Read in ----
 

## Ideally, you should be reading in the raw matrix, but you don't have that yet!
  obj <- list()
  
  for (j in samples) {
    cat("\n\n")
    obj[[j]] <- as.Seurat(file = paste0("../GIMR_GWCCG_200802_IRIVOI_10X/211215_A00152_0509_BHWJ5MDSX2_summary/", j,"/summary/sample_feature_bc_matrix"),
                         sample.id = j, h5 = FALSE,
                         min.cells = 3, min.features = 200)
  }
 
 
################################################################################################################################ #
## QC of counts matrices ----
 
## Basic qc
  for (j in samples) {
    logi <- TRUE
    if (j == "NHA_6") {
      logi <- FALSE
    } else {
      logi <- TRUE
    }
    obj[[j]] <- qc.features(obj[[j]], cc = logi)
  }
 
## Plot
  # dummy columns for NHA_6, in which cell cycle quantification fails
  obj$NHA_6@meta.data[,c("S.Score", "G2M.Score", "Phase")] <- NA
  
  # extract data for plotting
  p <- lapply(obj, function(x) x@meta.data)
  p <- do.call("rbind", p)
  
  p <- melt(p, id.vars = "orig.ident")
  p$Celltype <- splitter(as.character(p$orig.ident), "_", 1)  
  p$value <- as.numeric(p$value)
  levels(p$variable) <- c("UMIs Per Nucleus", "Genes Per Nucleus", 
                          "Mito %", "Small Ribo %", "Large Ribo %", "MALAT1 %", 
                          "S Score", "G2M Score", "Phase")
  
  # ggplot
  qc.plot <- function(vars, trans = "identity", lim = c(0,NA)) {
    x <- p[which(p$variable %in% vars),]
    
    ggplot(x, aes(x = orig.ident, y = value, fill = Celltype, colour = Celltype)) +
      geom_violin(colour = "black", scale = "width") +
      geom_boxplot(fill = "white", outlier.shape = NA, width = 0.20) +
      scale_fill_manual(values = ct.col) +
      scale_y_continuous(limits = lim, expand = c(0,0), labels = scales::comma, trans = trans) +
      scale_colour_manual(values = ct.col) +
      facet_wrap(~variable, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
            axis.title = element_blank(), panel.grid.major.x = element_blank(),
            legend.position = "bottom")  
  }
  
  pdf("QC/Expression Per Nucleus.pdf", height = 6, width = maxw)
  qc.plot(c("UMIs Per Nucleus", "Genes Per Nucleus"), trans = "log10", lim = c(100, NA)) 
  qc.plot(c("UMIs Per Nucleus", "Genes Per Nucleus"), trans = "identity", lim = c(100, NA)) 
  dev.off()
  
  pdf("QC/Other Transcripts Per Nucleus.pdf", height = maxh, width = maxw)
  qc.plot(c( "Mito %", "Small Ribo %", "Large Ribo %", "MALAT1 %")) + geom_hline(yintercept = 10, linetype = 2)
  dev.off()
  
  
## tSNE visualisation for clustering
  x <- lapply(obj, norm.bySeurat, min.nUMI = 1000)
  x <- lapply(x, RunPCA)
  x <- lapply(x, clust.bySeurat, dim = 35)
  
  DimPlot(x$NHA_1, group.by = "Phase")
  
  
  ## CODE WAS LOST DURING CRASH! These generated the files "tSNE -*")
  

  
# ################################################################################################################################ #
# ## Gene-level QC ----
#   
#   
# ## Distribution of gene expression
#   x <- rowSums(obj$NHA_1@assays$RNA@counts)
#   y <- x/ncol(obj$NHA_1)
#   
#   p <- data.frame(TotalCounts = x, AverageCounts = y, LogAverage = )
#   
#   ggplot(p, aes(y = y)) +
#     geom_boxplot() +
#     scale_y_continuous(trans = pseudo_log_trans(sigma = 10, base = 1.5))
#   
#   
# ## Fraction zeros
#   x <- apply(obj$NHA_1@assays$RNA@counts, 1, function(x) length(which(x == 0)))
#   y <- x/ncol(obj$NHA_1)
#   
#   p <- data.frame(Fraction = y,
#                   Bin = cut(y, breaks = c(0,1,2,3,4,5,6,7,8,9,9.5, 9.9, 9.99, 10)/10))
#   
#   lab()
#   
#   ggplot(p, aes(x = y)) +
#     geom_density() +
#     scale_y_continuous()
#   
#   ggplot(p, aes(x = Bin)) +
#     geom_bar()
#   
#   z <- lapply(obj, function(x) {
#     x <- apply(x@assays$RNA@counts, 1, function(x) length(which(x == 0)))
#     y <- x/ncol(x)
#     
#     
#   })
#   
#   zz <- do.call("rbind", z)
#   
#   
#   p <- data.frame(Fraction = y,
#                     Bin = cut(y, breaks = c(0,1,2,3,4,5,6,7,8,9,9.5, 9.9, 9.99, 10)/10))
#   
  
 
################################################################################################################################ #
## Barcode assignment ----
  
## Read in the whitelists
  bc.whitelist <- read.table("../Whitelists/3M-february-2018.txt")
  bc.whitelist <- as.character(bc.whitelist$V1)
  
  # protospacers: the guide RNAs cloned in
  x <- read.csv("../Whitelists/Neg Library.csv")
  x <- data.frame(GuideID = paste0("Neg_", 1:nrow(x)),
                  GuideSequence = x$GuideSequence,
                  GuideCoord = NA,
                  TargetCat = "Negative",
                  TargetID  = x$TargetID,
                  TargetCoord = NA)
  x$TargetID[which(is.na(x$TargetID))] <- "NoTarget"
  
  y <- read.csv("../Whitelists/Pos Library.csv")
  y <- data.frame(GuideID = y$GuideID,
                  GuideSequence = y$GuideSequence,
                  GuideCoord = ".",
                  TargetCat = "Promoter",
                  TargetID  = y$TargetID,
                  TargetCoord = ".")
  
  ps.whitelist <- list()
  ps.whitelist$NHA <- read.csv("../Whitelists/NHA Enh Library.csv")
  ps.whitelist$SY5Y <- read.csv("../Whitelists/SY5Y Enh Library.csv")
  ps.whitelist <- lapply(ps.whitelist, function(z) rbind(z, y, x))
  ps.whitelist <- do.call("rbind", ps.whitelist)
  ps.whitelist$Celltype <- splitter(rownames(ps.whitelist), "\\.", 1)
  ps.whitelist$GuideSequence <- substr(ps.whitelist$GuideSequence, 
                                       start = 1, 
                                       stop = 16) # keeps only the first sixteen nt, as sixteen is the length of protospacer sequence that get_barcodes.py called
  ps.whitelist <- unique(ps.whitelist)
  
  write.csv(ps.whitelist, file = "../Whitelists/Protospacer Whitelist.csv")

## Read in assignments
  g <- list()
  for (j in samples) {
    # read
    k <- gsub("_", "-", j)
    
    if (j %in% c("NHA_6", "SY5Y_5", "SY5Y_6")) {
      x  <- read.table(paste0("../../Results/1_GuideAssignment/", j, "_screen.txt"), sep = "\t", header = TRUE)  # use assignments from the screen's sequencing
    } else {
      x  <- read.table(paste0("../../Results/1_GuideAssignment/", k, "_hnPCR2.txt"), sep = "\t", header = TRUE) # use assignments from hnPCR2
    }
    
  
    y <- ps.whitelist[which(substr(ps.whitelist$Celltype, 1, 3) == substr(j, 1, 3)),]
    
    # is it a called barcode from cellranger?
    x$cr_call <- x$cell %in% colnames(obj[[j]])
    
    # what protospacer is it targeting?
    # m <- match(x$barcode, ps.whitelist$GuideSequence)
    # x$GuideID <- ps.whitelist$GuideID[m]  
    # x$TargetCat <- ps.whitelist$TargetCat[m]  
    # x$TargetID <- ps.whitelist$TargetID[m]  
    # x$TargetCt <- ps.whitelist$Celltype[m]  

    
    # err... repeat the above, but removing the first nucleotide (which is nearly always G. might be a cloning artifact?) 
    m <- match(substr(x$barcode, 2, 1000), substr(y$GuideSequence, 1, 15))
    x$GuideID <- y$GuideID[m]  
    x$TargetCat <- y$TargetCat[m]  
    x$TargetID <- y$TargetID[m]  
    x$TargetCt <- y$Celltype[m]  
    
    # output
    g[[j]] <- x
    
  }

  
## Plot library-level analysis
  # wrangle
  p <- do.call("rbind", g)
  p$Library <- splitter(rownames(p), "\\.", 1)
  
  p <- p[,c("Library", "cr_call", "TargetCat", "umi_count")]
  # p <- melt(p, id.vars = "Library")
  
  pdf(file = "../../Results/1_GuideAssignment/Library-level Statistics.pdf", height = 4, width = 6)
  ggplot(p, aes(x = Library, fill = TargetCat)) +
    geom_bar() +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Number of distinct barcode-guide assignments") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom", panel.border = element_blank(), axis.line.y = element_line())  
  
  ggplot(p, aes(x = umi_count, fill = Library)) +
    geom_bar(colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Number of UMIs supporting barcode-guide assignments") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          legend.position = "right", panel.border = element_blank(), axis.line.y = element_line())  
  
  ggplot(p, aes(fill = cr_call, x = Library)) +
    geom_bar(colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Library") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom", panel.border = element_blank(), axis.line.y = element_line())  
  dev.off()  
  
  pdf(file = "../../Results/1_GuideAssignment/Library-level Statistics (3+ UMIs).pdf", height = 4, width = 6)
  p2 <- p[which(p$umi_count > 2),]
  ggplot(p2, aes(x = Library, fill = TargetCat)) +
    geom_bar() +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Number of distinct barcode-guide assignments") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom", panel.border = element_blank(), axis.line.y = element_line())  
  
  ggplot(p2, aes(x = umi_count, fill = Library)) +
    geom_bar(colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Number of UMIs supporting barcode-guide assignments") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          legend.position = "right", panel.border = element_blank(), axis.line.y = element_line())  
  
  ggplot(p2, aes(fill = cr_call, x = Library)) +
    geom_bar(colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Library") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom", panel.border = element_blank(), axis.line.y = element_line())  
  dev.off()  
  
  
## Plot (non-negative) guide-level analysis
  pdf(file = "../../Results/1_GuideAssignment/Guide-level Statistics (3+ UMIs).pdf", height = 4, width = maxw)
  
  p <- do.call("rbind", g)
  p <- p[which(p$umi_count > 2),]
  p$Library <- splitter(rownames(p), "\\.", 1)
  
  p <- p[,c("cell", "GuideID", "TargetID", "TargetCat", "TargetCt", "Library")]
  
  # p <- melt(p, id.vars = "Library")
  p$ID <- paste0(p$TargetID, "_", p$TargetCt)
  
  p$ID <- factor(p$ID, levels = sort(unique(paste0(ps.whitelist$TargetID, "_", ps.whitelist$Celltype)))) # ensures that missing guides are annotated with zero
  
  ## On guides per cell
    y <- as.data.frame(table(p$cell, p$TargetCt))
    y <- y[-which(y$Freq == 0),]
    y$Bin <- cut(y$Freq, c(0, 1, 5, 10, 1000))
    levels(y$Bin) <- c("1", "2-5", "6-10", ">11")
    
    ggplot(y, aes(x = Freq, fill = Var2)) +
      geom_bar() +
      scale_fill_manual(values = ct.col) +
      facet_wrap(~Var2) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      labs(y = "Count", x = "Number of Guides Per Cell") +
      theme(panel.border = element_blank(), legend.position = c(0.85,085),
            legend.title = element_blank(), axis.line.y = element_line(), panel.grid.major.x = element_blank())  
    
    ggplot(y, aes(x = Bin, fill = Var2)) +
      geom_bar() +
      facet_wrap(~Var2) +
      scale_fill_manual(values = ct.col) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      labs(y = "Count", x = "Number of Guides Per Cell") +
      theme(panel.border = element_blank(), legend.position = c(0.85,085),
            legend.title = element_blank(), axis.line.y = element_line(), panel.grid.major.x = element_blank())  
    
  ## On cells per guide
    # tabulate
    y <- as.data.frame(table(p$ID, p$TargetCt))
    y <- dcast(y, Var1~Var2)  
    
    y$Total <- apply(y[,-1], 1, max)
    y$Ct <- splitter(as.character(y$Var1), "_", 2)
    
   
    # plot density
    ggplot(y, aes(x = Total, colour = Ct, fill = Ct)) +
      geom_density(alpha = 0.3) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      scale_fill_manual(values = ct.col) +
      scale_colour_manual(values = ct.col) +
      labs(y = "Density", x = "Number of Cells Per Guide") +
      theme(panel.border = element_blank(), axis.line.y = element_line(), 
            legend.title = element_blank(), legend.position = c(0.85,085))  
    
    # plot bin
    y$Bin <- cut(y$Total, c(-1, 0, 20, 50, 100, 150, 200, 1000, 1000000))
    levels(y$Bin) <- c("0", "1-20", "21-50", "51-100", "101-150", "151-200", "200-1000", ">1000")
    ggplot(y, aes(x = Bin, fill = Ct)) +
      geom_bar(position = "dodge") +
      facet_wrap(~Ct) +
      scale_fill_manual(values = ct.col) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      labs(y = "Count", x = "Number of Cells Per Guide") +
      theme(panel.border = element_blank(), axis.line.y = element_line(), legend.position = c(0.85,085),
            panel.grid.major.x = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))  
    
  dev.off() 
    
  
## Logical annotation for whether each cell has a given guide 
  min.umi <- 3
  for (lib in samples) {
    print(lib)
    
    # collect variables
    x <- ps.whitelist[which(ps.whitelist$Celltype == splitter(lib, "_", 1)),]
    y <- g[[lib]]
    z <- obj[[lib]]@meta.data
    
    # threshold assignments to those with at least min.umi
    y <- y[which(y$umi_count >= min.umi),]
    
    # match
    z$AnyGuide <- rownames(z) %in% y$cell # annotate cells that have at least one guide
    z[,unique(x$TargetID)] <- FALSE # create columns for each of the guide targets (note: this therefore pools guides with the same target for power)
    
    use <- rownames(z)[which(z$AnyGuide)]
    for (j in use) {
      w <- y$TargetID[which(y$cell == j)]
      w <- unique(w)
      if (anyNA(w)) w <- w[-which(is.na(w))]
      
      z[j,w] <- TRUE # for the cell in the loop, sets TRUE to the columns which are hits
    }
    
    obj[[lib]]@meta.data <- z
  }
  
## Create pooled object
  ## NHA
    # merge
    nha <- merge(x = obj$NHA_1, y = obj[2:8], add.cell.ids = samples[1:8])
    
    # restrict to cells with a guide
    nha <- subset(nha, cells = which(nha$AnyGuide))
  
    # normalise expression
    nha <- norm.bySeurat(nha, min.nUMI = 1000)
    
    # cluster
    nha <- RunPCA(nha)
    ElbowPlot(nha, ndims = 50) 
    nha <- clust.bySeurat(nha, dim = 20) # 20 is a reasonable elbow
    
    pdf(file = "QC/Pooled Non-Integrated tSNE - NHA.pdf", height = maxw, width = maxw)
    DimPlot(nha, group.by = "Phase", pt.size = 0.2)
    DimPlot(nha, group.by = "orig.ident", pt.size = 0.2)
    dev.off()  
    
    # save
    save(nha, file = "NHA Pooled.rda")
    
  ## SY5Y
    # merge
    sy <- merge(x = obj$SY5Y_1, y = obj[10:18], add.cell.ids = samples[9:18])
    
    # restrict to cells with a guide
    sy <- subset(sy, cells = which(sy$AnyGuide))
    
    # normalise expression
    sy <- norm.bySeurat(sy, min.nUMI = 1000)
    
    # cluster
    sy <- RunPCA(sy)
    ElbowPlot(sy, ndims = 50) 
    sy <- clust.bySeurat(sy, dim = 20) # 20 is a reasonable elbow
    
    
    pdf(file = "QC/Pooled Non-Integrated tSNE - SY5Y.pdf", height = maxw, width = maxw)
    DimPlot(sy, group.by = "Phase")
    DimPlot(sy, group.by = "orig.ident")
    dev.off()  
    
    # save
    save(sy, file = "SY5Y Pooled.rda")
    