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
  
## A threshold
  min.umi <- 3

################################################################################################################################ #
## Load processed libraries ----
  
  
load("../../Data/Preprocessed/NHA Pooled (No Guide Annotation).rda")

  
################################################################################################################################ #
## Read in whitelists ----
  
  
## Protospacers
  x <- read.csv("../../Data/Whitelists/Neg Library.csv")
  x <- data.frame(GuideID = paste0("Neg_", 1:nrow(x)),
                  GuideSequence = x$GuideSequence,
                  GuideCoord = NA,
                  TargetCat = "Negative",
                  TargetID  = x$TargetID,
                  TargetCoord = NA)
  x$TargetID[which(is.na(x$TargetID))] <- "NoTarget"
  
  y <- read.csv("../../Data/Whitelists/Pos Library.csv")
  y <- data.frame(GuideID = y$GuideID,
                  GuideSequence = y$GuideSequence,
                  GuideCoord = ".",
                  TargetCat = "Promoter",
                  TargetID  = y$TargetID,
                  TargetCoord = ".")
  y$GuideID <- paste0("Pos_", y$GuideID)
  y$GuideID <- gsub("_A$", "_g1", y$GuideID)
  y$GuideID <- gsub("_B$", "_g2", y$GuideID)
  
  # z <- read.csv("../../Data/Whitelists/NHA Enh Library.csv") # note: this does not include the Negs that are in this library, but that is inconsequential as they appear in your variable x anyway
  z <- read.csv("../../../TransferFromRNA/FullLibrary_Selection/Results/Final_List/SR_final_library_design/NHA Enh Library_with Neg Cont.csv")
  g <- grep("NEG", z$GuideID)
  z$GuideID[g] <- paste0("Neg_E_", 1:length(g))
  z$TargetCat[g] <- "Negative"
  z <- z[which(z$TargetCat != "Promoter"),] # as these are also in the positive library
  
  guide.list <- rbind(z, y, x)
  guide.list$Celltype <- "NHA"
  
  # fix an artefact in $TargetCoord
  guide.list$TargetCoord <- splitter(guide.list$TargetCoord, "-exp-", 1)
  
  # save
  write.csv(guide.list, file = "../../Data/Whitelists/Protospacer Whitelist (NHA).csv")
  

################################################################################################################################ #
## Read in assignments ----
  
g <- list()
  
## Loop
  for (j in samples) {
    
    print(j)
    
    # read in
    x  <- read.table(paste0("GuideAssignment/", j, "_GuideAssignments.txt"), sep = "\t", header = TRUE) 
    
    # is it a called barcode from cellranger?
    y <- colnames(nha)[which(nha$Library == j)] # barcodes for a given library
    y <- gsub(paste0(j, "_"), "", y) # removes the library tag from the barcode
    x$cr_call <- x$cell %in% y
    
    # what protospacer is it targeting? old version based on match
    ps <- guide.list$GuideSequence
    short.guides <- which(nchar(ps) == 19) 
    ps[short.guides] <- paste0(ps[short.guides], "G") # our pos guides are 19nt instead of 20nt, so to ensure the match works I append the nucleotide that follows in the vector
    ps <- paste0("G", ps) # because the output of getprotospacers.py always leads with a G before the protospacer
    
    
    # m <- match(substr(x$barcode, 2, 1000), substr(guide.list$GuideSequence, 1, 15)) # the match removes the first nucleotide (which is nearly always G. might be a cloning artifact?)
    m <- match(x$barcode, ps)
    x$GuideID <- guide.list$GuideID[m]
    x$TargetCat <- guide.list$TargetCat[m]
    x$TargetID <- guide.list$TargetID[m]
    x$TargetCt <- guide.list$Celltype[m]

    # output
    g[[j]] <- x
    
  }

  save(g, file = "GuideAssignment/GuideAssignments.rda")
  
################################################################################################################################ #
## Plot library-level analysis ----
  
## Setup dataframe 
  # wrangle
  p <- do.call("rbind", g)
  p$Library <- splitter(rownames(p), "\\.", 1)

  
## Functions  
  # number of distinct guides in a library
  guidePlot_distinct <- function(x) {
    ggplot(x, aes(x = Library, fill = TargetCat)) +
    geom_bar() +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Number of distinct barcode-guide assignments") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom", panel.border = element_blank(), axis.line.y = element_line())  
  }
  
  # umis supporting each barcode
  guidePlot_nUMIs <- function(x) {
    ggplot(x, aes(x = umi_count, fill = Library)) +
      geom_bar(colour = "black") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      theme_bw() +
      labs(y = "Count", x = "Number of UMIs supporting barcode-guide assignments") +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
            panel.grid.major.x = element_blank(),
            legend.position = "right", panel.border = element_blank(), axis.line.y = element_line())  
  }
  
  # fraction of assignments in cells called by cr
  guidePlot_called <- function(x) {
    ggplot(x, aes(fill = cr_call, x = Library)) +
    geom_bar(colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Library") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom", panel.border = element_blank(), axis.line.y = element_line())  
  }
  
## Basic library-level statistics
  pdf(file = "GuideAssignment/Library-level Statistics (No UMI Threshold).pdf", height = 4, width = 6)
  guidePlot_distinct(p)
  guidePlot_nUMIs(p)
  guidePlot_called(p)
  dev.off()
  
## Plot effect of UMI threshold on stringency of results
  q <- data.frame(Count = p$umi_count,
                  Call = p$cr_call,
                  Whitelist = !(is.na(p$GuideID)))
  q <- q[-which(q$Count == 0),]
  
  # calculate the percentage of TRUE calls and whitelists as a function of count
  
  thresh.check <- function(thresh, x = q) {
    r <- x[which(x$Count >= thresh),]
    n <- nrow(r)
    res <- data.frame(Thresh = thresh,
                      Call = length(which(r$Call)) / n,
                      Whitelisted = length(which(r$Whitelist)) / n,
                      Both = length(which(r$Call & r$Whitelist)) / n,
                      Total = length(which(r$Call & r$Whitelist)))
    rownames(res) <- paste0("Thresh_", thresh)
    return(res)
  }
  
  q <- lapply(c(1,2,3,4,5,10,20), thresh.check)
  q <- do.call("rbind", q)
  q <- melt(q, id.vars = c("Thresh", "Total"))
  levels(q$variable) <- c("Valid Barcode", "Valid Protospacer", "Valid Barcode & Protospacer")
  
  pdf(file = "GuideAssignment/Threshold Justification.pdf", height = 4, width = maxw)
  ggplot(q, aes(fill = as.factor(Thresh), y = value, x = variable)) +
    geom_col(position = "dodge", colour = "black") +
    theme_bw() +
    scale_fill_carto_d(palette = "Magenta") +
    scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
    labs(y = "Fraction of Barcode-Protospacer Assignments") +
    guides(fill = guide_legend(title = "       Minimum UMI count for\nBarcode-Protospacer Assignment", nrow = 1)) +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(), axis.title.x = invis, legend.position = "bottom")
  dev.off()
  
  
## Filter to high-confidence calls with 3+ UMIs
  p <- p[which(p$umi_count >= min.umi),]
  
## Replot
  pdf(file = "GuideAssignment/Library-level Statistics.pdf", height = 4, width = 6)
  guidePlot_distinct(p)
  guidePlot_nUMIs(p)
  guidePlot_called(p)
  dev.off()
  
  
################################################################################################################################ #
## Guide-level analysis ----
  
  
## Reuse the object "p"

## Add an id, and ensure that missing guides are annotated with zero
  p$ID <- paste0(p$TargetID, "_", p$TargetCt)
  p$ID <- factor(p$ID, levels = sort(unique(paste0(guide.list$TargetID, "_", guide.list$Celltype)))) 
  
  
## On guides per cell
  # wrangle
  y <- as.data.frame(table(p$cell, p$TargetCt))
  y <- y[-which(y$Freq == 0),]
  y$Bin <- cut(y$Freq, c(0, 1, 5, 10, 1000))
  levels(y$Bin) <- c("1", "2-5", "6-10", ">11")
  
  pdf(file = "GuideAssignment/Guides Per Cell.pdf", height = 4, width = maxw)
  ggplot(y, aes(x = Freq)) +
    geom_bar() +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Number of Guides Per Cell") +
    theme(panel.border = element_blank(), legend.position = c(0.85,085),
          legend.title = element_blank(), axis.line.y = element_line(), panel.grid.major.x = element_blank())  
  
  ggplot(y, aes(x = Bin)) +
    geom_bar() +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Number of Guides Per Cell") +
    theme(panel.border = element_blank(), legend.position = c(0.85,085),
          legend.title = element_blank(), axis.line.y = element_line(), panel.grid.major.x = element_blank())  
  dev.off()  
  
## On cells per guide
  y <- as.data.frame(table(p$ID))
  
  pdf(file = "GuideAssignment/Cells Per Guide.pdf", height = 4, width = maxw)
  
  # plot density
  ggplot(y, aes(x = Freq)) +
    geom_density(alpha = 0.3) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Density", x = "Number of Cells Per Guide") +
    theme(panel.border = element_blank(), axis.line.y = element_line(), 
          legend.title = element_blank(), legend.position = c(0.85,085))  
  
  # plot bin
  y$Bin <- cut(y$Freq, c(-1, 0, 20, 50, 100, 150, 200, 1000, 1000000))
  levels(y$Bin) <- c("0", "1-20", "21-50", "51-100", "101-150", "151-200", "200-1000", ">1000")
  ggplot(y, aes(x = Bin)) +
    geom_bar(position = "dodge") +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Number of Cells Per Guide") +
    theme(panel.border = element_blank(), axis.line.y = element_line(), legend.position = c(0.85,085),
          panel.grid.major.x = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))  
  
  dev.off() 
  
  
################################################################################################################################ #
## gRNA expression matrix ----

## Process the tall guide expression matrix
  # filter to reduce burden
  tall <- lapply(g, function(x) {
    x <- x[which(x$cr_call & !(is.na(x$TargetID)) & x$umi_count > 0),]
    x <- x[,c("cell", "GuideID", "umi_count")]
  })
  
  # rename cells to account for library
  for (j in names(tall)) tall[[j]]$cell <- paste0(j, "_", tall[[j]]$cell)
  
  # bind 
  tall <- do.call("rbind", tall)
  
## Initialise a wide expression matrix
  guide.exp.raw <- matrix(ncol = ncol(nha), nrow = nrow(guide.list)) %>% as.data.frame()
  colnames(guide.exp.raw) <- colnames(nha)    
  rownames(guide.exp.raw) <- guide.list$GuideID
  guide.exp.raw[,] <- 0
  
  # in theory, this can be quickly achieved with dcast, but it cannot catch exceptions with zero data
  for (j in colnames(guide.exp.raw)) {
    print(j)
    
    if (j %in% tall$cell) {
      
      x <- tall[which(tall$cell == j),]
      guide.exp.raw[x$GuideID,j] <- x$umi_count 
      
    } else {
    
        guide.exp.raw[,j] <- 0
        
    }
    
    
  }
  
## Create a filtered version
  min.guide.umi <- 3
  
  guide.cpm <- guide.exp <- guide.frac <- list()
  
  for (j in colnames(guide.exp.raw)) {
    print(j)
    
    # get guide UMI, removing those belwo threshold
    x <- guide.exp.raw[,j]
    x[which(x < min.guide.umi)] <- 0
    
    if (sum(x) == 0) {
      
      guide.cpm[[j]] <- guide.exp[[j]] <- guide.frac[[j]] <- 0
      
    } else {
      
      guide.cpm[[j]] <- x / ((nha@meta.data[j,"UMI_Total"] + sum(x)) / 10^6) # cpm
      guide.exp[[j]] <- x # umi
      guide.frac[[j]] <- x / sum(x) # fraction
      
    }
    
    
    
  }

  guide.cpm <- do.call("cbind", guide.cpm)  %>% as.data.frame()  
  guide.exp <- do.call("cbind", guide.exp)  %>% as.data.frame()  
  guide.frac <- do.call("cbind", guide.frac)  %>% as.data.frame()  
  
  
  rownames(guide.cpm) <- rownames(guide.exp) <- rownames(guide.frac) <- rownames(guide.exp.raw)

    
################################################################################################################################ #
## Annotation ---- 
  
# ## Reuse "p" again
#   # quick check just in case the code was missed:
#   p <- p[which(p$umi_count >= min.umi),]
#   
#   # filter to remove na guides
#   p <- p[-which(is.na(p$GuideID)),]
#   
#   # rename cells to contain both barcode and library (i.e., like the nha seurat pooled object)
#   p$cell <- paste0(p$Library, "_", p$cell)
#   
#   # filter to only assignments in your called cells!
#   p <- p[which(p$cell %in% colnames(nha)),]
  

## Add guide assignments
  z <- nha@meta.data
  z$TransductionPool <- z$LogMOI <- z$MOI <- z$AnyGuide <- "." # for stats
  z <- z[,-which(colnames(z) %in% c("nCount_RNA", "nFeature_RNA", "RNA_snn_res.2", "seurat_clusters"))]
  z <- cbind(z, t(guide.frac > 0.02)) # adding TRUE/FALSE for guide assignments
  
## Add target assignments
  groups <- unique(guide.list$TargetID[which(guide.list$TargetCat %in% c("Enh", "Promoter"))])
  for (j in groups) {
    print(j)
    x <- guide.list$GuideID[which(guide.list$TargetID == j)] # get list of guides that constitute a target group
    y <- apply(z[,x], 1, any) # are any of them true?
    
    if (substr(j, 1, 3) == "Enh") {
      z[,j] <- y # output  
    } else {
      z[,paste0("Pos_", j)] <- y # output
    }
    
  }

  
## Resolve those recently added metadata columns
  # moi
  z$MOI <- rowSums(z[,guide.list$GuideID])
  z$LogMOI <- log(z$MOI + 1)
  z$AnyGuide <- z$MOI != 0
  
  ## Transduction pool
    ## For each cell, count its number of different positive, negative, and enhancer guides (pre-pooling)
    e <- guide.list$GuideID[c(which(guide.list$TargetCat == "Enh"), grep("Neg_E_", guide.list$GuideID))] # note that the first 20 pos guides also appear in this pool, but will be ignored for simplicity
    p <- guide.list$GuideID[which(guide.list$TargetCat == "Promoter")]
    n <- guide.list$GuideID[setdiff(which(guide.list$TargetCat == "Negative"), grep("Neg_E_", guide.list$GuideID))]
  
    x <- data.frame(Pri = rowSums(z[,e]),
                    Pos = rowSums(z[,p]),
                    Neg = rowSums(z[,n]))

    ## Crude categorisation
      cat <- rep("Ambiguous", nrow(x))
      cat[which(x$Pri > 0)] <- "Primary" # simply put: has any guide from the primary pool, excluding those 20 pos guides noted above
      cat[which(x$Pos > 0 & ((x$Pri + x$Neg) == 0))] <- "Positive" # simply put: has positive guides only
      cat[which(x$Neg > 0 & ((x$Pri + x$Pos) == 0))] <- "Negative" # simply put: has negative guides only
      cat[which(rowSums(x) == 0)] <- "None" # simply put: has no guides above threshold!
    
    ## Plot
      x$Cat <- factor(cat)
      tab <- table(x$Cat)
      levels(x$Cat) <- paste0(names(tab), " Transduction Pool (n=", tab, ")")
      
      p <- melt(x)
      p$Bin <- cut(p$value, c(-1:10, 15, 20, 25, 50, 1000 ))
      levels(p$Bin) <- as.character(c(0:10, "11-15", "16-20", "21-25", "26-50", "51+"))
      
      pdf(file = "GuideAssignment/Transduction Pool Assignment.pdf", height = 5, width = maxw)
      
      ggplot(p[which(p$value > 0),], aes(x = Bin, colour = variable, fill = variable)) +
        geom_bar(colour = "black") +
        facet_wrap(~Cat, drop = FALSE, scales = "free", ncol = 2) +
        theme_bw() +
        scale_y_continuous(expand = c(0,0)) +
        labs(x = "Number of Guides From a Given Pool Per Cell", y = "Number of Cells") +
        theme(panel.grid.major.x = invis, legend.position = "bottom", axis.text.x = element_text(hjust = 1, angle = 90))
      dev.off()
      
    ## Add
      z$TransductionPool <- cat

## Return metadata
  nha@meta.data <- z
 
  
################################################################################################################################ #
## Refilter cells ----  
  
  
## Here, to remove outliers and potential doublets, filter any cell with > 99th percentile of MOI=1 cells. On a per-library basis.
  outlier.umi.cells <- list()
  
  for (j in unique(z$Library)) {
    k <- z[which(z$Library == j),]
    m <- quantile(k$UMI_Total[k$MOI == 1], 0.99)
    outlier.umi.cells[[j]] <- rownames(k)[which(k$UMI_Total > (m))]
  }
  outlier.umi.cells <- do.call("c", outlier.umi.cells) # a total of 1316 cells
  
  table(substr(outlier.umi.cells, 1, 5)) # the number of cells called as outliers per library
  
  # NHA_1 NHA_2 NHA_3 NHA_4 NHA_5 NHA_7 NHA_8 
  #   173   209   229   213   191   172   130 
  
## Perform filtering
  keep <- which(!(rownames(z) %in% outlier.umi.cells))
  nha <- subset(nha, cells = keep) 
  
  
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
  meta <- as.data.frame(nha@meta.data)
  write.csv(meta, file = "../../Data/Preprocessed/Metadata.csv")
  save(meta, file = "../../Data/Preprocessed/Metadata.rda")
  
  
## Guide Expression Matrices
  save(guide.exp, guide.frac, guide.cpm, guide.exp.raw, file = "../../Data/Preprocessed/Guide Expression Matrices.rda")
  
## As Seurat object in rda
  save(nha, file = "../../Data/Preprocessed/NHA Pooled (Final).rda")
  
    