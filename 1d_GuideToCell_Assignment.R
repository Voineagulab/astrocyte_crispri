## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----


## Generic
  rm(list = ls()); gc()
  # setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Data/Preprocessed/")
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
  
  
## Barcodes
  # bc.whitelist <- read.table("../Whitelists/3M-february-2018.txt")
  # bc.whitelist <- as.character(bc.whitelist$V1)
  
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
  
  # guide.list <- list()
  # guide.list$NHA <- read.csv("../Whitelists/NHA Enh Library.csv")
  # guide.list$SY5Y <- read.csv("../Whitelists/SY5Y Enh Library.csv")
  # guide.list <- lapply(guide.list, function(z) rbind(z, y, x))
  # guide.list <- do.call("rbind", guide.list)
  guide.list <- rbind(z, y, x)
  guide.list$Celltype <- "NHA"
  guide.list$GuideSequence <- substr(guide.list$GuideSequence, 
                                       start = 1, 
                                       stop = 16) # keeps only the first sixteen nt, as sixteen is the length of protospacer sequence that get_barcodes.py called
  # guide.list <- unique(guide.list)
  
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
    
    # what protospacer is it targeting?
    m <- match(substr(x$barcode, 2, 1000), substr(guide.list$GuideSequence, 1, 15)) # the match removes the first nucleotide (which is nearly always G. might be a cloning artifact?) 
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
    # scale_fill_manual(values = ct.col) +
    # scale_colour_manual(values = ct.col) +
    labs(y = "Density", x = "Number of Cells Per Guide") +
    theme(panel.border = element_blank(), axis.line.y = element_line(), 
          legend.title = element_blank(), legend.position = c(0.85,085))  
  
  # plot bin
  y$Bin <- cut(y$Freq, c(-1, 0, 20, 50, 100, 150, 200, 1000, 1000000))
  levels(y$Bin) <- c("0", "1-20", "21-50", "51-100", "101-150", "151-200", "200-1000", ">1000")
  ggplot(y, aes(x = Bin)) +
    geom_bar(position = "dodge") +
    # facet_wrap(~Ct) +
    # scale_fill_manual(values = ct.col) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(y = "Count", x = "Number of Cells Per Guide") +
    theme(panel.border = element_blank(), axis.line.y = element_line(), legend.position = c(0.85,085),
          panel.grid.major.x = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))  
  
  dev.off() 
    
################################################################################################################################ #
## Annotation ---- 
  
## Create dataframe
  # of guides
  # guide.annot <- matrix(nrow = nrow(guide.list), ncol = ncol(nha)) %>% as.data.frame()
  # colnames(guide.annot) <- colnames(nha)
  # rownames(guide.annot) <- guide.list$GuideID
  
  # of targets
  # target.annot <- matrix(nrow = length(unique(guide.list$TargetID)), ncol = ncol(nha)) %>% as.data.frame()
  # colnames(target.annot) <- colnames(nha)
  # rownames(target.annot) <- unique(guide.list$TargetID)
  

## Reuse "p" again
  # quick check just in case the code was missed:
  p <- p[which(p$umi_count >= min.umi),]
  
  # filter to remove na guides
  p <- p[-which(is.na(p$GuideID)),]
  
  # rename cells to contain both barcode and library (i.e., like the nha seurat pooled object)
  p$cell <- paste0(p$Library, "_", p$cell)
  
  # filter to only assignments in your called cells!
  p <- p[which(p$cell %in% colnames(nha)),]
  
## Setup some new metadata columns in NHA (this initialisation is for better column ordering)
  nha$TransductionPool <- nha$LogGuideUMIs <- nha$GuideUMIs <- nha$LogMOI <- nha$MOI <- nha$AnyGuide <- "."

  nha$AnyGuide <- colnames(nha) %in% p$cell
  
  
## Annotate the cellranger object with true/false for each guide
  # store metadata in new object
  z <- nha@meta.data
  z[,guide.list$GuideID] <- FALSE # create columns for each of the guides
    
  # annotate each cell with guides
  use <- rownames(z)[which(z$AnyGuide)]
  for (j in use) {
    print(which(use == j))
    w <- p$GuideID[which(p$cell == j)]
    w <- unique(w) # not truly necessary for the current loop (using GuideID rather than TargetID)
    if (anyNA(w)) w <- w[-which(is.na(w))]
    
    z[j,w] <- TRUE # for the cell in the loop, sets TRUE to the columns which are hits
  }

  # annotate with MOI
  z$MOI <- rowSums(z[,which(colnames(z) %in% guide.list$GuideID)])
  z$LogMOI <- log(z$MOI + 1)
      
## Annotate the cellranger object with true/false for each target
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
  
  
## Return metadata
  nha@meta.data <- z
  
# save(obj, file = "NHA Individual Libraries.rda")



  
################################################################################################################################ #
## Categorise cells by guide transduction pool ---- 
  

## For each cell, count its number of different positive, negative, and enhancer guides (pre-pooling)
  e <- guide.list$GuideID[which(guide.list$TargetCat == "Enh")]
  p <- guide.list$GuideID[which(guide.list$TargetCat == "Promoter")]
  n <- guide.list$GuideID[which(guide.list$TargetCat == "Negative")]

  x <- data.frame(Enh = rowSums(nha@meta.data[,e]),
                  Pos = rowSums(nha@meta.data[,p]),
                  Neg = rowSums(nha@meta.data[,n]))
  
## Crude categorisation
  cat <- rep("Ambiguous", nrow(x))
  # cat[which(x$Enh > (x$Pos + x$Neg))] <- "Enh" # simply put: has a dominance of enhancer guides
  cat[which(x$Enh > 0)] <- "Primary" # simply put: has a dominance of enhancer guides
  cat[which(x$Pos > 0 & (x$Enh + x$Neg) == 0)] <- "Positive" # simply put: has positive guides only
  cat[which(x$Neg > 0 & (x$Enh + x$Pos) == 0)] <- "Negative" # simply put: has negative guides only
  cat[which(rowSums(x) == 0)] <- "None" # simply put: has negative guides only
  
  # table
  table(cat)
  
  # Ambiguous       Enh       Neg      None       Pos 
  #     132     36296      3554      4125      5139 
  
  # looking at the cells with no cat, seems that many are truly neg or pos but have a noisy assignment, but we'll maintain these definitions
  
  
## Add to metadata
  nha@meta.data$TransductionPool <- cat
  
  
## Plot distribution of counts
  x$Cat <- factor(cat)
  tab <- table(x$Cat)
  levels(x$Cat) <- paste0(names(tab), " Transduction Pool (n=", tab, ")")
  
  p <- melt(x)
  p$Bin <- cut(p$value, c(-1:10, 15, 20, 25, 50, 1000 ))
  levels(p$Bin) <- as.character(c(0:10, "11-15", "16-20", "21-25", "26-50", "51+"))
  
  pdf(file = "GuideAssignment/Transduction Pool Assignment.pdf", height = 5, width = maxw)
  ggplot(p, aes(x = value, colour = variable, fill = variable)) +
    geom_density() +
    facet_wrap(~Cat, scales = "free", ncol = 2) +
    scale_x_continuous(limits = c(0, 25))
  
  ggplot(p[which(p$value > 0),], aes(x = Bin, colour = variable, fill = variable)) +
    geom_bar(colour = "black") +
    facet_wrap(~Cat, drop = FALSE, scales = "free", ncol = 2) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid.major.x = invis, legend.position = "bottom", axis.text.x = element_text(hjust = 1, angle = 90))
  dev.off()
     
## UMAP
  # setup
   p <- data.frame(nha@meta.data, nha@reductions$umap@cell.embeddings)
   umap.plot.2 <- function(meta, facet = TRUE) {
     plot <- ggplot(p, aes_string(x = "UMAP_1", y = "UMAP_2", colour = meta)) +
       geom_point(size = 1) +
       theme_void() +
       scale_colour_carto_c(palette = "Geyser") +
       
       theme(legend.position = c(0.15, 0.15)) 
     
     if (facet) {
       plot + facet_wrap(~orig.ident) 
     } else {
       plot
     }
   }
  
  # transduction pool
  pdf(file = "Pooled/UMAP - Transduction Pool Assignment.pdf", height = maxw, width = maxw)
  umap.plot.2("TransductionPool", facet = FALSE) + scale_colour_lancet() 
  umap.plot.2("TransductionPool", facet = FALSE) + scale_colour_lancet() + facet_wrap(~TransductionPool) + NoLegend()
  dev.off()
    
  # moi
  p[,"log(MOI+1)"] <- log(p$MOI + 1)
  pdf(file = "Pooled/UMAP - MOI.pdf", height = maxw, width = maxw)
  umap.plot("MOI", facet = FALSE)
  umap.plot("log(MOI+1)", facet = FALSE)
  dev.off()
  
# ## In the enhancer transduction pool, evaluate the sensibility of its pos and neg guides
#   # that is, are these pos and neg guides in the enh library whitelist?
#   e <- which(nha$TransductionPool == "Enh")
#   e <- nha@meta.data[e,]
#   
#   e2 <- apply(e[grep("Neg", colnames(e))], 2, table)
#   e2 <- apply(e2, )

 
  
################################################################################################################################ #
## Save ---- 

save(nha, file = "../../Data/Preprocessed/NHA Pooled (Annotated With Guides).rda")
  

  
      
   
    