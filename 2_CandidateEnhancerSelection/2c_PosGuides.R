## This script selects 250 positive control guides (targeting 125 TSSs) from Gasperini et al., 2019
## These were evaluated in K562 cells in a similar CRISPRi screen
## We selected guides based on the expression of the TSS in NHA cells as well as SY5Y cells
## SY5Y is a neuroblastoma cell line in which we plan to carry out a CRISPRi screen using the same positive control library


setwd("/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/")

## Gasperini's positive controls ("At this threshold, 94% (357 of 381) of TSS-targeting positive controls repressed their associated genes")
gasp <- read.csv("Results/FromSefi/Gasperini Positive Control Guides from Table S1.csv")
gasp$Gene <- gsub("TSS_", "", gasp$Gene)


################################################################################################################################ #
## NHA ----

## Setup
mean.exp.nha <- read.csv("Results/SC_Ast_meanExp.csv")
exp.genes.nha <- mean.exp.nha$X[which(mean.exp.nha$x > 0.3)]
exp.genes.nha <- as.character(exp.genes.nha)


## I note that exp genes is sorted... So I will take sort the guide Matrix by this, and take the top 50 genes 
  # first, get common genes
  common.nha <- exp.genes.nha[which(exp.genes.nha %in% gasp$Gene)]
  
  # get the top 50 highest expressed of these
  # common <- common[1:125] # exploiting the fact that exp.genes.nha and therefore common is sorted high expression, descending
  
  gasp.nha <- gasp[which(gasp$Gene %in% common.nha),]
  
  # for each gene, choose one of the two guides
  # gasp.nha <- gasp.nha[seq(1, nrow(gasp.nha), 2),] # per advice from Sefi, this is being cut; we want at least two guides per gene, helps with the statistics!
  
  # relabel the genes
  gasp.nha$Gene <- paste0(gasp.nha$Gene, c("_A", "_B"))
  
  # save
  # write.csv(gasp.nha, file = "Results/FromSefi/Gasperini Selected Positive Controls (All).csv", row.names = FALSE)

################################################################################################################################ #
## SY5Y ----  
  
  
## Get expression SY5Y
  load("../Bulk_ATAC_RNAseq/GOK8505-GOK8837/GOK8505A2/sy5y.exp.noFilter.rda")
  exp.genes.sy5y <- sy5y.exp[which(sy5y.exp$RPKM > 6),]
  
## Convert to gene symbol
  load("../../BrainCellularComposition/Data/Preprocessed/geneInfo.rda") 
  load("../../BrainCellularComposition/Data/Preprocessed/exonicLength.rda") 
  source("../../BrainCellularComposition/Scripts/Fun_Preprocessing.R")
  
  rownames(exp.genes.sy5y) <-  exp.genes.sy5y$EnsID
  exp.genes.sy5y <- addSymbol(exp.genes.sy5y)
  
## Rank by expression
  exp.genes.sy5y <- exp.genes.sy5y[order(exp.genes.sy5y$RPKM, decreasing = TRUE),]
  
## Get common genes
  exp.genes.sy5y <- rownames(exp.genes.sy5y)
  
  common.sy5y <- exp.genes.sy5y[which(exp.genes.sy5y %in% gasp$Gene)]
  
## Get guide sequences
  gasp.sy5y <- gasp[which(gasp$Gene %in% common.sy5y),]
  

################################################################################################################################ #
## Get unified list ----  
  
## Get intersection of NHA and SY5Y
  common <- intersect(common.nha, common.sy5y)

  positive.controls <- gasp[which(gasp$Gene %in% common),]
  positive.controls$guideID <- paste0(positive.controls$Gene, c("_A", "_B"))

## Annotate with expression ranks
  # nha
  m <- match(positive.controls$Gene, exp.genes.nha)
  positive.controls$NHA.exp.rank <- m # as exp.genes.nha is ordered by highest expression, getting the index is equivalent to the rank
  positive.controls$NHA.exp.rank <- (rank(positive.controls$NHA.exp.rank) + 0.5) / 2 # gets the rank of expression within these common genes
  
  # sy5y
  m <- match(positive.controls$Gene, exp.genes.sy5y)
  positive.controls$SY5Y.exp.rank <- m # as exp.genes.nha is ordered by highest expression, getting the index is equivalent to the rank
  positive.controls$SY5Y.exp.rank <- (rank(positive.controls$SY5Y.exp.rank) + 0.5) / 2 # gets the rank of expression within these common genes

  # combined rank 
  positive.controls$Combined.rank <- positive.controls$NHA.exp.rank + positive.controls$SY5Y.exp.rank
  positive.controls$Combined.rank <- (rank(positive.controls$Combined.rank) + 0.5) / 2
  
## Output
  positive.controls <- positive.controls[order(positive.controls$Combined.rank),] # sorts by combined rank
  
  write.csv(positive.controls, "Results/Positive Controls From Gasperini - All Hits.csv")

