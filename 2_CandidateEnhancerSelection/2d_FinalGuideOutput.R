## This script simply saves and reformats the output of previous scripts in 2_CandidateEnhancerSelection

################################################################################################################################ #
## Setup ----

## Setup
setwd("/Volumes/Data1/PROJECTS/CROPseq/FullLibrary_Selection/")
rm(list = ls())
options(stringsAsFactors = FALSE)



################################################################################################################################ #
## +ve Control Library ----


## Read in
p.cont <- read.csv("Results/Positive Controls From Gasperini - All Hits.csv", row.names = 1)

## Get top 250 guides for top 125 genes!
  # the below rearrangement ensures that, in cases of tied combined rank, the one with the higher rank in Ast is prioritised
  p.cont <- p.cont[order(p.cont$NHA.exp.rank),]
  p.cont <- p.cont[order(p.cont$Combined.rank),]

  pos.lib <- p.cont[1:250,]
  
## Rearrange columns
  pos.lib <- pos.lib[,c(3,1,2,4,5,6)]
  colnames(pos.lib)[1:3] <- c("GuideID", "GuideSequence", "TargetID")
  
  
## Save
  write.csv(pos.lib, file = "Results/Final_List/Positive Control Library.csv", row.names = FALSE)


################################################################################################################################ #
## NHA Library ----

## Read in files
nha.long <- read.csv("Results/NHA/FromSefi/Final Guides, Long Regions.csv")
nha.short <- read.csv("Results/NHA/FromSefi/Final Guides, Short Regions.csv")
p.cont <- read.csv("Results/Positive Controls From Gasperini - All Hits.csv")
n.cont <- NA

## Step 1: combine the enhancer guides
  nha <- rbind(nha.long, nha.short)
  
## Step 2: create informative columns
  ## An enhancer id
    nha$TargetID <- factor(nha$Input)
    levels(nha$TargetID) <- paste0("Enh", 1:length(levels(nha$TargetID)))
    
  ## An enhancer region tag
    nha$TargetCoord <- nha$Target.Alias
    
  ## A guide category
    nha$TargetCat <- "Enh"
    # nha$TargetCat[grep("exp", nha$Input)] <- "EnhShort" # the string exp will find any of the short enhancers
    
  ## Guide coordinates
    x <- nha$Input
    x <- sapply(strsplit(x, "-") , "[", 1)
    y <- sapply(strsplit(x, ":") , "[", 1) # chromosome
    x <- sapply(strsplit(x, ":") , "[", 2)
    x <- as.numeric(x)
    x <- x + nha$sgRNA.Cut.Position..1.based. # start point of the sgRNA
    
    z <- rep(1, nrow(nha))
    z[nha$Orientation == "antisense"] <- -1
    z <- z * 20 # end point of the sgRNA, accounting for strand!
    
    nha$GuideCoord <- paste0(y, ":", x, "-", (x + z))
    
  ## Guide Sequence
    nha$GuideSequence <- nha$sgRNA.Sequence # just renaming
    
  ## Guide ID
    y <- x <- as.character(nha$TargetID) # get the target id of every guide
    for (j in 1:length(x)) {
      z <- x[1:j-1] 
      z <- length(which(z == x[j])) # counts how many guides earlier in the list have the same target
      z <- z + 1 # add 1 to this number
      y[j] <- paste0(y[j], "_g", z) # the guide id is defined as targetID_gZ
    }
    
    nha$GuideID <- paste0(y, "_", nha$GuideCoord)
    
  ## Get relevant columns
    nha <- nha[,c("GuideID", "GuideSequence", "GuideCoord", "TargetCat", "TargetID", "TargetCoord")]
    
    
## On positive controls
  enh.pos <- p.cont
  
  enh.pos <- p.cont[1:20,] # selecting the 20 best guides
  enh.pos <- enh.pos[,c(3,1)]
  colnames(enh.pos) <- c("GuideID", "GuideSequence")
  
  # add info columns
  enh.pos$GuideCoord <- "."
  enh.pos$TargetCat <- "Promoter"
  enh.pos$TargetID <- gsub("_A|_B", "", enh.pos$GuideID)
  enh.pos$TargetCoord <- "."
    
  
## Export
  nha <- rbind(nha, enh.pos)
  write.csv(nha, file = "Results/Final_List/NHA Enh Library.csv", row.names = FALSE)

  
  