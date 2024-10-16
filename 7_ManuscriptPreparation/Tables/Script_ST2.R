
library(tidyr)


## Define table number
  numb <- 2
    
## Set directory
  dirs <- list.files("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/") 
  dir <- dirs[grep(numb, dirs)]
  setwd(paste0("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/", dir))

## Load data
  df <- read.csv("../../../FullScale/Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  df <- df[,-which(colnames(df) == "Celltype")] # extraneous
  
## Sort the enhancer guides by ID
  w <- which(df$TargetCat == "Enh")
  s <- gsub("Enh", "", df$TargetID[w])
  s <- as.numeric(s)

  df[w,] <- df[order(s),]
  
## Add a column showing the region used as a search space for the guide
  df$GuideSearchSpace <- df$TargetCoord
  
  used <- read.csv("../../../FullLibrary_Selection/Results/FromSefi/Final Guides, Short Regions.csv")
  coord <- strsplit(used$Input, "-exp-")
  coord.expanded <- sapply(coord, "[", 2)
  coord.orig <- sapply(coord, "[", 1)
  m <- match(coord.orig, df$GuideSearchSpace)
  
  df$GuideSearchSpace[m] <- coord.expanded

## Add a column for the transduction pool
  df$LentiviralPool <- "."
  
  # the enhancer pool
  df$LentiviralPool[which(df$TargetCat == "Enh")] <- "Test Pool"
  df$LentiviralPool[grep("Neg_E", df$GuideID)] <- "Test Pool"
  
  # the positive control pool
  w <- which(df$TargetCat == "Promoter")
  df$LentiviralPool[w] <- "Positive Control"
  df$TargetCoord[w] <- paste0("TSS_", df$TargetID[w])
  df$GuideSearchSpace[w] <- df$TargetCoord[w]
  
  # the negative control pool
  df$LentiviralPool[setdiff(which(df$TargetCat == "Negative"),  grep("Neg_E", df$GuideID))] <- "Negative Control"
  
  
  # common to the enhancer and positive pools
  g <- grep("Promoter", df$TargetCat)
  g <- g[1:20]
  df$LentiviralPool[g] <- "Test Pool & Positive Control"
  

## Finalise column names
  colnames(df) <- c("Guide ID",
                    "Guide Sequence",
                    "Guide Binding Site",
                    "Target Classification",
                    "Target ID",
                    "Target Coordinates",
                    "Target Coordinates (Extended)",
                    "Lentiviral Pool")
  
## Write
  write.csv(df, paste0("ST", numb, ".csv"), row.names = FALSE)
  
## Legend
  nl <- "\n"
  leg <- paste0("Supplementary Table ", numb, ": List of sgRNAs.", nl, nl,
                "Guide ID: ID of guide.", nl,
                "Guide Sequence: 20bp protospacer of guide. Note that positive controls (from Gasperini et al., 2019) are 19bp.", nl,
                "Guide Binding Site: Coordinates of guide binding.", nl,
                "Target Classification: Category of genomic element targeted.", nl,
                "Target ID: ID of target genomic element.", nl,
                "Target Coordinates: Coordinates of target genomic element. Positive controls guides (from Gasperini et al., 2019) target the TSS of the gene, as described in the original publication.", nl,
                "Target Coordinates (Extended): Coordinates of extended target genomic element, where elements <200bp in length were expanded to 200bp when selecting guides.", nl,
                "Lentiviral Pool: Lentiviral transduction pool in which the guide was included.")
  write.table(leg, file = "Legend.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  