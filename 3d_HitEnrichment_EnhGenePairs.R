## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----

## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/EnhGenePairs/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(readxl)
  library(tidyr)
  library(liftOver)
  library(rtracklayer)
  library(rcartocolor)

  source("../../../Scripts/Functions.R")

## Load
  # expression
  # load("../../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  
  # results from earlier scripts
  res.final <- read.csv("../../2_DE/Enh/Results Final.csv")  
  annot.enh <- read.csv("../Chromatin/Final - Annotation Logical.csv")
  annot.gene <- read.csv("../Genes/Final - Annotation Logical.csv")
  

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  sig.colours <- c("black", "firebrick1")

## Setup dataframe
  ## Generic dataframe to store annotations pertinent at the pair level
    repEGP <- res.final[,c("Pair", "Enh", "Gene", "HitPermissive")]
    rownames(repEGP) <- repEGP$Pair
    colnames(repEGP)[4] <- "Hit"

  ## List of hit pairs
    hit.egp.vector <- repEGP$Pair[which(repEGP$Hit)]
    hit.egp <- repEGP[which(repEGP$Hit),]
  
## Functions
  # intersect with a “left outer join”. that is, for each feature in A report each overlap with B. if no overlaps are found, report a NULL feature for B
  # useful for binary calls
  loj <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-loj", 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
  
  # intersect with full report
  wawb <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-wa", "-wb", 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
 
  # read in a bed file
  read.bed <- function(dir) read.delim(dir, header = FALSE)
  
  # convert coordinate id to enhancer id
  convert.id <- function(x) {
    m <- match(x, guides$TargetCoord)
    y <- guides$TargetID[m]
    return(y)
  }
  
  
## Commonly used paths
  # candidate coordinates
  nha_dir_38 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed" 
  nha_dir_19 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed" 


################################################################################################################################ #
## Comparison to ABC ----
  
# ## Read inABC
#   abc <- read.delim("/mnt/Data1/PROJECTS/Farbod_ABC/ABC_results/predictions/EnhancerPredictions_NHAPeaks.bed", sep="\t")
#   colnames(abc) <- c("NHA_chr", "NHA_start", "NHA_end", "NHA_id", "NHA_score", "NHA_strand",
#                   "ABC_chr", "ABC_start", "ABC_end", "ABC_id", "Gene", "Score", "CellType", "ABC_Score")
#   
#   abc$ABC_id <- gsub("intergenic|", "", abc$ABC_id, fixed = TRUE)
#   abc$NHA_id <- paste0(abc$NHA_chr, ":", abc$NHA_start, "-", abc$NHA_end)
#   abc <- abc[which(abc$NHA_id %in% x$Enh.Pos), ]
# 
# ## Evaluate ABC pairs
#   # get enhancer and pair annotation for ABC
#   m <- match(abc$NHA_id, guides$TargetCoord)
#   abc$Enh <- guides$TargetID[m]
#   abc$Pair <- paste0(abc$Enh, "_", abc$Gene)  
#   
#   # match
#   abc <- abc[which(abc$Pair %in% x$Pair), c("Enh", "Gene", "Pair", "NHA_id", "ABC_id", "ABC_Score")]
# 
#   # save
#   write.csv(abc, file = "ABC Raw.csv")
#   
# ## Plot
#   m <- match(y$Pair, abc$Pair)
#   y$ABC.Score <- abc$ABC_Score[m]
#   y$ABC.Hit <- y$Pair %in% abc$Pair
#   
#   p <- table(y$HitCategory, y$ABC.Hit)
#   p <- p / rowSums(p)
#   p <- data.frame(ABC.Rate = p[,2])
#   
#   p$Cat <- factor(rownames(p), levels = levels(y$HitCategory))
#   sig.colours2 <- c(pal_lancet()(2), "grey80")
#   
#   pdf(file = "ABC Replication Rate.pdf", height = 3, width = 5)
#   ggplot(p, aes(x = Cat, y = ABC.Rate, fill = Cat)) +
#     geom_col(colour = "black") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis, legend.position = "none") +
#     scale_fill_manual(values = sig.colours2) +
#     labs(y = "Fraction of Pairs Called by ABC") +
#     scale_y_continuous(limits = c(0,1), expand = c(0,0))
#   dev.off()  
#   
    
 
  
################################################################################################################################ #
## Comparison to activity correlation-based enhancer-gene pairs  ----
  
  
## From Dong et al, 2022
  # used lasso regression in >1000 samples to link eRNA to mRNA, 500kb window
  
## Process input
  # read in
  dong_lasso <- readxl::read_xlsx("../../../../PublicData/Dong2022_eRNA_sortedBrain/41588_2022_1170_MOESM3_ESM.xlsx", sheet = "Supplementary Tables 5")
  
  # annotate enhancers, and filter to those active in NeuN-
  dong_ids <- readxl::read_xlsx("../../../../PublicData/Dong2022_eRNA_sortedBrain/41588_2022_1170_MOESM3_ESM.xlsx", sheet = "Supplementary Tables 2")
  m <- match(dong_lasso$enhancer_name, dong_ids$gene_name) # poorly-named column for the latter, but it is the enhancer id
  dong_lasso$Source <- dong_ids$source[m]
  dong_lasso <- dong_lasso[which(dong_lasso$Source != "neuron"),] # from ~35k to ~18k
  
  # make bed
  dong_bed <- gsub("glia_", "", dong_lasso$enhancer) %>% strsplit(., "_")
  dong_bed <- data.frame(chr = sapply(dong_bed, "[", 1),
                         start = sapply(dong_bed, "[", 2),
                         end = sapply(dong_bed, "[", 3),
                         id = dong_lasso$enhancer_name,
                         Gene = dong_lasso$gene_name,
                         Coef = dong_lasso$coef)  
  
  dong_in <- "../../../../PublicData/Dong2022_eRNA_sortedBrain/EnhGeneLasso.bed"
  dong_out <- "DongOverlap.bed"
  write.bed(dong_bed, dong_in)  

## Intersect
  wawb(nha_dir_38, dong_in, dong_out) # I believe it is hg38. my evidence is indirect: 1) the paper states it needed to liftover to hg19, and 2) the supplementary code uses FANTOM5 hg38 coordinates
  
## Read in
  
  
## Wrangle into a more suitable dataframe
  # ## First, subset to your peaks that are in dong
  #   y <- dong_ids[which(dong_ids$source != "neuron"),]
  #   write.bed(y, "../../../../PublicData/Dong2022_eRNA_sortedBrain/Peaks_NonNeuron.bed")
  #   wawb(nha_dir_38, "../../../../PublicData/Dong2022_eRNA_sortedBrain/Peaks_NonNeuron.bed", "../../../../PublicData/Dong2022_eRNA_sortedBrain/Peaks_NonNeuron_Intersect.bed")
  #   y <- read.bed("../../../../PublicData/Dong2022_eRNA_sortedBrain/Peaks_NonNeuron_Intersect.bed")
  #   y <- y[,1:4] %>% unique()
  #   y$Enh <- convert.id(y$V4)
  #   
  #   # add screen-derived hit genes
  #   y$ScreenHit <- y$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
  #   y$ScreenGene <- NA
  #   
  #   for (j in y$Enh) {
  #     if (j %in% hit.egp$Enh) {
  #       z <- hit.egp[which(hit.egp$Enh == j),]
  #       y$ScreenGene[which(y$Enh == j)] <- paste(z$Gene, collapse = "/")
  #     }
  #   }
  #   
  #   # filter columns
  #   y <- y[,5:7]
    
  ## Read in lasso-derived EGPs
    x <- read.delim(dong_out, header = FALSE)
    x <- unique(x)
    colnames(x) <- c("chr", "start", "end", "id", "dong.chr", "dong.start", "dong.end", "dong.id", "dong.gene", "dong.coef")
    x$Enh <- convert.id(x$id)
    x <- x[,c("Enh", "dong.id", "dong.gene", "dong.coef")]
    
    # # add to above dataframe "y"
    # y$LassoHit <- y$Enh %in% x$Enh
    # y$LassoGene <- NA
    # y$LassoCoef <- NA
    # 
    # for (j in y$Enh) {
    #   if (j %in% x$Enh) {
    #     z <- x[which(x$Enh == j),]
    #     y$LassoGene[which(y$Enh == j)] <- paste(z$dong.gene, collapse = "/")
    #     y$LassoCoef[which(y$Enh == j)] <- paste(z$dong.coef, collapse = "/")
    #   }
    # }
    
  ## What fraction of lasso-predicted links from our enhancers have experimental support?
    # z <- x
    
    # is the enhancer a screen hit?
    x$ScreenHit <- x$Enh %in% hit.egp$Enh
    
    # is the screen gene consistent?
    x$ScreenConsistent <- NA
    
    for (j in 1:nrow(x)) {
      if (x$ScreenHit[j]) {
        y <- hit.egp[which(hit.egp$Enh == x$Enh[j]),]
        x$ScreenConsistent[j] <- x$dong.gene[j] %in% y$Gene
      }
    }
    
    # categorise
    x$Support <- "."
    x$Support[!(x$ScreenHit)] <- "No Activity"
    x$Support[which(!(x$ScreenConsistent))] <- "Gene Target\nInconsistent"
    x$Support[which(x$ScreenConsistent)] <- "Gene Target\nConsistent"
    x$Support <- factor(x$Support, levels = levels(as.factor(x$Support))[3:1])
    
    
    # plot
    pdf(file = "Dong - Lasso-predicted Target Consistency.pdf", height = 3, width = 5)
    cols <- c("grey80", carto_pal(3, "ArmyRose")[c(1,3)])
    ggplot(x, aes(x = Support, fill = Support)) +
      geom_bar(width = 0.7, colour = "black") +
      scale_fill_manual(values = cols) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", 
            panel.grid.major.x = invis) +
      labs(x = "Activity of Overlapping Peak in Experimental Data", y = "Number of Lasso-predicted EGPs") +
      scale_y_continuous(expand = c(0,0))
    
    ggplot(x, aes(x = Support, fill = Support, y = dong.coef)) +
      geom_quasirandom(shape = 21, colour = "black", size = 2) +
      scale_fill_manual(values = cols) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none", 
            panel.grid.major.x = invis) +
      labs(x = "Activity of Overlapping Peak in Experimental Data", y = "Number of Lasso-predicted EGPs") +
      scale_y_continuous(expand = c(0,0))
    dev.off()
      
    
        

## Herring 2022
  

  
    
################################################################################################################################ #
## TADs ----
  
  
## Read in
  tad <- read.delim("../../../../FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38.bed", header = FALSE)
  
    
## Some characteristics
  ## Size distribution
    p <- data.frame(Size = tad$V3 - tad$V2)
    p$Size <- p$Size / 10^6
    
    pdf(file = "TAD/Size Distribution.pdf", height = 3, width = 3)
    ggplot(p, aes(x = ".", y = Size)) +
      geom_quasirandom(alpha = 0.5) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
      labs(y = "Size Distribution of 1131 TADs (MB)", x = "") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.text.x = invis, panel.grid.major.x = invis)
    dev.off()
  
## Do TADs overlap?
  ## Setup
    x <- tad
    x$ID <- paste0("TAD_", 1:nrow(x))
    write.bed(x, "TAD/SelfIntersection_In.bed")
    
  ## Intersect with self
    call <- paste("intersectBed",
                  "-a", "TAD/SelfIntersection_In.bed",
                  "-b", "TAD/SelfIntersection_In.bed",
                  "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                  ">", "TAD/SelfIntersection_Out.bed")
    system(call, intern = FALSE, wait = TRUE) 
    
  ## Read in
    self <- read.delim("TAD/SelfIntersection_Out.bed", header = FALSE) # the logic here is to intersect bed with itself
    colnames(self) <- c("A_chr", "A_start", "A_end", "A_ID", "B_chr", "B_start", "B_end", "B_ID", "Overlap")
  
    # remove duplicates
    self <- self[-which(self$A_ID == self$B_ID),] # where the a and b peaks are the same
    self <- self[-which(as.numeric(splitter(self$A_ID, "_", 2)) > as.numeric(splitter(self$B_ID, "_", 2))),] # where the a and b peak pair is the same as a b and a peak pair
  
  ## Calculate the fraction of overlap, where 1 means that a tad is entirely within another
    self$OverlapFraction <- apply(self[,-c(1,4,5,8)], 1, function (x) {
      sizeA <- x[2] - x[1]
      sizeB <- x[4] - x[3]
      frac <- x[5] / min(sizeA, sizeB)
      
    })
  
  ## Plot
    self$OverlapComplete <- factor(self$OverlapFraction == 1)
    levels(self$OverlapComplete) <- c("2 TAD Windows\nPartially Overlap", "Smaller TAD Wholly\nInside Larger TAD")
    pdf(file = "TAD/Self Overlap Distribution.pdf", height = 3, width = 3)
    ggplot(self, aes(x = OverlapComplete)) +
      geom_bar() +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "Count of TAD Intersections") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.title.x = invis)
    dev.off()
  
 
    
    
## Intersect EGPs with TADs
  ## Create a bed file whose coordinates are the intervening section between E and P
    x <- res.final[,c("Pair", "Enh.Pos", "Gene.TSS")]
    
    # get the chromosome
    x$chr <- splitter(x$Enh.Pos, ":", 1)
    
    # get the start, in this case the gene TSS
    x$Coord1 <- splitter(x$Gene.TSS, ":", 2) %>% as.numeric()
    
    # get the end, which is Enh's end or start (whichever is closer!)
    x$Coord2 <- apply(x, 1, function(x) {
      e.coord <- splitter(x[2], ":", 2)
      e.left <- splitter(e.coord, "-", 1) %>% as.numeric() 
      e.right <- splitter(e.coord, "-", 2) %>% as.numeric() 
      t.start <- as.numeric(x[5])
      use.left <- (abs(t.start - e.left)) < (abs(t.start - e.right))
      
      if (use.left) {
        return(e.left)
      } else {
        return(e.right)
      }
    })
    
    # the start coordinate is whichever is smaller of Coord1 and Coord2
    x$start <- apply(x[,c("Coord1", "Coord2")], 1, min)
    x$end <- apply(x[,c("Coord1", "Coord2")], 1, max) # and the opposite for the end
    
    # write
    x <- x[,c("chr", "start", "end", "Pair")]
    write.bed(x, "TAD/Pair_Intervening_Coord.bed")
    
    
  # here, reuse the file from the CTCF section, which contains the intervening section between E and P
  call <- paste("intersectBed",
                "-a", "TAD/Pair_Intervening_Coord.bed",
                "-b", "TAD/SelfIntersection_In.bed", # despite the name, this is simply TAD coordinates with an annotation
                "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                ">", "TAD/Tad_Vs_Pair_Overlap.bed")
  
  system(call, intern = FALSE, wait = TRUE) 
  
  ## Read in
    tad_overlap <- read.delim("TAD/Tad_Vs_Pair_Overlap.bed", header = FALSE)
    colnames(tad_overlap) <- c("EGP_Chr", "EGP_Start", "EGP_End", "EGP", "TAD_Chr", "TAD_Start", "TAD_End", "TAD_ID", "Overlap_bp")
    
  ## How much of the EGP is contained within a tad?
    tad_overlap$EGP_Distance <- (tad_overlap$EGP_End - tad_overlap$EGP_Start)
    tad_overlap$Overlap_Frac <- tad_overlap$Overlap_bp / tad_overlap$EGP_Distance
    
  
  ## Now categorise each EGP
    # first, get all EGPs
    pair2tad <- res.final[,c("Pair", "Gene.Distance.Bin", "HitPermissive")]
    pair2tad$Gene.Distance.Bin <- factor(pair2tad$Gene.Distance.Bin, levels = c("2-10kb", "10-50kb", "50-100kb", "100-500kb"))
    
    # define basic true/false statements to determine categories
    within <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac == 1)] 
    cross <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac < 1)] 
    
    # is the EGP fully contained within a single TAD?
    pair2tad$WithinTAD <- pair2tad$Pair %in% within
  
    # does the EGP cross a tad boundary?
    pair2tad$CrossTAD <- pair2tad$Pair %in% cross
    
    # combinations
    pair2tad$WithinAndCrossTAD <- pair2tad$WithinTAD & pair2tad$CrossTAD
    pair2tad$WithinNotCrossTAD <- pair2tad$Pair %in% setdiff(within, cross)
    pair2tad$CrossNotWithinTAD <- pair2tad$Pair %in% setdiff(cross, within)
    pair2tad$NoTAD <- rowSums(pair2tad[,4:8]) == 0 # fortunately 0, as was defined this way...
  
  ## Save
    write.csv(pair2tad, "TAD/Pair-TAD Annotation.csv")
    
    
    
################################################################################################################################ #
## Concordance in enhancer and gene annotation ----
          
## I have dropped this section from analyses
  # whilst there were many strong coenrichments, with strong OR (5+) and nominally-significant p-values, n is too low to be reliable.
    
  

  
  
  
################################################################################################################################ #
## Save ----
  
write.csv(egp_replicate, file = "Final - Pair Replication.csv")    
write.csv(egp_coannot, file = "Final - Pair Coannotation.csv")
write.csv(coannot_stats, file = "Final - Pair Coannotation Stats.csv")
  