################################################################################################################################ #
#Intersection all NHA intergenic Peaks detected with TAD boundaries, super enhancers and validated enhancers 
################################################################################################################################ #

remove(list=ls())
gc()
library(ggplot2)
library(scales)
library(cowplot)
library(reshape2)
library(readxl)
library(tidyr)
library(liftOver)
library(rtracklayer)
library(rcartocolor)

##Load predefined functions
#source("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Scripts/Functions.R")  
source("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/2.Scripts/Functions.R")

##Set local directory
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub")

#################################################################################################################################
## Intersect EGPs with TAD boundary (iPSC-derived astrocyte data from Rajarajan et al)
##################################################################################################################################

#load data
#All_EGPs <- read.csv("Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")
All_EGPs <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")
## Create a bed file whose coordinates are the intervening of EGPs
x_JP=All_EGPs[c("Pair", "Enh.Pos", "TSS", "Enh.start", "Enh.end")]
x_JP$chr= splitter(x_JP$Enh.Pos, ":", 1)
x_JP$Coord1 = x_JP$TSS

# get the end, which is Enh's end or start (whichever is closer!)
##Take coordinate 2 from my dataset
x_JP$use.left <- (abs(x_JP$TSS - x_JP$Enh.start)) < (abs(x_JP$TSS - x_JP$Enh.end))
x_JP$Coord2=ifelse(x_JP$use.left==TRUE, x_JP$Enh.start, x_JP$Enh.end)

# the start coordinate is whichever is smaller of Coord1 and Coord2
x_JP$start=apply(x_JP[,c("Coord1", "Coord2")], 1, min)
x_JP$end=apply(x_JP[,c("Coord1", "Coord2")], 1, max)

# write corresponding bed file 
x_JP <- x_JP[,c("chr", "start", "end", "Pair")]
write.bed(x_JP, "Scripts/6.Predictive_models/Sanity_checks/Allintergenic_pair_Intervening_Coord.bed")

###run beedtols
call=paste("intersectBed",
      "-a", "Scripts/6.Predictive_models/Sanity_checks/Allintergenic_pair_Intervening_Coord.bed",
      "-b", "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/TAD_SelfIntersection_In.bed", #  TAD coordinates with an annotation
      "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
      ">", "Scripts/6.Predictive_models/Sanity_checks/Tad_Vs_Allintergenic_Overlap.bed")
system(call, intern = FALSE, wait = TRUE)
## Read in
tad_overlap <- read.delim("Scripts/6.Predictive_models/Sanity_checks/Tad_Vs_Allintergenic_Overlap.bed", header = FALSE)
colnames(tad_overlap) <- c("EGP_Chr", "EGP_Start", "EGP_End", "EGP", "TAD_Chr", "TAD_Start", "TAD_End", "TAD_ID", "Overlap_bp")

## Compute how much of the EGP is contained within a tad
tad_overlap$EGP_Distance <- (tad_overlap$EGP_End - tad_overlap$EGP_Start)
tad_overlap$Overlap_Frac <- tad_overlap$Overlap_bp / tad_overlap$EGP_Distance

## Now categorise each EGP
pair2tad_JP <- All_EGPs[,c("Pair","pass_rf", "Enh.Pos")]

# define basic true/false statements to determine categories
within <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac == 1)] 

cross <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac < 1)] 

# is the EGP fully contained within a single TAD?
pair2tad_JP$WithinTAD <- pair2tad_JP$Pair %in% within

# does the EGP cross a tad boundary?
pair2tad_JP$CrossTAD <- pair2tad_JP$Pair %in% cross

# Define combinations
pair2tad_JP$WithinAndCrossTAD <- pair2tad_JP$WithinTAD & pair2tad_JP$CrossTAD
pair2tad_JP$WithinNotCrossTAD <- pair2tad_JP$Pair %in% setdiff(within, cross)
pair2tad_JP$CrossNotWithinTAD <- pair2tad_JP$Pair %in% setdiff(cross, within)
pair2tad_JP$NoTAD <- rowSums(pair2tad_JP[,4:8]) == 0 

## Save data
write.csv(pair2tad_JP, "Scripts/6.Predictive_models/Sanity_checks/Pair-TAD Annotation_allintergenic.csv")


##################################################################################################################################
# Prepare bed file on all NHA intergenic peaks detected using now only Peak coordinates (presumably enhancers)
##################################################################################################################################

## Create a bed file with enhancers coordinates; enh start point and enh end point
####Preparing bed file with all  intergenic ATAC peaks
All_EGPs <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")
All_EGPs=All_EGPs[c("Enh", "Enh.Pos", "TSS", "Enh.start", "Enh.end")]
All_EGPs$chr= splitter(All_EGPs$Enh.Pos, ":", 1)

# write bed file
All_EGPs <- unique(All_EGPs[,c("chr", "Enh.start", "Enh.end", "Enh")])
write.bed(All_EGPs, "Scripts/6.Predictive_models/Sanity_checks/UniqueIntergenicPeaks.bed")


################################################################################################################################################################################################
# Intersect intergenic peaks with super-enhancers (NHA data from Hnisz et al 2013
###############################################################################################################################################################################################
## Packages, functions, and libraries
library(Rsamtools)
library(rtracklayer)

## Defining intersect function with a “left outer join”. that is, for each feature in A report each overlap with B. if no overlaps are found, report a NULL feature for B
# useful for binary calls
loj <- function(a, b, out, fun = "-loj") {
  call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                fun, 
                ">", out)
  
  system(call, intern = FALSE, wait = TRUE)    
  print("Complete!")
}


## Load bed superenhancer data
hnisz_dir_in <- "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/PROCESSED/liftOver_hg38/Astrocytes.hg38.super.bed" # note: this has super

## Run bedtools and intersect NHA intergenic peaks with super enhancers 
UniqueIntergenicPeaks.bed <- "Scripts/6.Predictive_models/Sanity_checks/UniqueIntergenicPeaks.bed"
call=paste("intersectBed",
           "-a", UniqueIntergenicPeaks.bed,
           "-b", hnisz_dir_in,
           "-loj", 
           ">", "Scripts/6.Predictive_models/Sanity_checks/SuperEnh_Vs_UniqueIntergenicPeaks.bed")

system(call, intern = FALSE, wait = TRUE)

#Load intersect file
read.bed <- function(dir) read.delim(dir, header = FALSE)
hnisz_overlap <- read.bed("Scripts/6.Predictive_models/Sanity_checks/SuperEnh_Vs_UniqueIntergenicPeaks.bed")

# filter to overlapping intergenic peaks with superenhancers
hnisz_overlap <- hnisz_overlap[which(hnisz_overlap$V9 != -1),]

# Save data
write.csv(hnisz_overlap, "Scripts/6.Predictive_models/Sanity_checks/SuperEnh_Vs_UniqueIntergenicPeaks.csv")


################################################################################################################################################################################################
# Intersect intergenic peaks with experimentally-validated enhancers in K562 cells from Yao 2022
###############################################################################################################################################################################################

## Read enhancers data from Yao and preparing corresponding bed file
#dir.create("Enh_Yao")
known_data <- read_xlsx("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/Yao_NatBiotech2022/41587_2022_1211_MOESM3_ESM.xlsx",
                        sheet = "SupplementaryTable2",
                        skip = 1)
write.bed(known_data, "Scripts/6.Predictive_models/Sanity_checks/ST2_KnownEnhancers.bed")

## Intersect NHA Intergenic peaks with Yao Superenhancers
# note that the above has no window to either side; this consistent with other analyses in this study
call=paste("intersectBed",
           "-a", UniqueIntergenicPeaks.bed,
           "-b", "Scripts/6.Predictive_models/Sanity_checks/ST2_KnownEnhancers.bed",
           "-wa", "-wb", 
           ">", "Scripts/6.Predictive_models/Sanity_checks/UniqueIntergenicPeaks_vs_YaoEnhancers.bed")

system(call, intern = FALSE, wait = TRUE)


## Save overlappin NHA with Yao superenhancers as csv file
known <- read.bed("Scripts/6.Predictive_models/Sanity_checks/UniqueIntergenicPeaks_vs_YaoEnhancers.bed")
write.csv(known, "Scripts/6.Predictive_models/Sanity_checks/UniqueIntergenicPeaks_vs_YaoEnhancers.csv")

