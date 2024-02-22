setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable8_Beluga/")
options(stringsAsFactors = FALSE)
source("../../../FullScale/Scripts/Functions.R")
# res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
guides <- read.csv("../../../FullScale/Data/Whitelists/NHA Enh Library.csv")

## 8A: known variants
  ## Read in
    kn <- read.csv("../../../EnhancerPredictionModels/Results/Beluga/KnownVariants/BelugaDiseaseScores.csv")

  ## Filter columns
    # add coordinates of enhancer
    m <- match(kn$Enh, guides$TargetID)
    kn$Enh_Coord <- guides$TargetCoord[m]
    kn$Enh_Coord <- splitter(kn$Enh_Coord, "-exp-", 1)
    
    # rearrange
    kn <- kn[,c("Enh", "Enh_Coord", "Hit",  "SNP", "SNP.pos", "KnownVariants_DIS_max_score", "KnownVariants_Escore")]
    
    # colnames
    colnames(kn)[3:7] <- c("Hit_Enhancer", "rsID", "SNP_position", "DIS_max_score", "Escore_recalculated") 
    
    # filter column
    kn <- kn[,-which(colnames(kn) == "Escore_recalculated")]

  ## Output
    write.csv(kn, file = "8A_SNPs.csv", row.names = FALSE)
    
## 8B: ISM
  ## Read in
    bel <- read.csv("../../../EnhancerPredictionModels/Results/Beluga/BelugaVariants/BelugaVariants.csv")
  
  ## Filter columns
    # add coordinates of variant
    bel$Variant_Position <- paste0("chr", bel$Chr, ":", bel$end)  # using "pos" rather than "end", as this matches to our SNP annotation in other dataframes
    
    # add coordinates of enhancer
    m <- match(bel$Enh, guides$TargetID)
    bel$Enh_Coord <- guides$TargetCoord[m]
    bel$Enh_Coord <- splitter(bel$Enh_Coord, "-exp-", 1)
    
    # remove known variants, as those appear in the previous table
    bel <- bel[,-grep("Known|SNP", colnames(bel))]
    
    # remove escore
    bel <- bel[,-grep("Escore", colnames(bel))]
    
    # remove extraneous annotation
    bel <- bel[,-which(colnames(bel) %in% c("Enh.pos", "Enh.size", "chrom", "end", "pos"))]
    
    # remove disease p-value (as it's the -log10 of the score)
    bel <- bel[,-which(colnames(bel) == "DIS_max_score_p")]
    
  ## Adjust column names
    colnames(bel) <- c("Enh", "Variant_PositionRelative", "ReferenceAllele", "MaxDiseaseScore", "MaxDiseaseScore_Significant", "Variant_Position", "Enh_Coord")

 ## Reorder
    bel <- bel[,c(1,7,6,2,3,4,5)]
    
  ## Save
    write.csv(bel, "8B_ISM.csv", row.names = FALSE)
  

  
## Read in Z scores
  z <- read.delim("../../../EnhancerPredictionModels/Results/Beluga/BelugaVariants/Webtool_Variants/1f9426dd-3bd4-48a8-b44a-84fbe98ef4fa_8_posplus1_FEATURE_zscore.tsv")
  z <- z[,c(1:8, grep("Astrocytes", colnames(z)))]
  