##############################################################################
#This script Create Suplementary Table 9E
##############################################################################
remove(list=ls()) ; gc()


#Load EGPs data
All_EGPs <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")

#Load CRISPRi screen data
res.final <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv") #Load Crispri hits
Supplementary_Table_9=All_EGPs

##############################################################################
#Sort ttseq data
##############################################################################
#Load intergenic data
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/4_EnhancerTranscription/TTseq/FeatureCounts_intergenic.rda")
intergenic_ttseq=cbind(counts_intergenic$Pos$counts, counts_intergenic$Neg$counts)
colnames(intergenic_ttseq)=c("Pos_ttseq","Pos_rnaseq","Neg_ttseq","Neg_rnaseq")
intergenic_ttseq=as.data.frame(intergenic_ttseq)
intergenic_ttseq$TTseq_Total <- intergenic_ttseq$Pos_ttseq + intergenic_ttseq$Neg_ttseq
intergenic_ttseq$RNAseq_Total <- intergenic_ttseq$Pos_rnaseq + intergenic_ttseq$Neg_rnaseq
intergenic_ttseq$NascentEnriched <- intergenic_ttseq$TTseq_Total > intergenic_ttseq$RNAseq_Total

################################################################################################################################ 
## Binary classification of peaks to be expressing eRNA 
################################################################################################################################
threshes <- c(2, 3, 5, 10)

for (use.thresh in threshes) {
  u <- paste0("Unidirectional_", use.thresh) 
  b <- paste0("Bidirectional_", use.thresh)
  c <- paste0("Category_", use.thresh)
  
  intergenic_ttseq[,u] <- (rowSums(intergenic_ttseq[,c("Pos_ttseq", "Neg_ttseq")] >= use.thresh) >= 1) & intergenic_ttseq$NascentEnriched
  intergenic_ttseq[,b] <- (rowSums(intergenic_ttseq[,c("Pos_ttseq", "Neg_ttseq")] >= use.thresh) == 2) & intergenic_ttseq$NascentEnriched
  
  intergenic_ttseq[,c] <- "Not transcribed"
  intergenic_ttseq[which(intergenic_ttseq[,u]),c] <- "Unidirectional"
  intergenic_ttseq[which(intergenic_ttseq[,b]),c] <- "Bidirectional"
  intergenic_ttseq[,c] <- factor(intergenic_ttseq[,c], levels = c("Not transcribed", "Unidirectional", "Bidirectional"))
}

#Sort columns name
colnames(Supplementary_Table_9)[1]="EnhID"
Supplementary_Table_9$Enh=Supplementary_Table_9$Enh.Pos

###########################
#Provide ttseq information 
############################
Supplementary_Table_9$Enh.TTseqCount_NegStrand=intergenic_ttseq[match(Supplementary_Table_9$EnhID, rownames(intergenic_ttseq)),]$Neg_ttseq
Supplementary_Table_9$Enh.TTseqCount_PosStrand=intergenic_ttseq[match(Supplementary_Table_9$EnhID, rownames(intergenic_ttseq)),]$Pos_ttseq
Supplementary_Table_9$Enh.RNAseqCount_NegStrand=intergenic_ttseq[match(Supplementary_Table_9$EnhID, rownames(intergenic_ttseq)),]$Neg_rnaseq
Supplementary_Table_9$Enh.RNAseqCount_PosStrand=intergenic_ttseq[match(Supplementary_Table_9$EnhID, rownames(intergenic_ttseq)),]$Pos_rnaseq
Supplementary_Table_9$Enh.eRNA=intergenic_ttseq[match(Supplementary_Table_9$EnhID, rownames(intergenic_ttseq)),]$Category_3

#Define Nascen rna as follow:
Supplementary_Table_9$Enh.NascentEnriched <- ifelse(
  rowSums(Supplementary_Table_9[c("Enh.TTseqCount_NegStrand", "Enh.TTseqCount_PosStrand")]) >= 3 &
    (rowSums(Supplementary_Table_9[c("Enh.TTseqCount_NegStrand", "Enh.TTseqCount_PosStrand")]) >
       rowSums(Supplementary_Table_9[c("Enh.RNAseqCount_NegStrand", "Enh.RNAseqCount_PosStrand")])),
  TRUE, 
  FALSE
)

#Set order of columns as follow:
Column_order <- c("Enh", "EnhID", "Gene", "Pair", "EnsID", "Distance", "TSS", "Gene.Nearest", 
                  "Gene.RNAseq_Counts", "Gene.RNAseq_CPM", "Gene.RNAseq_RPKM", "Gene.Housekeeping", 
                  "Gene.StabilityIndex", "Enh.chr", "Enh.start", "Enh.end", "Enh.Pos", 
                  "ATACseq.Pileup", "Enh.Midpoint", "Chip.H3K27ac", "Chip.H3K4me3", 
                  "Enh.TTseqCount_NegStrand", "Enh.TTseqCount_PosStrand", 
                  "Enh.RNAseqCount_NegStrand", "Enh.RNAseqCount_PosStrand", "Enh.NascentEnriched", 
                  "Enh.eRNA", "numTSSEnhGene", "rf_pred", "Gene.Tested", "pass_rf")

Supplementary_Table_9=Supplementary_Table_9[Column_order]

#Define pair column asenh pos + gene
Supplementary_Table_9$Pair=paste0(Supplementary_Table_9$Enh.Pos,"_", Supplementary_Table_9$Gene)

##Subset EGPs by those predicted by EGRF and  > 0.5PRKM (this give Stable 9E)
Supplementary_Table_9E=Supplementary_Table_9[Supplementary_Table_9$pass_rf==TRUE & Supplementary_Table_9$Gene.RNAseq_RPKM>0.5  & !All_EGPs$Enh.Pos %in% res.final$Enh.Pos,]
dim(Supplementary_Table_9E)

#Export Tables
write.csv(Supplementary_Table_9E, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable9_RandomForests/9E_RandomForestPredictedHits.csv")

#Subset By RNA-seq expression > 0.5 RPKM and enhancers not tested in our CRISPRi screening. This generates Stable9G (EGrf EGP all predictions)
Supplementary_Table_AllPreds=Supplementary_Table_9[Supplementary_Table_9$Gene.RNAseq_RPKM>0.5 & !All_EGPs$Enh.Pos %in% res.final$Enh.Pos,]

write.csv(Supplementary_Table_AllPreds, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable9_RandomForests/9_AllRandomForestPredictedHits.csv")
