#################################################
##This script applies EGrf to predict EGPs for intergenic ATAC-seq peaks not included in the CRISPRi screen
##############################################

#Load data
library(tictoc)
library(ggpointdensity)
library(cowplot)

# Set the working directory to the specified path.
path <-"/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/"
setwd(paste0(path, "/1.Data/InputData"))

# Source additional R scripts containing custom functions.
source(paste0(path, "/2.Scripts/Functions.R"))

# Connect to the Ensembl database for gene information.
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

# Load CRISPRi screen data
all_rf_factors <- read.csv("Astrocytes/CRISPRi_Results Final.csv")
# Load  ATAC-seq data including all peaks called 
allpeaks <- read.csv("Astrocytes/NHA_All_ATACseq_Peaks_Annotated.csv")
# Read gene annotation information and rename coloumns for consistency.
geneinfo <- read.table("GeneInfo.txt")
geneinfo$TSS <- sub(".*:","",geneinfo$TSS)
colnames(geneinfo)[colnames(geneinfo) == "Symbol"] <- "Gene"

# Extract position components for intergenic peaks and rename columns for consistency.
allpeaks$Enh.size <- allpeaks$end - allpeaks$start
allpeaks$Enh.Midpoint <- allpeaks$start + round(allpeaks$Enh.size / 2)
colnames(allpeaks)[1:3] <- c("Enh.chr", "Enh.start", "Enh.end")
colnames(allpeaks)[colnames(allpeaks) == "name"] <- "Enh"
colnames(allpeaks)[colnames(allpeaks) == "id"] <- "Enh.Pos"
allpeaks <- allpeaks[allpeaks$Enh.chr %in% unique(geneinfo$Chr),]
allpeaks$Enh.chr <- sub("chr", "", allpeaks$Enh.chr)

#Save Set/location/bed files
intergenic_peaks <- allpeaks[allpeaks$GENCODE32 == F,]
# write.table(intergenic_peaks[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")], 
#              file = "Data/PeakData/Intergenic_Peaks.bed",
#              row.names = F, col.names = F, quote = F, sep = "\t")
intergenic_loc <- intergenic_peaks[,c("Enh", "Enh.chr", "Enh.start", "Enh.end")]
colnames(intergenic_loc) <- c("GENE", "CHR", "START", "END")
intergenic_loc$CHR <- paste0("chr",  intergenic_loc$CHR)
# write.table(intergenic_loc, "Data/PeakData/Intergenic_Peaks.loc", row.names = F, 
#              col.names = T, quote = F, sep = "\t")
# write.table(intergenic_loc[,c( "CHR", "START", "END","GENE")], 
#              "Data/PeakData/Intergenic_Peaks_chr.bed", row.names = F, 
#              col.names = F, quote = F, sep = "\t")
# write.table(intergenic_peaks$Enh, "Data/PeakData/Intergenic_Peaks.set"
#             , row.names = F, col.names = F, quote = F, sep = "\t")
#Complete peak annotations
intergenic_peaks$tested_enh <- intergenic_peaks$Enh.Pos %in% all_rf_factors$Enh.Pos
intergenic_peaks$HitPermissive <- intergenic_peaks$Enh.Pos %in% all_rf_factors[all_rf_factors$HitPermissive > 0, "Enh.Pos"]
intergenic_peaks$HitPermissive_NegZ <- intergenic_peaks$Enh.Pos %in% all_rf_factors[all_rf_factors$HitPermissive_NegZ > 0, "Enh.Pos"]

#Define variables for the model
selected_vars_exp <- c("Chip.H3K4me3", "Chip.H3K27ac","Gene.Nearest", 
                       "Distance","ATACseq.Pileup", "numTSSEnhGene", "Gene.StabilityIndex")

# Add Astrocyte ChIP-seqdata from ENCODE NHA samples.
folder <- "Astrocytes/Processed_ChIP_ENCODE/Intergenic/" # remove un-necessary files
filenames <- list.files("Astrocytes/Processed_ChIP_ENCODE/Intergenic/", pattern = ".tab")
filenames <- filenames[str_detect(filenames, pattern = "Intergenic")]
intergenic_peaks <- addChipBigWigs(intergenic_peaks, filenames,"Astrocytes/Processed_ChIP_ENCODE/Intergenic/")

# #Add TF Footprinting data generated with Tobias
# intergenic_peaks <- addTobiasBound(intergenic_peaks, merge_col = "Enh",bound_file = "Astrocytes/bound_overlaps_intergenic.bed",
#                unbound_file = "Astrocytes/unbound_overlaps_intergenic.bed")

# # Add Enhancer TTseq data, i.e eRNA expression
# ttseq <- read.csv("Astrocytes/Results_Table_intergenic.csv")
# colnames(ttseq) <- paste0("TTseq.", colnames(ttseq))
# colnames(ttseq)[1] <- "Enh"
# intergenic_peaks <- merge(intergenic_peaks, ttseq[,c("Enh", "TTseq.TTseq_Total", "TTseq.TT_Enrich", "TTseq.RNAseq_Total", "TTseq.Ratio_TTversusRNA")])

#Update column names to match RF data
colnames(intergenic_peaks)[colnames(intergenic_peaks) == "fold_enrichment"] <- "Peaks.fold_enrichment"
colnames(intergenic_peaks)[colnames(intergenic_peaks) == "pileup"] <- "ATACseq.Pileup"

# Save intergenic peaks
write.csv(intergenic_peaks, "Results/Get_All_Peak_Preds/Annotated_intergenic.csv", row.names = F)

##Create EGPs for all genes in 500kb
EGPs <- data.frame(EnsID = as.character(), Gene = as.character(),Gene.Distance = as.numeric(), Enh = as.character())
# Loop through each unique enhancer peak in the 'intergenic_peaks' data frame.
for (peak in unique(intergenic_peaks$Enh)) {
  # Extract the midpoint of the current enhancer peak.
  midpoint <- intergenic_peaks[intergenic_peaks$Enh == peak,"Enh.Midpoint"]
  # Construct the chromosome identifier by concatenating 'chr' with the chromosome number of the current enhancer peak.
  chr <- paste0("chr",intergenic_peaks[intergenic_peaks$Enh == peak,"Enh.chr"])
  # Filter 'geneinfo' data frame to get genes on the same chromosome as the current enhancer peak.
  chr_genes <- geneinfo[geneinfo$Chr == chr,]
  # Calculate the absolute distance of each gene's TSS from the enhancer peak's midpoint.
  chr_genes$Gene.Distance <- abs(as.numeric(chr_genes$TSS) - midpoint) 
  # Select genes that are within 500kb of the enhancer peak.
  genes <- chr_genes[chr_genes$Gene.Distance < 500000,c("EnsID","Gene", "Gene.Distance", "TSS")] 
  # If there are any genes close enough to the enhancer peak, add them to the 'EGPs' data frame.
  if (nrow(genes) > 0) {
    EGPs <- rbind(EGPs, cbind(genes,peak))
  }
}

#Loop over each unique enhancer and calculate which gene is the nearest
EGPs$Gene.Nearest <- FALSE
for (enh in unique(EGPs$peak)) {
  EGPs[EGPs$Enh == enh & EGPs$Gene.Distance == min(EGPs[EGPs$Enh == enh,"Gene.Distance"]),"Gene.Nearest"] <- TRUE
}

#Add bulk Expression
colnames(EGPs)[colnames(EGPs) == "peak"] <- "Enh"
feature_counts <- read.csv("Astrocytes/RNAseq_FeatCounts.csv")
colnames(feature_counts) <- paste0("Gene.", colnames(feature_counts))
colnames(feature_counts)[1] <- "Gene"
EGPs <- merge(EGPs,feature_counts, all.x = T)

#Add gene stability index
EGPs <- addGeneStability(EGPs, geneinfo = geneinfo, ensembl)

#Combine feature data with EG pairs dataframe
intergenic_predvals <- intergenic_peaks[,colnames(intergenic_peaks) %in% 
                                          c(selected_vars_exp, "Enh", "Enh.Pos","Enh.Midpoint", "Enh.chr", "Enh.start", "Enh.end")]
EGPs$Pair <- paste0(EGPs$Enh, "_", EGPs$Gene)
EGPs <- merge(EGPs,intergenic_predvals, all.x = T)

#Add Hits / Tested annotation
EGPs$Enh.Pos_Gene <- paste0(EGPs$Enh.Pos, "_", EGPs$Gene)
EGPs$EGP.Tested <- EGPs$Enh.Pos_Gene %in% all_rf_factors$Enh.Pos_Gene
EGPs$EGP.HitPermissive <- EGPs$Enh.Pos_Gene %in% all_rf_factors[all_rf_factors$HitPermissive,]$Enh.Pos_Gene

##Add number of Transcription Start Sites (TSS) between Enhancers and genes, and filtering out va
EGPs$Enh.chr=paste0("chr", sep="",EGPs$Enh.chr)
EGPs$numTSSEnhGene <- addnumTSSbetween(EGPs, geneinfo)
colnames(EGPs)[colnames(EGPs) == "Gene.Distance"] <- "Distance"
selected_vars_exp[!selected_vars_exp %in% colnames(EGPs)]

#Get Predictions from RF model
model_EGrf <- readRDS("../../3.Predictions/RF_Results/Astrocytes/EGPModelEGrf_Power_subset_HitPermissive_NegZ_cv.rds") 

#Checking and Handling missing values
apply(EGPs, 2, function(x){ sum(is.na(x))})
EGPs[is.na(EGPs$Gene.StabilityIndex),"Gene.StabilityIndex"] <- 0

#Predicting enh-gene interactions by using ranger package
preds <- predict(model_EGrf, EGPs) 
EGPs$rf_pred <- preds$predictions[,2]
EGPs$Gene.Tested <- EGPs$Gene %in% all_rf_factors$Gene

##Filtering outputs by a 80% precision threshold  
EGrf80precisionCutoff <- read.table("../../3.Predictions/RF_Results/Astrocytes/EGrf80ptPrecisionCutoff.txt")$x
EGPs$pass_rf <- EGPs$rf_pred > EGrf80precisionCutoff & ! EGPs$EGP.Tested

# Supp table - filter for  expressed genes with a permissive threshold of RPKM >0.1; 
# remove TT-seq data which is not used and thus confusing.
# Also remove EGP.HitPermissive and EGP.Tested which are false for all entries, because these are EGPs not tested in the screen.

supp_table<-EGPs[((EGPs$pass_rf)&(EGPs$Gene.RNAseq_RPKM > 0.1)), -grep("TTseq", colnames(EGPs)) ]
supp_table<-supp_table[ , -grep("EGP", colnames(supp_table))] 

#Save data
setwd("../../3.Predictions/RF_Results/")
write.csv(EGPs, "Astrocytes_All_Intergenic_Predictions/All_EGPs.csv", row.names = F)
write.table(EGPs[EGPs$pass_rf,c("Pair","Gene.RNAseq_RPKM")], "Astrocytes_All_Intergenic_Predictions/EGPs_pred_80Prec_fc.txt", row.names = F, quote = F, col.names = F)
write.csv(supp_table, "Astrocytes_All_Intergenic_Predictions/STable_AllIntergenicPreds.csv", row.names = F)
