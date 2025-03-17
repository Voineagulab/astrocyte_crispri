###This script generates the training variables for random forest prediction models applied to Astrocyte data 

#04-01-23
#Sam Bagot
#Modified by IV and JP: 1.1_CreateEGPdataframe_JP_IV
rm(list=ls())

# Set the working directory to the specified path.
path <-"/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/"
setwd(paste0(path, "/1.Data/InputData"))
# Load necessary libraries.
library("stringr")
library("seqinr")
library("tidyr")
library("PRROC")
library("biomaRt")
library("ggplot2")
library("readxl")
library("rapportools")
# Connect to the Ensembl database for gene information.
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
# Source additional R scripts containing custom functions.
source(paste0(path, "/2.Scripts/Functions.R"))

# Read in the list of variables used for EGrf and TAP-seq
final.vars<-read.csv("Astrocytes/AstroVars.csv")
############################## Load and format CRISPRi screening data  #######################################################
##############################################################################################################################
# Load CRISPRi screen data
all_rf_factors <- read.csv("Astrocytes/CRISPRi_Results Final.csv")

# Create a logical column indicating negative log fold change and hit permissive, i.e. positive EG pairs
all_rf_factors$HitPermissive_NegZ <- all_rf_factors$HitPermissive & all_rf_factors$logfc.vst < 0

##Loading power calculation data
power <- read.csv("Astrocytes/CRISPRi_Power Simulation - Power per EGP.csv")

##Retain positive EG pairs, and negative EG pairs with >= 80% power to detect a 15% downregulation in gene expression
all_rf_factors<- all_rf_factors[all_rf_factors$Pair %in% power[power$WellPowered015,]$Pair |
                                           all_rf_factors$HitPermissive_NegZ,]

# Create a new column by combining Enhancer Position and Gene.
all_rf_factors$Enh.Pos_Gene <- paste0(all_rf_factors$Enh.Pos,"_",all_rf_factors$Gene)
#write.table(unique(all_rf_factors$Gene), file = "Data/GeneData/Genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Read gene annotation information and merge it with the main data frame.
geneinfo <- read.table("GeneInfo.txt")
# Remove duplicate ENSIDs to avoid redundancy.
geneinfo <- geneinfo[match(unique(geneinfo$Symbol), geneinfo$Symbol),]
# Merge gene information with all_rf_factors based on Gene and Symbol columns.
all_rf_factors <- merge(all_rf_factors, geneinfo[,c("EnsID", "Symbol")], by.x = "Gene", by.y = "Symbol")

# Extract position components for enhancer and gene positions.
all_rf_factors <- extract_pos(all_rf_factors,"Enh.Pos", out_name = "Enh.")
all_rf_factors <- extract_pos(all_rf_factors,"Gene.Pos", out_name = "Gene.")

# Extract unique enhancer positions for candidates and hit enhancers
unique_enh <- unique(all_rf_factors[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")])
hit_enhs <- unique(all_rf_factors[all_rf_factors$HitPermissive,c("Enh.chr", "Enh.start", "Enh.end", "Enh")])
nonhit_enhs <- unique(all_rf_factors[!all_rf_factors$HitPermissive,c("Enh.chr", "Enh.start", "Enh.end", "Enh")])

# Write bed files for enhancer data.
#write.bed(unique_enh, file = "Data/PeakData/Tested_Enh_nochr.bed")
#write.bed(hit_enhs, "Data/PeakData/Hit_Enhancers.bed")
#write.bed(nonhit_enhs, "Data/PeakData/NonHit_Enhancers.bed")

# Add "chr" prefix to chromosome names and write to a BED file.
unique_enh$Enh.chr <- paste0("chr", unique_enh$Enh.chr)
#write.table(unique_enh, file = "Data/PeakData/Tested_Enh.bed", row.names = FALSE, quote = FALSE, col.names = FALSE , sep = "\t")

# Generate a bed file specifying the Enhancer Midpoint and Gene TSS. This will be used to obtain contact data from Hi-C.
all_rf_factors$Enh.size <- all_rf_factors$Enh.end - all_rf_factors$Enh.start

all_rf_factors$Gene.TSS <- sub(".*:","",all_rf_factors$Gene.TSS)
all_rf_factors$Enh.Midpoint <- all_rf_factors$Enh.end - round(all_rf_factors$Enh.size / 2)
#write.table(all_rf_factors[,c("Enh.chr", "Enh.Midpoint", "Gene.TSS", "Pair")], file = "Data/GeneData/Astrocytes_HG38_Enh_Midpoint_Gene_TSS.txt", row.names = F, quote = F, sep = "\t")

# Define gene promoters as TSS +/- 500bp
gene_promoters <- unique(all_rf_factors[,c('Gene', "Gene.chr", "Gene.TSS", "Gene.TSS")])
gene_promoters$Gene.TSS <- as.numeric(gene_promoters$Gene.TSS) - 500
gene_promoters$Gene.TSS.1 <- as.numeric(gene_promoters$Gene.TSS.1) + 500
# Write the promoter data to a BED file for further analysis.
#write.bed(gene_promoters[,c(2,3,4,1)], file = "Data/GeneData/Gene_Promoters.bed")

############################## Add training variable data to the CRISPRi screen EG pairs  #######################################################
##############################################################################################################################

##########################
# Add Enhancer TTseq data, i.e eRNA expression
##########################
ttseq <- read.csv("Astrocytes/TTseq_Results Table.csv")
# Rename columns for consistency and clarity.
colnames(ttseq) <- paste0("TTseq.", colnames(ttseq))
colnames(ttseq)[1] <- "Enh"
# Merge TTseq data with the main dataset.
all_rf_factors <- merge(all_rf_factors, ttseq[,c("Enh", "TTseq.TTseq_Total", "TTseq.TPM_TTseq_Total","TTseq.TT_Enrich", "TTseq.RNAseq_Total", "TTseq.Ratio_TTversusRNA")])
# Compute the difference between TTseq and RNAseq total counts.
all_rf_factors$TTseq.Diff_TTRNA <- all_rf_factors$TTseq.TTseq_Total - all_rf_factors$TTseq.RNAseq_Total


##########################
# Add ATAC seq data
##########################
# Load ATAC-seq data
nha_peaks <- read.csv("Astrocytes/NHA_Peaks_Annotated.csv") 

# Select specific columns from nha_peaks for further analysis.
nha_peaks <- nha_peaks[, c("pileup", "id")]

# Clean up variable names  for better readability.
colnames(nha_peaks) <- c("ATACseq.Pileup", "Enh.Pos")

# Merge the ATAC seq data (nha_peaks and narrowPeaks) with the all_rf_factors dataframe.
all_rf_factors <- merge(all_rf_factors, nha_peaks)

##############
# Add TF Footprinting data generated with Tobias
##############
# Load Tobias TF count data
Tobias_TF <- read.csv("Astrocytes/TOBIAS_TF_counts.csv")
# Select relevant columns for TF counts.
Tobias_TF <- Tobias_TF[, c("Enh", "Bound_TF_counts", "Unbound_TF_counts", "Exp_Bound_TF_counts", "Exp_Unbound_TF_counts")]
# Rename columns in Tobias_TF for clarity
colnames(Tobias_TF)[colnames(Tobias_TF) != "Enh"] <- paste0("Tobias.", colnames(Tobias_TF)[colnames(Tobias_TF) != "Enh"])
# Merge Tobias data with the all_rf_factors dataframe.
all_rf_factors <- merge(all_rf_factors, Tobias_TF)
# Calculate the ratio of bound TF counts to the total of bound and unbound TF counts.
all_rf_factors$Tobias.Bound_ratio <- all_rf_factors$Tobias.Bound_TF_counts / (all_rf_factors$Tobias.Bound_TF_counts + all_rf_factors$Tobias.Unbound_TF_counts)

# Replace any NA (not available) values with 0 in the all_rf_factors dataset.
all_rf_factors[is.na(all_rf_factors)] <- 0

########################
# Add Astrocyte ChIP-seq and open chromaatin state data from ENCODE NHA samples.
#######################
# Define the folder containing chip data.
folder <- "Astrocytes/Processed_ChIP_ENCODE/"
filenames <- paste0(final.vars$Accession[!is.na(final.vars$Accession)], ".tab")
#add ChIP-seq data
all_rf_factors <- addChipBigWigs(all_rf_factors, filenames, folder)

########################
# Add Hi-C data
########################
# Load Hi-C interaction data 
hicContacts <- read.csv("Astrocytes/Astrocyte_Enh_HiC.csv")
hicContacts <- hicContacts[,4:8]
# Rename columns for clarity.
colnames(hicContacts) <- c("Pair", "HiC_interaction_astrocyte_cerebellum", 
                           "HiC_interaction_astrocyte_spinal", "HiC_interaction_NeuN_pos", "HiC_interaction_NeuN_neg")
# Merge Hi-C data with the all_rf_factors dataframe 
all_rf_factors <- merge(all_rf_factors, hicContacts, by = "Pair", all.x = T)

##################
## Add Variant Disease Scores obtained with the Webtool implementation of the Beluga Deeplearning Model
##################
# Load Beluga disease score predictions
beluga_DIS_res <- read.csv("Astrocytes/MaxBelugaDiseaseScores.csv")
# Merge the Beluga predictions with the all_rf_factors dataframe.
all_rf_factors <- merge(all_rf_factors, beluga_DIS_res, all.x = T)
# Replace NA values with 0.
all_rf_factors[is.na(all_rf_factors$Beluga.MaxDisScore), "Beluga.MaxDisScore"] <- 0

#################
## Add Housekeeping information and gene stability index
################
colnames(geneinfo)[colnames(geneinfo) == "Symbol"] <- "Gene"
all_rf_factors <- addGeneStability(all_rf_factors, geneinfo = geneinfo, ensembl)

############################
## Add ENCODE rE2G training variable data as well as ENCODE rE2G & ABC predictions for NHA Astrocytes from ENCFF440FMQ.bed
############################
# Define column names for the ENCODE and ABC prediction data.
ENCODE_colnames <- c("Encode.Enh.chr", "Encode.Enh.start", "Encode.Enh.end", "Encode.Enh",
                     "class", "Gene", "EnsID", "TargetGeneTSS",
                     "isSelfPromoter", "CellType", 'numTSSEnhGene.Feature',
                     "distanceToTSS.Feature", "normalizedDNase_enh.Feature",
                     "normalizedDNase_prom.Feature", "numNearbyEnhancers.Feature",
                     "sumNearbyEnhancers.Feature", "ubiquitousExpressedGene.Feature",
                     "numCandidateEnhGene.Feature", "3DContactAvgHicTrack2.Feature",
                     "3DContactAvgHicTrack2_squared.Feature", "activityEnhDNaseOnlyAvgHicTrack2_squared.Feature",
                     "activityPromDNaseOnlyAvgHicTrack2.Feature", "ABCScore_ENCODE",
                     "DNAseOnly_Score")

# Read in a bed file containing the overlap between the CRISPRi screen enhancers and ENCODE rE2G predictions for astrocytes
# This bed file is generated by the following script: OverlapENCODE_rE2G.sh
# The first 4 columns correspond to the astrocyte CRISPRi screen peaks, followed by the overlapping ENCODE data
pred.file <- "Astrocytes/ENCFF440FMQ.bed"
# Process the bed file to extract and format the necessary data.
# "DNAseOnly_Score" refers to the ENCODE rE2G DNase only model predictions.

  res <- read.table(pred.file)
  colnames(res) <- c("Enh.chr", "Enh.start", "Enh.end", "Enh", ENCODE_colnames)
  res <- res[,c("Enh", "EnsID", "CellType", "ABCScore_ENCODE", "DNAseOnly_Score")]
  # Append the cell type to the score column names for differentiation.
  ct <- unique(res$CellType)
  colnames(res)[4:5] <- paste0(colnames(res)[4:5], ".", ct)
  # Filter out rows where EnsID is missing and remove the 'CellType' column.
  res <- res[!is.na(res$EnsID), !colnames(res) %in% c("CellType")]
  preds <- res

# Aggregate the prediction data for non-unique  Enhancer Gene pairs, calculating the mean of scores.
# This would occur if multiple ENCODE NHA enhancers overlapped the same peak in our ATAC-seq data.
astro_preds <- aggregate(.~Enh+EnsID,preds, mean, na.rm = T,na.action = NULL)

all_rf_factors <- merge(all_rf_factors, astro_preds, all.x = T)
#Keep the variable names for NHA astrocytes short for plotting purposes
colnames(all_rf_factors)[colnames(all_rf_factors) %in% c("ABCScore_ENCODE.astrocyte", "DNAseOnly_Score.astrocyte")] <- c("ABCScore_ENCODE", "rE2G.DNAseOnly")

############################## Final formatting and save  #######################################################
#################################################################################################################

# Modify chromosome notation in geneinfo for consistency.
geneinfo$Chr <- sub("chr", "", geneinfo$Chr)
# Include TSS information in the all_rf_factors dataframe.
all_rf_factors$TSS <- all_rf_factors$Gene.TSS 
# Calculate and add the number of TSS between each gene and its associated enhancer.
all_rf_factors$numTSSEnhGene <- addnumTSSbetween(all_rf_factors, geneinfo)

# Add a column indicating the presence of overlapping enhancer-gene pairs in the ENCODE data.
all_rf_factors$ENCODE_overlap.Exists <- !is.na(all_rf_factors$rE2G.DNAseOnly)
# Order the dataframe by the 'Pair' column
all_rf_factors <- all_rf_factors[order(all_rf_factors$Pair),]

#Considering these as missed predictions but when evaluating performance will use intersecting EG pairs only
all_rf_factors[!all_rf_factors$ENCODE_overlap.Exists, c("rE2G.DNAseOnly", "ABCScore_ENCODE")] <- 0

all_rf_factors[is.na(all_rf_factors$Gene.StabilityIndex), "Gene.StabilityIndex"] <- 0 #Assume if not in results 0 stability

colnames(all_rf_factors)[colnames(all_rf_factors) == "Gene.Distance"] <- "Distance"

# Save data
write.csv(all_rf_factors, "../TrainingData/TrainingDataframe_Astrocytes.csv", row.names = FALSE) 
