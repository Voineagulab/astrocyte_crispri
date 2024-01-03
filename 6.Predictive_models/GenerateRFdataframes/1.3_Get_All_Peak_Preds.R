################
##
##
##Predicting on all peaks
##
##
###############
library(tictoc)
library(ggpointdensity)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/") #/Volumes/share/
source("Scripts/Header_functions.R")
source("Scripts/Annotation_Functions.R")
library(cowplot)
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
all_rf_factors <- read.csv("Data/all_rf_factors.csv")
allpeaks <- read.csv("../FullLibrary_Selection/Results/Peaks_Annotated.csv")
geneinfo <- read.table("../FullScale/Data/Whitelists/GeneInfo.txt")
geneinfo$TSS <- sub(".*:","",geneinfo$TSS)
colnames(geneinfo)[colnames(geneinfo) == "Symbol"] <- "Gene"

#######
##
##Annotate peaks with nearby gene expression data
##
#######
allpeaks$Enh.size <- allpeaks$end - allpeaks$start
allpeaks$Enh.Midpoint <- allpeaks$start + round(allpeaks$Enh.size / 2)
colnames(allpeaks)[1:3] <- c("Enh.chr", "Enh.start", "Enh.end")
colnames(allpeaks)[colnames(allpeaks) == "name"] <- "Enh"
colnames(allpeaks)[colnames(allpeaks) == "id"] <- "Enh.Pos"
allpeaks <- allpeaks[allpeaks$Enh.chr %in% unique(geneinfo$Chr),]
allpeaks$Enh.chr <- sub("chr", "", allpeaks$Enh.chr)
#############
#Save Set/location/bed files
intergenic_peaks <- allpeaks[allpeaks$GENCODE32 == F,]
write.table(intergenic_peaks[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")], 
            file = "Data/PeakData/Intergenic_Peaks.bed",
            row.names = F, col.names = F, quote = F, sep = "\t")
intergenic_loc <- intergenic_peaks[,c("Enh", "Enh.chr", "Enh.start", "Enh.end")]
colnames(intergenic_loc) <- c("GENE", "CHR", "START", "END")
intergenic_loc$CHR <- paste0("chr",  intergenic_loc$CHR)
write.table(intergenic_loc, "Data/PeakData/Intergenic_Peaks.loc", row.names = F, 
            col.names = T, quote = F, sep = "\t")
write.table(intergenic_loc[,c( "CHR", "START", "END","GENE")], 
            "Data/PeakData/Intergenic_Peaks_chr.bed", row.names = F, 
            col.names = F, quote = F, sep = "\t")
write.table(intergenic_peaks$Enh, "Data/PeakData/Intergenic_Peaks.set"
            , row.names = F, col.names = F, quote = F, sep = "\t")


selected_vars_exp <- c("Chip.H3K4me3", "Chip.H3K27ac","Gene.Nearest", 
                       "Distance","ATACseq.Pileup", "numTSSEnhGene", "Gene.StabilityIndex")
##########
#Complete peak annotations
intergenic_peaks$tested_enh <- intergenic_peaks$Enh.Pos %in% all_rf_factors$Enh.Pos
intergenic_peaks$HitPermissive <- intergenic_peaks$Enh.Pos %in% all_rf_factors[all_rf_factors$HitPermissive > 0, "Enh.Pos"]
intergenic_peaks$HitPermissive_NegZ <- intergenic_peaks$Enh.Pos %in% all_rf_factors[all_rf_factors$HitPermissive_NegZ > 0, "Enh.Pos"]
#Add Chip
folder <- "Results/Process_Chip/"
filenames <- list.files("Results/Process_Chip/", pattern = ".tab")
filenames <- filenames[str_detect(filenames, pattern = "Intergenic")]
intergenic_peaks <- addChipBigWigs(intergenic_peaks, filenames,"Results/Process_Chip/")

#Tobias Bound 
intergenic_peaks <- addTobiasBound(intergenic_peaks, merge_col = "Enh",bound_file = "Results/Tobias/Footprint/BINDetect/bound_overlaps_intergenic.bed",
               unbound_file = "Results/Tobias/Footprint/BINDetect/unbound_overlaps_intergenic.bed")
#
#Add TTseq/RNAseq info
ttseq <- read.csv("../FullScale/Results/4_EnhancerTranscription/TTseq/Results_Table_intergenic.csv")
colnames(ttseq) <- paste0("TTseq.", colnames(ttseq))
colnames(ttseq)[1] <- "Enh"
intergenic_peaks <- merge(intergenic_peaks, ttseq[,c("Enh", "TTseq.TTseq_Total", "TTseq.TT_Enrich", "TTseq.RNAseq_Total", "TTseq.Ratio_TTversusRNA")])

#update column names to match RF
colnames(intergenic_peaks)[colnames(intergenic_peaks) == "fold_enrichment"] <- "Peaks.fold_enrichment"
colnames(intergenic_peaks)[colnames(intergenic_peaks) == "pileup"] <- "ATACseq.Pileup"

# bedtools intersect -wb -a /mnt/Data0/PROJECTS/CROPSeq/PublicData/Yao_NatBiotech2022/ST2_KnownEnhancers.bed \
# -b /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Intergenic_Peaks_chr.bed > Intergenic_KnownEnhancers.bed
knownenhancers <- read.table("Results/Get_All_Peak_Preds/Intergenic_KnownEnhancers.bed")
intergenic_peaks$knownEnhancer <- intergenic_peaks$Enh %in% knownenhancers$V9

write.csv(intergenic_peaks, "Results/Get_All_Peak_Preds/Annotated_intergenic.csv", row.names = F)
#intergenic_peaks <- read.csv("Results/Get_All_Peak_Preds/Annotated_intergenic.csv")

###############
##
##Create EGPs
##
##
###############

#Create candidate EGPs for all genes in 500kb
EGPs <- data.frame(EnsID = as.character(), Gene = as.character(),Gene.Distance = as.numeric(), Enh = as.character())
for (peak in unique(intergenic_peaks$Enh)) {
  midpoint <- intergenic_peaks[intergenic_peaks$Enh == peak,"Enh.Midpoint"]
  chr <- paste0("chr",intergenic_peaks[intergenic_peaks$Enh == peak,"Enh.chr"])
  chr_genes <- geneinfo[geneinfo$Chr == chr,]
  chr_genes$Gene.Distance <- abs(as.numeric(chr_genes$TSS) - midpoint) 
  genes <- chr_genes[chr_genes$Gene.Distance < 500000,c("EnsID","Gene", "Gene.Distance", "TSS")] 
  if (nrow(genes) > 0) {
    EGPs <- rbind(EGPs, cbind(genes,peak))
  }
}
#Gene nearest
EGPs$Gene.Nearest <- FALSE
for (enh in unique(EGPs$Enh)) {
  EGPs[EGPs$Enh == enh & EGPs$Gene.Distance == min(EGPs[EGPs$Enh == enh,"Gene.Distance"]),"Gene.Nearest"] <- TRUE
}

#Add bulk Expression
colnames(EGPs)[colnames(EGPs) == "peak"] <- "Enh"
feature_counts <- read.csv("Data/GeneData/RNAseq_FeatCounts.csv")
colnames(feature_counts) <- paste0("Gene.", colnames(feature_counts))
colnames(feature_counts)[1] <- "Gene"
EGPs <- merge(EGPs,feature_counts, all.x = T)

#EGPs <- read.csv("tmp.csv")
#Add gene stability 
warning("Check this")
EGPs <- addGeneStability(EGPs, geneinfo = geneinfo, ensembl)

#Add peak features (see above annotations)
intergenic_predvals <- intergenic_peaks[,colnames(intergenic_peaks) %in% 
                                          c(selected_vars_exp, "Enh", "Enh.Pos","Enh.Midpoint", "Enh.chr", "Enh.start", "Enh.end")]
EGPs$Pair <- paste0(EGPs$Enh, "_", EGPs$Gene)
EGPs <- merge(EGPs,intergenic_predvals, all.x = T)



#Add Hits / Tested annotation
EGPs$Enh.Pos_Gene <- paste0(EGPs$Enh.Pos, "_", EGPs$Gene)
EGPs$EGP.Tested <- EGPs$Enh.Pos_Gene %in% all_rf_factors$Enh.Pos_Gene
EGPs$EGP.HitPermissive <- EGPs$Enh.Pos_Gene %in% all_rf_factors[all_rf_factors$HitPermissive,]$Enh.Pos_Gene

# 
# #Add Tads  
# tads <- read.delim("../FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38.bed", header = FALSE)
# colnames(tads) <- c("chr", "start", "end")
# tads$chr<- sub("chr","", tads$chr)
# same_tad <- sameTADs(EGPs, tads)
# EGPs <- merge(EGPs,same_tad, all.x = T)
# EGPs[is.na(EGPs$EGP.Shared_TADs),"EGP.Shared_TADs"] <- 0
# EGPs$TAD <- EGPs$EGP.Shared_TADs > 0 

EGPs$numTSSEnhGene <- addnumTSSbetween(EGPs, geneinfo)
colnames(EGPs)[colnames(EGPs) == "Gene.Distance"] <- "Distance"
selected_vars_exp[!selected_vars_exp %in% colnames(EGPs)]
selected_vars[!selected_vars %in% colnames(EGPs)]
#Get Predictions from RF

model_EGrf <- readRDS("Results/ComparisonPlots/SavedModels/EGPModelEGrf_Power_subset_HitPermissive_NegZ_cv.rds") #TODO
apply(EGPs, 2, function(x){ sum(is.na(x))})

EGPs[is.na(EGPs$Gene.StabilityIndex),"Gene.StabilityIndex"] <- 0
preds <- predict(model_EGrf, EGPs) 
EGPs$rf_pred <- preds$predictions[,2]
EGPs$Gene.Tested <- EGPs$Gene %in% all_rf_factors$Gene


EGrf80precisionCutoff <- read.table("Results/ComparisonPlots/80ptPrecisionCutoff.txt")$x
EGPs$pass_rf <- EGPs$rf_pred > EGrf80precisionCutoff & ! EGPs$EGP.Tested

write.table(EGPs[EGPs$pass_rf,c("Pair","Gene.RNAseq_RPKM")], "Results/Get_All_Peak_Preds/EGPs_pred_80Prec_fc.txt", row.names = F, quote = F, col.names = F)
write.csv(EGPs, "Results/Get_All_Peak_Preds/All_EGPs.csv", row.names = F)
#EGPs <- read.csv("Results/Get_All_Peak_Preds/All_EGPs.csv")



