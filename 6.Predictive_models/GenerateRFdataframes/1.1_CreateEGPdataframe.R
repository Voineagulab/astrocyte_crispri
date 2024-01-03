#Inital data exploration
#04-01-23
#Sam Bagot
#enhancer prediciton project
#HG38?? 
#packages
#library(kmer)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/")
library("stringr")
library(seqinr)
library("tidyr")
source("Scripts/Header_functions.R")
source("Scripts/Annotation_Functions.R")
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

#Load
all_rf_factors <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")  #temp_data/Results_Final_COPY.csv
all_rf_factors$Enh.Pos_Gene <- paste0(all_rf_factors$Enh.Pos,"_",all_rf_factors$Gene)
all_rf_factors$HitPermissive_NegZ <- all_rf_factors$HitPermissive & all_rf_factors$logfc.vst < 0
write.table(unique(all_rf_factors$Gene), file = "Data/GeneData/Genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)



#Add EnsID
geneinfo <- read.table("../FullScale/Data/Whitelists/GeneInfo.txt")
geneinfo <- geneinfo[match(unique(geneinfo$Symbol), geneinfo$Symbol),]#Remove duplicate ENSIDs (this may be a problem later but don't want duplicates)
all_rf_factors <- merge(all_rf_factors, geneinfo[,c("EnsID", "Symbol")], by.x = "Gene", by.y = "Symbol")


#Add Fasta sequence originally this was to make kmers but this has been removed (In github history)
NHA_fasta <- read.table("Data/PeakData/PeaksFasta.tab")
colnames(NHA_fasta) <- c("Enh", "Fasta")
NHA_fasta$Fasta <- toupper(NHA_fasta$Fasta)

all_rf_factors <- merge(all_rf_factors,NHA_fasta)
hit_enhancers <- unique(all_rf_factors[which(all_rf_factors$HitCore == TRUE),]$Enh)
##current assumptions -> 

extract_pos <- function(data, pos_col, out_name = "") {
  data[,paste0(out_name, "chr")] <-  sub("chr(.*):.*","\\1",data[,pos_col])
  data[,paste0(out_name, "start")] <- as.numeric(sub(".*:(.*)-.*","\\1",data[,pos_col]))
  data[,paste0(out_name, "end")] <- as.numeric(sub(".*:.*-(.*)","\\1",data[,pos_col]))
  data[,paste0(out_name, "size")] <- data[,paste0(out_name, "end")] - data[,paste0(out_name, "start")]
  return(data)
}

#Distance vs. Hit comparisons 
hist(log10(all_rf_factors$Gene.Distance))
hist(log10(all_rf_factors[all_rf_factors$HitCore == TRUE, ]$Gene.Distance))
t.test(all_rf_factors$Gene.Distance, all_rf_factors[all_rf_factors$HitCore == TRUE, ]$Gene.Distance)

#separate out components of position
all_rf_factors <- extract_pos(all_rf_factors,"Enh.Pos", out_name = "Enh.")
all_rf_factors <- extract_pos(all_rf_factors,"Gene.Pos", out_name = "Gene.")
unique_enh <- unique(all_rf_factors[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")])
hit_enhs <- unique(all_rf_factors[all_rf_factors$HitPermissive,c("Enh.chr", "Enh.start", "Enh.end", "Enh")])
nonhit_enhs <- unique(all_rf_factors[!all_rf_factors$HitPermissive,c("Enh.chr", "Enh.start", "Enh.end", "Enh")])




write.bed(unique_enh, file = "Data/PeakData/Tested_Enh_nochr.bed")
write.bed(hit_enhs, "Data/PeakData/Hit_Enhancers.bed")
write.bed(nonhit_enhs, "Data/PeakData/NonHit_Enhancers.bed")

unique_enh$Enh.chr <- paste0("chr", unique_enh$Enh.chr)
write.table(unique_enh, file = "Data/PeakData/Tested_Enh.bed", row.names = FALSE, quote = FALSE, col.names = FALSE , sep = "\t")

#For creating Hi-C contacts (Hamid Rokny)
all_rf_factors$Enh.size <- all_rf_factors$Enh.end - all_rf_factors$Enh.start
all_rf_factors$Gene.TSS <- sub(".*:","",all_rf_factors$Gene.TSS)
all_rf_factors$Enh.Midpoint <- all_rf_factors$Enh.end - round(all_rf_factors$Enh.size / 2)
write.table(all_rf_factors[,c("Enh.chr", "Enh.Midpoint", "Gene.TSS", "Pair")],file = "Data/GeneData/Astrocytes_HG38_Enh_Midpoint_Gene_TSS.txt", row.names = F, quote = F, sep = "\t")

#makeWindowBed(unique(all_rf_factors[,c('Gene', "Gene.chr", "TSS"]), )
gene_promoters <- unique(all_rf_factors[,c('Gene', "Gene.chr", "Gene.TSS", "Gene.TSS")])
gene_promoters$Gene.TSS <-  as.numeric(gene_promoters$Gene.TSS) - 500
gene_promoters$Gene.TSS.1 <-  as.numeric(gene_promoters$Gene.TSS.1) + 500
write.bed(gene_promoters[,c(2,3,4,1)], file = "Data/GeneData/Gene_Promoters.bed")



#Add TTseq data
ttseq <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/4_EnhancerTranscription/TTseq/Results Table.csv")
colnames(ttseq) <- paste0("TTseq.", colnames(ttseq))
colnames(ttseq)[1] <- "Enh"
all_rf_factors <- merge(all_rf_factors, ttseq[,c("Enh", "TTseq.TTseq_Total", "TTseq.TPM_TTseq_Total","TTseq.TT_Enrich", "TTseq.RNAseq_Total", "TTseq.Ratio_TTversusRNA")])
all_rf_factors$TTseq.Diff_TTRNA <- all_rf_factors$TTseq.TTseq_Total - all_rf_factors$TTseq.RNAseq_Total


#################
#ABC scores 
#CHANGED FROM CHIP_predictions/ to Old_Peaks_NHA_peaks (currently using ABC_Predictions which are the old ones that are wrong)
#The expression table is wrong plus so we load all predictions and subset.
#################
abc_exp <- read_pred("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Farbod_ABC/enhancer-prediction-models-main/ABC/outputs/AllPeaks_NHAChip_predictions/EnhancerPredictionsAllPutative_NHAPeaks.filtered.bed", pred_non_cols) #pred_cols
abc_non <- read_pred("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Farbod_ABC/enhancer-prediction-models-main/ABC/outputs/AllPeaks_NHAChip_predictions/EnhancerPredictionsNonExpressedGenes_NHAPeaks.filtered.bed", pred_non_cols)

length(unique(all_rf_factors$Enh.Pos))

#show all EGPs & Enhancers are in ABC data
sum(all_rf_factors$Enh.Pos_Gene %in% abc_exp$Enh.Pos_Gene | all_rf_factors$Enh.Pos_Gene %in% abc_non$Enh.Pos_Gene)
sum(unique(all_rf_factors$Enh.Pos) %in% unique(abc_exp$Enh.Pos) | unique(all_rf_factors$Enh.Pos) %in% unique(abc_non$Enh.Pos))
unique(all_rf_factors[! (all_rf_factors$Enh.Pos_Gene %in% abc_exp$Enh.Pos_Gene | all_rf_factors$Enh.Pos_Gene %in% abc_non$Enh.Pos_Gene), "Enh.Pos"])
all_rf_factors[which(! all_rf_factors$Enh.Pos_Gene %in% abc_non$Enh.Pos_Gene & ! all_rf_factors$Enh.Pos_Gene %in% abc_exp$Enh.Pos_Gene),]

#currently just taking the first ABC score that so we don't double up
abc_non <- abc_non[match(unique(abc_non$Enh.Pos_Gene), abc_non$Enh.Pos_Gene),c("Enh.Pos_Gene", "ABC_Score", "hic_contact")]
abc_exp <- abc_exp[match(unique(abc_exp$Enh.Pos_Gene), abc_exp$Enh.Pos_Gene),c("Enh.Pos_Gene", "ABC_Score", "hic_contact")]
joined_abc <- rbind(abc_non, abc_exp)
all_rf_factors <- merge( all_rf_factors, joined_abc, by = "Enh.Pos_Gene", all.x = TRUE)
all_rf_factors$ABC.Exists <- all_rf_factors$Enh.Pos_Gene %in% abc_exp$Enh.Pos_Gene | all_rf_factors$Enh.Pos_Gene %in% abc_non$Enh.Pos_Gene

#ABC performance
aucPlots(all_rf_factors$ABC_Score,all_rf_factors$HitPermissive )
aucPlots(all_rf_factors[all_rf_factors$Gene.Exp > 0.1,]$ABC_Score,all_rf_factors[all_rf_factors$Gene.Exp > 0.1,]$HitPermissive )


#############
#Add ATAC seq data
#############
nha_peaks <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Final_List/NHA_Peaks_Annotated.csv") 
str(nha_peaks)
narrowPeaks <- read.table("../PreprocessedData/Bulk_ATAC_RNAseq/GOK8505-GOK8837/GOK8837A1/Mapping/MACS2/NHA_ATAC_S3.filtered.BAM_peaks.narrowPeak")
colnames(narrowPeaks) <- NarrowPeaksHeader
narrowPeaks <- narrowPeaks[narrowPeaks$name %in% nha_peaks$name, c("name","Peakscore", "signalValue", "pValue")]
narrowPeaks <- merge(narrowPeaks, nha_peaks[,c("id", "name")])
nha_peaks <- nha_peaks[, c("pileup", "id", "X.LOG10.pvalue.", "fold_enrichment", "nRep")]

#Clean up variable names
colnames(nha_peaks) <- c("ATACseq.Pileup", "Enh.Pos", "Peaks.LOG10P", "Peaks.fold_enrichment", "Peaks.nRep")
colnames(narrowPeaks) <- c("name", "NarrowPeaks.Peakscore", "NarrowPeaks.signalValue", "NarrowPeaks.pValue", "id")

#all_rf_factors <- all_rf_factors[,! colnames(all_rf_factors) %in% c("NarrowPeaks.Peakscore", "NarrowPeaks.signalValue", "NarrowPeaks.pValue", "Peaks.Pileup", "Peaks.LOG10P", "Peaks.fold_enrichment", "Peaks.nRep")]
all_rf_factors <- merge(all_rf_factors, nha_peaks)
all_rf_factors <- merge(all_rf_factors, narrowPeaks[,c("id","NarrowPeaks.Peakscore", "NarrowPeaks.signalValue", "NarrowPeaks.pValue")], by.x = "Enh.Pos", by.y = "id")


##############
#Add Tobias Data
##############
#Total Bound, Unbound & Expressed 
Tobias_TF <- read.csv("Results/Tobias/Summaries/TOBIAS_TF_counts.csv")
Tobias_TF <- Tobias_TF[, c("Enh", "Bound_TF_counts", "Unbound_TF_counts", "Exp_Bound_TF_counts", "Exp_Unbound_TF_counts")]
colnames(Tobias_TF)[colnames(Tobias_TF) != "Enh"] <-  paste0("Tobias.", colnames(Tobias_TF)[colnames(Tobias_TF) != "Enh"])
all_rf_factors <- merge(all_rf_factors, Tobias_TF)
all_rf_factors$Tobias.Bound_ratio <- all_rf_factors$Tobias.Bound_TF_counts/ (all_rf_factors$Tobias.Bound_TF_counts + all_rf_factors$Tobias.Unbound_TF_counts) 

#Proportion of TFs shared with target Gene
# More_TFs <- read.csv("Results/Tobias/Summaries/Gene_Enh_shared_TFs.csv")[,1:5]
# More_TFs_exp <- read.csv("Results/Tobias/Summaries/Exp_Gene_Enh_shared_TFs.csv")[,1:5]
# colnames(More_TFs)[colnames(More_TFs) != "Pair"] <- paste0("Tobias.", colnames(More_TFs)[colnames(More_TFs) != "Pair"])
# colnames(More_TFs_exp)[colnames(More_TFs_exp) != "Pair"] <- paste0("Tobias.", colnames(More_TFs_exp)[colnames(More_TFs_exp) != "Pair"])
# all_rf_factors <- merge(all_rf_factors, More_TFs) #, "Intersecting_TFs"
# all_rf_factors <- merge(all_rf_factors, More_TFs_exp)

#Add matrix of tobias TFs (This automatically selected expressed)
all_rf_factors<- addTobiasBound(all_rf_factors, bound_file = "Results/Tobias/Footprint/BINDetect/bound_overlaps.bed",
                                unbound_file = "Results/Tobias/Footprint/BINDetect/unbound_overlaps.bed") 
all_rf_factors[is.na(all_rf_factors)] <- 0  


############
##super enhancers annotation 
###########
super_enh<- read.csv("../FullScale/Results/3_HitEnrichment/Chromatin/Final - Annotation Logical.csv")
super_enh <- super_enh[,c("Enh", "Superenhancer")]
colnames(super_enh) <- c("Enh", "Hnisz_SuperEnhancers")
all_rf_factors <- merge(all_rf_factors, super_enh)

#TADs
TADs <- read.csv("../FullScale/Results/3_HitEnrichment/TAD - EGP Annotation.csv")
all_rf_factors <- merge(all_rf_factors, TADs[,c("Pair","WithinTAD")], by = "Pair", all.x = TRUE)

############
#Add chip astrocyte datasets used by Enformer
###########
folder <- "Results/Process_Chip/"
filenames <- list.files(folder, pattern = ".tab")
filenames <- filenames[str_detect(filenames, pattern = "Fulco|Gasp|Intergenic", negate = T)]
all_rf_factors <- addChipBigWigs(all_rf_factors, filenames, folder )
all_rf_factors$ABC_Approx <- all_rf_factors$Chip.H3K27ac / all_rf_factors$Gene.Distance
#Performs better than actual ABC
aucPlots(all_rf_factors$ABC_Approx, all_rf_factors$HitPermissive, "ABC Approx")


############
#Add Hi-C  from Hamid
############
hicContacts <- read.csv("Data/Hamid_HiC/Astrocyte_Enh_HiC.csv")
hicContacts <- hicContacts[,4:8]
colnames(hicContacts) <- c("Pair", "HiC_interaction_astrocyte_cerebellum", 
                           "HiC_interaction_astrocyte_spinal", "HiC_interaction_NeuN_pos", "HiC_interaction_NeuN_neg")
all_rf_factors <- merge(all_rf_factors, hicContacts, by = "Pair", all.x = T)
all_rf_factors$hic_contact <- all_rf_factors$HiC_interaction_NeuN_neg #overwrite HiC from ABC with NeuN_minus

########
#bulk Gene expression (and TT.seq)
########
feature_counts <- read.csv("Data/GeneData/RNAseq_FeatCounts.csv")
colnames(feature_counts) <- paste0("Gene.", colnames(feature_counts))
colnames(feature_counts)[1] <- "Gene"
all_rf_factors <- merge(all_rf_factors,feature_counts, all.x = T)
ref_feature_counts <- read.csv("Data/GeneData/Ref_RNAseq_FeatCounts.csv")
colnames(ref_feature_counts) <- paste0("Gene.", colnames(ref_feature_counts))
colnames(ref_feature_counts)[1] <- "Gene"
all_rf_factors <- all_rf_factors[,! colnames(all_rf_factors) %in% c("Gene.Ref_RNAseq_Counts", "Gene.Ref_RNAseq_CPM", "Gene.Ref_RNAseq_RPKM")]
all_rf_factors <- merge(all_rf_factors,ref_feature_counts, all.x = T)



##############
## Add New ENCODE & ABC predictions
##############
ENCODE_colnames <-c("Encode.Enh.chr",	"Encode.Enh.start",	"Encode.Enh.end", "Encode.Enh",
                    "class",	"Gene",	"EnsID",	"TargetGeneTSS",
                    "isSelfPromoter",	"CellType",	'numTSSEnhGene.Feature',
                    "distanceToTSS.Feature",	"normalizedDNase_enh.Feature",
                    "normalizedDNase_prom.Feature",	"numNearbyEnhancers.Feature",
                    "sumNearbyEnhancers.Feature",	"ubiquitousExpressedGene.Feature",
                    "numCandidateEnhGene.Feature",	"3DContactAvgHicTrack2.Feature",
                    "3DContactAvgHicTrack2_squared.Feature",	"activityEnhDNaseOnlyAvgHicTrack2_squared.Feature",
                    "activityPromDNaseOnlyAvgHicTrack2.Feature",	"ABCScore_ENCODE",
                    "DNAseOnly_Score")
pred.files <- list.files("Results/ENCODE_rE2GPredictions/",pattern = "ENCFF.*.bed", full.names = T)
astro_preds <- lapply(pred.files, function (file) {
  res <- read.table(file)
  colnames(res) <- c("Enh.chr", "Enh.start", "Enh.end", "Enh", ENCODE_colnames)
  res <- res[,c("Enh","EnsID", "CellType", "ABCScore_ENCODE","DNAseOnly_Score")]
  ct <- unique(res$CellType)
  colnames(res)[4:5] <- paste0(colnames(res)[4:5], ".", ct)
  res <- res[! is.na(res$EnsID),! colnames(res) %in% c("CellType")]
  return(res)
})
#Note we do not have predictions for all EGPs so will have to conisder how to plot them
#Their LOCs/LINCs have no ENSIDs and don't use the same naming convention
preds  <- data.frame(Enh = character(), EnsID = character())
for (df in astro_preds) {
  preds <- merge(preds, df, all = T)
}
astro_preds <- aggregate(.~Enh+EnsID,preds, mean, na.rm = T,na.action = NULL)
#all_rf_factors <- all_rf_factors[,! str_detect(colnames(all_rf_factors) ,"ABCScore_ENCODE|DNAseOnly_Score")]
all_rf_factors <- merge(all_rf_factors, astro_preds, all.x = T)
colnames(all_rf_factors)[colnames(all_rf_factors) %in% c("ABCScore_ENCODE.astrocyte", "DNAseOnly_Score.astrocyte")] <- c("ABCScore_ENCODE", "rE2G.DNAseOnly")

##################
##Add beluga predictions (currently not used)
##################
beluga_DIS_res <- read.csv("Results/Beluga/KnownVariants/MaxBelugaDiseaseScores.csv")
all_rf_factors <- merge(all_rf_factors, beluga_DIS_res, all.x = T)
all_rf_factors[is.na(all_rf_factors$Beluga.MaxDisScore),"Beluga.MaxDisScore"] <- 0

#################
##Add House keeping/sc stability index
################
colnames(geneinfo)[colnames(geneinfo) == "Symbol"] <- "Gene"
#Updating this function added another more housekeeping genes
#all_rf_factors <- all_rf_factors[,! str_detect(colnames(all_rf_factors) ,"Gene.Housekeeping|Gene.StabilityIndex")]
all_rf_factors <- addGeneStability(all_rf_factors, geneinfo = geneinfo, ensembl)

#Add num TSS between Gene and enhancer
geneinfo$Chr <- sub("chr", "", geneinfo$Chr)
all_rf_factors$TSS <- all_rf_factors$Gene.TSS 
all_rf_factors$numTSSEnhGene <- addnumTSSbetween(all_rf_factors, geneinfo)



#Save file
all_rf_factors$ENCODE_overlap.Exists <- ! is.na(all_rf_factors$rE2G.DNAseOnly)
all_rf_factors <- all_rf_factors[order(all_rf_factors$Pair),] #stop changes in predictions due to seed
write.csv(all_rf_factors, "Data/all_rf_factors.csv", row.names = FALSE) 
#all_rf_factors <- read.csv("Data/all_rf_factors.csv")
#Get test stats
removed_variables <- removed_variables[removed_variables != "HitCore"]
random_forest_data <- all_rf_factors[, ! colnames(all_rf_factors) %in% removed_variables]
getTestStats(random_forest_data)


