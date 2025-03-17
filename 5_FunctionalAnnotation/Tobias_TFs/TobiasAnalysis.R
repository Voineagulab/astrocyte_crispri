#This script analyzes Tobias outputs
#
#23-03-03
#
#Sam Bagot
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPseq/EnhancerPredictionModels/")
source("Scripts/Header_functions.R")
source("Scripts/Tobias_TFs/TobiasRFunctions.R")
source("Scripts/Annotation_Functions.R")
library(tidyr)
library(gprofiler2)

###############
## load data
############### 
homo_TFs <- read.csv("Results/Tobias/ExpressedTFs/Expressed_TFs.csv")
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/Bulk_ATAC_RNAseq/RESULTS/RNAseq/geneData.rda")
RPKM <- geneData$geneRPKM
########
#Load bound & unbound overlaps 
final_results <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")  #temp_data/Results_Final_COPY.csv
gene_info <- read.table("../FullScale/Data/Whitelists/GeneInfo.txt")
#There are some duplicates 
gene_info <- unique(gene_info[,c("Symbol", "EnsID")])
gene_info <- gene_info[match(unique(gene_info$Symbol), gene_info$Symbol),]
final_results <- merge(final_results, unique(gene_info[,c("Symbol", "EnsID")]), by.x = "Gene", by.y = "Symbol")

#Load and format overlaps data
unbound_overlaps <- read.table("Results/Tobias/Footprint/BINDetect/unbound_overlaps.bed")
bound_overlaps <- read.table("Results/Tobias/Footprint/BINDetect/bound_overlaps.bed")
unbound_overlaps <- update_overlaps(unbound_overlaps, final_results )
bound_overlaps <- update_overlaps(bound_overlaps, final_results)
exp_bound_overlaps <- bound_overlaps[bound_overlaps$TF_Name %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names,]
exp_unbound_overlaps <- unbound_overlaps[unbound_overlaps$TF_Name %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names,]


##########
##
##This section compares the total number of TF Bound and Unbound between Hits and non-hits 
##
#########
#aggregate number of overlaps and merge objects
bound_counts <- setnames(aggregate(by = list(bound_overlaps$Enh.Pos),bound_overlaps$TF_Name, length),c("Enh.Pos", "Bound_TF_counts"))
unbound_counts  <- setnames(aggregate(by = list(unbound_overlaps$Enh.Pos),unbound_overlaps$TF_Name, length),c("Enh.Pos", "Unbound_TF_counts"))
exp_bound_counts <- setnames(aggregate(by = list(exp_bound_overlaps$Enh.Pos),exp_bound_overlaps$TF_Name, length),c("Enh.Pos", "Exp_Bound_TF_counts"))
exp_unbound_counts  <- setnames(aggregate(by = list(exp_unbound_overlaps$Enh.Pos),exp_unbound_overlaps$TF_Name, length),c("Enh.Pos", "Exp_Unbound_TF_counts"))


unbound_bound <- merge(bound_counts, unbound_counts, all = T)
exp_unbound_bound <- merge(exp_bound_counts, exp_unbound_counts, all = T)
unbound_bound <- merge(unbound_bound, exp_unbound_bound, all = T)
unbound_bound$Hit <- unbound_bound$Enh.Pos %in% unique(final_results[final_results$HitPermissive == TRUE,]$Enh.Pos)
unbound_bound <- merge(unbound_bound, unique(final_results[,c("Enh.Pos", "Enh")]), all.y = T)
unbound_bound[is.na(unbound_bound)] <- 0
t.test(unbound_bound[!unbound_bound$Hit,]$Exp_Bound_TF_counts, unbound_bound[unbound_bound$Hit,]$Exp_Bound_TF_counts)
t.test(unbound_bound[!unbound_bound$Hit,]$Exp_Unbound_TF_counts, unbound_bound[unbound_bound$Hit,]$Exp_Unbound_TF_counts)
write.csv(unbound_bound, "Results/Tobias/Summaries/TOBIAS_TF_counts.csv", row.names = FALSE)
#unbound_bound <- read.csv("Results/Tobias/Summaries/TOBIAS_TF_counts.csv")


#########
#Create tracks for UCSC browser with header
#########

write.table(exp_bound_overlaps[,c("TF_chr", "TF_start","TF_end","TF_Name")], file = "Results/Tobias/Summaries/expressed_combined_bound_overlaps.bed",row.names = F, col.names = F, quote = F, sep = "\t")
system("sed -i '1s/^/track name=\"Bound exp TFs\"\n/' Results/Tobias/Summaries/expressed_combined_bound_overlaps.bed")
bound_combined <- fread("Results/Tobias/Footprint/BINDetect/combined.bed")
bound_combined$V4 <- sub(".*\\.(.*)_.*","\\1",bound_combined$V4)
bound_combined <- bound_combined[bound_combined$V4 %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names,]
write.table(bound_combined[,1:4], file = "Results/Tobias/Summaries/expressed_combined_bound.bed",row.names = F, col.names = F, quote = F, sep = "\t")
system("sed -i '1s/^/track name=\"Bound exp TFs\"\n/' Results/Tobias/Summaries/expressed_combined_bound.bed")


####################
##
##Individual TF analysis
##
###################
#Are any TFs never bound or Never Unbound
unique(unbound_overlaps$TF_Name)[! unique(unbound_overlaps$TF_Name) %in% unique(bound_overlaps$TF_Name)]
sum(unbound_overlaps$TF_Name == "ONECUT2")
unique(bound_overlaps$TF_Name)[ ! unique(bound_overlaps$TF_Name) %in% unique(unbound_overlaps$TF_Name)]
sum(bound_overlaps$TF_Name == "ZBED1")

#Sanity check for Bound TFs being more highly expressed (~%5 more of Bound TFs are expressed)
summary <- compare_overlaps(unbound_overlaps, bound_overlaps, homo_TFs, RPKM)
write.csv(summary, "Results/Tobias/Summaries/TF_Expressed.csv")

union(unbound_overlaps[! unbound_overlaps$TF_Name %in%  homo_TFs[homo_TFs$V1_Exists & homo_TFs$V2_Exists,]$concat_names,"TF_Name"], 
      bound_overlaps[! bound_overlaps$TF_Name %in%  homo_TFs[homo_TFs$V1_Exists & homo_TFs$V2_Exists,]$concat_names,"TF_Name"])

#Get bound and unbound Overlaps for each enhancer
enhancers <- read.csv("Data/Enhancer_factors.csv")
tmp <- enhancers[,c("Enh", "Enh.Pos", "HitPermissive", "HitPermissive_NegZ")]
bound_matrix <- addTobiasBound(tmp, bound_file = "Results/Tobias/Footprint/BINDetect/bound_overlaps.bed",
               unbound_file = "Results/Tobias/Footprint/BINDetect/unbound_overlaps.bed")
bound_matrix_unexp <- addTobiasBound(tmp, bound_file = "Results/Tobias/Footprint/BINDetect/bound_overlaps.bed",
                               unbound_file = "Results/Tobias/Footprint/BINDetect/unbound_overlaps.bed", expressed = F)
write.csv(bound_matrix, "Results/Tobias/Summaries/Bound_Matrix_Expressed.csv", row.names = F)
write.csv(bound_matrix_unexp, "Results/Tobias/Summaries/Bound_Matrix.csv", row.names = F)

#Get Boolean Fisher Test Results for expressed TFs only
bound_matrix_TrueFalse <- bound_matrix[,! colnames(bound_matrix) %in% c("Enh", "Enh.Pos")] > 0 

test_boundMatrix <- do.call("rbind",getTestStats(bound_matrix_TrueFalse, FT.only = T))
write.csv(test_boundMatrix[order(test_boundMatrix$p),], "Results/Tobias/Summaries/FisherTestsExpressedTFMatrix.csv", row.names = F)

# test_boundMatrix <- getTestStats(bound_matrix_TrueFalse)
# colnames(test_boundMatrix)[colnames(test_boundMatrix) == "T.Stat_Odds_ratios"] <- "Odds_ratio" #OR "T-Stat_Odds_ratios"
# test_boundMatrix <- test_boundMatrix[,colnames(test_boundMatrix) != "Method"]
# test_boundMatrix$Variable <- sub("_counts","",test_boundMatrix$Variable)
# write.csv(test_boundMatrix, "Results/Tobias/Summaries/FisherTestsExpressedTFMatrix_oldformat.csv", row.names = F)
#test_boundMatrix <- read.csv("Results/Tobias/Summaries/FisherTestsExpressedTFMatrix.csv")
