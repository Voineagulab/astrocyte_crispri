#################
##
##Tobias Expressed genes. 
## @author: Sam Bagot
##This script processes TF binding motifs downloaded from JASPAR 
##
#################
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPseq/EnhancerPredictionModels/")
source("Scripts/Header_functions.R")
source("Scripts/Tobias_TFs/TobiasRFunctions.R")

#These need to be generated BEFORE running the TOBIAS scripts. 
K562_exp <- read.csv("Data/K562Data/RNAseq/Ref_RNAseq_FeatCounts.csv")
K562_homo_TFs <- gethomo_TFs(K562_exp[,c("EnsID", "Ref_RNAseq_RPKM")], K562_exp$EnsID)
write.table(paste0("\\.",K562_homo_TFs[K562_homo_TFs$TF.1_Expressed & K562_homo_TFs$TF.2_Expressed,"concat_names"], "_"), file = "Results/Tobias/ExpressedTFs/K562_Expressed_TFs_regex.txt", quote = FALSE, row.names =  FALSE, col.names = FALSE)
write.csv(K562_homo_TFs, file = "Results/Tobias/ExpressedTFs/K562_Expressed_TFs.csv", row.names = F)

load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/Bulk_ATAC_RNAseq/RESULTS/RNAseq/geneData.rda")
RPKM <- geneData$geneRPKM
homo_TFs <- gethomo_TFs(RPKM[,c("EnsID", "NHA")], geneData$geneCounts$EnsID)
write.table(paste0("\\.",homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,"concat_names"], "_"), file = "Results/Tobias/ExpressedTFs/Expressed_TFs_regex.txt", quote = FALSE, row.names =  FALSE, col.names = FALSE)
write.csv(homo_TFs, file = "Results/Tobias/ExpressedTFs/Expressed_TFs.csv", row.names = F)

RPKM_noNA <- RPKM[!is.na(RPKM$Symbol),]

#679/680 the fused protein is omitted Dimers TF combos where both genes are in data 
sum(homo_TFs$TF.1_Exists & homo_TFs$TF.2_Exist)
homo_TFs[! homo_TFs$TF.1_Exists & ! homo_TFs$TF.2_Exist,]
#630/631 unique TFs exist
length(union(homo_TFs[homo_TFs$TF.2_Exists,"TF.2"], homo_TFs[homo_TFs$TF.1_Exists,"TF.1"]))
#631 unique TFs in data
union(homo_TFs$TF.2, homo_TFs$TF.1)
#387 genes are in the RPKM data 
sum(RPKM$EnsID %in% homo_TFs$ENS.1 | RPKM$EnsID %in% homo_TFs$ENS.2)
#331/630 are >= 0.5 RPKM (there is an NA)
length(union(homo_TFs[homo_TFs$TF.1_Expressed,"TF.1"],homo_TFs[homo_TFs$TF.2_Expressed,"TF.2"]))
#281
length(union(K562_homo_TFs[K562_homo_TFs$TF.1_Expressed,"TF.1"],K562_homo_TFs[K562_homo_TFs$TF.2_Expressed,"TF.2"])) #K562
#416 
length(union(homo_TFs[homo_TFs$TF.1_RPKM_Exists,"TF.1"],homo_TFs[homo_TFs$TF.2_RPKM_Exists,"TF.2"]))
