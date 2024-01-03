#Compare Families of TFs 

setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPseq/EnhancerPredictionModels/")
library("ggplot2")
library(stringr)
library(ggbeeswarm)

#Investigating Hits without binding and Non-hits with
#enhancers[enhancers$Tobias.Unbound_TF_counts == 0 & enhancers$HitPermissive == 0, "Enh"]
#enhancers[enhancers$Tobias.Bound_TF_counts == 0 & enhancers$HitPermissive != 0, c("Enh","Z", "Tobias.Bound_TF_counts", "Tobias.Unbound_TF_counts" )]

#BiocManager::install("JASPAR2022")
#Both methods results in non-0 exit status
#BiocManager::install("TFBSTools", ) 
#install.packages("/Users/z5205171/Downloads/CNEr_1.38.0.tar.gz", repos = NULL, type="source")
#Use Jaspar Matrices
library(JASPAR2022)
library(TFBSTools)
opts<-list()
opts[["all_versions"]] <- TRUE
opts[["species"]] = "Homo sapiens" #Homo Sapiens
opts[["collection"]] ="CORE"
JASPAR_PFMatrixList = getMatrixSet(JASPAR2022, opts) 
saveRDS(JASPAR_PFMatrixList, file = "Results/Tobias/JASPAR_Matrices.rds")

families  <- lapply(JASPAR_PFMatrixList, function(x) {
  name <- x@name
  matrixClass <- x@matrixClass
  family <- x@tags$family
  #profileMatrix <- x@profileMatrix
  return(cbind(name, family, matrixClass)) #, profileMatrix
})
families <- as.data.frame(do.call("rbind",families ))
families$concat_names <- sub("::","",families$name)
families <- unique(families)
write.csv(families, "Results/Tobias/Summaries/JASPAR_TF_Families.csv", row.names = F)

#ALL FAMILIES ARE PRESENT IF ANYTHING APPEARS THEN THEY ARE MISSING
test_boundMatrix <- read.csv("Results/Tobias/Summaries/FisherTestsExpressedTFMatrix.csv")
test_boundMatrix[! test_boundMatrix$concat_names %in% families$concat_names,]
TF_Families <- read.delim("../PublicData/AnimalTFDB/Homo_sapiens_TF_AnimalTFDB.txt", sep="\t", header=TRUE)
homo_TFs <- read.csv( "Results/Tobias/ExpressedTFs/Expressed_TFs.csv")


test_boundMatrix <- test_boundMatrix[str_detect(test_boundMatrix$Variable, "Tobias.Bound_"),]
test_boundMatrix$concat_names <- sub("Tobias.Bound_","",test_boundMatrix$Variable)
test_boundMatrix<- merge(test_boundMatrix, homo_TFs[,c("TF.1", "TF.2", "concat_names")], by = "concat_names")
test_boundMatrix$LOG10FDR <- -log10(test_boundMatrix$FDR)
test_boundMatrix$JASPAR_Family <- families[match(test_boundMatrix$concat_names, families$concat_names),"family"]
test_boundMatrix$JASPAR_Class <- families[match(test_boundMatrix$concat_names, families$concat_names),"matrixClass"]

#Update names so they fit on plots
unique(test_boundMatrix$JASPAR_Class)
test_boundMatrix$JASPAR_Class_original <- test_boundMatrix$JASPAR_Class
test_boundMatrix$JASPAR_Class <- sub("factors", "",test_boundMatrix$JASPAR_Class)
test_boundMatrix$JASPAR_Class <- sub("domain", "",test_boundMatrix$JASPAR_Class)
test_boundMatrix$JASPAR_Class <- sub("box", "",test_boundMatrix$JASPAR_Class)
test_boundMatrix$JASPAR_Class <- sub("zinc finger", "ZF",test_boundMatrix$JASPAR_Class)
test_boundMatrix$JASPAR_Class <- sub(".*\\((.*)\\)", "\\1",test_boundMatrix$JASPAR_Class)
test_boundMatrix$JASPAR_Class <- str_trim(test_boundMatrix$JASPAR_Class)
test_boundMatrix$JASPAR_Class <- as.factor(test_boundMatrix$JASPAR_Class)
unique(test_boundMatrix$JASPAR_Class)
colnames(test_boundMatrix)[colnames(test_boundMatrix) == "T.Stat_Odds_ratios"] <- "Odds Ratio"

write.csv(test_boundMatrix, "Results/Tobias/Summaries/JASPAR_group_AnnotatedTFs.csv",row.names = F)
# 
# load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/Bulk_ATAC_RNAseq/RESULTS/RNAseq/geneData.rda")
# RPKM <- geneData$geneRPKM
# RPKM_noNA <- RPKM[!is.na(RPKM$Symbol),]
# #test_boundMatrix <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Summaries/JASPAR_group_AnnotatedTFs.csv")
# #Quick sanity check for expression and pvalue
# exp_vs_sig <- merge(test_boundMatrix, RPKM_noNA[,c("Symbol", "NHA")], by.x = "TF.1",by.y =  "Symbol", all.x = T)
# exp_vs_sig <- merge(exp_vs_sig, RPKM_noNA[,c("Symbol", "NHA")], by.x = "TF.2",by.y =  "Symbol", all.x = T)
# # exp_vs_sig[is.na(exp_vs_sig$NHA.x),"NHA.x"] <- median(exp_vs_sig$NHA.x)
# # exp_vs_sig[is.na(exp_vs_sig$NHA.y),"NHA.y"] <- median(exp_vs_sig$NHA.y)
# exp_vs_sig$Exp <- sqrt(exp_vs_sig$NHA.x^2 + exp_vs_sig$NHA.y^2)
# cor(exp_vs_sig$NHA.x, exp_vs_sig$LOG10FDR)
# cor(exp_vs_sig$NHA.y, exp_vs_sig$LOG10FDR)
# 
# plot_details <- geom_text() + facet_wrap(vars(JASPAR_Class)) + scale_color_discrete(guide = guide_legend(ncol = 1) )
# pdf("Results/Tobias/Plots/Expression_vs_significance.pdf")
#   ggplot(exp_vs_sig[!is.na(exp_vs_sig$Exp),], aes(log2(Exp), LOG10FDR, colour = JASPAR_Class, label = concat_names)) +
#     geom_text(data = exp_vs_sig[!is.na(exp_vs_sig$Exp) & exp_vs_sig$LOG10FDR > -log10(0.05),]) + 
#     geom_point(data = exp_vs_sig[!is.na(exp_vs_sig$Exp) & exp_vs_sig$LOG10FDR <= -log10(0.05),]) + 
#     facet_wrap(vars(JASPAR_Class)) + scale_color_discrete(guide = guide_legend(ncol = 1) )
#   ggplot(exp_vs_sig[!is.na(exp_vs_sig$Exp),], aes(log2(Exp), LOG10FDR, colour = JASPAR_Class, label = concat_names)) +
#     geom_point() +facet_wrap(vars(JASPAR_Class)) + stat_smooth(method = "lm", se = F)+ scale_color_discrete(guide = guide_legend(ncol = 1) )
#   ggplot(exp_vs_sig[!is.na(exp_vs_sig$Exp),], aes(log2(Exp), LOG10FDR)) +
#     geom_point(aes(colour = JASPAR_Class,)) + stat_smooth()+ scale_color_discrete(guide = guide_legend(ncol = 1) )
# dev.off()
# 
# 
# #Now in Manuscript TranscriptionFactors
# familyTFPlots <- function (test_boundMatrix, x = "JASPAR_Family", y = "LOG10FDR") {
#   maxTFs <- setNames(aggregate(by = list(test_boundMatrix[,x]), test_boundMatrix[,y], max), c(x, y))
#   max_TFs<- merge(maxTFs, test_boundMatrix[,c(x, y, "concat_names")])
#   p1 <- ggplot(test_boundMatrix, aes(x = .data[[x]], y=  .data[[y]])) +
#     geom_beeswarm() + theme_bw() + theme(axis.text.x = element_text(face = "bold",vjust = 0.5, hjust=1,  angle = 90)) + #geom_beeswarm()
#     geom_text(data=max_TFs, size = 2.5,nudge_y=0.2, aes(x = .data[[x]], y=  .data[[y]], label=concat_names)) +
#     labs(title = paste(x,"annotation"))
#   if (y == "T.Stat_Odds_ratios") {
#     p1 <- p1 + geom_hline(yintercept=1, colour="grey" )
#   } else {
#     p1 <- p1 + geom_hline(yintercept=-log10(0.05), colour="lightgrey", linetype = "dashed" )
#   }
#   return(p1)
# }
# 
# pdf("Results/Tobias/Plots/FamilyFTPlots.pdf", width = 10)
# familyTFPlots(test_boundMatrix, x = "JASPAR_Class", y = "LOG10FDR")
# familyTFPlots(test_boundMatrix, x = "JASPAR_Class", y = "T.Stat_Odds_ratios")
# familyTFPlots(test_boundMatrix, x = "JASPAR_Family", y = "LOG10FDR")
# familyTFPlots(test_boundMatrix, x = "JASPAR_Family", y = "T.Stat_Odds_ratios")
# dev.off()
# 
# test_bound_sig <- test_boundMatrix[test_boundMatrix$FDR < 0.05,] 
# pdf("Results/Tobias/Plots/FamilyFTPlots.pdf", width = 10)
# familyTFPlots(test_bound_sig, x = "JASPAR_Class", y = "LOG10FDR")
# familyTFPlots(test_bound_sig, x = "JASPAR_Class", y = "T.Stat_Odds_ratios")
# familyTFPlots(test_bound_sig, x = "JASPAR_Family", y = "LOG10FDR")
# familyTFPlots(test_bound_sig, x = "JASPAR_Family", y = "T.Stat_Odds_ratios")
# dev.off()
# 
# 
# test_boundMatrix[! test_boundMatrix$TF.1 %in%  TF_Families$Symbol,]
# homo_TFs$concat_names
# test_boundMatrix
# 
# 
# ftb.sorted=ftb[order(ftb$ftest_pval.adj), ]
# 
# #CRC Significance test
# crcs <- read.csv("../PublicData/SaintAndre_CRC/CRC_Astro.csv")
# test_boundMatrix$crc <- test_boundMatrix$concat_names %in% crcs$Astrocytes
# fisher.test(test_boundMatrix$crc, test_boundMatrix$FDR<0.05)
# 
# 
# test_boundMatrix[order(test_boundMatrix$T.Stat_Odds_ratios, decreasing = T),"concat_names"]
# test_boundMatrix[order(test_boundMatrix$LOG10FDR, decreasing = T),"concat_names"]
# 
# #Correlate pvalue & odds ratio with expression
# feature_counts <- read.csv("Data/GeneData/RNAseq_FeatCounts.csv")
# colnames(feature_counts) <- paste0("Gene.", colnames(feature_counts))
# colnames(feature_counts)[1] <- "Gene"
# test_boundMatrix<- merge(test_boundMatrix, feature_counts[,c("Gene","Gene.RNAseq_RPKM")], all.x = T, by.x = "TF.1", by.y = "Gene")
# test_boundMatrix <- merge(test_boundMatrix, feature_counts[,c("Gene","Gene.RNAseq_RPKM")], all.x = T, by.x = "TF.2", by.y = "Gene")
# test_boundMatrix$Gene.RNAseq_RPKM_geometric_mean <- sqrt((test_boundMatrix$Gene.RNAseq_RPKM.x^2+ test_boundMatrix$Gene.RNAseq_RPKM.x^2))
# 
# cor(test_boundMatrix[!is.na(test_boundMatrix$Gene.RNAseq_RPKM_geometric_mean),]$Gene.RNAseq_RPKM_geometric_mean, 
#     test_boundMatrix[!is.na(test_boundMatrix$Gene.RNAseq_RPKM_geometric_mean),]$T.Stat_Odds_ratios)
# cor(test_boundMatrix[!is.na(test_boundMatrix$Gene.RNAseq_RPKM_geometric_mean),]$Gene.RNAseq_RPKM_geometric_mean, 
#     test_boundMatrix[!is.na(test_boundMatrix$Gene.RNAseq_RPKM_geometric_mean),]$LOG10FDR)
# 
# 
