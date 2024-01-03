#######################
##
## Enhancer predictions & technical models for subsetting EGPs
## @author: sam bagot
## @data: 15-11-23
##
####################

setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPseq/EnhancerPredictionModels/")
source("Scripts/Header_functions.R")
#source("Scripts/ScoreGeneration/Comparison_Functions.R")
#dir.create("Results/ComparisonPlots/Enhancermodels/")

#Read in Enhancer dataframes
enhancers <- read.csv("Data/Enhancer_factors.csv")
non_technical <- c( "Hnisz_SuperEnhancers", "Tobias.Exp_Bound_TF_counts", "CHIP.H3K4me3.astrocyte_mean", #"CHIP.H3K4me1_astrocyte_mean","Peaks.fold_enrichment",
                    "Chip.H3K27ac",  "ATACseq.Pileup", "TTseq.TTseq_Total", "TTseq.Ratio_TTversusRNA") #"Tobias.Bound_FOSL2JUND",
#, "CHIP.H3K4me1_astrocyte_mean", "Peaks.fold_enrichment""NarrowPeaks.signalValue", "Chip.H3K4me1",
enh_allgenes_tech <- c("AllGenes.Max_RNAseq_RPKMLT500kb","AllGenes.Max_RNAseq_RPKMLT50kb", 
                       "AllGenes.ClosestGeneExp", "AllGenes.LT500KB", "AllGenes.ClosestGeneDist", 
                       "AllGenes.LT50KB")

#This is main function to get enhancers that have a strong change of detection 
getEnhancerBelowPrecision <- function(enhancers, suffix = "", hits = "HitPermissive_NegZ", technical, train_var,probability) {
  model <- getRFmodel(enhancers, technical,train_var = train_var, importance =  paste0('Results/ComparisonPlots/Enhancermodels/Enh_importance_allgenes_tech',suffix,'.csv'))
  saveRDS(model, paste0("Results/ComparisonPlots/Enhancermodels/EnhancerTechModel",suffix,".rds"))
  enhancers <- getRFmodelPreds(enhancers, model, colname = "all_genes_tech_rf_ref",
                               pred.method = "oob", train_var = train_var, probability = probability) 
  if (probability) {
    pr_techrf_allgenes <- pr.curve(scores.class0 = enhancers$all_genes_tech_rf_ref, weights.class0 = enhancers[,hits] > 0, curve = T, rand.compute = T)
    prec_50_score_allgenes <- pr_techrf_allgenes$curve[which.min(abs(pr_techrf_allgenes$curve[,2] - 0.50)), 3] 
    below_tech_enhancers <- enhancers[enhancers$all_genes_tech_rf_ref > prec_50_score_allgenes,]
  } else {
    pr_techrf_allgenes <- pr.curve(scores.class0 = -enhancers$all_genes_tech_rf_ref, weights.class0 = enhancers[,hits] > 0, curve = T, rand.compute = T)
    prec_50_score_allgenes <- - pr_techrf_allgenes$curve[which.min(abs(pr_techrf_allgenes$curve[,2] - 0.50)), 3] 
    below_tech_enhancers <- enhancers[enhancers$all_genes_tech_rf_ref < prec_50_score_allgenes,]
  }
  pdf(paste0("Results/ComparisonPlots/Enhancermodels/aucTechnical",suffix,".pdf"))
    print(plot(pr_techrf_allgenes, rand.plot = TRUE))
  dev.off()
  saveRDS(prec_50_score_allgenes, paste0("Results/ComparisonPlots/Enhancermodels/EnhancerTechcutoff",suffix,".rds"))
  write.table(below_tech_enhancers$Enh, file = paste0("Results/ComparisonPlots/Enhancermodels/Enhancers_below_tech",suffix,".txt"), quote = F, row.names = F, col.names = F)
  print(paste("Enhancers below technical", nrow(below_tech_enhancers)))
  return(enhancers)
}

getEnhancerNontechnical <- function (below_tech_enhancers,suffix = "", hits = "HitPermissive_NegZ",  non_technical,train_var, probability = F) {
  model <- getRFmodel(below_tech_enhancers, non_technical, train_var = train_var, 
                      importance = paste0("Results/ComparisonPlots/Enhancermodels/Allgenes_belowtech_prec",suffix,".csv") )
  saveRDS(model, paste0("Results/ComparisonPlots/Enhancermodels/EnhancerBiologicalModel",suffix,".rds"))
  below_tech_enhancers <- getRFmodelPreds(below_tech_enhancers, model, colname =  "rf_above_tech_prec", 
                                          pred.method = "oob",train_var = train_var, probability = probability)
  pdf(paste0("Results/ComparisonPlots/Enhancermodels/aucNontech",suffix,".pdf"))
   getRFperformance(below_tech_enhancers[,train_var], below_tech_enhancers$rf_above_tech_prec, 
                    below_tech_enhancers[,hits] > 0, train_var = train_var, 
                    plot.title = suffix)
  dev.off()
  return(below_tech_enhancers)
}

getEnhancerModels <- function(enhancers, suffix = "", hits = "HitPermissive_NegZ", train_var = "Z", technical, non_technical, probability = F) {
  enhancers <- getEnhancerBelowPrecision(enhancers = enhancers,suffix =  suffix, 
                                                    hits = hits, technical = technical, 
                                                    train_var = train_var,probability = probability)
  # below_tech_enhancers <- getEnhancerNontechnical(below_tech_enhancers = below_tech_enhancers,
  #                                                 suffix = suffix, hits = hits, 
  #                                                 non_technical = non_technical, 
  #                                                 train_var = train_var,probability = probability)
  return(enhancers)
}


enhancers <- getEnhancerModels(enhancers, technical = enh_allgenes_tech, non_technical = non_technical )
write.csv(enhancers, "Data/Enhancer_factors_pluspredictions.csv", row.names = F)
