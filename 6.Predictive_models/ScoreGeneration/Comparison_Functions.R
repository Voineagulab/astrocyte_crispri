#ComparisonPlots and RF building functions

library("biomaRt")
library("cowplot")
library("stringr")
library("seqinr")

###################
#Variables used in models
selected_vars_exp <- c("Chip.H3K4me3", "Chip.H3K27ac","Gene.Nearest", 
                   "Distance","ATACseq.Pileup", "numTSSEnhGene", "Gene.StabilityIndex")
# selected_vars_extended <- c("Chip.H3K4me3", "TTseq.TTseq_Total", "Chip.H3K27ac",
#                             "Gene.Nearest", "Distance", "Tobias.Exp_Bound_TF_counts",  
#                             "ATACseq.Pileup", "numTSSEnhGene", 
#                             "Gene.StabilityIndex", "Beluga.MaxDisScore")

selected_vars_exp_replace_ATAC <- c(selected_vars_exp[selected_vars_exp != "ATACseq.Pileup"], "DNAse") #Note this is not quite the same I use bigwig here
plus_hic_vars <- c(selected_vars_exp, "HiC_interaction_NeuN_neg") #May need to replace this with actual contacts not from ABC pipeline #HiC_interaction_NeuN_neg
all_vars_astro <- c(selected_vars_exp, "Tobias.Exp_Bound_TF_counts","TTseq.TTseq_Total", 
                    "TTseq.Ratio_TTversusRNA", "Beluga.MaxDisScore","HiC_interaction_NeuN_neg")
tap_seq_vars <- c("Distance","HiC_interaction_NeuN_neg", 
                  "Chip.H3K27ac", "Chip.H3K4me3","DNAse", 
                  "Chip.POLR2A", "Chip.H3K4me1", "Chip.H3K27me3")

#Variables used in encode models
selected_vars_exp_encode <- c(selected_vars_exp[selected_vars_exp != c( "ATACseq.Pileup") ], 
                          "DNAse.Pileup")

tap_seq_vars_encode<- c(tap_seq_vars[!tap_seq_vars %in% c("DNAse","HiC_interaction_NeuN_neg")],
                        "DNAse.Pileup", "hic_contact")
combined_model_encode <- c(selected_vars_exp_encode,"hic_contact",
                           "normalizedEP300_enhActivity", "activity_prom")  
all_vars_encode <- c(combined_model_encode,"Tobias.Exp_Unbound_TF_counts", "Tobias.Exp_Bound_TF_counts",
                     "TTseq.TTseq_Total", "TTseq.Ratio_TTversusRNA")
                            #OR "Gene.TTseq_Counts", "Gene.RNAseq_RPKM"
replaceEnhactivitywithExp <- c(combined_model_encode[combined_model_encode != "activity_prom"],"Gene.TTseq_Counts", "Gene.RNAseq_RPKM")
encode_ext_model <- read.table("Results/ENCODE_rE2GPredictions/ENCODE_Predictors.txt")$V1

negativePredictors <- c("ClosestGeneDist","Distance", "Gene.StabilityIndex", "numTSSEnhGene", 
                        "Tobias.Exp_Unbound_TF_counts", "activity_prom")

plot_cols <- c(EGrf = pals$Primary[8],EGrf.Extended = pals$Primary[7],
               TAPseq.rf = pals$Primary[4], 
               Distance = "#000000",
               ABC_Score = pals$Primary[5],rE2G.DNAseOnly = pals$Primary[2],
               rE2G.Extended = pals$Primary[6], #ENCODE.
               rE2G.Extended.rf = "grey"
) 

#######################
#gets pr & roc AUC curves for a list of variables
getAUCLists <- function (df, aucvars, Hit, probability) {
  prList <- list()
  rocList <- list()
  if (sum(aucvars %in% colnames(df)) < length(aucvars)) warning("missing columns in vars list")
  aucvars <- aucvars[aucvars %in% colnames(df)]
  for (var in aucvars) {
    #negative predictors flip the variables
    
    if (var %in%  negativePredictors) {
      score <- 1 / df[,var]
    } else if (str_detect(var, "rf") & probability != T) {
      score <- - df[,var]
    } else {
      score <- df[,var]
    }
    prList <- append(prList, list(pr.curve(scores.class0 = score, weights.class0 = df[,Hit] > 0, curve = TRUE, rand.compute = TRUE)))
    rocList <- append(rocList, list(roc.curve(scores.class0 = score, weights.class0 = df[,Hit] > 0, curve = TRUE, rand.compute = TRUE)))
  }
  names(prList)  <- aucvars
  names(rocList)  <- aucvars
  return(list(prList, rocList))
}

#This function trains all models in a list
#Variable list is a named list of variables
trainAllModels <- function(all_rf_factors, variables.list,noexp_name = "noexp_rf",probability = F,
                           suffix = "", train_var = "logfc.vst", save_models = F, pred_method = "cv", folder = "Results/ComparisonPlots/") {
  #OOB with CV splitting unique enhancers and genes for training so test is unique from training (10x10)
  final_suffix <- paste0(suffix, "_",train_var,"_",pred_method )
  
  for (name in names(variables.list)) {
    message(name)
    message(paste("Vars:", paste0(variables.list[[name]], collapse = ",")))
    model<- getRFmodel(all_rf_factors, subset = variables.list[[name]], train_var = train_var, 
                       importance = paste0(folder,"/EG_Variable_Importance_",name,final_suffix ,".csv")) #
    all_rf_factors <- getRFmodelPreds(all_rf_factors, model , colname = name, pred.method = pred_method, 
                                      train_var = train_var, probability = probability, 
                                      save_models =  paste0(folder,"/SavedModels/EGPModel",name,final_suffix,".rds"))
  }
  return(all_rf_factors)
}


updateVarNames <- function(names) {
  names <- sub("Chip.","",names)
  names[names == "activity_prom"] <- "Promoter.Activity"
  names[names == "hic_contact"] <- "HiC.Contact"
  names[names =="normalizedEP300_enhActivity"] <- "EP300"
  names[names =="DNAse.Pileup"] <- "DNAse"
  names[names =="Tobias.Exp_Bound_TF_counts"] <- "Bound.TFfootprints"
  names[names =="Tobias.Exp_Unbound_TF_counts"] <- "Unbound.TFfootprints"
  names[names =="TTseq.TTseq_Total"] <- "TTseq.Counts"
  names[names =="TTseq.Ratio_TTversusRNA"] <- "TTseq.RNAseqRatio"
  return(names)
}
getVarsExcel<- function (varslist) {
  
  variables <- data.frame(variables = sort(unique(unlist(varslist))))
  for (model  in names(varslist)) {
    variables[,model] <- variables$variables %in% varslist[[model]]
  }
  variables$variables<- updateVarNames(variables$variables)
  return(variables)
}
