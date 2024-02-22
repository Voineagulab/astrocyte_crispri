##################
## This script carries out all the random forest model training and prediction presented in Green et al 2024.
## @author: Sam Bagot
## @date: 14-11-23
## Modified by IV and JP 2_TrainEGPModels_IVtest
#################
rm(list=ls())
path <-"/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/"
setwd(paste0(path, "/1.Data/"))
################### Source functions
source(paste0(path,"/2.Scripts/Functions.R"))

#########################################################
#Define variables for each model
#########################################################

#########Astrocytes
#EGrf variables
selected_vars_exp <- c("Chip.H3K4me3", "Chip.H3K27ac","Gene.Nearest", 
                       "Distance","ATACseq.Pileup", "numTSSEnhGene", "Gene.StabilityIndex")
#EGrf models with one variable removed at a time
removed.vars_astro <- list()
for (var in selected_vars_exp) {
  removed.vars_astro[[paste0(var, ".removed")]] <- selected_vars_exp[selected_vars_exp != var]
}

#TAPseq.rf variables
tap_seq_vars <- c("Distance","HiC_interaction_NeuN_neg", 
                  "Chip.H3K27ac", "Chip.H3K4me3","DNAse", 
                  "Chip.POLR2A", "Chip.H3K4me1", "Chip.H3K27me3")

# Preparing a list of variables for all astrocyte models
vars.list <- append(list(TAPseq.rf = tap_seq_vars,
                         EGrf = selected_vars_exp,
                         EGrf.PlusHiC = c(selected_vars_exp, "HiC_interaction_NeuN_neg"),
                         EGrf.PlusTobias = c(selected_vars_exp,"Tobias.Exp_Unbound_TF_counts", "Tobias.Exp_Bound_TF_counts"), 
                         EGrf.PlusTTseq = c(selected_vars_exp,"TTseq.TTseq_Total", "TTseq.Ratio_TTversusRNA"),
                         EGrf.PlusBelugaDIS = c(selected_vars_exp,"Beluga.MaxDisScore"),
                         EGrf.ReplaceATACwithDNAse = c(selected_vars_exp[selected_vars_exp != "ATACseq.Pileup"], "DNAse")), 
                    removed.vars_astro)

#########K562
# EGrf variables (using DNase-seq rather than ATAC-seq for open chromatin state)
selected_vars_exp_encode <- c(selected_vars_exp[selected_vars_exp != c( "ATACseq.Pileup") ], 
                              "DNAse.Pileup")
# TAPseq.rf variables (using the variable names form the K562 dataframe)
tap_seq_vars_encode<- c(tap_seq_vars[!tap_seq_vars %in% c("DNAse","HiC_interaction_NeuN_neg")],
                        "DNAse.Pileup", "hic_contact")

# EGrf-extended variables
combined_model_encode <- c(selected_vars_exp_encode,"hic_contact",
                           "normalizedEP300_enhActivity", "activity_prom")  
# rE2G-extended.rf variables
encode_ext_model <- read.csv("InputData/K562/K562Vars.csv")
encode_ext_model <- encode_ext_model$Variable[encode_ext_model$rE2G_extended_rf]

# EGrf-extended models with one variable removed at a time
removed.vars <- list()
for (var in combined_model_encode) {
  removed.vars[[paste0(var, ".removed")]] <- combined_model_encode[combined_model_encode != var]
}
# Preparing a list of variables for all K562 models
ENCODE.vars.list <- append(list(EGrf = selected_vars_exp_encode,
                                TAPseq.rf = tap_seq_vars_encode,
                                EGrf.Extended = combined_model_encode,
                                rE2G.Extended.rf = encode_ext_model),
                                removed.vars)
#########################################################
# Running Astrocyte models
#########################################################
#Load Astrocytes Datasets (hg38)
all_rf_factors_powered <- read.csv("TrainingData/TrainingDataframe_Astrocytes.csv")

message("Running Astrocyte Models")

###all models are trained here and saved. 
###Javier Comment - random forest models are trained using the ranger package with hyperparameters tuned from the caret model training (functions defined in Comparison_Functions.R file)
all_rf_factors_powered <- trainAllModels(all_rf_factors_powered, save_models = T,
                                         variables.list = vars.list,suffix = "_Power_subset",
                                         probability = T, train_var = "HitPermissive_NegZ",
                                         folder="../3.Predictions/RF_Results/Astrocytes/")

write.csv(all_rf_factors_powered, "../3.Predictions/RF_Results/Astrocytes/Astrocyte_trainingData_pluspredictions.csv", row.names = F)

#Generating a precision-recall (PR) curve and finding the score (probability threshold) at which the precision is closest to 0.8
points<- pr.curve(scores.class0 = all_rf_factors_powered$EGrf,weights.class0 = all_rf_factors_powered$HitPermissive_NegZ, curve = T)
points <- points$curve

#0.8 precision = 0.33 recall and 0.506 probability
score_80ptprecision<- points[which.min(abs(points[,2] - 0.8)),3]
write.table(score_80ptprecision, file = "../3.Predictions/RF_Results/Astrocytes/EGrf80ptPrecisionCutoff.txt")

#########################################################
# Running K562 models
#########################################################
#K562 hg38
K562_EGPs <- read.csv("TrainingData/TrainingDataframe_K562.csv")
colnames(K562_EGPs)[colnames(K562_EGPs) == "Gene.Distance"] <- "Distance"

#Handling missing values
#The 37 missing EGPs are nearly all hits unfortunately (from Gasperini) 
#they are also removed by their code so I will have to drop them
#There are also missing values for 937 EGPs which have everything except distance set to 0 (or a default)
K562_EGPs <- K562_EGPs[K562_EGPs$ENCODE_overlap.Exists ,] #& ! K562_EGPs$Missing.Values
#Only one gene without expression so just removing it 
#K562_EGPs <- K562_EGPs[! is.na(K562_EGPs$Gene.RNAseq_RPKM),]
K562_EGPs[is.na(K562_EGPs$Gene.StabilityIndex), "Gene.StabilityIndex"] <- 0 #Assume if not in results 0 stability
#apply(K562_EGPs, 2, function(x) {(sum(is.na(x)))})

message("Running K562 models")

K562_EGPs <- trainAllModels(all_rf_factors = K562_EGPs, train_var = "Significant_NegEffectSize",
                            variables.list = ENCODE.vars.list, suffix = "_K562_EGPs",
                            probability = T,
                            folder = "../3.Predictions/RF_Results/K562")
write.csv(K562_EGPs, "../3.Predictions/RF_Results/K562/K562_trainingData_pluspredictions.csv", row.names = F)
