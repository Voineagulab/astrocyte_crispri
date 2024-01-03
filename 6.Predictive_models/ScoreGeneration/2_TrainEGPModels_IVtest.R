##################
##
## Prediction Comparison plots
## @author: Sam Bagot
## @date: 14-11-23
##
#################
rm(list=ls())
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/")
source("Scripts/Header_functions.R")
source("../Manuscript/Figs/FinalFigureFunctions.R")
source("Scripts/ScoreGeneration/Comparison_Functions.R")


#Load Astrocytes Datasets (hg38)
all_rf_factors <- read.csv("Data/all_rf_factors.csv")
#If not in the data currently setting to 0 but need another solution
#Considering these as missed predictions but also subset plots to intersecting only
all_rf_factors[! all_rf_factors$ENCODE_overlap.Exists, c("rE2G.DNAseOnly", "ABCScore_ENCODE")] <- 0
all_rf_factors[is.na(all_rf_factors$Gene.StabilityIndex), "Gene.StabilityIndex"] <- 0 #Assume if not in results 0 stability
all_rf_factors$ABC_Score <- all_rf_factors$ABCScore_ENCODE #Replace my run of ABC on our ATAC seq with the new ENCODE one
colnames(all_rf_factors)[colnames(all_rf_factors) == "Gene.Distance"] <- "Distance"

message("Running Voineagu Models")
removed.vars_astro <- list()
for (var in selected_vars_exp) {
  removed.vars_astro[[paste0(var, ".removed")]] <- selected_vars_exp[selected_vars_exp != var]
}
# for (var in all_vars_astro) {
#   removed.vars[[paste0("All.",var, ".removed")]] <- all_vars_astro[all_vars_astro != var]
# }

vars.list <- append(list(TAPseq.rf = tap_seq_vars,
                  EGrf = selected_vars_exp,
                  EGrf.PlusHiC = plus_hic_vars, 
                  EGrf.PlusTobias = c(selected_vars_exp,"Tobias.Exp_Unbound_TF_counts", "Tobias.Exp_Bound_TF_counts"), 
                  EGrf.PlusTTseq = c(selected_vars_exp,"TTseq.TTseq_Total", "TTseq.Ratio_TTversusRNA"),
                  EGrf.PlusBelugaDIS = c(selected_vars_exp,"Beluga.MaxDisScore"),
                  #EGrf.All = all_vars_astro,
                  EGrf.ReplaceATACwithDNAse = selected_vars_exp_replace_ATAC
                  ), removed.vars_astro)

power <- read.csv("../FullScale/Results/2_DE/Power/Power Simulation - Power per EGP.csv")
all_rf_factors_powered <- all_rf_factors[all_rf_factors$Pair %in% power[power$WellPowered015,]$Pair |
                                           all_rf_factors$HitPermissive_NegZ,]
all_rf_factors_powered <- trainAllModels(all_rf_factors_powered, save_models = T,
                                         variables.list = vars.list,suffix = "_Power_subset",
                                 probability = T, train_var = "HitPermissive_NegZ")
print(head(all_rf_factors_powered$EGrf, n = 20))
write.csv(all_rf_factors_powered, "Data/all_rf_factors_subset_pluspredictions.csv", row.names = F)


points<- pr.curve(scores.class0 = all_rf_factors_powered$EGrf,weights.class0 = all_rf_factors_powered$HitPermissive_NegZ, curve = T)
points <- points$curve
#0.8 precision = 0.33 recall and 0.506 probability
score_80ptprecision<- points[which.min(abs(points[,2] - 0.8)),3]
write.table(score_80ptprecision, file = "Results/ComparisonPlots/EGrf80ptPrecisionCutoff.txt")

#K562 hg38
K562_EGPs <- read.csv("Results/K562/ENCODE_EGP.csv")
colnames(K562_EGPs)[colnames(K562_EGPs) == "Gene.Distance"] <- "Distance"
#The 37 missing EGPs are nearly all hits unfortunately (from Gasperini) 
#they are also removed by their code so I will have to drop them
#There are also missing values for 937 EGPs which have everything except distance set to 0 (or a default)
K562_EGPs <- K562_EGPs[K562_EGPs$ENCODE_overlap.Exists ,] #& ! K562_EGPs$Missing.Values
#Only one gene without expression so just removing it 
K562_EGPs <- K562_EGPs[! is.na(K562_EGPs$Gene.RNAseq_RPKM),]
K562_EGPs[is.na(K562_EGPs$Gene.StabilityIndex), "Gene.StabilityIndex"] <- 0#Assume if not in results 0 stability
apply(K562_EGPs, 2, function(x) {(sum(is.na(x)))})

message("Running ENCODE models")
removed.vars <- list()
for (var in combined_model_encode) {
  removed.vars[[paste0(var, ".removed")]] <- combined_model_encode[combined_model_encode != var]
}
# for (var in all_vars_encode) {
#   removed.vars[[paste0("All.",var, ".removed")]] <- all_vars_encode[all_vars_encode != var]
# }
ENCODE.vars.list <- append(list(EGrf = selected_vars_exp_encode,
                         TAPseq.rf = tap_seq_vars_encode,
                         EGrf.Extended = combined_model_encode,
                         EGrf.Extended.PlusTobias = c(combined_model_encode,"Tobias.Exp_Unbound_TF_counts", "Tobias.Exp_Bound_TF_counts"), 
                         EGrf.Extended.PlusTTseq = c(combined_model_encode,"TTseq.TTseq_Total", "TTseq.Ratio_TTversusRNA"),
                         EGrf.Extended.ReplaceEnhActivitywithExp = replaceEnhactivitywithExp,
                         #EGrf.Extended.All = all_vars_encode,
                         rE2G.Extended.rf = encode_ext_model),removed.vars)
K562_EGPs <- trainAllModels(all_rf_factors = K562_EGPs, train_var = "Significant_NegEffectSize",
                            variables.list = ENCODE.vars.list, suffix = "_K562_EGPs",
                            probability = T )
write.csv(K562_EGPs, "Results/K562/ENCODE_EGP_plusPredictions.csv", row.names = F)



warning("This doesn not work on linux")
library("xlsx")
variables_astro <- getVarsExcel(vars.list[str_detect(names(vars.list), "EGrf|E2G|TAP")])
write.xlsx(variables_astro, sheetName = "AstroVars", "Results/ComparisonPlots/RFvars.xlsx", row.names = F)
variables_K562 <- getVarsExcel(ENCODE.vars.list[str_detect(names(ENCODE.vars.list), "EGrf|E2G|TAP")])
write.xlsx(variables_K562, sheetName = "K562Vars","Results/ComparisonPlots/RFvars.xlsx", append = T, row.names = F)
#variable descriptions are manually annotated to create the final file