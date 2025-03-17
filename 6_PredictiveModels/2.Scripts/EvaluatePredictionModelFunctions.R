#This script define the functions sourced in EvaluatePredictionModels.R
library(ggplot2)
library(ranger)
library(caret)
library(data.table)
library(PRROC)
library(stringr)
library(biomaRt)
library(xgboost)
library(tidyr)
library(rapportools)
library(parallel)
library(readxl)
library(biomaRt)
library(cowplot)
library(stringr)
library(seqinr)

##Set seed
seed = 135643
set.seed(seed)


# Function to create a confusion matrix and its visualization.
confusionM <- function(r, #Dataframe containing predictions and actual Hits.
                       abc = FALSE, # Boolean indicating if ABC method is used
                       thresholded = FALSE, # Boolean indicating if the data is already thresholded.
                       title = "", # Title for the confusion matrix plot
                       useZ = T) { # Boolean indicating if Z-score should be used
  #  # Applying different thresholds
  if (thresholded == FALSE & abc == TRUE) {
    r$above_threshold <- as.numeric(r$ABC.Exists)
    title <- "ABC Confusion Matrix"
  } else if (thresholded == FALSE & useZ == F) {
    r$above_threshold <- as.numeric(abs(r$prediction) > 0.4)
    title <- "RF Confusion Matrix"
  } else if (thresholded == FALSE){
    r$above_threshold <- as.numeric(abs(r$prediction) > 2)
    title <- "RF Confusion Matrix"
  } 
  # Converting the Hits column to numeric.
  r$Hits <- as.numeric(r$Hits)
  # Creating the confusion matrix 
  cm <- confusionMatrix(factor(r$above_threshold),factor(r$Hits), positive = "1")
  # Plotting the confusion matrix.
  plotCaretCM(cm, title)
  return(cf)
}

###############################################################
###Defining Color Palettes and Plot Styling  for Data Visualization
###############################################################
pals <- list()

# hits
pals$Hits <- c("#3A5F9A", "#F29D6C")
pals$Hits_Darker <- c("#2d4976", "#e47636")

# continuous
pals$grn2orng <- c("#b51a00", "#ee5900", "#ff9d68", "#feceb8", "#ffe2d6", "#d9daab", "#84c053", "#4c852e", "#003f01")

# sets of neutral colour
pals$One <- "#354a69" # dark-grey-blue
pals$Two <- c("#6F083D", "#8DCACE")

# main discrete palette
pals$Primary <- c("#fc9d77", "#3d8f80", "#2c2c54", "#f3ccde", "#8dcace", "#3a5f9a", "#6f083d", "#dd5d2b")
pals$Primary_Darker <- c("#E4784E", "#317266", "#1C1C35", "#DA97B6", "#43979D", "#2D4976", "#39041F", "#B64519")

## Themes
invis <- element_blank()
ts <- function(x) element_text(size = x)
theme_resizeText <- theme_bw() + theme(axis.text = ts(6), legend.text = ts(6),
                                       axis.title = ts(7))
standard_theme <- theme_bw() + theme(panel.grid = element_blank(), text = ts(6), strip.text = ts(7),
                                     axis.title = ts(7), axis.line = element_line(), legend.position = "bottom")
## Sizing
# a4
h_a4 <- 29.7 / 2.54
w_a4 <- 21.0 / 2.54

h_margin <- (29.7 - (2.54)) / 2.54
w_margin <- (21 - (2.54)) / 2.54

## Text
text90 <- element_text(angle = 90, hjust = 0.5, vjust = 1)


# Negative predictors used to get PR and ROCK AUC curves
negativePredictors <- c("ClosestGeneDist","Distance", "Gene.StabilityIndex", "numTSSEnhGene",
                        "Tobias.Exp_Unbound_TF_counts", "activity_prom")


##Setting colors to different variables
plot_cols <- c(EGrf = pals$Primary[8],EGrf.Extended = pals$Primary[7],
               TAPseq.rf = pals$Primary[4], 
               Distance = "#000000",
               ABC_Score = pals$Primary[5],rE2G.DNAseOnly = pals$Primary[2],
               rE2G.Extended = pals$Primary[6], #ENCODE.
               rE2G.Extended.rf = "grey"
) 


#Get summary of performance metrics (i.e. AUPRC, AUROC, PR70 confusion matrix statistics)
getAUCLists <- function (df, #Dataframe containing the dataset.
                         aucvars, #Variables for which AUC metrics are to be calculated.
                         Hit, #Column name indicating the target variable.
                         probability) { # Flag to indicate whether to use probabilities in calculations.
  #Initializing lists to store various metrics.
  prList <- list()
  rocList <- list()
  cms <- list()
  pr70s <- list()
  #Check for missing columns in the dataset and warn if any.
  if (sum(aucvars %in% colnames(df)) < length(aucvars)) warning("missing columns in vars list")
  aucvars <- aucvars[aucvars %in% colnames(df)]
  # Iterating over each variable to compute the AUC metrics.
  for (var in aucvars) {
    #Handling negative predictors by flipping the score.
    if (var %in%  negativePredictors) {
      score <- 1 / df[,var]
      ## Handling random forest scores when probability is not TRUE.
    } else if (str_detect(var, "rf") & probability != T) {
      score <- - df[,var]
    } else {
      score <- df[,var]
    }
    #Computing the Precision-Recall curve.
    pr <- pr.curve(scores.class0 = score, weights.class0 = df[,Hit] > 0, curve = TRUE, rand.compute = TRUE)
    # Finding the PR curve point at 70% recall.
    pr_70 <- pr$curve[which.min(abs(pr$curve[,1] - 0.7)),]
    # Creating a confusion matrix
    cms <- append(cms, list(confusionMatrix(as.factor(score > pr_70[3]), as.factor(df[,Hit] > 0), positive = "TRUE")))
    # Appending the PR curve data at 70% recall.
    pr70s <- append(pr70s, list(pr_70))
    #Storing the PR and ROC curves in their respective lists.
    prList <- append(prList, list(pr))
    rocList <- append(rocList, list(roc.curve(scores.class0 = score, weights.class0 = df[,Hit] > 0, curve = TRUE, rand.compute = TRUE)))
  }
  names(prList)  <- aucvars
  names(rocList)  <- aucvars
  names(pr70s) <- aucvars
  names(cms) <- aucvars
  return(list(AUPRC = prList, AUROC = rocList, PR70 = pr70s,CMs = cms))
}

###Rename variables
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


################
#Functions for generating plots
################
# Function to create and plot AUC (Area Under Curve) for PR (Precision-Recall) and ROC (Receiver Operating Characteristic) curves.
aucPlots <- function (prediction, #Predictions from the model.
                      Hits, #Actual hits data.
                      title = "", #Plots title
                      doplot = T, #Boolean indicating whether to plot
                      addRecallLine = F) { #Whether to add a recall line,
  # Check for all zero predictions.
  if (all(prediction == 0)) { warning("prediction values all 0")}
  # Calculating PR and ROC curves.
  pr <- pr.curve(scores.class0 = prediction, weights.class0 = Hits, curve = TRUE, rand.compute = TRUE)
  roc <- roc.curve(scores.class0 = prediction, weights.class0 = Hits, curve = TRUE, rand.compute = TRUE)
  # Creating titles for PR and ROC plots.
  pr.title = paste("PR:", signif(pr$auc.integral, 3),  title)
  roc.title = paste("ROC:", signif(roc$auc, 3),  title)
  # Plotting PR curve 
  p1 <- ggplot(as.data.frame(pr$curve), aes(V1, V2, colour = V3)) +
    geom_line(linewidth = 1) + labs(x = "Recall", y = "Precision",colour = "Predictor", title = pr.title) + 
    scale_y_continuous( limits = c(0,1)) + theme_bw() +
    scale_x_continuous( limits = c(0,1)) + geom_hline(yintercept = sum(Hits > 0) / length(Hits), colour = "black")
  if (addRecallLine) {
    p1 <- p1 + geom_vline(linetype = "dashed", xintercept = 0.7, colour = "darkred")
  }
  # Plotting ROC curve 
  p2 <- ggplot(as.data.frame(roc$curve), aes(V1, V2, colour = V3)) +
    geom_line(linewidth = 1) + labs(x = "FPR", y = "Recall",colour = "Predictor", title = roc.title) + 
    scale_y_continuous( limits = c(0,1)) + theme_bw() +
    scale_x_continuous( limits = c(0,1)) + geom_abline(slope = 1, intercept = 0, colour = "black")
  if (addRecallLine) {
    p2 <- p2 + geom_hline(linetype = "dashed",yintercept = 0.7, colour = "darkred")
  }
  print(p1)
  print(p2)
  aucs <- as.data.frame(rbind(pr$auc.integral, roc$auc) )
  rownames(aucs) <- c("PR_AUC", "ROC_AUC")
  aucs <- signif(aucs, digits = 3)
  return(aucs)
}

##Defining thresholds
DNAse.Threshold = 0.202140768
E2G.Extended.Threshold = 0.319162213


####Bind AUClist into a single dataframe
bindAUCsList <- function (aucsList, #List of AUC data
                          column= "curve") { #The specific column in the AUC data to be extracted
  aucsList <- do.call("rbind",lapply(names(aucsList),function(model) {
    return(data.frame(Model = model,aucsList[[model]][[column]] ))}   
  ))
  return(aucsList)
}

###Plot Curves using metrics stored in getAUCLists
allAUCurves <- function (aucsList,  #List of AUC data
                         plot_cols, #Colors to use for different models in the plot
                         x.lab, #label for the x axis
                         y.lab, # label for the y axis
                         add.ENCODErecall.point = F){ ## Whether to add ENCODE recall points to the plot.
  # Binding AUC data into a single data frame.
  aucsList <- bindAUCsList(aucsList)
  #Creating the AUC plot
  aucPlots <- ggplot(aucsList, aes(X1, X2, colour = Model)) +  
    geom_line(aes(linewidth = Model,alpha = Model )) + 
    standard_theme + labs(x = x.lab, y = y.lab) + 
    scale_y_continuous( limits = c(0,1)) + 
    scale_linewidth_manual(values = c("EGrf" = 1, "EGrf.Extended" = 1), 
                           na.value = 0.5, guide = "none") +
    scale_alpha_manual(values = c("EGrf" = 1, "EGrf.Extended" = 1), 
                       na.value = 0.5, guide = "none") +
    scale_x_continuous( limits = c(0,1)) +
    scale_colour_manual(values = plot_cols) #, guide = guide_legend(ncol =3)
  # Adding a vertical line for recall if x-axis is "Recall".
  if (x.lab == "Recall") {
    aucPlots <- aucPlots + geom_vline(xintercept = 0.7, linetype = "dashed", colour = pals$One)
  }
  #adding ENCODE recall points to the plot.
  if (add.ENCODErecall.point) {   
    DNAse.only <- aucsList[aucsList$Model == "rE2G.DNAseOnly", ]
    DNAse.point <- DNAse.only[which.min(abs(DNAse.only[, "X3"]- DNAse.Threshold)), ]
    aucPlots <- aucPlots + geom_point(data = DNAse.point, aes(X1, X2, colour = Model), size = 2.5, show.legend=FALSE)
    if ("rE2G.Extended" %in% aucsList$Model) {
      DNAse.only <- aucsList[aucsList$Model == "rE2G.Extended", ]
      DNAse.point <- DNAse.only[which.min(abs(DNAse.only[, "X3"]- E2G.Extended.Threshold)), ]
      aucPlots <- aucPlots + geom_point(data = DNAse.point, aes(X1, X2, colour = Model), size = 2.5, show.legend=FALSE)
    }
  }
  return(aucPlots)
}

###Create and output precision-recall (PR) and a receiver operating characteristic (ROC) plots
# outputPlot <- function (prPlot, #The precision-recall plot
#                         rocPlot, #the ROC plot
#                         plot.type) { # The type of plot to output ('combined', 'pr', or 'roc').
#   if (plot.type == "combined") {
#     # Combining PR and ROC plots with legends.
#     le <- get_legend(prPlot)
#     prPlot <- prPlot + theme(legend.position = "none")
#     rocPlot <- rocPlot + theme(legend.position = "none")
#     plots <- plot_grid(prPlot, rocPlot)
#     out <- plot_grid(plots, le, rel_heights = c(1,0.2), ncol = 1)
#   } else if (plot.type == "pr") {
#     out <- prPlot
#   } else if (plot.type == "roc") {
#     out <- rocPlot
#   }
#   print(out)
#   return(out)
# }


####
outputPlot <- function (prPlot, #The precision-recall plot
                        rocPlot, #the ROC plot
                        plot.type) { # The type of plot to output ('combined', 'pr', or 'roc').

  # Function to modify plot theme for increased label size
  modify_plot_theme <- function(plot) {
    plot + theme(
      axis.title = element_text(size = rel(1.5)),  # Adjust axis title size
      axis.text = element_text(size = rel(1.5)),   # Adjust axis text/labels size
      legend.title = element_text(size = rel(1)),# Adjust legend title size
      legend.text = element_text(size = rel(1))  # Adjust legend text size
    )
  }

  if (plot.type == "combined") {
    # Apply theme modifications to both plots
    prPlot <- modify_plot_theme(prPlot)
    rocPlot <- modify_plot_theme(rocPlot)

    # Combining PR and ROC plots with legends.
    le <- get_legend(prPlot)
    prPlot <- prPlot + theme(legend.position = "none")
    rocPlot <- rocPlot + theme(legend.position = "none")
    plots <- plot_grid(prPlot, rocPlot)
    out <- plot_grid(plots, le, rel_heights = c(1,0.2), ncol = 1)
  } else if (plot.type == "pr") {
    out <- modify_plot_theme(prPlot)
  } else if (plot.type == "roc") {
    out <- modify_plot_theme(rocPlot)
  }

  print(out)
  return(out)
}


#### Performs comparison between multiple models and generates corresponding plots
aucComparison <- function(df, ###Dataframe containing the data
                          Hit = "HitPermissive_NegZ", ##Column indicating hits
                          plot_cols, #Color mapping for the plot
                          probability = T, #Whether to use probability
                          plot.type = "combined", # The type of plot ('combined', 'roc', 'pr').
                          add.recall.point = F) { #Whether to add recall points.
  # Validating plot type.
  plottypes <- c("combined", "roc", "pr")
  if (! plot.type %in% plottypes) {stop(paste(c("select plot type from: ", plottypes),collapse = ","))}
  #Filtering plot colors and variables based on the dataframe.
  plot_cols <- plot_cols[names(plot_cols) %in% colnames(df)] 
  aucvars <- names(plot_cols)[names(plot_cols)  %in% colnames(df)]
  #Get summary of performance metrics (i.e. AUPRC, AUROC, PR70 confusion matrix statistics) for the different models 
  aucLists <- getAUCLists(df, aucvars, Hit, probability)
  #Creating PR and ROC plots.
  prPlot <- allAUCurves(aucLists[[1]], plot_cols,"Recall","Precision", add.recall.point) 
  rocPlot <- allAUCurves(aucLists[[2]], plot_cols,"FPR","Recall")  
  #Outputting the plot based on specified type.
  plot<- outputPlot(prPlot, rocPlot, plot.type)
  if(plot.type != "combined") {
    return(plot)
  } else {
    return(aucLists)
  }
}

##Plot subsets of Precision-Recall plots for the different models 
plotSubsetsAUCs <- function (df, #dataframe with data
                             hit = "HitPermissive_NegZ", #column containing hits 
                             plot_cols, #color mapping for plots
                             probability= T, #Whether to use probability
                             dataset.name = "") { #Name of the dataset
  # Creating a directory for the dataset.
  dir.create(dataset.name, showWarnings = F)
  # Saving the overall performance plot
  pdf(paste0(dataset.name,"/Performance_",dataset.name,".pdf"), height = w_margin / 2, width = w_margin)
  aucLists <-  aucComparison(df,Hit = hit,  plot_cols = plot_cols, probability = probability)
  dev.off()
  # Generating and saving plots for subsets of data based on model categories.
  # Excluding EGrf models.
  pdf(paste0(dataset.name,"/Performance_",dataset.name,"_NOEGRF.pdf"), height = w_margin / 2, width = w_margin)
  aucComparison(df, Hit = hit,
                plot_cols = plot_cols[str_detect(names(plot_cols),"EGrf" ,negate = T)], 
                probability = probability)
  dev.off()
  # Including only EGrf models.
  pdf(paste0(dataset.name,"/Performance_",dataset.name,"_EGRF_only.pdf"), height = w_margin / 2, width = w_margin)
  aucComparison(df, Hit = hit,
                plot_cols = plot_cols[str_detect(names(plot_cols),"EGrf")], 
                probability = probability)
  dev.off()
  # Including only ENCODE models.
  pdf(paste0(dataset.name,"/Performance_",dataset.name,"_ENCODE_only.pdf"), height = w_margin / 2, width = w_margin)
  aucComparison(df, Hit = hit,
                plot_cols = plot_cols[str_detect(names(plot_cols),"rE2G|ABC")], 
                probability = probability)
  dev.off()
  # Aggregating performance metrics at 70% recall level and saving.
  PR70 <- as.data.frame(do.call("rbind",aucLists[[3]]))
  colnames(PR70) <- c("Recall_70", "Precision_70Recall", "ModelThreshold_70Recall")
  PR70$F1Score_70Recall <- (2 * PR70$Recall_70 * PR70$Precision_70Recall) / (PR70$Recall_70 + PR70$Precision_70Recall)
  PR70$AUPRC <- do.call("rbind",lapply(aucLists[[1]], function(x) {return(x$auc.integral) }))
  PR70$AUROC <- do.call("rbind",lapply(aucLists[[2]], function(x) {return(x$auc) }))
  write.csv(PR70, paste0(dataset.name,"/Recall70.csv"))
  saveRDS(aucLists[[4]], paste0(dataset.name,"/CMsRecall70.rds"))
  return(aucLists)
}


###################
#Bootstrapping plots
###################
#Generate bar plots for PR curves
barplot.plot <- function(prcs, #Dataframe containing precision-recall curve statistics.
                         plot_cols = plot_cols, ##Color mapping for the plot
                         x.title = "AUPRC", #Title for x axis
                         y.column= "name") { #Column name to be used for y-axis

  # Adjusting p-values for multiple testing using FDR and marking significant results.
  prcs$p.value <-  p.adjust(prcs$p.value,"fdr")
  prcs$significant <- ""
  if (nrow(prcs[prcs$p.value < 0.05,]) > 0) {
    prcs[prcs$p.value < 0.05,]$significant <- " *"
  }
  # Constructing names for the bar plot.
  #prcs$name <- paste0(prcs$model, "\n","p=",signif(prcs$p.value,3),"", prcs$significant) # display or not pvalue in y axis
  prcs$name <- paste0(prcs$model)
  prcs$name <- factor(prcs$name,levels = unique(prcs[order(prcs$means),"name"]))
  prcs$model <- factor(prcs$model,levels = unique(prcs[order(prcs$means),"model"]))
  # Adding additional colors for specific models.
  addcols <- prcs$model[str_detect(prcs$model, "Replace|Plus")]
  newcols <- rep(pals$Primary[1], length(addcols))
  names(newcols) <- addcols
  # Creating the bar plot using ggplot2.
  plot <- ggplot(prcs, aes(means, .data[[y.column]], fill = model)) +
    geom_bar(stat = "identity")  + labs(x = x.title, y = "Model") +
    scale_fill_manual(values = c(plot_cols,newcols), guide = "none") +
    geom_errorbar(aes(xmin = lower.bound, xmax = upper.bound), width = 0.2) +
    standard_theme + theme(legend.position = "none")# +
  return(plot)
}

##calculates p-values and delta values for bootstrapped results
getP_and_delta <- function (bindRes, #Bound results from bootstrap analysis.
                            means, #Mean values of the bootstrap results.
                            comparison.col) { #The column used for comparison in the bootstrap analysis.
  # Calculating p-values for greater and lesser comparisons.
  p.greater <- apply(bindRes, 2, function(x) { (1+sum(x >= bindRes[,comparison.col])) / (length(x)+1) })
  p.lesser <- apply(bindRes, 2, function(x) { (1+sum(x <= bindRes[,comparison.col])) / (length(x)+1) })
  # Calculating delta values.
  delta <-  means - median(bindRes[,comparison.col])
  # Forming a dataframe with calculated values.
  df <- data.frame(p.greater,  p.lesser)
  df$p.twosided <- apply(df,1, min) *2 #*2 since a two sided test?
  df$delta <- delta
  return(df)
} 

## Function to summarize bootstrap results.
bootstrap.sum <- function(bootstrapAUCs, #List of bootstrap AUCs.
                          variables, #Variables involved in the bootstrap analysis.
                          p.value_col= "EGrf") { #The column used for p-value calculation.
  #Binding bootstrap AUC results into a single dataframe.
  bindRes <- do.call("rbind",bootstrapAUCs)
  # Calculating median, lower, and upper bounds.
  means <- apply(bindRes, 2, median)
  lower <- apply(bindRes, 2, function(x){sort(x)[ceiling(length(x)* 0.025)]})
  upper <- apply(bindRes, 2, function(x){sort(x)[floor(length(x) * 0.975)]})
  # Getting p-values and delta values.
  pd <- getP_and_delta(bindRes,means,  p.value_col) 
  # Combining results into a dataframe.
  out <- data.frame(model = names(means), means, lower.bound = lower, #means - margin, 
                    upper.bound = upper, #means + margin, 
                    delta = pd$delta,
                    p.value = pd$p.twosided)
  return(out)
}

#Function to perform bootstrap analysis for AUC metrics.
bootstrap.AUCs <- function (df, #Dataframe containing the dataset.
                            variables, # Variables used in the AUC calculation.
                            bootstraps = 1000, #Number of bootstrap iterations
                            hit = "HitPermissive_NegZ", #Column name indicating the hit 
                            probability = T, #whether to use probability
                            p.value_comp) { #Column used for p-value comparison.
  set.seed(seed) ## Setting seed for reproducibility
  bootstrapPRAUCs <- list()
  bootstrapROCs <- list()
  # Looping through the number of bootstraps to calculate AUCs.
  for(i in 1:bootstraps) {
    smpl <- sample(1:nrow(df), size = nrow(df), replace = T)
    aucLists <- getAUCLists(df[smpl,], variables, Hit = hit, probability = probability)
    bootstrapPRAUCs[[i]] <- unlist(lapply(aucLists[[1]], function(x) {return(x$auc.integral) }))
    bootstrapROCs[[i]] <- unlist(lapply(aucLists[[2]], function(x) {return(x$auc) }))
  }
  # Summarizing the bootstrap results for PR and ROC AUCs.
  bootstrapPRAUCs <- bootstrap.sum(bootstrapPRAUCs, variables,  p.value_comp)
  bootstrapROCs <- bootstrap.sum(bootstrapROCs, variables,  p.value_comp)
  
  return(list(PR = bootstrapPRAUCs,ROC = bootstrapROCs))
}


# Function to create barplots showing AUPRC in rf models adding variables one at a time, 
comparison.barplot <- function(df, #Dataframe containing the dataset
                               variables, # Variables used in the comparison.
                               bootstraps = 1000, #Number of bootstrap iterations
                               p.value_comp = "EGrf", #Column used for p-value comparison
                               hit = "HitPermissive_NegZ", #Column name indicating the hit
                               probability = T, #Boolean indicating whether to use probability
                               plot_cols, #Color mapping for plots.
                               boot.summary = NULL, #Precomputed bootstrap summary
                               dataset.name) { # Name of the dataset for file naming.
  # Calculating bootstrap summary if not provided.
  if (is.null(boot.summary)){
    boot.summary <- bootstrap.AUCs(df, variables, bootstraps = bootstraps, 
                                   hit, probability, p.value_comp = p.value_comp)
  }
  # Filtering out results from the summaries.
  PR.sum <- boot.summary[["PR"]][str_detect(boot.summary[["PR"]]$model, "All", negate = T),]
  ROC.sum <- boot.summary[["ROC"]][str_detect(boot.summary[["ROC"]]$model, "All", negate = T),]
  # Preparing model colors for the plots.
  model.cols <- plot_cols[names(plot_cols) %in% c("rE2G.DNAseOnly","rf.TAPseq", "rf.rE2G.EXT", 
                                                  "TAPseq.rf", "rE2G.Extended.rf",
                                                  "rE2G.Extended", "EGrf", "EGrf.Extended")]
  # Generating bar plots for AUROC and AUPRC.
  roc.p1 <- barplot.plot(ROC.sum[ROC.sum$model %in% names(model.cols),],model.cols,  x.title = "AUROC")
  pr.p1 <- barplot.plot(PR.sum[PR.sum$model %in% names(model.cols),], model.cols)
  pdf(paste0(dataset.name,"/Bootstrapping_",dataset.name,"_AUCs.pdf"), width = w_margin / 3, height = h_margin *1 / 6)
  print(pr.p1)
  dev.off()
  pdf(paste0(dataset.name,"/Bootstrapping_",dataset.name,".pdf"), width = w_margin / 3, height = h_margin *2 / 10)
  print(barplot.plot(PR.sum[! PR.sum$model %in% names(plot_cols[str_detect(names(plot_cols),"EGrf", negate = T)]),], 
                     plot_cols))
  dev.off()
  pdf(paste0(dataset.name,"/Bootstrapping_",dataset.name,"Plusonly.pdf"), width = w_margin / 3, height = h_margin *1 / 6)
  print(barplot.plot(PR.sum[str_detect(PR.sum$model,"EGrf"),], plot_cols))
  dev.off()
  pdf(paste0(dataset.name,"/Bootstrapping_",dataset.name,"ROCs.pdf"), width = w_margin / 3, height = w_margin *4 / 10)
  print(barplot.plot(ROC.sum[! ROC.sum$model %in% names(plot_cols[str_detect(names(plot_cols),"EGrf", negate = T)]),], 
                     plot_cols,x.title = "ROC"))
  dev.off()
  write.csv(boot.summary, paste0(dataset.name, "/BootstrappingResults.csv"), row.names = F)
  return(boot.summary)
}


# Function to plot barplots of median AUPRC comparing performance between models in Astrocytes and K562 cell
benchmarkK562Astro <- function (bootstrap_Astro, #Bootstrap results for Astrocytes
                                bootstrap_K562, #Bootstrap results for K562
                                cols, #Color mapping for plots
                                y.column = "model") { #column name for  y-axis
  
  # Filtering models based on the provided colors.
  bootstrap_Astro <- bootstrap_Astro[bootstrap_Astro$model %in% names(cols),]
  bootstrap_K562 <- bootstrap_K562[bootstrap_K562$model %in% names(cols),]
  # Adding cell type information.
  bootstrap_Astro$CellType <- "Astro"
  bootstrap_K562$CellType <- "K562"
  #  # Combining data from both cell types.
  combined <- rbind(bootstrap_Astro, bootstrap_K562)
  #Generating the bar plot with facet for each cell type.
  plot <- barplot.plot(combined, plot_cols = cols, y.column = y.column) +
    scale_x_continuous(limits = c(0,0.85))+ 
    facet_wrap(vars(CellType), scales = "free",ncol = 1) +
    theme(
      axis.title = element_text(size = rel(1.5)),  # Adjust axis title size
      axis.text.x =  element_text(size = rel(1.5)),   # Adjust axis text/labels size
      axis.text.y =  element_text(size = rel(1.5)),
      strip.text  =  element_text(size = rel(1.5)),
      legend.title = element_text(size = rel(1)),# Adjust legend title size
      legend.text = element_text(size = rel(1))  # Adjust legend text size
    )
  return(plot)
}


# Function to plotting how median AUPRC values change when individual variables are removed from a model, in relation to their permutation-based variable importance scores

ImportancePlots <- function(boostrap_PR, # Bootstrap results for precision-recall.
                            colours) { #Color mapping for different variables.
  #Generating the plot using ggplot2
  ggplot(boostrap_PR, aes(`Permutation Variable Importance`, means, colour= Variable)) +
    geom_errorbar(aes(ymin = lower.bound, ymax = upper.bound)) +      
    geom_point() + scale_shape_manual(values = c(16,15)) +
    scale_colour_manual(values = colours) +
    geom_hline(aes(yintercept = prline), colour = pals$One, linetype= "dashed") +
    standard_theme + labs(y = "AUPRC removing variable") + 
    theme(legend.position = "right",
          axis.title = element_text(size = rel(2)),  # Adjust axis title size
          axis.text.x =  element_text(size = rel(2)),   # Adjust axis text/labels size
          axis.text.y =  element_text(size = rel(2)),
          #strip.text  =  element_text(size = rel(1.5)),
          legend.title = element_text(size = rel(2)),# Adjust legend title size
          legend.text = element_text(size = rel(2))) 
}


### Function to combine bootstrapped PR data with variable importance scores.
combineImportance <- function(boostrap_PR, ##Dataframe containing bootstrapped PR data
                              importance, #Dataframe containing variable importance scores.
                              colours, #Colour mapping for the plot.
                              line.var = "EGrf", #The baseline model for comparison
                              CellType) { #The cell type to be considered 
  colnames(importance) <- c("Variable", "Permutation Variable Importance")
  boostrap_PR <- boostrap_PR[boostrap_PR$model != "Distance", ]
  boostrap_PR$Variable <- sub(".removed","", boostrap_PR$model)
  boostrap_PR$Basemodel <- line.var
  boostrap_PR$prline <- boostrap_PR[boostrap_PR$model == line.var,"means"]
  boostrap_PR <- merge(boostrap_PR,importance) #removes rows that aren't variables
  boostrap_PR$CellType <- CellType
  boostrap_PR$Variable <- updateVarNames(boostrap_PR$Variable)
  # Plotting the combined data using ImportancePlots function
  print(ImportancePlots(boostrap_PR, colours))
  return(boostrap_PR)
}


# Function to plot the AUPRC for individual variables to evaluate performance improvement.
aucIndividualVariables <- function(df, #Dataframe containing the dataset.
                                   aucvars, #: Variables to be used for AUC calculation.
                                   cols, #Color mapping for the plots.
                                   Hit, #Column name indicating the hit status.
                                   probability = T) { #Boolean indicating whether to use probability
  #Get summary of performance metrics (i.e. AUPRC, AUROC, PR70 confusion matrix statistics) for the different models 
  aucLists <- getAUCLists(df, aucvars, Hit, probability = probability)
  # Binding AUC lists and selecting the 'auc.integral' column.
  AUCs <- bindAUCsList(aucLists[[1]], "auc.integral")
  colnames(AUCs)[2] <- "AUPRC"
  AUCs$Model <- updateVarNames(AUCs$Model)
  AUCs$Model <- factor(AUCs$Model,AUCs[order(AUCs$AUPRC),"Model"])
  ggplot(AUCs, aes(AUPRC, Model, fill = Model)) +
    geom_bar(stat = "identity")  + 
    scale_fill_manual(values = cols, na.value = pals$One) +
    standard_theme + theme(legend.position = "none") +
    scale_x_continuous(limits = c(0,0.5)) +
    geom_vline(xintercept = sum(df[,Hit] > 0) / nrow(df)) +
    theme(
          axis.title = element_text(size = rel(1.5)),  # Adjust axis title size
          axis.text.x =  element_text(size = rel(1.5)),   # Adjust axis text/labels size
          axis.text.y =  element_text(size = rel(1.5)),
          #strip.text  =  element_text(size = rel(1.5)),
          legend.title = element_text(size = rel(1.5)),# Adjust legend title size
          legend.text = element_text(size = rel(1.5)))
}
