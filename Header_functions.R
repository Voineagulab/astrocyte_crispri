#Header File 
library(ggplot2)
library(ranger)
library(caret)
library(data.table)
#library(pROC)
library(PRROC)
library(stringr)
library("biomaRt")
library("xgboost")
library("tidyr")
library("rapportools")
library("parallel")
library(readxl)
#install.packages("xgboost")

seed = 135643
set.seed(seed)

#This is for running TT.tests & Fisher tests on the entire data.frame
removed_variables <- c("X", "Enh_Gene","ENSG.targetgene","HitCore", "HitCategory", "Fasta", "Enh",  "Enh.Pos", "Pair", "Gene", "Gene.Pos", "Gene.TSS",
                       "FDR.N50", "FDR.SCEPTRE", "FDR.N250", "P.N250", "P.N50", "P.NB", "P.SCEPTRE", "X", "X.1",# 
                       "nCells", "logfc", "logfc.vst", "Enh.chr", "Enh.start", "Enh.end", "Gene.chr", "Gene.start", "Tobias.Intersecting_TFs",
                       "Gene.end", "Enh.Pos_Gene", "PeakId", "Gene.Distance.Bin", "Exp_Intersecting_TFs", "Intersecting_TFs", "start", "end", "chr", 'TSS', 'EnsID', "Hnisz_SuperEnhancers_OLD")
#This is for making an Enhancer level data.frame (remove these columns then take unique)
gene_dependant <- c("X", "Pair", "Gene", "Gene.Pos", "Gene.TSS", "Gene.Distance","Distance", "Gene.Distance.Bin", "Gene.Nearest", "Gene.Upstream", "Gene.Exp", "Z", "logfc", "logfc.vst",
                    "P.SCEPTRE", "P.NB", "P.N50", "P.N250", "FDR.SCEPTRE", "FDR.N50", "FDR.N250", "HitPermissive","HitPermissive_NegZ", "HitCore", "HitCategory", "ABC_Score", "ABC.Exists", 
                    "Tobias.Pair_TF_overlaps","Tobias.Prop_shared_TFs", "Tobias.Intersecting_TFs_Exp_Enh_TF_in_Gene", "Tobias.Exp_Sum_Gene_TFs_in_Enh_TFs", "Tobias.Exp_Pair_TF_overlaps", "Tobias.Exp_Prop_shared_TFs", 
                    "Tobias.Exp_Intersecting_TFs","Tobias.Intersecting_TFs", "Tobias.Exp_Sum_Enh_TFs_in_Gene_TFs", "Tobias.Sum_Gene_TFs_in_Enh_TFs", "Tobias.Sum_Enh_TFs_in_Gene_TFs",  "Tobias.Bound_ratio",
                    "Gene.chr", "Gene.start", "Gene.end", "Gene.size", "Enh.Pos_Gene", "WithinTAD", "ABC_Score_OldChip", "ABC_OldChip.Exists",
                    "ClosestGeneExp", "MaxExpression", "MaxExpressionLT50kb", "NearestTested","ClosestGeneDist", "NGenesTested", "LT50KB",
                    "ABC_Score", "ABC.Exist", "ABC_Score_ATAC", "ABC_ATAC.Exists", "ABC_Score_CHIP", "ABC_CHIP.Exists", "ABC_Score_Nasser", "ABC_Nasser.Exist", "ABC_Approx", 
                    "HiC_interaction_astrocyte_cerebellum", "HiC_interaction_astrocyte_spinal", "HiC_interaction_NeuN_pos", "HiC_interaction_NeuN_neg",
                    "Gene.RNAseq_RPKM","Gene.TTseq_Counts", "Gene.Ref_RNAseq_RPKM", "Gene.RNAseq_CPM", "Gene.Ref_RNAseq_CPM", "Gene.RNAseq_Counts", "Gene.Ref_RNAseq_Counts",
                    "Gene.TTseq_RPKM", "Gene.TTseq_CPM")

#Original abc variables
ABC_vars <- c("ABC.Exists", "ABC_Score", "NHA_score", "ABC_Expressed")

enhancer_removed <- c("Enh", "Enh.Pos", "Enh.chr", "Enh.start", "PeakId", "Enh.end")
norm <- function(x){(x-min(x))/(max(x)-min(x))}

#Unused function to get training weights (increase sampling of hits)
getweights <- function(t, Hits = NA) {
  if (is.na(Hits)) {
    stop("You must set Hits") 
  }
  t <- as.data.frame(t)
  if(is.null(t[,])) { stop("Hits Not present") }
  hit_ratio <- sum(t[,Hits] > 0) / nrow(t)
  t[,"weights"] <- NA
  t[t[,Hits] == 0,"weights"] <- hit_ratio
  t[t[,Hits] > 0,"weights"] <- 1 - hit_ratio
  weights <- t[,"weights"]
  if(is.na(sum(weights))){stop("invalid values in weights")}
  return(weights)
}

#This is only used in archived files
randomiseSets <- function(rfd, train_size = .8, Hits = NA, create = F) {
  if (is.na(Hits)) {
    stop("You must set Hits") 
  }
  no_hits <<- length(which(rfd[,Hits] > 0))
  no_nonhits <<- length(which(rfd[,Hits] == 0))
  train_hits_rows <<- sample(which(rfd[,Hits] > 0), no_hits * train_size)
  train_nonhits <<- sample(which(rfd[,Hits] == 0), no_nonhits * train_size)
  if (create == T) {
    train_set <<- rfd[c(train_hits_rows,train_nonhits ),]
    test_set <<- rfd[! 1:nrow(rfd) %in% c(train_hits_rows,train_nonhits ),]
  }
}

#This is only used in archived files
plotCaretCM <- function (cm, title = "") {
  title <- paste(title,"Prec:", signif(cm$byClass["Precision" ],3),"Recall:", signif(cm$byClass["Recall"],3), "Balanced Acc:", signif(cm$byClass["Balanced Accuracy"], 3) )
  plot <- ggplot(as.data.frame(cm$table), aes(Reference, Prediction,  fill = Freq, label = Freq)) +
          geom_tile() +
          geom_label() +
          theme_bw() + 
          scale_fill_gradient(low="white", high="darkgreen") +
          scale_y_discrete(limits = rev) +
          theme(panel.border = element_blank(), 
                panel.grid.major = element_blank()) +
          labs(title = title)
  print(plot)
  return(plot)
}

#This is only used in comparison plots and archived files
confusionM <- function(r, abc = FALSE, thresholded = FALSE, title = "", useZ = T) {
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
  r$Hits <- as.numeric(r$Hits)
  cm <- confusionMatrix(factor(r$above_threshold),factor(r$Hits), positive = "1")
  plotCaretCM(cm, title)
  return(cf)
}


#This is only used in archived files
#R is Data.frame with cols $Hits, $prediction, $
statsRF <- function (r, pdf.title = NULL, plot.title = "Ranger") {
  if (! "Z" %in% colnames(r)) {
    warning("Skipping StatsRF (no Z value in data)")
    return(r)
  }
  r$Pair <- rownames(r)
  wt <- wilcox.test(r[r$Hits == TRUE, "prediction"], r[r$Hits == FALSE, "prediction"])
  my_lm <- lm(r$prediction ~ r$Z)
  my_auc <- auc( as.factor(r$Hits), as.numeric(abs(r$prediction)) )
  if (! is.null(pdf.title)) {
    pdf(file = paste0("results/Ranger_Implementation/", pdf.title ,".pdf"))
  }
  plot.title <- paste("AUC:", my_auc,"Wilcox",signif(wt$p.value, 3), "lm:grad",signif(my_lm$coefficients[2], 3), plot.title)
  print(ggplot(r, aes(Z, prediction,  colour = as.factor(Hits))) + # , label = Pair , label = Pair
          geom_point() + 
          stat_smooth(method = "lm", formula = y~ x)  +
          labs(title = plot.title)) +
          theme_classic()
  if (! is.null(pdf.title))  {
    dev.off()
  } 
  return(r)
}


#points(pr$curve[signif(pr$curve[,1],2) == 0.75,], col = "black", pch = 19, lwd = 5 ),
#text(1,1,txt, col = "black")
aucPlots <- function (prediction, Hits, title = "", doplot = T, addRecallLine = F) {
  if (all(prediction == 0)) { warning("prediction values all 0")}
  pr <- pr.curve(scores.class0 = prediction, weights.class0 = Hits, curve = TRUE, rand.compute = TRUE)
  roc <- roc.curve(scores.class0 = prediction, weights.class0 = Hits, curve = TRUE, rand.compute = TRUE)
  pr.title = paste("PR:", signif(pr$auc.integral, 3),  title)
  roc.title = paste("ROC:", signif(roc$auc, 3),  title)
  
  p1 <- ggplot(as.data.frame(pr$curve), aes(V1, V2, colour = V3)) +
          geom_line(linewidth = 1) + labs(x = "Recall", y = "Precision",colour = "Predictor", title = pr.title) + 
          scale_y_continuous( limits = c(0,1)) + theme_bw() +
          scale_x_continuous( limits = c(0,1)) + geom_hline(yintercept = sum(Hits > 0) / length(Hits), colour = "black")
  if (addRecallLine) {
    p1 <- p1 + geom_vline(linetype = "dashed", xintercept = 0.7, colour = "darkred")
  }
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


pred_cols <- c("NHA_chr", "NHA_start", "NHA_end", "NHA_id", "NHA_score", "NHA_strand","ABC_chr", "ABC_start", "ABC_end", "ABC_id", "Gene", "Score", "CellType", "ABC_Score")
pred_non_cols <- c("NHA_chr", "NHA_start", "NHA_end", "NHA_id", "NHA_score", "NHA_strand", "ABC_chr", "ABC_start", "ABC_end", "ABC_id",	"class", "activity_base",	"Gene",	"TargetGeneTSS",	"TargetGeneExpression",	
                   "TargetGenePromoterActivityQuantile",	"TargetGeneIsExpressed",	"distance",	"isSelfPromoter",	"powerlaw_contact",	"powerlaw_contact_reference",	"hic_contact",	"hic_contact_pl_scaled",	"hic_pseudocount",	"hic_contact_pl_scaled_adj",	"ABC.Score.Numerator",	"ABC_Score",	"powerlaw.Score.Numerator",	"powerlaw.Score",	"CellType")
#read prediction files 
read_pred <- function (file, colnames) {
  pred <- fread(file)
  colnames(pred) <- colnames
  pred$Enh.Pos_Gene <- paste0( pred$NHA_chr, ":",  pred$NHA_start, "-",  pred$NHA_end, "_" , pred$Gene)
  pred$Enh.Pos <- sub("_.*", "", pred$Enh.Pos_Gene)
  return(pred)
}

#I currently use this but maybe shouldn't
getENSID <- function(symbols, attributes = c("ensembl_gene_id","hgnc_symbol"), filters = "hgnc_symbol", aggregate = TRUE) {
  mart <- useMart("ensembl", "hsapiens_gene_ensembl")
  #datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  ENSB <- getBM(attributes=attributes, filters=filters,values = symbols,
                mart=mart, uniqueRows=T)
  #one has no ENSB ID because it is a fusion protein
  if (aggregate == TRUE) {
  ENSB <- aggregate(ENSB$ensembl_gene_id, by = list(ENSB[,filters]), FUN = function (x) { paste0(x, collapse = ",") } )
  }
  
  return(ENSB)
}

#Gavin's method
run.EnhancerFisher <- function(data, data.column, Hit = "HitPermissive") {
  total <- nrow(data)
  logical_data <- data[,data.column] %>% as.logical()
  Hits <- data[,Hit]
  # run stats
  ft <- table(logical_data, Hits) %>% fisher.test()
  # output stats
  out <- data.frame(Variable = data.column,
                    Total_TRUE = sum(Hits),
                    fraction_total_TRUE = sum(Hits) / total,
                    Total_Hit_TRUE = sum(Hits & logical_data),
                    Fraction_Hit_TRUE = sum(Hits & logical_data) / sum(logical_data),
                    p = ft$p.value,
                    OR = ft$estimate,
                    Lower = ft$conf.int[1],
                    Upper = ft$conf.int[2],
                    row.names = data.column)
  return(out)
}

#My old method
old_testStats <- function(data, print = F, Hit = "HitPermissive") {
  test_stats <- data.frame(col = character(), p.value = numeric(), odds.ratio = numeric(), method = character())
  for (col in colnames(data)) {
    if (str_detect(col,"ABC") | col %in% c("Enh")) { #Hit|
      message(paste("Skipping", col))
      next
    } else if (is.boolean(data[1,col])) {
      if (length(unique(data[,c(col)])) == 1) { message(paste(col," is constant")); next}
      ft <- fisher.test(table(data[,c(col)], data[,Hit] > 0))
      test_stats <- rbind(test_stats, c(col, ft$p.value, ft$estimate, "FT (Odds ratio)"))
      pval <- ft$p.value
    } else {
      tt <- t.test(data[data[,Hit] > 0,col], data[data[,Hit] == 0,col])
      test_stats <- rbind(test_stats, c(col, tt$p.value, tt$statistic, "TT (T stat)"))
      pval <- tt$p.value
    }
    if (print == T) {
      print(ggplot(data, aes(data[,col], colour = data[,Hit] > 0)) + geom_density() + labs(title = paste(col, "pval =", signif(pval, 3)), colour = Hit ))
    }
  }
  colnames(test_stats) <- c("Variable", "P", "T-Stat_or_Odds_ratios", "Method")
  test_stats$P <- as.numeric(test_stats$P)
  test_stats$FDR <- p.adjust(test_stats$P, method = "fdr")
  test_stats$BON <- p.adjust(test_stats$P, method = "bonferroni")
  test_stats <- test_stats[order(test_stats$P),]
  return(test_stats)

}
#get test stats for all cols vs hit permissive
getTestStats <- function(data, print = F, Hit = "HitPermissive", FT.only = F) {
  if (FT.only == T) {
    #returning different object
    test_stats <- lapply(colnames(data),function(column) {run.EnhancerFisher(data, column, Hit) })
  } else {
    test_stats <- old_testStats(data, print, Hit)
  }
  return(test_stats)
}

#change names from bigwigs 
accession2Name <- function (mylist) {
  accessionNames <- list(
    "ENCFF499UDS" = "DNASE.astrocyte_hippocampus", 
   "ENCFF901UBX" = "DNASE.astrocyte_spinal_cord",
  "ENCFF382FZE" = "DNASE.astrocyte_cerebellum",
  "ENCFF320CHE" = "DNAse",
  "ENCFF033TTE" = "CHIP.CTCF.astrocyte",
  "ENCFF058GEM" = "Chip.H3K27ac",
  "ENCFF184NZS" =  "Chip.H3K4me3", 
  "ENCFF656QVF" = "Chip.H3K9ac",
  "ENCFF736AZD" = "Chip.H4K20me1",
  "ENCFF424JNY" = "CHIP.CTCF.astrocyte_spinal_cord",
  "ENCFF510YMH" = "CHIP.H3K4me3.astrocyte_spinal_cord",
  "ENCFF072YMW" = "CHIP.H3K4me3.astrocyte_cerebellum",
  "ENCFF757YRI" =  "CHIP.CTCF.astrocyte_cerebellum",
  "ENCFF476HBN" = "Chip.H3K27me3",
  "ENCFF242JDN" = "Chip.H3K4me1",
  "ENCFF702AYE" = "Chip.POLR2A",
  "ENCFF640EWX" = "CHIP.H3K4me3.K562"
  )
  outList <- mylist
  names(outList) <- mylist
  outList[outList[outList %in% names(accessionNames)]] <- accessionNames[outList[outList %in% names(accessionNames)]]
  return(outList)
}

#Tobias function
TFmatrix <- function (overlaps, outnumeric = FALSE, expressed_genes = NULL, TF_col = "V8", target_col = "V4") {
  warning("This function was recently changed and results may be bugged")
  overlaps[,TF_col] <- sub(".*\\.(.*)_.*","\\1",overlaps[,TF_col])
  agg <- aggregate(by = list(overlaps[,target_col], overlaps[,TF_col]), overlaps[,TF_col], length)
  matrix<- spread(agg, Group.2, x)
  names <- matrix[,"Group.1"]
  matrix[is.na(matrix)] <- 0
  if (outnumeric == FALSE) {
    matrix[matrix > 0] <- TRUE
    matrix[,2:ncol(matrix)] <- apply(matrix[,2:ncol(matrix)], MARGIN = 2, as.logical)
    matrix[,"Group.1"] <- names
  } 
  if (! is.null(expressed_genes)) {
    expressed_genes <- c(expressed_genes, "Group.1")
    matrix <- matrix[,colnames(matrix) %in% expressed_genes]
  }
  return(matrix)
}

#Tobias function
getTFPairSummary <- function (final_results, gene_bound_overlaps, bound_overlaps, title ="", TF_cols, ensID = "EnsID", enhID = "Enh.Pos") {
  gene_data <- aggregate(gene_bound_overlaps$TF_Name, by = list(gene_bound_overlaps[,ensID]), function(x) {paste0(x, collapse = ",")})
  enh_data <- aggregate(bound_overlaps$TF_Name, by = list(bound_overlaps[,enhID]), function(x) {paste0(x, collapse = ",")})
  tf_df <- data.frame( ) #list()
  for (pair in unique(final_results$Pair)) { 
    tmp <- final_results[final_results$Pair == pair,][1,] #Only get first instance?
    enh_TFs <- unlist(str_split(enh_data[enh_data$Group.1 == tmp[,enhID],"x"], ",")) 
    gene_TFs  <- unlist(str_split(gene_data[gene_data$Group.1 == tmp[,ensID],"x"], ","))
    tf_df <- rbind(tf_df, cbind(pair,  
                               sum(enh_TFs %in% gene_TFs), 
                               sum(gene_TFs %in% enh_TFs),
                               sum(enh_TFs %in% gene_TFs) / length(enh_TFs),
                              length(intersect(enh_TFs, gene_TFs)) / length(union(enh_TFs, gene_TFs)),
                              paste0(intersect(enh_TFs, gene_TFs), collapse = ",")
                              ))
  }
  
  colnames(tf_df) <- c("Pair" , paste0(title, TF_cols))
  # tf_df$Gene <- sub(".*_", "", tf_df$Pair)
  # tf_df$Enh <- sub("_.*", "", tf_df$Pair)
  tf_df[is.na(tf_df)] <- 0
  #removed section binding to data.frame
  return(tf_df)
}
TF_cols <- c("Sum_Enh_TFs_in_Gene_TFs",  "Sum_Gene_TFs_in_Enh_TFs" , "Pair_TF_overlaps", "Prop_shared_TFs", "Intersecting_TFs")

#Tobias function
update_overlaps <- function (overlaps, final_results, gene = FALSE, Hitcol = "HitPermissive", factor = NULL, subset = T) {
  if (gene == TRUE) {
    if (subset == T) {
      overlaps <- overlaps[,c(1:4, 11:ncol(overlaps))]
    }
    if (is.null(factor))
      factor <- "EnsID"
    a <- c("Peak_chr", "Peak_start", "Peak_end", "Peak_ID", "Gene_chr", "Gene_start", "Gene_end", factor) #,"Gene"
    
  } else {
    if (is.null(factor))
      factor <- "Enh.Pos"
    a <- c("Enh_chr", "Enh_start", "Enh_end", factor)
  }
  colnames(overlaps) <- c(a, "TF_chr", "TF_start", "TF_end", "TF_Name")
  overlaps$TF_Name <- sub("(.*)_.*", "\\1", overlaps$TF_Name)
  overlaps$TF_Name <- sub(".*\\.(.*)", "\\1", overlaps$TF_Name)
  if (!is.null(Hitcol)) {
    overlaps$Hits <- overlaps[,factor] %in% final_results[final_results[,Hitcol],factor]
  }
  return(overlaps)
}


#Enhancer annotation function
addEGPdatatoEnh <- function (EGP, Enh, doHits = F, geneCol, bulkExp = F ) {
  Enh$Enh.Midpoint <- Enh$Enh.start + round(Enh$Enh.size / 2)
  if ("Gene.Nearest" %in% colnames(EGP))
    Enh <-  merge(Enh, setNames(aggregate(EGP$Gene.Nearest , by = list(EGP$Enh), FUN = sum), c("Enh", "NearestTested")), all.x = T)
  if ("Gene.Exp" %in% colnames(EGP))
    Enh <-merge(Enh, setNames(aggregate(EGP$Gene.Exp, by = list(EGP$Enh), FUN = max), c("Enh", "MaxExpression")), all.x = T)
  Enh <- merge(Enh, setNames(aggregate(EGP$Gene.Distance , by = list(EGP$Enh), FUN = min), c("Enh", "ClosestGeneDist")), all.x = T)
  
  if (bulkExp == T) {
    expcols <- c("Gene.TTseq_Counts","Gene.RNAseq_RPKM","Gene.Ref_RNAseq_RPKM")
    expcols <- expcols[expcols %in% colnames(EGP)]
    for (col in expcols) {
      outname <- sub("Gene.","",col)
      #print(col)
      Enh <- merge(Enh, setNames(aggregate(EGP[,col], by = list(EGP$Enh), FUN = max), c("Enh", paste0("Max",outname))), all.x = T)
      Enh <-  merge(Enh, setnames(aggregate(EGP[EGP$Gene.Distance < 50000,col], by = list(EGP[EGP$Gene.Distance < 50000,]$Enh), FUN = max),  c("Enh", paste0("Max",outname,"LT50kb"))), all.x  = T)
      Enh <- merge(Enh, unique(EGP[,c("Gene.Distance", "Enh", col)]), by.x = c("ClosestGeneDist", "Enh"), by.y = c("Gene.Distance", "Enh"), all.x = T)
      colnames(Enh)[colnames(Enh) == col] <- paste0("Closest",col)
      Enh[is.na(Enh[,paste0("Max",outname)]) | is.infinite(Enh[,paste0("Max",outname)]), paste0("Max",outname)] <- 0
      Enh[is.na(Enh[,paste0("Max",outname,"LT50kb")]) | is.infinite(Enh[,paste0("Max",outname,"LT50kb")]), paste0("Max",outname,"LT50kb")] <- 0
      Enh[is.na(Enh[,paste0("Closest",col)]), paste0("Closest",col)] <- 0
      #Enh[is.na(Enh$NGenesTested),c("NearestTested",paste0("Closest",col))] <- 0
    }
  }
  
  if (doHits) {
    if ("Z" %in% colnames(EGP)) {
      Enh <-merge(Enh, setNames(aggregate(EGP$Z , by = list(EGP$Enh), FUN = min), c("Enh", "Z")), all.x = T)
      Enh <-merge(Enh, setNames(aggregate(EGP$Z , by = list(EGP$Enh), FUN = max), c("Enh", "MaxZ")), all.x = T)
    }
    if ("HitCore" %in% colnames(EGP))
      Enh <- merge(Enh, setNames(aggregate(EGP$HitCore , by = list(EGP$Enh), FUN = sum), c("Enh", "HitCore")), all.x = T)
    if ("HitPermissive" %in% colnames(EGP))
      Enh <- merge(Enh, setNames(aggregate(EGP$HitPermissive , by = list(EGP$Enh), FUN = sum), c("Enh", "HitPermissive")), all.x = T)
    if ("HitPermissive_NegZ" %in% colnames(EGP))
      Enh <- merge(Enh, setNames(aggregate(EGP$HitPermissive_NegZ, by = list(EGP$Enh), FUN = sum), c("Enh", "HitPermissive_NegZ")), all.x = T)
  }
  if ("logfc.vst" %in% colnames(EGP))
    Enh <- merge(Enh, setNames(aggregate(EGP$logfc.vst, by = list(EGP$Enh), FUN = min), c("Enh", "logfc.vst")), all.x = T)
  
  Enh <- merge(Enh, setNames(aggregate(EGP[,geneCol], by = list(EGP$Enh), FUN = length), c("Enh", "NGenesTested")), all.x = T)
  Enh <- merge(Enh, setNames(aggregate(EGP$Gene.Distance < 50000, by = list(EGP$Enh), FUN = sum), c("Enh", "LT50KB")), all.x = T)
  Enh <- merge(Enh, setNames(aggregate(EGP$Gene.Distance < 100000, by = list(EGP$Enh), FUN = sum), c("Enh", "LT100KB")), all.x = T)
  Enh <- merge(Enh, setNames(aggregate(EGP$Gene.Distance > 100000, by = list(EGP$Enh), FUN = sum), c("Enh", "GT100KB")), all.x = T)
  Enh <- merge(Enh, setNames(aggregate(EGP$Gene.Distance < 200000, by = list(EGP$Enh), FUN = sum), c("Enh", "LT200KB")), all.x = T)
  Enh <- merge(Enh, setNames(aggregate(EGP$ABC_Approx, by = list(EGP$Enh), FUN = max), c("Enh", "ABC_Approx.Max")), all.x = T)
  Enh <- merge(Enh, setnames(aggregate(EGP$ABC_Score, by = list(EGP$Enh), FUN = max), c("Enh", "ABC_Score.Max")), all.x = T)
  
  if ("Gene.Exp" %in% colnames(EGP)) {
    Enh <- merge(Enh, unique(EGP[,c("Gene.Distance", "Enh", "Gene.Exp")]), by.x = c("ClosestGeneDist", "Enh"), by.y = c("Gene.Distance", "Enh"), all.x = T)
    colnames(Enh)[colnames(Enh) == "Gene.Exp"] <- "ClosestGeneExp"
    Enh <-  merge(Enh, setnames(aggregate(EGP[EGP$Gene.Distance < 50000,]$Gene.Exp, by = list(EGP[EGP$Gene.Distance < 50000,]$Enh), FUN = max),  c("Enh", "MaxExpressionLT50kb")), all.x  = T)
    Enh[is.na(Enh$MaxExpression), "MaxExpression"] <- 0
    
    Enh[is.na(Enh$MaxExpressionLT50kb), "MaxExpressionLT50kb"] <- 0
    Enh[is.infinite(Enh$MaxExpressionLT50kb), "MaxExpressionLT50kb"] <- 0
    Enh[is.infinite(Enh$MaxExpression),"MaxExpression"] <- 0
    Enh[is.na(Enh$ClosestGeneExp), "ClosestGeneExp"] <- 0
    Enh[is.na(Enh$NGenesTested),c("NearestTested","ClosestGeneExp")] <- 0
  }
  Enh[is.na(Enh$ClosestGeneDist),"ClosestGeneDist"] <- 1000000
  Enh[is.na(Enh$ABC_Score.Max),"ABC_Score.Max"] <- 0
  Enh[is.na(Enh$ABC_Approx.Max),"ABC_Approx.Max"] <- 0
  Enh[is.na(Enh$NGenesTested),c("NGenesTested", "LT50KB", "LT100KB", "LT200KB", "GT100KB")] <- 0
  return(Enh)
}


#Also Enhancer annotation function
add_allgene_tech <- function (enhancers, counts.file, RPKM.col = "Ref_RNAseq_RPKM") {
  feature_counts <- read.csv(counts.file)
  colnames(feature_counts)[1] <- "Symbol"
  geneinfo <- read.table("../FullScale/Data/Whitelists/GeneInfo.txt")
  geneinfo$Chr <- sub("chr", "",geneinfo$Chr)
  geneinfo$TSS <- as.numeric(sub(".*:", "",geneinfo$TSS))
  geneinfo <- merge(geneinfo, feature_counts)
  #Initialise variables
  enhancers[,c("AllGenes.ClosestExpressedGene","AllGenes.ClosestGeneDist")] <- 1000000
  enhancers[,c("AllGenes.LT50KB","AllGenes.LT500KB", 
               "AllGenes.Max_RNAseq_RPKMLT500kb", "AllGenes.Max_RNAseq_RPKMLT50kb",
               "AllGenes.ClosestGeneExp")] <- 0
  enh_list <- mclapply(as.list(unique(enhancers$Enh)), mc.cores = 10, function(enh) {
    enh_row <- enhancers[enhancers$Enh == enh,]
    chr_genes <- geneinfo[geneinfo$Chr  == enh_row[,"Enh.chr"],]
    LT50kb <- abs(chr_genes$TSS - enh_row$Enh.Midpoint) < 50000
    LT500kb <- abs(chr_genes$TSS - enh_row$Enh.Midpoint) < 500000
    enh_row$AllGenes.ClosestGeneDist <- min(abs(chr_genes$TSS - enh_row$Enh.Midpoint))
    enh_row$AllGenes.ClosestExpressedGene <- min(abs(chr_genes[chr_genes[,RPKM.col] > 1,]$TSS - enh_row$Enh.Midpoint))
    enh_row$AllGenes.LT50KB <- sum(LT50kb)
    enh_row$AllGenes.Max_RNAseq_RPKMLT50kb <- max(chr_genes[LT50kb,RPKM.col])
    enh_row$AllGenes.LT500KB <- sum(LT500kb)
    enh_row$AllGenes.Max_RNAseq_RPKMLT500kb <- max(chr_genes[LT500kb,RPKM.col])
    enh_row$AllGenes.ClosestGeneExp <- chr_genes[which.min(abs(chr_genes$TSS - enh_row$Enh.Midpoint)),RPKM.col]
    return(enh_row)
  })
  enh_list <- do.call("rbind",enh_list)
  enh_list[is.infinite(enh_list$AllGenes.Max_RNAseq_RPKMLT500kb), "AllGenes.Max_RNAseq_RPKMLT500kb"] <- 0
  enh_list[is.infinite(enh_list$AllGenes.Max_RNAseq_RPKMLT50kb), "AllGenes.Max_RNAseq_RPKMLT50kb"] <- 0
  return(enh_list)
}






#Random forest function
getTrainControl <- function(method = "cv", number = 10, modelsTested = 6, rf.method) {
  # if (rf.method == "ranger") {
  #   seeds <- readRDS("Data/rf_seeds.RData")
  # } 
  set.seed(seed)
  seeds <- vector(mode = "list",length = (number + 1))
  for(i in 1:number) seeds[[i]] <-  sample.int(n=1000, modelsTested )
  seeds[[number+1]] <-  sample.int(n=1000, 1)
  tc <-  trainControl(method= method, number= number, savePredictions = TRUE, seeds = seeds)
  return(tc)
}
#subset = Variable subset, importance = file to output importance to, train_var = variable used to train, 
#Method has several meanings (cv vs. oob and "ranger" vs. xgboost, please check) This returns a caret model
getRFmodel <- function (all_rf_factors, subset, train_var = "Z", importance = NULL, train.method = "ranger", weights = NULL) { #, probability = F
  random_forest_data <- all_rf_factors[,c(train_var, subset)]
  colnames(random_forest_data)[colnames(random_forest_data) == train_var] <- "train_var"
  if(str_detect(train_var, "Hit|Significant")) {
    random_forest_data$train_var <- as.factor(random_forest_data$train_var > 0)
  }
  set.seed(seed)
  train_control <- getTrainControl(rf.method = train.method)
  #caret train module gets hyper parameters
  model <- train(train_var ~., data=random_forest_data, trControl=train_control, 
                 method= train.method,weights = weights, seed = seed) # , probability = probability
  if (! is.null(importance) ) {
    set.seed(seed)
    rf <- ranger(train_var ~ ., random_forest_data, seed = seed, mtry = model$bestTune$mtry, 
                 splitrule = model$bestTune$splitrule, num.trees = 500, 
                 min.node.size = model$bestTune$min.node.size, importance = "permutation")
    write.csv(rf$variable.importance[order(rf$variable.importance, decreasing = T)], file = importance)
  }
  return(model)
}
#vars,
getRFmodelPreds <- function (all_rf_factors, model, colname = "rf_pred", pred.method = "oob", train_var = NULL, probability = F, save_models= NULL) { #vars = NULL, 
  all_rf_factors[,colname] <- 0 
  vars <- colnames(model$trainingData)[colnames(model$trainingData) != ".outcome"]
  if (is.null(train_var)) {
    stop("Must set target var train_var")
  }
  set.seed(seed)
  rf <- ranger(all_rf_factors[,train_var] ~ ., data = all_rf_factors[,vars],num.trees = 1000,
               seed = seed,mtry = model$bestTune$mtry, 
               splitrule = model$bestTune$splitrule, 
               min.node.size = model$bestTune$min.node.size, 
               probability = probability)
  #Use for Enhancer - Gene Pairs
  if (pred.method == "cv") {
    #Actual cross folding with unique genes and enhancers in each fold
    all_rf_factors <- crossFoldGeneLevel(all_rf_factors, model, vars = vars,
                                         train_var = train_var, colname  = colname, 
                                         probability = probability )
  } else if (pred.method == "oob") {
    #Predictions for Enhancers
    if (probability == F) {
      all_rf_factors[,colname] <- rf$predictions
    } else {
      all_rf_factors[,colname] <- rf$predictions[,"TRUE"]
    }
  }
  if (!is.null(save_models)) {
    saveRDS(rf, file =save_models)
  }
  
  return(all_rf_factors)
}


plotZvPred <- function (Z, Pred, Hits,train_var,title = "", doplot = T, facet = F) {
  correlation <- signif(cor(Z,Pred),3)
  if (doplot == T) {
    df <- data.frame(Z, Pred, Hits)
    plot.title <- paste("Cor=",correlation, title)
    my_plot<- ggplot(df, aes(Z, Pred,  colour = as.factor(Hits))) +
            geom_point() + 
            stat_smooth(method = "lm", formula = y~ x)  +
            labs(title = plot.title, x = train_var) +
            theme_bw() +
            scale_y_reverse() + 
            scale_x_reverse()
    print(my_plot)
    if (facet == T) {
      print(my_plot +facet_wrap(vars(Hits)))
    }
  }
  return(correlation)
}

getRFperformance <- function (train_col, pred_col, hits, train_var, plot.title, doplot = T, probability = F) {
  if (T == probability) {
    perf <- aucPlots(pred_col, hits, title = plot.title, doplot = doplot)
  } else {
    perf <- aucPlots(- pred_col, hits, title = plot.title, doplot = doplot)
  }
  cor <- plotZvPred(Z = train_col, Pred = pred_col, Hits = hits, train_var = train_var ,title = plot.title, doplot = doplot)
  perf <- rbind(perf, cor)
  return(perf)
}

write.bed <- function (data,file) {
  write.table(data, file, row.names = F, quote = F, sep = "\t", col.names = F)
}

#
#########Cross fold method
getFolds <- function(data, nfolds = 10) {
  data <- sample(data, size = length(data)) #randomises data order
  foldsize <- length(data) / nfolds
  folds <- list()
  for (i in 1:nfolds) {
    lower <- ceiling((foldsize * (i-1)) + 1)
    upper <- ceiling(foldsize * i)
    folds[[i]] <- data[lower:upper]
  }
  return(folds)
}
#model is a caret CV train model with best tune information. 
crossFoldGeneLevel <- function (df, model ,vars,train_var = "logfc.vst", nfolds = 10, colname = "rf_pred", probability = F) {
  df_order <- df[,"Pair"]
  set.seed(seed)
  if (probability == T ) {
    df[,train_var] <- as.factor(df[,train_var])
  }
  genefolds <- getFolds(unique(df[,"Gene"]),nfolds = 10)
  enhfolds <- getFolds(unique(df[,"Enh"]),nfolds = 10)
  out <- data.frame()
  message(paste0("Running CV performance\nRunning:",nfolds, "x", nfolds,"folds"))
  for (genefold in 1:nfolds) {
    for (enhfold in 1:nfolds) {
      train_set <- df[! df[,"Gene"] %in% genefolds[[genefold]] & ! df[,"Enh"] %in% enhfolds[[enhfold]],]
      test_set <- df[df[,"Gene"] %in% genefolds[[genefold]] & df[,"Enh"] %in% enhfolds[[enhfold]],]
      #message(paste("Test Set Size",nrow(test_set)))
      if (nrow(test_set) > 0) {
        set.seed(seed)
        rf <- ranger(train_set[,train_var] ~ ., data = train_set[,vars],num.trees = 1000,
                     seed = seed,mtry = model$bestTune$mtry, 
                     splitrule = model$bestTune$splitrule, 
                     min.node.size = model$bestTune$min.node.size, 
                     probability = probability)
        pred <- predict(rf, test_set)
        if (probability == F) {
          test_set[,colname] <- pred$predictions
        } else {
          test_set[,colname] <- pred$predictions[,"TRUE"]
        }
        out <- rbind(out, test_set)
      }
    }
  }
  #The variable has to remain a non-factor for some functions
  if (probability == T ) {
    out[,train_var] <- as.logical(out[,train_var])
  }
  #ensure input order is the same as output order (stops bugs with consistency)
  out <- out[match( df_order, out[,"Pair"]),]
  return(out)
}

