################
#
#This script define R Functions used for TOBIAS analyses.
#@author Sam Bagot
#@Date: 23-03-03
#
############

##WARNING THIS MAKES START POSITION WRONG (CURRENTLY OUT BY 1BP)
process_counts <- function (counts, final_results, TF_counts = "TF_Counts") {
  colnames(counts) <- c("Chr", "Start", "End", "Peak", "1", "strand", "2", "3", "4", "5",  TF_counts )
  counts$Enh.Pos <- paste0(counts$Chr, ":", counts$Start + 1, "-", counts$End)
  counts$Tested <- counts$Enh.Pos %in% final_results$Enh.Pos
  return(counts)
}

#Takes lsit of ensembl IDs (that can be ,comma seperated for each row) and checks if they are in a background
ens_Exists <- function (ens, bkg_ens) {
  split_string <- unlist(str_split(ens, ","))
  exists <- sum(split_string %in% bkg_ens) > 0  
  return(exists) 
}

summary_overlaps <- function (overlaps, homo_TFs, RPKM) {
  values <- nrow(overlaps)
  values <- append(values, sum(overlaps$TF_Name %in%  homo_TFs$concat_names))
  
  Exists <- overlaps$TF_Name %in%  homo_TFs[homo_TFs$TF.1_Exists & homo_TFs$TF.2_Exists,]$concat_names
  values <- append(values, sum(Exists))
  values <- append(values, mean(RPKM[RPKM$Symbol %in% overlaps[Exists,]$TF_Name,"NHA"], na.rm = TRUE))
  
  Expressed <- overlaps$TF_Name %in%  homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names
  values <- append(values, sum(Expressed))
  values <- append(values, mean(RPKM[RPKM$Symbol %in% overlaps[Expressed,]$TF_Name,"NHA"],  na.rm = TRUE))
  
  RPKM_Exists <- overlaps$TF_Name %in%  homo_TFs[homo_TFs$TF.1_RPKM_Exists & homo_TFs$TF.2_RPKM_Exists,]$concat_names
  values <- append(values, sum(RPKM_Exists))
  #values <- append(values, mean(RPKM[RPKM$Symbol %in% overlaps[Expressed,]$TF_Name,"NHA"],  na.rm = TRUE))
  
  names(values) <- c("nrow", "TF_in_Homo_TFs", "TF_Exist", "Mean_Expression_Exists", "RPKM_gt_0.5", "Mean_Expression_RPKM_gt_0.5", "RPKM_Exists")
  values[['Fraction_Exist']] <- values[['TF_Exist']] / values[['TF_in_Homo_TFs']] 
  values[['Fraction_RPKM_0.5']] <- values[['RPKM_gt_0.5']] / values[['TF_in_Homo_TFs']]
  values[['Fraction_RPKM_Exists']] <- values[['RPKM_Exists']] / values[['TF_in_Homo_TFs']]
  return(values)
}

#This function is not reusable and would need to be rewritten if being used again
getProportion <- function (unbound_overlaps, bound_overlaps, homo_TFs) {
  tf_sums <- merge(as.data.frame(table(unbound_overlaps$TF_Name)), as.data.frame(table(bound_overlaps$TF_Name)), by = "Var1", all = T)
  colnames(tf_sums) <- c("TF", "unbound_TFs", "bound_TFs")
  #tf_sums[is.na(tf_sums$unbound_TFs),]$unbound_TFs <- 0
  #tf_sums[is.na(tf_sums$bound_TFs),]$bound_TFs <- 0
  
  tf_sums <- merge(tf_sums, homo_TFs, by.x = "TF", by.y = "concat_names")
  tf_sums$pct_bound <- tf_sums$bound_TFs / (tf_sums$bound_TFs + tf_sums$unbound_TFs)
  write.csv(tf_sums, file = "Results/Tobias/Summaries/TF_sums_per_TF.csv", row.names = FALSE)
  
  
  expressed_TFs <- mean(tf_sums[tf_sums$TF.1_Expressed & tf_sums$TF.2_Expressed, "pct_bound"], na.rm = TRUE)
  expressed_TFs <- append(expressed_TFs, quantile(tf_sums[tf_sums$TF.1_Expressed & tf_sums$TF.2_Expressed, "pct_bound"], na.rm = T))
  unexpressed_TFs <- mean(tf_sums[! tf_sums$TF.1_Expressed & ! tf_sums$TF.2_Expressed, "pct_bound"], na.rm = T)
  unexpressed_TFs <- append(unexpressed_TFs, quantile(tf_sums[! tf_sums$TF.1_Expressed & ! tf_sums$TF.2_Expressed, "pct_bound"], na.rm = T))
  
  prop <- rbind(expressed_TFs, unexpressed_TFs)
  colnames(prop)[1] <- "Percent_Bound"
  write.csv(prop, file = "Results/Tobias/Summaries/Percent_expressed_TFs_bound.csv")
  return(prop)
  #bound_overlaps$TF_Name %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,"concat_names"]
}

#This function is not reusable (getProportion is too rigidly written)
compare_overlaps <- function (unbound_overlaps, bound_overlaps, homo_TFs, RPKM) {
  prop <- getProportion(unbound_overlaps, bound_overlaps, homo_TFs)
  summary <- as.data.frame(t(summary_overlaps(unbound_overlaps, homo_TFs, RPKM)))
  summary[2,] <- summary_overlaps(bound_overlaps, homo_TFs, RPKM)
  
  rownames(summary) <- c("Unbound_TF_counts", "Bound_TF_counts")
  return(summary)
}

TF_facts <- function(unbound_bound, factor = "TF_Counts") {
  Untested <- unbound_bound[unbound_bound$Tested == FALSE, factor]
  Tested <- unbound_bound[unbound_bound$Tested == TRUE, factor]
  values  <- mean(Untested)
  values <- append(values, mean(Tested))
  values <- append(values, wilcox.test(Tested, Untested)$p.value)
  
  Tested_Non_Hits <- unbound_bound[unbound_bound$Tested == TRUE & unbound_bound$Hit == FALSE, factor] 
  Hits <-unbound_bound[unbound_bound$Hit == TRUE, factor]
  
  values <- append(values, mean(Tested_Non_Hits))
  values <- append(values, mean(Hits))
  values <- append(values,wilcox.test(Tested_Non_Hits, Hits)$p.value)
  values <- append(values, quantile(unbound_bound[unbound_bound$Tested == TRUE & unbound_bound$Hit == FALSE, factor]))
  values <- append(values, quantile(unbound_bound[unbound_bound$Hit == TRUE , factor]))
  names(values) <- c("Mean_Untested", "Mean_Tested", "WT_Tested_v_Untested", "Mean_Tested_Non_Hits", "Mean_Hits", "WT_Hits_v_Non-hits", paste0("Quantile_Non_Hits_", names(values[7:11])) , paste0("Quantile_Hits_", names(values[12:16])) )
  return(values)
}

#RPKM is data.frame of format EnsID, RPKM
gethomo_TFs <- function (RPKM, all_EnsIDs) {
  homo_TFs <- read.table("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/TF_BindingSites/JASPAR2022_human_TF.txt")
  homo_TFs$V1 <- sub(".*\\.","",homo_TFs$V1)
  homo_TFs$V2 <- sub("(.*)::(.*)","\\2",homo_TFs$V1)
  homo_TFs$V1 <- sub("(.*)::(.*)","\\1",homo_TFs$V1)
  homo_TFs <- unique(homo_TFs)
  
  ENSB <- getENSID(union(homo_TFs$V1, homo_TFs$V2 ))
  
  homo_TFs <- merge(homo_TFs, ENSB, all.x = TRUE, by.x = "V1", by.y = "Group.1")
  homo_TFs <- merge(homo_TFs, ENSB, all.x = TRUE, by.x = "V2", by.y = "Group.1")
  colnames(homo_TFs) <- c("TF.1", "TF.2", "ENS.1", "ENS.2")
  homo_TFs$TF.1_Exists <- unlist(lapply(homo_TFs$ENS.1, ens_Exists, bkg_ens = all_EnsIDs))
  homo_TFs$TF.2_Exists <- unlist(lapply(homo_TFs$ENS.2, ens_Exists, bkg_ens = all_EnsIDs))
  
  homo_TFs$concat_names <- homo_TFs$TF.2
  homo_TFs[homo_TFs$TF.2 != homo_TFs$TF.1, ]$concat_names <- paste0(homo_TFs[homo_TFs$TF.2 != homo_TFs$TF.1, ]$TF.2, homo_TFs[homo_TFs$TF.2 != homo_TFs$TF.1, ]$TF.1)
  
  #homo_TFs$ENS.1 %in% RPKM[RPKM$NHA >= 0.5, "EnsID"]
  homo_TFs$TF.1_Expressed <- unlist(lapply(homo_TFs$ENS.1, ens_Exists, bkg_ens = RPKM[RPKM[,2] >= 0.5, "EnsID"]))
  homo_TFs$TF.2_Expressed <- unlist(lapply(homo_TFs$ENS.2, ens_Exists, bkg_ens = RPKM[RPKM[,2] >= 0.5, "EnsID"]))
  
  homo_TFs$TF.1_RPKM_Exists <- unlist(lapply(homo_TFs$ENS.1, ens_Exists, bkg_ens = RPKM$EnsID))
  homo_TFs$TF.2_RPKM_Exists <- unlist(lapply(homo_TFs$ENS.2, ens_Exists, bkg_ens = RPKM$EnsID))
  return(homo_TFs)
  #length(unique(RPKM[(RPKM$EnsID %in% homo_TFs$ENS.1 | RPKM$EnsID %in% homo_TFs$ENS.2) & RPKM$NHA >= 0.5,"Symbol"]))
}
