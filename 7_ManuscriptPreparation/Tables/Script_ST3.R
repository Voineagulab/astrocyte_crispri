setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable3_DE/")
options(stringsAsFactors = FALSE)
source("../../../FullScale/Scripts/Functions.R")
guides <- read.csv("../../../FullScale/Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)


## 3A: CRISPRi DE results on enhancers

  ## Read in
    x <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
  
  ## Add EnsID
    m <- match(x$Gene, geneInfo$Symbol)
    x$GeneID <- geneInfo$EnsID[m]
      
  ## Convert logfc to log2fc
    x$log2FC <- log2(exp(x$logfc.vst))
    
  ## Trim columns and reorder
    x <- x[,c("Pair", "Enh", "Gene", "GeneID", "log2FC", "Z", "P.SCEPTRE",  "FDR.SCEPTRE", "P.N50", "FDR.N50",  "HitPermissive", "Enh.Pos", "Gene.Pos", "Gene.TSS", "Gene.Distance", "Gene.Exp", "nCells"),]  
    colnames(x) <- c("Pair", "Enhancer", "GeneSymbol", "GeneID", "log2FC", "Z_Sceptre", "P_Sceptre", "FDR_Sceptre", "P_Empirical",  "FDR_Empirical", "Hit", "EnhancerCoord", "GeneCoord", "GeneTSS", "PairDistance", "GeneExpression", "nCells")
    
  ## Sort rows
    # by hits, enhancer, gene
    ord <- order(-as.numeric(x$Hit), as.numeric(gsub("Enh", "", x$Enhancer)), x$GeneSymbol)
    x <- x[ord,]
  
  
    
  ## Output
    write.csv(x, file = "3A_EnhancerGenePairDE.csv", row.names = FALSE)

## 3B: negative controls
    ## Read in
      load("../../../FullScale/Results/2_DE/Neg/SCEPTRE Output (Neg, Guide-level).rda")
      load("../../../FullScale/Results/2_DE/Neg/SCEPTRE Output (NegE, Guide-level).rda")
      
      neg <- list()
      neg$N250 <- de.neg.guidelevel
      neg$N50 <- de.negE
      
    ## Process
      neg <- lapply(neg, function(x) {
        colnames(x) <- c("Gene", "Guide", "Type", "SCEPTRE_P", "SCEPTRE_Z")
        x <- x[,c("Guide", "Gene", "SCEPTRE_P", "SCEPTRE_Z")]
      })
      
      neg <- do.call("rbind", neg)
      neg <- neg[order(as.numeric(gsub("Neg_|Neg_E_", "", neg$Guide)), neg$Gene),]
  
  ## Round to save disk space (lops off two dig)
      neg$SCEPTRE_Z <- round(neg$SCEPTRE_Z, 4)
      
  ## Save
      write.csv(neg, file = "3B_NegativeControls.csv", row.names = FALSE)
    
## 3C: positive controls
    
  ## Read in
    pos <- read.csv("../../../FullScale/Results/2_DE/Pos/Results Matrix.csv", row.names = 1)
    
  ## Convert logfc to log2fc
    pos$Log2FC <- log2(exp(pos$LogFC))
    
  ## Add note where the guide name and gene name are in disagreement
    pos$Notes <- "."
    pos$Notes[which(pos$Guides == "FAM96B")] <- "The tested gene (CIA02B) and guide id (FAM96B) are synonyms."
    pos$Notes[which(pos$Guides == "SEPT11")] <- "The tested gene (SEPTIN11) and guide id (SEPT11) are synonyms."
    
  ## Rename guide ids
    pos$Guides <- paste0("Pos_", pos$Guides)
    
    
  ## Clean columns
    
    pos <- pos[,c("Guides", "Gene", "Log2FC", "P", "FDR", "Hit", "Notes")]
    colnames(pos)[4:5] <- c("SCEPTRE_P", "SCEPTRE_FDR")
    
    
  ## Sort rows
    pos <- pos[order(pos$Gene),]
    
  ## Redo hits
    pos$Hit <- pos$SCEPTRE_FDR < 0.05
    
  # ## Add nCells
  #   load("../../../FullScale/Results/2_DE/Sceptre Input Files.rda", verbose = TRUE)
  #   nCellsPerTarget <- rowSums(sceptre.guide.pooled)
  #   m <- match(pos$Gene, splitter(names(nCellsPerTarget), "Pos_", 2))
  #   pos$nCells <- nCellsPerTarget[m]
    
  ## Save
    write.csv(pos, file = "3C_PositiveControls.csv", row.names = FALSE)
    

    
    
## 3E: Nanostring replication
  ## Read in  metadata
    load("../../../Validation_RNAseq/Nanostring/Results/Final/ProcessedData.rda", verbose = TRUE)  
    res_nano <- read.csv("../../../Validation_RNAseq/Nanostring/Results/Final/LinearModels.csv", row.names = 1)
    
  ## Metadata
    rownames(meta_nano) <- splitter(rownames(meta_nano), "\\.", 2)
  
    # add metadata column on whether the sample was tested
    meta_nano$AnalysisBatch <- "."
    meta_nano$AnalysisBatch[meta_nano$Batch == "Run_1"] <- "b1"
    meta_nano$AnalysisBatch[meta_nano$Batch %in% c("Run_2", "Run_3")] <- "b23"
    meta_nano$AnalysisBatch[meta_nano$Batch %in% c("Run_4", "Run_5")] <- "b45"
    
    meta_nano$AnalysisIncluded <- "Y"
    meta_nano$AnalysisIncluded[meta_nano$Input_ng == 300] <- "N"
    # meta_nano$AnalysisIncluded[meta_nano$Batch == "Run_1"] <- "N"
    meta_nano$AnalysisIncluded[meta_nano$Batch == "Run_1" & meta_nano$Enh == "Enh854_HSPB1" & meta_nano$Input_ng == 150] <- "Y" # the exception to Run1, as being tested
    meta_nano$AnalysisIncluded[meta_nano$Batch == "Run_1" & (is.na(meta_nano$Enh)) & meta_nano$Input_ng == 150] <- "N (but used as outgroup)" # background samples for the above
    
    # reorder
    meta_nano <- relocate(meta_nano, c(colnames(meta_nano)[1:6], "AnalysisBatch", "AnalysisIncluded"))
    
    # save
    # write.csv(meta_nano, "../../Tables/STable_Nanostring/SuppTable_Nanostring_Meta.csv")
  
  ## Process results table
    res_nano <- res_nano[,-9] # bonferroni correction, not used
    m <- match(res_nano$Group, meta_nano$Group)
    res_nano$Guide <- meta_nano$Guide[m]
    res_nano <- relocate(res_nano, "Guide")
    
    # add a note on FTH1 in batch1
    res_nano$Notes <- "."
    g <- which(res_nano$Batches == "b1" & res_nano$Group == "FTH1")
    res_nano$Notes[g] <- "A technical replicate of this was performed in batches 2+3, with higher n. For all analyses and plots in the manuscript (e.g. overall replication rate), we ignore this result"
    
    # annotate the pooled results
    w <- grep("Pool", res_nano$Group)
    res_nano$Notes[w] <- "A pool of guides was transduced into this group"
    res_nano$Guide[w] <- paste0("All guides for ", splitter(res_nano$Group[w], "_", 1))
    
    # annotate ANKRD1
    w <- which(res_nano$Guide == "Enh168_g3_chr10:90931016-90930996")
    res_nano$Notes[w] <- "In the enhancer screen, this guide silenced both ANKRD1 and PCGF5, thus replication was assessed for both"
    
    # filter columns and rename
    res_nano$Group <- splitter(res_nano$Guide, "_", 1)
    colnames(res_nano)[c(2,4)] <- c("Enhancer_targeted", "Gene_tested")
    res_nano$Enhancer_targeted[grepl("All guides", res_nano$Enhancer_targeted) & grepl("HSPB1", res_nano$Enhancer_targeted)] <- paste(res_nano$Enhancer_targeted[grepl("Enh", res_nano$Enhancer_targeted) & grepl("HSPB1", res_nano$Gene_tested)], collapse = "+")
    res_nano$Enhancer_targeted[grepl("All guides", res_nano$Enhancer_targeted) & grepl("LGALS3", res_nano$Enhancer_targeted)] <- paste(res_nano$Enhancer_targeted[grepl("Enh", res_nano$Enhancer_targeted) & grepl("LGALS3", res_nano$Gene_tested)], collapse = "+")
    res_nano$Enhancer_targeted[grepl("All guides", res_nano$Enhancer_targeted) & grepl("PTMA", res_nano$Enhancer_targeted)] <- paste(res_nano$Enhancer_targeted[grepl("Enh", res_nano$Enhancer_targeted) & grepl("PTMA", res_nano$Gene_tested)], collapse = "+")
    res_nano$Enhancer_targeted[grepl("All guides", res_nano$Enhancer_targeted) & grepl("NEAT1", res_nano$Enhancer_targeted)] <- paste(res_nano$Enhancer_targeted[grepl("Enh", res_nano$Enhancer_targeted) & grepl("NEAT1", res_nano$Gene_tested)], collapse = "+")
    
    # rename fold-change column
    colnames(res_nano) <- gsub("log2fc", "Log2FC", colnames(res_nano))
    
    # save
    write.csv(res_nano, "3E_NanostringResults.csv", row.names = FALSE)
    
    
    
## 3E: gene skipping
  ## Load
    skip <- read.csv("../../../FullScale/Results/3_HitEnrichment/EnhGenePairs/Intervening Gene Classification Between EGPs.csv", row.names = 1)
    
  ## Process
    colnames(skip)[4] <- "Hit"
    
  ## Save
    write.csv(skip, "3F_SkippingGenes.csv", row.names = FALSE)
    
    
