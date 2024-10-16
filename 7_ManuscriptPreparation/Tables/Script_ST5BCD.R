setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable5_TFs_RegulatoryCircuitry/")
options(stringsAsFactors = FALSE)
source("../../../FullScale/Scripts/Functions.R")
res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")

## 5A appears in a separate script

## 5B: binding of TFs to enhancers
  ## Load
    tfbs <- read.csv("../../../EnhancerPredictionModels/Results/Tobias/Summaries/Bound_Matrix_Expressed.csv")
    
  ## Clean columns
    tfbs <- tfbs[,-grep("Unbound", colnames(tfbs))]
    tfbs <- tfbs[,-grep("Tobias.Exp_Bound_TF_counts", colnames(tfbs))]
    tfbs <- tfbs[,-c(1,4)] # hit calls and enhancer coordinates
    tfbs$HitPermissive <- tfbs$HitPermissive > 0
    
    
    colnames(tfbs) <- gsub("Tobias.Bound_", "", colnames(tfbs))
  
  ## Add hit gene
    tfbs$LinkedGene <- NA
    h <- res.final[which(res.final$HitPermissive),]
    for (j in which(tfbs$Hit)) {
      h1 <- h[which(h$Enh == tfbs$Enh[j]),]
      tfbs$LinkedGene[j] <- paste(h1$Gene, collapse = "/")
    }
    
  ## Reorder columns
    tfbs <- relocate(tfbs, c("Enh", "HitPermissive", "LinkedGene"))
    
    colnames(tfbs)[2] <- "Hit"
    
  ## Save
    write.csv(tfbs, "5B_TobiasBound.csv", row.names = FALSE)

# ## 5C: enrichments old version
#   ## Load
#     tfEnr <- read.csv("../../../EnhancerPredictionModels/Results/Tobias/Summaries/FisherTestsExpressedTFMatrix.csv")
# 
#   ## Filter rows
#     tfEnr <- tfEnr[-grep("Unbound", tfEnr$Variable),]
#     tfEnr$Variable <- gsub("Tobias.Bound_", "", tfEnr$Variable)
#     tfEnr <- tfEnr[which(tfEnr$Variable %in% colnames(tfbs)),]
# 
#   ## Filter columns
#     tfEnr <- tfEnr[,c("Variable", "OR", "p", "Lower", "Upper")]
#     colnames(tfEnr)[1] <- "TF"
# 
#   ## Save
#     write.csv(tfEnr, "5C_TobiasBoundEnrich.csv", row.names = FALSE)
    
## 5C new version which only considers well-powered non-hit enhancers
#Tobias Class plots
  test_boundMatrix <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Summaries/JASPAR_group_AnnotatedTFs.csv")
  colnames(test_boundMatrix)[colnames(test_boundMatrix) == "Odds.Ratio"] <- "Odds Ratio"
  #Alter names so they fit on the plot
  test_boundMatrix$JASPAR_Class <- gsub("Nuclear receptors", "Nuc. recep.", test_boundMatrix$JASPAR_Class) # %>%
    # gsub("binding", "", .)# test_boundMatrix$JASPAR_Class <- gsub("Nuclear receptors", "Nuc. recep.", test_boundMatrix$JASPAR_Class) %>%
    # gsub("binding", "", .)
  
  # replace p-values with those from the well-powered subset of enhancers
  x <- read.csv("../../../EnhancerPredictionModels/Results/Tobias/Summaries/FisherTestsExpressedTFMatrix_wellpowered.csv")
  # y <- read.csv("../../../EnhancerPredictionModels/Results/Tobias/Summaries/FisherTestsExpressedTFMatrix.csv") # confirms that the existing p-values are from here
  m <- match(test_boundMatrix$Variable, x$Variable)
  test_boundMatrix$P <- x$p[m]
  test_boundMatrix$`Odds Ratio` <- x$OR[m]
  test_boundMatrix$FDR <- p.adjust(test_boundMatrix$P, method = "fdr")
  test_boundMatrix$BON <- p.adjust(test_boundMatrix$P, method = "bonf")
  
  # clean
  test_boundMatrix <- test_boundMatrix[,c(1,3,4,6,7,8,9,11,12)]
  notDimer <- test_boundMatrix$TF.1 == test_boundMatrix$TF.2
  test_boundMatrix$TF.2[notDimer] <- NA
  test_boundMatrix <- test_boundMatrix[-which(is.na(test_boundMatrix$P)),]
  colnames(test_boundMatrix)[1] <- "Name"
  test_boundMatrix <- test_boundMatrix[order(test_boundMatrix$P),]
  
  write.csv(test_boundMatrix, file = "5C_TobiasBoundEnrich_wellpowered.csv", row.names = FALSE)
    
## 5D: astrocyte-specificity of genes, enhancers, and TFs
  ## Read in the elements showing specificity to check
    net_ge <- read.csv("../../../IV/RESULTS/Publication/EData.EGheatmap_v2.csv", row.names = 1)
    net_tf <- read.csv("../../../IV/RESULTS/Publication/EData.TFheatmap_v2.csv", row.names = 1)
    
    to_check <- c(rownames(net_ge), rownames(net_tf))
    if (any(c("Enh8301", "Enh1621") %in% to_check)) {
      to_check <- to_check[-which(to_check %in% c("Enh8301", "Enh1621"))]
    }
    
  ## Get values
    # gene DE
    astMod_rna <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Astrocyte Markers - Herring 2022.csv", row.names = 1)
  
    # chromatin de  
    astMod_atac <- read.csv("../../../FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Cell-type Specificity Models.csv")
 
  ## Create table
    dat <- data.frame(Node = to_check,
                      Class = ".")
    
    dat$Class[which(dat$Node %in% rownames(net_tf))] <- "TF"
    dat$Class[grep("^Enh", dat$Node)] <- "FunctionalEnh"
    dat$Class[!(grepl("^Enh", dat$Node)) & dat$Node %in% rownames(net_ge)] <- "RegulatedGene"
    
    dat <- dat[order(dat$Class),]
    
  ## Add statistics for genes
    m <- match(dat$Node, astMod_rna$Gene)
    
    dat$Log2FC <- astMod_rna$log2fc[m]
    dat$P <- astMod_rna$P[m]
    dat$FDR <- astMod_rna$FDR[m]
    
  ## Add statistics for enhancers
    w <- which(dat$Class == "FunctionalEnh")
    m <- match(dat$Node[w], astMod_atac$X)
    dat$Log2FC[w] <- astMod_atac$Astro.EffectSize[m]
    dat$P[w] <- astMod_atac$Astro.P[m]
    dat$FDR[w] <- astMod_atac$Astro.FDR[m]
    
  ## Load family information
    fam <- read.csv("../../../EnhancerPredictionModels/Results/Tobias/Summaries/JASPAR_TF_Families.csv")
    fam <- fam[which(fam$name %in% rownames(net_tf)),]
    length(unique(fam$family))
    
  ## Output
    write.csv(dat, "5D_AstrocyteSpecificNodes.csv", row.names = FALSE)
    