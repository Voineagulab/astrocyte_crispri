################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace//")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Seurat)
  library(sceptre)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(rcartocolor)
  library(tidyverse)
  library(gprofiler2)
  library(disgenet2r)
  library(VennDiagram)
  library(stringr)
  library(SPARQL)
  library(RCurl)
  library(igraph)  
  library(readxl)
  

## Load
  source("../../Scripts/Functions.R")
  # load("../../Data/Preprocessed/NHA Pooled (Final).rda")
  # guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  # guides <- guides[which(guides$Celltype == "NHA"),]

## Plotting
  invis <- element_blank()
  sig.colours <- c("black", "firebrick1")

## Load results
  res.final <- read.csv("../2_DE/Enh/Results Final.csv")
  
## Set up gene lists
  # the hit and background genes, useful for enrichment
  hit.genes <- unique(res.final$Gene[which(res.final$HitPermissive)]) # hit genes
  bg <- unique(res.final$Gene) # highly-expressed genes as background
  
  # a dataframe to store gene-level logical annotations
  annotation.logi <- data.frame(Gene = bg, Hit = bg %in% hit.genes)
  annotation.logi <- annotation.logi[order(!(annotation.logi$Hit), annotation.logi$Gene),]
  
  # a dataframe to store fisher test results for overenrichments
  # (actually, a function to fill initialise/fill this)
  run.fisher.vsHits <- function(data = annotation.logi, data.column, rbind = FALSE, rbind.to = enrichments) {
    # get information
    total <- nrow(data)
    x <- data[,data.column] %>% as.logical()
    y <- data$Hit
    
    # run stats
    f <- table(x, y) %>% fisher.test()
    
    # output stats
    out <- data.frame(Total_TRUE = sum(x),
                      Fraction_Bg_TRUE = sum(x) / total,
                      Total_Hit_TRUE = sum(x & y),
                      Fraction_Hit_TRUE = sum(x & y) / sum(y),
                      p = f$p.value,
                      OR = f$estimate,
                      Lower = f$conf.int[1],
                      Upper = f$conf.int[2],
                      row.names = data.column)  
    
    if (rbind) {
      return(rbind(enrichments, out))
    } else {
      return(out)  
    }
    
  }
    
    
################################################################################################################################ #
## Temporal trends within Astrocytes ----
    
## Definitions:
  # maturation takes immature to mature neurons
  # ageing is mature neurons to death
       
## Maturation: Herring
  # includes snRNA-seq from foetal to adult, where adult is 20-40.
    
  ## Read in
    # file is called herring_trends
    load("../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/DE_Trends.rda")
    herring_trends <- herring_trends[-which(is.na(herring_trends$Astro)),]
    
  ## Add to annotation dataframes
    annotation.logi$AstMaturation_Herring2022 <- annotation.logi$Gene %in% herring_trends$Symbol
    annotation.logi$AstMaturation_Up_Herring2022 <- annotation.logi$Gene %in% herring_trends$Symbol[grep("up", herring_trends$Astro)]
    annotation.logi$AstMaturation_Dn_Herring2022 <- annotation.logi$Gene %in% herring_trends$Symbol[grep("down", herring_trends$Astro)]
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Herring2022", rbind = FALSE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Up_Herring2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Dn_Herring2022", rbind = TRUE)
  
## Maturation: Zhang  
  ## A seminal work which pioneered the immunopanning of astrocytes from the human brain. Many comparisons within:
    # mouse vs human
    # foetal vs adult
    # glioblastoma vs healthy
    # epilepsy vs healthy
    
  ## Read in (actually, results from all comparisons in the study)
    zhang_path <- "../../../PublicData/Zhang2016_ImmunpanningHumanAst/1-s2.0-S0896627315010193-mmc5.xlsx"
    zhang_sheets <- excel_sheets(zhang_path)
    
    zhang_data <- list()
    for (j in zhang_sheets) zhang_data[[j]] <- read_xlsx(zhang_path, sheet = j, skip = 1)
      
  ## Collect the genes
    zhang_data <- lapply(zhang_data, function(x) x$Gene) # column 1 always has gene symbol, but its first entry is always "Gene"
    
  ## Convert mouse symbol to human (where relevant)
    m <- grep("mouse", zhang_sheets) 
    zhang_data[m] <- lapply(zhang_data[m], convert_mouse2human, path.fix = FALSE)
    
  ## Rename levels
    names(zhang_data) <- gsub(" ", "_", names(zhang_data))
    
  ## Add to annotation (foetal vs adult only)
    annotation.logi$AstMaturation_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$up_fetal_vs_adult, zhang_data$down_fetal_vs_adult)
    annotation.logi$AstMaturation_Up_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$up_fetal_vs_adult)
    annotation.logi$AstMaturation_Dn_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$down_fetal_vs_adult)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Zhang2016", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Up_Zhang2016", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Dn_Zhang2016", rbind = TRUE)

  
## Maturation: Krawczyk  
  # similar to the above study, but with a much large scope, with ~4x as many samples
    
  # again, read in all relevant data here, but only output to the annotation dataframes the relevant maturation work
    
  ## Read in all
    krawczyk_data <- list()
    krawczyk_data$Peritumour <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-3_peritumour.xlsx", skip = 1)
    krawczyk_data$FocalCorticalDysplasia <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-9_FCD.xlsx", skip = 1)
    krawczyk_data$Maturation_Up <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Up")
    krawczyk_data$Maturation_Dn <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Down")
    krawczyk_data$Ageing <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-12_Ageing.xlsx", skip = 1)
  
  ## Extract gene symbol
    # krawczyk_data <- lapply(krawczyk_data, function(x) x$`Gene Name`)
    
  ## Add to annotation (maturation only)
    annotation.logi$AstMaturation_Krawczyk2022 <- annotation.logi$Gene %in% c(krawczyk_data$Maturation_Up$`Gene Name`, krawczyk_data$Maturation_Dn$`Gene Name`)
    annotation.logi$AstMaturation_Up_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Maturation_Up$`Gene Name`
    annotation.logi$AstMaturation_Dn_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Maturation_Dn$`Gene Name`
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Up_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Dn_Krawczyk2022", rbind = TRUE)
    
## Ageing: Krawczyk
  # already processed above
    
  ## Add to annotation (ageing only)
    annotation.logi$AstAgeing_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Ageing$`Gene Name`
    annotation.logi$AstAgeing_Up_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Ageing$`Gene Name`[which(krawczyk_data$Ageing$`log2(Fold Change)` > 0)]
    annotation.logi$AstAgeing_Dn_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Ageing$`Gene Name`[which(krawczyk_data$Ageing$`log2(Fold Change)` < 0)]
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Up_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Dn_Krawczyk2022", rbind = TRUE)
    
## Ageing: Palmer 2021
  # snRNA-seq of the ageing brain
    
  ## Read in
    ageing_palmer <- read_xlsx("../../../PublicData/snRNAseq/Palmer2021_Aging/pnas.2114326118.sd06.xlsx", sheet = "Ast")
    colnames(ageing_palmer)[1] <- "Gene"
    ageing_palmer <- ageing_palmer[which(ageing_palmer$p_val_adj < 0.05),]
    
  # note that a positive fold-change indicates higher expression in old vs. young brains.
  
  ## Add to annotation
    annotation.logi$AstAgeing_Palmer2021 <- annotation.logi$Gene %in% ageing_palmer$Gene
    annotation.logi$AstAgeing_Up_Palmer2021 <- annotation.logi$Gene %in% ageing_palmer$Gene[which(ageing_palmer$avg_logFC > 0)]
    annotation.logi$AstAgeing_Dn_Palmer2021 <- annotation.logi$Gene %in% ageing_palmer$Gene[which(ageing_palmer$avg_logFC < 0)]
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Palmer2021", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Up_Palmer2021", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Dn_Palmer2021", rbind = TRUE)
  

################################################################################################################################ #
## Astrocyte biology: Disease DEGs ----
    
# ## ASD snRNA-seq (Velmeshev et al, 2019)
#   # load in
#   asd <- readxl::read_xls("../../../PublicData/snRNAseq/Velmeshev2019_ASD/NIHMS1053005-supplement-Data_S4.xls", sheet = 1)
#   
#   # filter to ast DEGs
#   asd <- asd[grep("AST", asd$`Cell type`),]
#   
#   # add to annotation
#   annotation.logi$AstDisease_ASD_Velmeshev2019 <- annotation.logi$Gene %in% asd$`Gene name`
#   enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Velmeshev2019", rbind = TRUE)
  

## ASD snRNA-seq (Gandal et al, 2022)
  # read in
  asd2022 <- readxl::read_xlsx("../../../PublicData/snRNAseq/Gandal2022_ASD/41586_2022_5377_MOESM10_ESM.xlsx", sheet = "DEA_ASDvCTL_sumstats", skip = 1)
  colnames(asd2022) <- c("Celltype", "Region", "Gene", "P", "logFC", "FDR")  
  
  # filter to ast degs
  asd2022 <- asd2022[grep("ASTRO", asd2022$Celltype),] 
  
  # filter to PFC
  asd2022 <- asd2022[asd2022$Region == "PFC",] 
  
  # add to annotation
  annotation.logi$AstDisease_ASD_Gandal2022 <- annotation.logi$Gene %in% asd2022$Gene
  annotation.logi$AstDisease_ASD_Up_Gandal2022 <- annotation.logi$Gene %in% asd2022$Gene[which(asd2022$logFC > 0)]
  annotation.logi$AstDisease_ASD_Dn_Gandal2022 <- annotation.logi$Gene %in% asd2022$Gene[which(asd2022$logFC < 0)]
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Gandal2022", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Up_Gandal2022", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Dn_Gandal2022", rbind = TRUE)
  
## Multiple sclerosis snRNA-seq (Jakel et al, 2019) 
  # read in
  ms <- list(Ast1 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes1")),
             Ast2 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes2")))
  
  ms <- do.call("rbind", ms)
  ms$Cluster <- splitter(rownames(ms), "\\.", 1)
  
  # add to annotation
  annotation.logi$AstDisease_MS_Jakel2019 <- annotation.logi$Gene %in% ms$gene
  annotation.logi$AstDisease_MS_Up_Jakel2019 <- annotation.logi$Gene %in% ms$gene[which(ms$avg_logFC > 0)]
  annotation.logi$AstDisease_MS_Dn_Jakel2019 <- annotation.logi$Gene %in% ms$gene[which(ms$avg_logFC < 0)]
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_MS_Jakel2019", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_MS_Up_Jakel2019", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_MS_Dn_Jakel2019", rbind = TRUE)
  

## AD snRNA-seq (Sadick 2022)
  # read in
  ad <- list(Up = read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST6_AD_DEGs_Astro.xlsx", sheet = "Astro_upregulated_DEGs_dis"),
             Dn = read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST6_AD_DEGs_Astro.xlsx", sheet = "Astro_downregulated_DEGs_dis"))
  
  ad <- lapply(ad, function(x) {
    colnames(x)[1] <- "Gene"
    return(x)
  })
  
  ad <- do.call("rbind", ad)
  
  # add to annotation
  annotation.logi$AstDisease_AD_Sadick2022 <- annotation.logi$Gene %in% ad$Gene
  annotation.logi$AstDisease_AD_Up_Sadick2022 <- annotation.logi$Gene %in% ad$Gene[which(ad$log2FC > 0)]
  annotation.logi$AstDisease_AD_Dn_Sadick2022 <- annotation.logi$Gene %in% ad$Gene[which(ad$log2FC < 0)]
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Sadick2022", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Up_Sadick2022", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Dn_Sadick2022", rbind = TRUE)
  
  
## Glioblastoma (Zhang 2016)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_GBM_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$up_GBM_vs_Healthy, zhang_data$down_GBM_vs_Healthy)
    annotation.logi$AstDisease_GBM_Up_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$up_GBM_vs_Healthy)
    annotation.logi$AstDisease_GBM_Dn_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$down_GBM_vs_Healthy)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_GBM_Zhang2016", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_GBM_Up_Zhang2016", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_GBM_Dn_Zhang2016", rbind = TRUE)

## Peritumour (Krawczyk 2022)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_Peritumour_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Peritumour$`Gene Name`
    annotation.logi$AstDisease_Peritumour_Up_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Peritumour$`Gene Name`[which(krawczyk_data$Peritumour$`log2(Fold Change)` > 0)]
    annotation.logi$AstDisease_Peritumour_Dn_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Peritumour$`Gene Name`[which(krawczyk_data$Peritumour$`log2(Fold Change)` < 0)]
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Peritumour_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Peritumour_Up_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Peritumour_Dn_Krawczyk2022", rbind = TRUE)
    
## Epilepsy (Zhang 2016)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_Epilepsy_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$up_Epilepsy_vs_Healthy, zhang_data$down_Epilepsy_vs_Healthy)
    annotation.logi$AstDisease_Up_Epilepsy_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$up_Epilepsy_vs_Healthy)
    annotation.logi$AstDisease_Dn_Epilepsy_Zhang2016 <- annotation.logi$Gene %in% c(zhang_data$down_Epilepsy_vs_Healthy)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Epilepsy_Zhang2016", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Up_Epilepsy_Zhang2016", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Dn_Epilepsy_Zhang2016", rbind = TRUE)

## Focal Cortical Dysplasia (Krawczyk 2022)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_FocalCorticalDysplasia_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$FocalCorticalDysplasia$`Gene Name`
    annotation.logi$AstDisease_FocalCorticalDysplasia_Up_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$FocalCorticalDysplasia$`Gene Name`[which(krawczyk_data$FocalCorticalDysplasia$`log2(Fold Change)` > 0)]
    annotation.logi$AstDisease_FocalCorticalDysplasia_Dn_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$FocalCorticalDysplasia$`Gene Name`[which(krawczyk_data$FocalCorticalDysplasia$`log2(Fold Change)` < 0)]
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_FocalCorticalDysplasia_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_FocalCorticalDysplasia_Up_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_FocalCorticalDysplasia_Dn_Krawczyk2022", rbind = TRUE)
    

################################################################################################################################ #
## Astrocyte biology: Activation ----
    
 ## Results from IL-1α+TNF+C1q stimulation of iAstrocytes (Leng 2022)
  ## Read in
    leng_data <- list()
    leng_path <- "../../../PublicData/Leng2022_iAstrocyteScreen/ST1_BulkActivationDEGs.xlsx"
    
    # now: this supplementary table contains activation results for 4 iPSC-derived astrocyte lines, as well as organoid data
    leng_data$TIC_hiPSC_Leng <- read_xlsx(leng_path, sheet = "iAstrocytes")
    leng_data$TIC_hiPSC_TCW <- read_xlsx(leng_path, sheet = "TCW et al. astrocytes")
    leng_data$TIC_hiPSC_Li <- read_xlsx(leng_path, sheet = "Li et al. astrocytes")
    leng_data$TIC_hiPSC_Krencik <- read_xlsx(leng_path, sheet = "Krencik et al. astrocytes")
    # leng_data$TIC_Organoid_Barbar <- read_xlsx(leng_path, sheet = "Barbar et al.")
    
  ## Wrangle
    leng_data <- lapply(leng_data, function(x) {
      x <- x[,-c(7:8)] # removes two extraneous columns
      pCol <- grep("adj", colnames(x))
      # fcCol <- grep("FC|FoldChange", colnames(x))
      x$Hit <- (x[,pCol] < 0.05) # & (x[,fcCol] > 1)
      x <- x[which(x$Hit),]
      return(x)
    })
  
    names(leng_data) <- paste0("AstActivation_", names(leng_data))
    
    # enrichments <- enrichments[-grep("Activation", rownames(enrichments)),]
    # annotation.logi <- annotation.logi[,-grep("Activation", colnames(annotation.logi))]
    
  ## Add each to annotation dataframes
    for (j in names(leng_data)) {
      annotation.logi[,j] <- annotation.logi$Gene %in% leng_data[[j]]$gene
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE)
      
      lab_up <- gsub("hiPSC", "hiPSC_up", j)
      annotation.logi[,lab_up] <- annotation.logi$Gene %in% leng_data[[j]]$gene[which(leng_data[[j]]$log2FoldChange > 0)]
      enrichments <- run.fisher.vsHits(data.column = lab_up, rbind = TRUE)
      
      lab_dn <- gsub("hiPSC", "hiPSC_dn", j)
      annotation.logi[,lab_dn] <- annotation.logi$Gene %in% leng_data[[j]]$gene[which(leng_data[[j]]$log2FoldChange < 0)]
      enrichments <- run.fisher.vsHits(data.column = lab_dn, rbind = TRUE)
      
    }
    
  # ## An aside: how replicable are these calls across studies?
  #   leng_overlap <- annotation.logi[,c("Gene", "Hit", names(leng_data))]
  #   leng_overlap$Activation_Sum <- rowSums(leng_overlap[,names(leng_data)])
  #   
  #   # sort
  #   leng_overlap <- leng_overlap[order(-leng_overlap[,"Hit"], -leng_overlap[,"Activation_Sum"]), ]
  #   
  #   write.csv(leng_overlap, file = "Genes/Activation - Leng2022.csv")  
  #   
  # ## Replicability of activation status as a function of screen hit
  #   # tab <- table(leng_overlap$Activation_Sum, leng_overlap$Screen_Hit)
  #   # tab <- as.data.frame(tab) %>% dcast(Var1~Var2)
  #   # colnames(tab) <- c("nDatasets_SigActivation", "Screen_ns", "Screen_Hit")
  #   
  #   leng_replicability <- list()
  #   for (j in 0:4) {
  #     
  #     if (j == 0) {
  #       x <- leng_overlap$Activation_Sum == j
  #       label <- "0"
  #     } else {
  #       x <- leng_overlap$Activation_Sum >= j  
  #       label <- paste0(j, "+")
  #     }
  #     
  #     
  #     y <- leng_overlap$Hit
  #     
  #     fish <- table(x, y) %>% fisher.test()
  #     
  #     leng_replicability[[as.character(j)]] <- data.frame(GeneActivatedIn = label,
  #                                                         nGenes = sum(x),
  #                                                         Also_Screen_Hit = sum(x & y),
  #                                                         p = fish$p.value,
  #                                                         OR = fish$estimate,
  #                                                         Lower = fish$conf.int[1],
  #                                                         Upper = fish$conf.int[2])
  #     
  #     
  #   }
  #   
  #   leng_replicability <- do.call("rbind", leng_replicability)
  #   write.csv(leng_replicability, file = "Genes/Activation - Leng2022 Replicability Fisher.csv") 
    


# ## LPS activation timecourse of mouse astrocytes in vivo
#   ## From Hasel 2021
#     
#   ## Read in
#     hasel_genes <- list()
#     
#     # bulk
#     hasel_bulk <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST1_LPS_BulkTimecourse.xlsx", sheet = 2)
#     
#     get_hasel_bulkDEGs <- function(y) {
#       # genes passing padjusted threshold
#       padj <- hasel_bulk[[paste0("padj_", y, "_LPS")]]
#       padj <- gsub("E", "e", padj) 
#       padj <- as.numeric(padj)
#       pass.padj <- padj < 0.05
#       
#       # genes passing log2fc threshold
#       fc <- hasel_bulk[[paste0("l2f_", y, "_LPS")]]
#       pass.fc <- fc > 1
#       
#       # get degs
#       z <- hasel_bulk$`Gene Name`[pass.padj & pass.fc]
#       
#       # convert to human symbol
#       z <- convert_mouse2human(genes = z, path.fix = FALSE, return.vector = TRUE)
#       z <- z[-which(is.na(z))]
#       z <- unique(z)
#       
#       # return
#       return(z)
#     }
#     
#     hasel_genes$MouseLPS_3h <- get_hasel_bulkDEGs("3h")
#     hasel_genes$MouseLPS_24h <- union(get_hasel_bulkDEGs("24h_Male"), get_hasel_bulkDEGs("24h_Female"))
#     hasel_genes$MouseLPS_72h <- get_hasel_bulkDEGs("72h")
#     hasel_genes$MouseLPS_Anytime <- do.call("c", hasel_genes[grep("^MouseLPS_", names(hasel_genes))]) %>% unique()
#     hasel_genes$MouseLPS_Everytime <- intersect(hasel_genes$MouseLPS_3h, hasel_genes$MouseLPS_24h) %>% intersect(hasel_genes$MouseLPS_72h)
#     
#     
#     # rat
#     hasel_rat <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST7_Rat.xlsx", sheet = 2)
#     hasel_rat$Human <- convert_rat2human(hasel_rat$SYMBOL)
#     
#     get_hasel_ratDEGs <- function(y) {
#       # genes passing padjusted threshold
#       padj <- hasel_rat[[paste0(y, "_padj")]]
#       padj <- gsub("E", "e", padj) 
#       padj <- as.numeric(padj)
#       pass.padj <- padj < 0.05
#       
#       # genes passing log2fc threshold
#       fc <- hasel_rat[[paste0(y, "_l2f")]]
#       pass.fc <- fc > 1
#       
#       # get degs
#       z <- hasel_rat$Human[pass.padj & pass.fc]
#       z <- unique(z)
#       
#       # return
#       return(z)
#     }
#     
#     hasel_genes$Rat_IFNgamma <- get_hasel_ratDEGs("Ifng_vs_Cnt") # interferon gamma
#     hasel_genes$Rat_IFNbeta <- get_hasel_ratDEGs("Ifnb_vs_Cnt") # interferon beta
#     hasel_genes$Rat_TIC <- get_hasel_ratDEGs("TIC_vs_Cnt") # tnfa, ilf1 alpha, c1q
#     
#     # single cell
#     hasel_sc <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST6_scRNAseq_DEGs.xlsx", skip = 1)
#     hasel_genes$MouseLPS_SC <- convert_mouse2human(hasel_sc$Gene_name) %>% unique()
#     
#     # rename
#     names(hasel_genes) <- paste0("AstActivation_", names(hasel_genes), "_Hasel2021")
#     
#   ## Add each list to the dataframe
#     for (j in names(hasel_genes)) {
#       annotation.logi[,j] <- annotation.logi$Gene %in% hasel_genes[[j]]
#       enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE)
#     }