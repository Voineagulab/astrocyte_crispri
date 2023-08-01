################################################################################################################################ #
## Setup ----

## Directories
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(tidyverse)
  library(rcartocolor)
  invis <- element_blank()

## Load NHA data
  source("../../Scripts/Functions.R")
  load("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Data/Preprocessed/NHA Pooled (Final).rda")
  load("../../../PreprocessedData/Bulk_ATAC_RNAseq/RESULTS/RNAseq/geneData.rda")

## Get expression: sc
  n <- list()
  n$SC_Norm <- rowMeans(nha@assays$RNA@data[,which(nha$AnyGuide)]) # normalised by Seurat
  
  n$SC_PB <- rowSums(nha@assays$RNA@counts[,which(nha$AnyGuide)]) # pseudobulk 
  n$SC_PB <- n$SC_PB / (sum(n$SC_PB) / 10^6) # pseudobulk cpm
  n$SC_PB <- log2(n$SC_PB + 1) # pseudobulk log2 cpm
  
  # filter  
  keep <- which(n$SC_Norm > 2^-6)
  n$SC_Norm <- n$SC_Norm[keep] # for our data, same threshold as always
  n$SC_PB <- n$SC_PB[keep]
  
## Get expression: bulk
  n$Bulk <- geneData$geneCPM
  n$Bulk <- n$Bulk[which(n$Bulk$NHA > 1),] # threshold
  n$Bulk <- setNames(n$Bulk$NHA, n$Bulk$Symbol)
  n$Bulk <- log2(n$Bulk + 1)
  
## Common genes
  common <- intersect(names(n$Bulk), names(n$SC_PB))
  n <- lapply(n, function(x) x[common])
  n <- do.call("cbind", n) %>% as.data.frame()
  
  
################################################################################################################################ #
## Comparison to Herring 2022 ----

## Rename for ease
  e <- n
  
## Load Herring 2022 data
  load("/mnt/Data0/PROJECTS/GWAS_Enrichment/SAM/data/cpm_all_df.rda")
  samps <- splitter(colnames(cpm_all_df), "\\.\\.", 1) %>% unique()
  herring <- list()
  for (j in samps) {
    print(j)
    g <- grep(j, colnames(cpm_all_df))
    if (length(g) == 1) {
      herring[[j]] <- cpm_all_df[,g]
    } else {
      herring[[j]] <- rowMeans(cpm_all_df[,grep(j, colnames(cpm_all_df))])
    }
  }
  herring <- do.call("cbind", herring) %>% as.data.frame()

#   herring <- cpm_all_df
#   colnames(herring) <- splitter(colnames(herring), "\\.\\.", 1)

## Rename columns
  x <- str_count(colnames(herring), "_")
  colnames(herring)[(x == 2)] <- sub("_", "-", colnames(herring)[(x == 2)])
  
## Expression filter
  herring <- herring[which(rowSums(herring > 1) > (ncol(herring) / 5)),] # above 1 cpm in at least 20% of samples

## Bind to our data
  common <- intersect(rownames(herring), rownames(e))
  herring <- herring[common,]
  e <- e[common,]
  
## Subset Herring to samples with at least 90% of these genes expressed
  herring <- herring[-grep("Vas", colnames(herring))]
  herring <- herring[,-rem]

## Log normalise
  herring <- log2(herring+1)
  

  

## Correlate
  cor <- as.data.frame(cor(herring, e), method = "s")
  # cor$ID <- splitter(rownames(cor), "\\.", 1)
  cor$ID <- rownames(cor)
  cor$Celltype <- splitter(cor$ID, "_", 1)
  cor$Stage <- splitter(splitter(cor$ID, "_", 2), "\\.\\.", 1)
  cor$Stage <- factor(cor$Stage, levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
  # cor$Celltype[grep("L2|L3|L4|L5", cor$Celltype)] <- "Exc"
  # cor$Celltype[grep("ID2|LAMP|PV|SST|VIP", cor$Celltype)] <- "Inh"  
  


## Plot
  library(ggsci)
  pdf(file = "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Cell-type Validation - Stats To Herring.pdf", height = 4, width = 8)
  ggplot(cor, aes(x = Celltype, fill = Stage, y = SC_PB)) +
    # geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    # geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
    scale_fill_lancet() +
    geom_col(position = "dodge") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = invis, panel.grid.major.x = invis, axis.line.y = element_line(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "Spearman Correlation to NHA Pseudobulk")
  
   ggplot(cor, aes(x = Celltype, fill = Stage, y = Bulk)) +
    # geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    # geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
    scale_fill_lancet() +
    geom_col(position = "dodge") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = invis, panel.grid.major.x = invis, axis.line.y = element_line(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "Spearman Correlation to NHA Bulk")
  dev.off()



## Integrative UMAP


## SingleR  
  library(SingleR)
  library(scuttle)
  
  # make the test data
  x <- normalizeCounts(nha@assays$RNA@counts[common,which(nha$AnyGuide)], log = TRUE)
  x <- as.data.frame(x)
  x <- apply(x, 2, as.numeric)
  rownames(x) <- common
  
  # make a reference
  y <- apply(herring, 2, as.numeric)
  rownames(y) <- common
  
  # run
  # start at 8:55, finish at 9:01
  sr <- SingleR(test = x,
                ref = y,
                labels = colnames(y),
                clusters = NULL)

  sr2 <- SingleR(test = x,
                ref = y,
                genes = "all",
                labels = colnames(y),
                clusters = NULL)
  
  sr3 <- SingleR(test = x,
                ref = y,
                tune.thresh = 0.01,
                labels = colnames(y),
                clusters = NULL)
    
  table(sr3$first.labels)
  table(sr3$pruned.labels)
  
################################################################################################################################ #
## Comparison to iAstrocytes (Leng 2022) ----
