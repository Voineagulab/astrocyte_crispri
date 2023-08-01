setwd("~/")
library(tidyverse)
library(readxl)

## Gene annotation
  BiocManager::install("EnsDb.Hsapiens.v86")
  library(EnsDb.Hsapiens.v86)
  
  # get list
  edb <- EnsDb.Hsapiens.v86 
  biotypes <- genes(edb, columns = c("gene_name", "gene_biotype"))
  biotypes <- biotypes@elementMetadata %>% as.data.frame()

  # filter to main types
  # keep.types <- c("lincRNA", "protein_coding")
  # biotypes <- biotypes[which(biotypes$gene_biotype %in% keep.types),]


## Sadick brain subtype markers
  sadick <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST3_ctMarkers.xlsx", sheet = "LEN_so_astro_r2_DEGs")
  sadick <- sadick[which(sadick$p_val_adj < 0.05),]
  m <- match(sadick$gene, biotypes$gene_name)
  sadick$Biotype <- biotypes$gene_biotype[m]
  
  sadick_lnc <- sadick[which(sadick$Biotype == "lincRNA"),]
  
  # can you visualise these genes?
  
  # For example, astrocyte clusters 0, 4, and 8 express unique sets of transcripts involved in synapse assembly, organization, and transmission (cluster 0: EGFR, LRRC4C, and EPHB1; cluster 4: DCLK1, NTNG1, and several semaphorins; cluster 8: EPHA4, AKAP12, and NLGN4X). 
  # Astrocyte clusters 4 and 6 highly express transcripts involved in glutamate signaling (GRIA1, GRIK4, and SHISA6). Clusters 2 and 5 express transcripts involved in extracellular matrix organization (cluster 2: ADAMTSL3 and L3MBTL4; cluster 5: ADAMTSL3 and FBN1), 
  # whereas cluster 5 also expresses transcripts involved in actin cytoskeletal organization (SORBS1 and SPIRE1). 
  # The inclusion of ADAMTSL3 in clusters 2 and 5 may point to a protective role of (some) AD or aged astrocytes, as their protective role has been reported in ischemia and cerebrovascular integrity in the APP/PS1 mouse model of AD (Cao et al., 2019). 
  # In contrast, cluster 3 astrocytes express transcripts involved in acute inflammatory responses (e.g., SERPINA3, C3, and OSMR) that we have reported on previously in both mice (Hasel et al., 2021; Liddelow et al., 2017) and humans (Barbar et al., 2020). 
  # Astrocyte cluster 1 is highly enriched for transcripts involved in oxidative stress (PSAP, COX1, and ND1/3) and associated with Aβ trafficking (e.g., APOE and CLU) and processing (e.g., ITM2B/2C). The inclusion of AD-risk genes, APOE and CLU, with integral membrane protein (ITM2B/2C) genes associated with cerebral amyloid angiopathy (Nelson et al., 2013; Vidal et al., 1999) in the same astrocyte cluster suggests a putative interaction. 
  # Astrocyte clusters 1 and 6 are both enriched in a number of metallothioneins and other transcripts involved in response to metal ions. 
  # Finally, cluster 7 expresses transcripts associated with apoptotic signaling and response to DNA damage.
  
  # Thus:
    # Protective: 2 and 5
    # Activated: 3
    # Stress: 1 

  # Per GO
    # 0: synapse regulations and organisation
    # 1: Respose to metal iions
    # 2: None
    # 3: vasculate, inflammation, signal transduction
    # 4: cell growth, regulation of ion transporters
    # 5: ECM, neuroinflammation, microglial cell actication, response to axon injury, glial activation
    # 6: metal ions
    # 7: apoptosis
    # 8: synpase regulation

## Leng activation DEGs
  leng <- list()
  leng_path <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Leng2022_iAstrocyteScreen/ST1_BulkActivationDEGs.xlsx"
  
  # now: this supplementary table contains activation results for 4 iPSC-derived astrocyte lines, as well as organoid data
  leng$Leng <- read_xlsx(leng_path, sheet = "iAstrocytes")
  leng$TCW <- read_xlsx(leng_path, sheet = "TCW et al. astrocytes")
  leng$Li <- read_xlsx(leng_path, sheet = "Li et al. astrocytes")
  leng$Krencik <- read_xlsx(leng_path, sheet = "Krencik et al. astrocytes")
  # leng$TIC_Organoid_Barbar <- read_xlsx(leng_path, sheet = "Barbar et al.")
  
  # annotate
  leng <- lapply(leng, function(x) {
    # hit
    x$SigP <- x$padj < 0.05
    x$SigFC <- abs(x$log2FoldChange) > 1
    x$Sig <- x$SigP & x$SigFC
    
    # biotype
    m <- match(x$gene, biotypes$gene_name)
    x$Biotype <- biotypes$gene_biotype[m]
    
    # return
    return(x[,-8])
  })


  ## Percentage hits
    sapply(leng, function(x) sum(x$Sig) / nrow(x))
    # Leng       TCW        Li   Krencik  
    # 0.09212063          0.09241109          0.18684952          0.27055926           
    
    sapply(leng, function(x) sum(x$SigP) / nrow(x))
    # Leng       TCW        Li   Krencik  
    #       0.2808958           0.2784934           0.3492102           0.4328708            
    
  ## Ratios
    table(leng$Leng$Sig, leng$Leng$Biotype == "lincRNA") %>% fisher.test() # OR = 0.83 ns
    table(leng$TCW$Sig, leng$TCW$Biotype == "lincRNA") %>% fisher.test() # OR = 1 ns
    table(leng$Li$Sig, leng$Li$Biotype == "lincRNA") %>% fisher.test() # OR = 0.38 e-7
    table(leng$Krencik$Sig, leng$Krencik$Biotype == "lincRNA") %>% fisher.test() # OR = 0.59 e-4
    
    
  ## Get lincRNA
    leng_lnc <- lapply(leng, function(x) x[which(x$Biotype == "lincRNA"),])
    leng_lnc <- do.call("rbind", leng_lnc)
    leng_lnc$Protocol <- strsplit(rownames(leng_lnc), "\\.") %>% sapply("[", 1)
    leng_lnc <- leng_lnc[,c("Protocol", "gene", "baseMean", "log2FoldChange", "padj", "Sig", "Biotype")]
    
    
## Leng subtypes
  iras <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/Leng2022_iAstrocyteScreen/ST4_ClusterMarkers.xlsx", sheet = "DEGs between IRAS2 vs. IRAS1")
  iras <- as.data.frame(iras)
  m <- match(iras$gene, biotypes$gene_name)
  iras$Biotype <- biotypes$gene_biotype[m]
  iras <- iras[which(iras$Biotype == "lincRNA"),]

  
## Astrocyte maturation
  mat <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/SuppTable3_DE.xlsx", skip = 2, sheet = "A Astro")
  mat <- mat[,c("gene_name", "trend_class", "F.p.value")]
  m <- match(mat$gene_name, biotypes$gene_name)
  mat$Biotype <- biotypes$gene_biotype[m]
  mat <- mat[which(mat$Biotype == "lincRNA"),]
  
  
## Foetal astrocyte markers
  foetal <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/snRNAseq/Zhong2018_FoetalSmartSeq/41586_2018_BFnature25980_MOESM2_ESM.xlsx", skip = 4, sheet = 1)
  foetal <- foetal[which(foetal$cluster == "Astrocytes"),] # 353
  foetal <- foetal[which(foetal$avg_diff > 0),] # 326
  m <- match(foetal$gene, biotypes$gene_name)
  foetal$Biotype <- biotypes$gene_biotype[m]
  foetal <- foetal[which(foetal$Biotype == "lincRNA"),]
  
  
## Astrocyte markers
  adult <- read.csv(file = "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Genes/Astrocyte Markers - Herring 2022.csv")
  adult <- adult[which(adult$AstMarker),]
  m <- match(adult$Gene, biotypes$gene_name)
  adult$Biotype <- biotypes$gene_biotype[m]
  adult <- adult[which(adult$Biotype == "lincRNA"),]
  

## What percentage are expressed in our scRNA-seq and bulk RNA-seq?