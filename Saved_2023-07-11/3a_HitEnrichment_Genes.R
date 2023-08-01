## This script annotates hit genes according to various public resources, including:

  # ontology enrichments using gprofiler2, and disease annotations from disgenet

  # phenotypes from disgenet, (mouse) ko phenotypes, human constraint databases

  # cell-type- and tissue-regulated genes  

  # astrocyte-specific differential expression in a range of disorders and phenotypes
  

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/")
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
    
  
## For hit genes only, store detailed annotation information in this dataframe
  annotation.detailed <- data.frame(Gene = sort(hit.genes))
  
## Quick load
  # load("Genes/Final.rda")

################################################################################################################################ #
## GO ----
  
## Run
  go <- gprofiler2::gost(query = hit.genes, custom_bg = bg, significant = TRUE, evcodes = TRUE, organism = "hsapiens", user_threshold = 0.05, correction_method = "fdr")
  go <- as.data.frame(go$result)
  go$parents <- sapply(go$parents, function(x) paste(x, collapse = "_"))

  
## Save
  g <- apply(go, 2, function(x) { # this code makes the dataframe suitable for csv format
    
    if (class(x) != "character") {
      return(x)
    } else {
      return(gsub(",", "_", x))
    }
    
  })
  
  write.csv(g, file = "Genes/GO.csv", quote = FALSE, row.names = FALSE)  
  

  
  
################################################################################################################################ #
## Disgenet ----
  

## Once-off authorisation
  disgenet_api_key <- get_disgenet_api_key(
    email = "i.voineagu@unsw.edu.au", 
    password = "UNSWlab123" )
  Sys.setenv(DISGENET_API_KEY = disgenet_api_key)
  
  
## Run
  # disease associations for each gene
  disgenet.g2d <- gene2disease(gene = hit.genes, 
                      database = "ALL", 
                      score = c(0.2, 1), 
                      verbose = TRUE)
  
  x <- disgenet.g2d@qresult
  x <- x[order(x$gene_symbol),]
  x <- relocate(x, "gene_symbol", "disease_name")
  write.csv(x, file = "Genes/Disgenet - Gene2Disease.csv", row.names = FALSE)
  
  
  # disease set enrichments for the list of hits
  disgenet.enrich <- disease_enrichment(entities = hit.genes, 
                                        database = "ALL", 
                                        universe = "CUSTOM", 
                                        custom_universe = bg, 
                                        verbose = TRUE)
  
  x <- disgenet.enrich@qresult
  write.csv(x, file = "Genes/Disgenet - Enrichment.csv", row.names = FALSE)
  

## Save
  save(disgenet.enrich, disgenet.g2d, file = "Genes/Disgenet.rda")
  
  
## Plot enrichment
  # function
  denrich <- function(p) {
    ggplot(p, aes(x = Description, y = -log10(FDR), fill = Ratio, size = Count)) +
      geom_point(shape = 21, colour = "black") +
      # geom_col(width = 0.1, colour = "black", fill = "black") +
      coord_flip() +
      scale_fill_viridis_c(option = "A") +
      theme_bw() +
      scale_y_continuous(limits = c(0, NA)) +
      theme(axis.title.y = invis) +
      labs()
  }
  
  # "Mental Disorders"
  p <- disgenet.enrich@qresult[grepl("Mental Disorders", disgenet.enrich@qresult$disease_class_name),]
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / length(hit.genes)
  p <- p[which(p$FDR < 0.05),]
  p <- p[1:10,] # top 10
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Genes/Disgenet Enrichment - Mental Disorders.pdf", height = 4, width = 8)
  denrich(p)
  dev.off()
  
  # "Nervous System Diseases"
  p <- disgenet.enrich@qresult[grepl("Nervous System Diseases", disgenet.enrich@qresult$disease_class_name),]
  p <- p[!(grepl("Neoplasms|Mental|Cardiovascular|Musculo", p$disease_class_name)),]
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / length(hit.genes)
  p <- p[which(p$FDR < 0.05),]
  p <- p[1:10,]
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Genes/Disgenet Enrichment - Nervous System Diseases.pdf", height = 3.5, width = 8)
  denrich(p)
  dev.off()
  
  # "Neoplasms"
  p <- disgenet.enrich@qresult[grepl("Neoplasms", disgenet.enrich@qresult$disease_class_name),]
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / length(hit.genes)
  p <- p[1:15,]
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Genes/Disgenet Enrichment - Neoplasms.pdf", height = 4, width = 8)
  denrich(p)
  dev.off()
  
## Add disgenet annotations to each gene
  annotation.detailed$Disgenet <- "."
  for (j in annotation.detailed$Gene) {
    x <- disgenet.g2d@qresult[which(disgenet.g2d@qresult$gene_symbol == j),]
    annotation.detailed$Disgenet[which(annotation.detailed$Gene == j)] <- paste(x$disease_name, collapse = "; ")
  }
  
  annotation.detailed$Disgenet <- gsub(",", "", annotation.detailed$Disgenet)
  
  
################################################################################################################################ #
## Miscellaneous gene lists ----
  
  
## SFARI
  sfari <- read.csv("../../../PublicData/SFARI-Gene_genes_10-12-2022release_10-13-2022export.csv")
  m <- match(annotation.detailed$Gene, sfari$gene.symbol)
  
  annotation.detailed$SFARI_Score <- sfari$gene.score[m]
  annotation.detailed$SFARI_Syndromic <- sfari$syndromic[m]
  
  
## Phenotype in mouse knockouts
  # access from https://www.informatics.jax.org/batch
  write.table(bg, file = "../../../PublicData/MGIPhenotypes/BackgroundGenes.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  # read in
  mgi <- read.delim("../../../PublicData/MGIPhenotypes/MGIBatchReport_20230309_050302.txt")
  mgi <- mgi[which(mgi$Term != ""),]

  # extract annotation
  annotation.detailed$MGI <- NA
  
    for (j in 1:nrow(annotation.detailed)) {
      print(j)
      g <- annotation.detailed$Gene[j]
      
      if (g %in% mgi$Input) {
        
        x <- mgi[which(mgi$Input == g),]
        annotation.detailed$MGI[j] <- paste(x$Term, collapse = "; ")
        
      } else {
        
        next
        
      }
    }
  
  annotation.detailed$MGI.Lethal <- grepl("lethal", annotation.detailed$MGI, ignore.case = TRUE)
  
  
## gnomAD
  gnomad <- read.delim("../../../PublicData/gnomAD_full_constraint_metrics.tsv")
  gnomad <- gnomad[which(gnomad$canonical == "true"),] # canonical transcript only
  m <- match(annotation.detailed$Gene, gnomad$gene)
  
  # pLI: probability of loss of function intolerance, where >0.9 is the accepted cut-off for being intolerant
  annotation.detailed$pLI <- round(gnomad$pLI[m], 3)
  
  # the new gnomad metric is observed/expected loss of function mutations. they recommend the upper bound of 90% confidence interval, so called: LOEUF (0.35 is the accepted cut-off)
  annotation.detailed$LOEUF <- gnomad$oe_lof_upper[m]  
  

################################################################################################################################ #
## Cell-type and tissue expression patterns ----
  

## Housekeeping genes, i.e. those with little variation across cells and tissues
  ## From Lin 2019: Evaluating stably expressed genes in single cells
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748759
    # the Supp table in the paper is not pre-filtered, so download from here: http://www.maths.usyd.edu.au/u/pengyi/software/scHK/scHK_human.xlsx
    
    ## Read in and process
      hkg <- readxl::read_xlsx("../../../PublicData/HousekeepingGenes/scHK_human.xlsx", sheet = 1) # this is the human list
      hkg <- hkg$`Gene Symbol`
      
      
    ## Add to universal annotations
      annotation.logi$Housekeeping_Lin2019 <- annotation.logi$Gene %in% hkg
      enrichments <- run.fisher.vsHits(data.column = "Housekeeping_Lin2019")
      
  ## From Eisenberg 2013:
    # https://www.tau.ac.il/~elieis/HKG/
    
    ## Read in
      eisenberg <- read.table("../../../PublicData/HousekeepingGenes/Eisenberg2013_Genes.txt")
      eisenberg <- eisenberg$V1
        
    ## Add to universal annotation
      annotation.logi$Housekeeping_Eisenberg2013 <- annotation.logi$Gene %in% eisenberg
      enrichments <- run.fisher.vsHits(data.column = "Housekeeping_Eisenberg2013", rbind = TRUE)
    
    
    
## Brain-specific genes
  ## Use GTEx V8, pre-transformed to median-per-tissue
    
  ## Load
    # # read in
    # gtex <- read.table("../../../PublicData/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep = "\t", skip = 2, header = TRUE)
    # 
    # # subset to used genes
    # gtex <- gtex[which(gtex$Description %in% res.final$Gene),]
    # 
    # # resolve the duplicate
    # elfn2 <- geneInfo$EnsID[which(geneInfo$Symbol == "ELFN2")]
    # gtex <- gtex[-which(gtex$Description == "ELFN2" & grepl(elfn2, gtex$Name)),] # this symbol has two rows for separate EnsIDs, so remove the EnsID which does not correspond to your records
    # 
    # # set rownames
    # rownames(gtex) <- gtex$Description
    # gtex <- gtex[,-c(1:2)]
    # 
    # # save
    # save(gtex, file = "../../../PublicData/GTEx/GTEx_V8_AllSamps.rda")
    
    # load
    load(file = "../../../PublicData/GTEx/GTEx_V8_AllSamps.rda")
    
    # metadata
    gtex_meta <- read.delim(file = "../../../PublicData/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
    rownames(gtex_meta) <- gsub("-", "\\.", gtex_meta$SAMPID)
    gtex_meta <- gtex_meta[colnames(gtex),] # same order as the expression data
    
    # remove cell-line-derived samples
    remove <- grep("Cells - ", gtex_meta$SMTSD)
    gtex_meta <- gtex_meta[-remove,]
    gtex <- gtex[,-remove]
  
  ## Compare BA9 to all non-brain tissues
    # sample filtering
    remove <- which((gtex_meta$SMTS == "Brain") & (gtex_meta$SMTSD != "Brain - Frontal Cortex (BA9)"))
    m <- gtex_meta[-remove,]
    e <- gtex[,-remove]
  
    # metadata
    isBrain <- m$SMTS == "Brain"
    
    # transform expression data
    e <- log2(e + 0.5)
    e <- t(e)
    
    # lm
    brainSpecific <- lm(e~isBrain)
    brainSpecific <- summary(brainSpecific)
    
    brainSpecific <- lapply(brainSpecific, function(z) {
      data.frame(P = as.numeric(z$coefficients["isBrainTRUE", "Pr(>|t|)"]),
                 log2fc = as.numeric(z$coefficients["isBrainTRUE", "Estimate"]))
    })
    
    brainSpecific <- do.call("rbind", brainSpecific)
    rownames(brainSpecific) <- splitter(rownames(brainSpecific), " ", 2)
    brainSpecific$FDR <- p.adjust(brainSpecific$P, method = "fdr")
    brainSpecific$BrainMarker <- brainSpecific$FDR < 0.05 & brainSpecific$log2fc > 1
    brainSpecific$Gene <- rownames(brainSpecific)
    brainSpecific <- brainSpecific[,c(5,2,1,3,4)]  
    
    write.csv(brainSpecific, "Genes/Tissue-Specificity - GTEx PFC vs non-brain.csv")
    
  ## Add to annotation
    annotation.logi$BrainSpecific_GTEx <- annotation.logi$Gene %in% brainSpecific$Gene[which(brainSpecific$BrainMarker)] 
    enrichments <- run.fisher.vsHits(data.column = "BrainSpecific_GTEx", rbind = TRUE)
    
    
## Astrocyte-specific genes, regardless of age
  ## Load pseudobulk data from Herring 2022
    load("../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/Pseudobulk_byGJS.rda")
    
  ## Filter
    # use.samps <- which(pb$Meta$nCells >= 10 & pb$Meta$Stage %in% c("Adolescence", "Adult"))
    use.samps <- which(pb$Meta$nCells >= 10)
    x <- pb$Exp[,use.samps]
    y <- pb$Meta[use.samps,]
      
  ## Process expression data
    # cpm
    x <- apply(x, 2, function(x) x / (sum(x) / 10^6))
    
    # log
    x <- log2(x + 0.5)
    
    # the rest
    x <- t(x)
  
  ## Process metadata
    y$IsAstrocyte <- y$Celltype == "Astro"
    # y$Age <- splitter(y$OriginalID, "_", 2) %>% sub("yr", "", .) %>% as.numeric()
    y$log10UMI <- log10(y$nUMI)
    
  ## Linear model
    # run
    # z <- lm(x~IsAstrocyte+Individual+Age+log10UMI, data = y)
    herring_astSpecific <- lm(x~IsAstrocyte+Individual+Stage+log10UMI, data = y)
    herring_astSpecific <- summary(herring_astSpecific)
    
    # extract effect of astrocytes
    herring_astSpecific <- lapply(herring_astSpecific, function(z) {
      data.frame(P = as.numeric(z$coefficients["IsAstrocyteTRUE", "Pr(>|t|)"]),
                 log2fc = as.numeric(z$coefficients["IsAstrocyteTRUE", "Estimate"]))
    })
  
    # return
    herring_astSpecific <- do.call("rbind", herring_astSpecific)
    rownames(herring_astSpecific) <- splitter(rownames(herring_astSpecific), " ", 2)
    herring_astSpecific$FDR <- p.adjust(herring_astSpecific$P, method = "fdr")
    herring_astSpecific$AstMarker <- herring_astSpecific$FDR < 0.05 & herring_astSpecific$log2fc > 1
    herring_astSpecific$Gene <- rownames(herring_astSpecific)
    herring_astSpecific <- herring_astSpecific[,c(5,2,1,3,4)]  
    write.csv(herring_astSpecific, file = "Genes/Astrocyte Markers - Herring 2022.csv")
    
  ## Add to annotation dataframes
    annotation.logi$AstMarkers_Herring2022 <- annotation.logi$Gene %in% herring_astSpecific$Gene[which(herring_astSpecific$AstMarker)]
    enrichments <- run.fisher.vsHits(data.column = "AstMarkers_Herring2022", rbind = TRUE)
   
     
## Astrocyte markers (foetal)
  ## Read in
    zhong_data <- read_xlsx("../../../PublicData/snRNAseq/Zhong2018_FoetalSmartSeq/41586_2018_BFnature25980_MOESM2_ESM.xlsx", skip = 4, sheet = 1)
    zhong_data <- zhong_data[which(zhong_data$cluster == "Astrocytes"),] # 353
    zhong_data <- zhong_data[which(zhong_data$avg_diff > 0),] # 326
  
  ## Add to annotation dataframes
    annotation.logi$AstMarkers_Foetal_Zhong2018 <- annotation.logi$Gene %in% zhong_data$gene
    enrichments <- run.fisher.vsHits(data.column = "AstMarkers_Foetal_Zhong2018", rbind = TRUE)
    
    

    

## Astrocyte subtypes (Sadick 2022)
  ## Read in
    subtypeMarkers <- read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST3_ctMarkers.xlsx", sheet = "LEN_so_astro_r2_DEGs")
    subtypeMarkers <- subtypeMarkers[which(subtypeMarkers$p_val_adj < 0.05),] # removes 3? How curious...
  
  ## Add to annotation dataframes
    annotation.logi$AstSubtypes_Sadick2022 <- annotation.logi$Gene %in% subtypeMarkers$gene
    enrichments <- run.fisher.vsHits(data.column = "AstSubtypes_Sadick2022", rbind = TRUE)
  
  ## Detailed annotation too
    annotation.detailed$AstSubtypes_Sadick2022 <- NA
  
    for (j in 1:nrow(annotation.detailed)) {
      print(j)
      g <- annotation.detailed$Gene[j]
      
      if (g %in% subtypeMarkers$gene) {
        
        x <- subtypeMarkers[which(subtypeMarkers$gene == g),]
        annotation.detailed$AstSubtypes_Sadick2022[j] <- paste(x$LEN_so_astro_r2_cluster, collapse = "; ")
        
      } else {
        
        next
        
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
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Herring2022", rbind = TRUE)
    # enrichments["AstMaturation_Herring2022",] <- run.fisher.vsHits(data.column = "AstMaturation_Herring2022", rbind = FALSE)
  
  ## Detailed annotation too
    annotation.detailed$AstMaturation_Herring2022 <- NA
  
    for (j in 1:nrow(annotation.detailed)) {
      print(j)
      g <- annotation.detailed$Gene[j]
      
      if (g %in% herring_trends$Symbol) {
        m <- match(g, herring_trends$Symbol)
        annotation.detailed$AstMaturation_Herring2022[j] <- herring_trends$Astro[m]
      } else {
        next
      }
     
    }
    

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
    
  ## Combine up and down
    m1 <- grep("mouse", names(zhang_data))
    zhang_data$human_vs_mouse <- do.call("c", (zhang_data[m1]))
    
    m2 <- grep("fetal", names(zhang_data))
    zhang_data$fetal_vs_adult <- do.call("c", (zhang_data[m2]))
    
    m3 <- grep("GBM", names(zhang_data))
    zhang_data$GBM_vs_Healthy <- do.call("c", (zhang_data[m3]))
    
    m4 <- grep("Epilepsy", names(zhang_data))
    zhang_data$Epilepsy_vs_Healthy <- do.call("c", (zhang_data[m4]))
      
  ## Add to annotation (foetal vs adult only)
    annotation.logi$AstMaturation_Zhang2016 <- annotation.logi$Gene %in% zhang_data$fetal_vs_adult
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Zhang2016", rbind = TRUE)

  
## Maturation: Krawczyk  
  # similar to the above study, but with a much large scope, with ~4x as many samples
    
  # again, read in all relevant data here, but only output to the annotation dataframes the relevant maturation work
    
  ## Read in all
    krawczyk_data <- list()
    krawczyk_data$Peritumour <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-3_peritumour.xlsx", skip = 1)
    krawczyk_data$FocalCorticalDysplasia <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-9_FCD.xlsx", skip = 1)
    krawczyk_data$Maturation <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Up")
    krawczyk_data$Maturation <- rbind(krawczyk_data$Maturation, read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Down"))
    krawczyk_data$Ageing <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-12_Ageing.xlsx", skip = 1)
  
  ## Extract gene symbol
    krawczyk_data <- lapply(krawczyk_data, function(x) x$`Gene Name`)
    
  ## Add to annotation (maturation only)
    annotation.logi$AstMaturation_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Maturation
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Krawczyk2022", rbind = TRUE)
    
## Ageing: Krawczyk
  # already processed above
    
  ## Add to annotation (ageing only)
    annotation.logi$AstAgeing_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Ageing
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Krawczyk2022", rbind = TRUE)
    
## Ageing: Palmer 2021
  # snRNA-seq of the ageing brain
    
  ## Read in
    ageing_palmer <- read_xlsx("../../../PublicData/snRNAseq/Palmer2021_Aging/pnas.2114326118.sd06.xlsx", sheet = "Ast")
    colnames(ageing_palmer)[1] <- "Gene"
    ageing_palmer <- ageing_palmer[which(ageing_palmer$p_val_adj < 0.05),]
    
  # note that a positive fold-change indicates higher expression in old vs. young brains.
  
  ## Add to annotation
    annotation.logi$AstAgeing_Palmer2021 <- annotation.logi$Gene %in% ageing_palmer$Gene
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Palmer2021", rbind = TRUE)
  

################################################################################################################################ #
## Astrocyte biology: Disease DEGs ----
    
## ASD snRNA-seq (Velmeshev et al, 2019)
  # load in
  asd <- readxl::read_xls("../../../PublicData/snRNAseq/Velmeshev2019_ASD/NIHMS1053005-supplement-Data_S4.xls", sheet = 1)
  
  # filter to ast DEGs
  asd <- asd[grep("AST", asd$`Cell type`),]
  
  # add to annotation
  annotation.logi$AstDisease_ASD_Velmeshev2019 <- annotation.logi$Gene %in% asd$`Gene name`
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Velmeshev2019", rbind = TRUE)
  

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
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Gandal2022", rbind = TRUE)
  
## Multiple sclerosis snRNA-seq (Jakel et al, 2019) 
  # read in
  ms <- list(Ast1 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes1")),
             Ast2 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes2")))
  
  ms <- do.call("rbind", ms)
  ms$Cluster <- splitter(rownames(ms), "\\.", 1)
  
  # add to annotation
  annotation.logi$AstDisease_MS_Jakel2019 <- annotation.logi$Gene %in% ms$gene
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_MS_Jakel2019", rbind = TRUE)
  

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
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Sadick2022", rbind = TRUE)
  
  
## Glioblastoma (Zhang 2016)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_GBM_Zhang2016 <- annotation.logi$Gene %in% zhang_data$GBM_vs_Healthy
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_GBM_Zhang2016", rbind = TRUE)

## Peritumour (Krawczyk 2022)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_Peritumour_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$Peritumour
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Peritumour_Krawczyk2022", rbind = TRUE)
    
## Epilepsy (Zhang 2016)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_Epilepsy_Zhang2016 <- annotation.logi$Gene %in% zhang_data$Epilepsy_vs_Healthy
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_Epilepsy_Zhang2016", rbind = TRUE)

## Focal Cortical Dysplasia (Krawczyk 2022)
  # data already processed above
  
  # add to annotation (GBM results only)
    annotation.logi$AstDisease_FocalCorticalDysplasia_Krawczyk2022 <- annotation.logi$Gene %in% krawczyk_data$FocalCorticalDysplasia
    enrichments <- run.fisher.vsHits(data.column = "AstDisease_FocalCorticalDysplasia_Krawczyk2022", rbind = TRUE)
    

    
    
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
    leng_data$TIC_Organoid_Barbar <- read_xlsx("leng_path", sheet = "Barbar et al.")
    
  ## Wrangle
    leng_data <- lapply(leng_data, function(x) {
      x <- x[,-c(7:8)] # removes two extraneous columns
      pCol <- grep("adj", colnames(x))
      # fcCol <- grep("FC|FoldChange", colnames(x))
      x$Hit <- (x[,pCol] < 0.05) # & (x[,fcCol] > 1)
      x <- x$gene[which(x$Hit)]
      return(x)
    })
  
    names(leng_data) <- paste0("AstActivation_", names(leng_data))
    
    # enrichments <- enrichments[-grep("Activation", rownames(enrichments)),]
    # annotation.logi <- annotation.logi[,-grep("Activation", colnames(annotation.logi))]
    
  ## Add each to annotation dataframes
    for (j in names(leng_data)) {
      annotation.logi[,j] <- annotation.logi$Gene %in% leng_data[[j]]
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE)
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
    


## LPS activation timecourse of mouse astrocytes in vivo
  ## From Hasel 2021
    
  ## Read in
    hasel_genes <- list()
    
    # bulk
    hasel_bulk <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST1_LPS_BulkTimecourse.xlsx", sheet = 2)
    
    get_hasel_bulkDEGs <- function(y) {
      # genes passing padjusted threshold
      padj <- hasel_bulk[[paste0("padj_", y, "_LPS")]]
      padj <- gsub("E", "e", padj) 
      padj <- as.numeric(padj)
      pass.padj <- padj < 0.05
      
      # genes passing log2fc threshold
      fc <- hasel_bulk[[paste0("l2f_", y, "_LPS")]]
      pass.fc <- fc > 1
      
      # get degs
      z <- hasel_bulk$`Gene Name`[pass.padj & pass.fc]
      
      # convert to human symbol
      z <- convert_mouse2human(genes = z, path.fix = FALSE, return.vector = TRUE)
      z <- z[-which(is.na(z))]
      z <- unique(z)
      
      # return
      return(z)
    }
    
    hasel_genes$MouseLPS_3h <- get_hasel_bulkDEGs("3h")
    hasel_genes$MouseLPS_24h <- union(get_hasel_bulkDEGs("24h_Male"), get_hasel_bulkDEGs("24h_Female"))
    hasel_genes$MouseLPS_72h <- get_hasel_bulkDEGs("72h")
    hasel_genes$MouseLPS_Anytime <- do.call("c", hasel_genes[grep("^MouseLPS_", names(hasel_genes))]) %>% unique()
    hasel_genes$MouseLPS_Everytime <- intersect(hasel_genes$MouseLPS_3h, hasel_genes$MouseLPS_24h) %>% intersect(hasel_genes$MouseLPS_72h)
    
    
    # rat
    hasel_rat <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST7_Rat.xlsx", sheet = 2)
    hasel_rat$Human <- convert_rat2human(hasel_rat$SYMBOL)
    
    get_hasel_ratDEGs <- function(y) {
      # genes passing padjusted threshold
      padj <- hasel_rat[[paste0(y, "_padj")]]
      padj <- gsub("E", "e", padj) 
      padj <- as.numeric(padj)
      pass.padj <- padj < 0.05
      
      # genes passing log2fc threshold
      fc <- hasel_rat[[paste0(y, "_l2f")]]
      pass.fc <- fc > 1
      
      # get degs
      z <- hasel_rat$Human[pass.padj & pass.fc]
      z <- unique(z)
      
      # return
      return(z)
    }
    
    hasel_genes$Rat_IFNgamma <- get_hasel_ratDEGs("Ifng_vs_Cnt") # interferon gamma
    hasel_genes$Rat_IFNbeta <- get_hasel_ratDEGs("Ifnb_vs_Cnt") # interferon beta
    hasel_genes$Rat_TIC <- get_hasel_ratDEGs("TIC_vs_Cnt") # tnfa, ilf1 alpha, c1q
    
    # single cell
    hasel_sc <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST6_scRNAseq_DEGs.xlsx", skip = 1)
    hasel_genes$MouseLPS_SC <- convert_mouse2human(hasel_sc$Gene_name) %>% unique()
    
    # rename
    names(hasel_genes) <- paste0("AstActivation_", names(hasel_genes), "_Hasel2021")
    
  ## Add each list to the dataframe
    for (j in names(hasel_genes)) {
      annotation.logi[,j] <- annotation.logi$Gene %in% hasel_genes[[j]]
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE)
    }
  
################################################################################################################################ #
## Save ----
  
    
write.csv(enrichments, file = "Genes/Final - Enrichments.csv")
write.csv(annotation.logi, file = "Genes/Final - Annotation Logical.csv", row.names = FALSE)
write.csv(annotation.detailed, file = "Genes/Final - Annotation Detailed.csv", row.names = FALSE)

save(enrichments, annotation.logi, annotation.detailed, file = "Genes/Final.rda")


################################################################################################################################ #
## Visualise ----


## A heatmap
  x <- annotation.logi
  rownames(x) <- x$Gene
  
  # filter columns
  keep.cols <- c("Hit",
                 "AstMarkers_Herring2022",
                 "AstMarkers_Foetal_Zhong2018",
                 "AstSubtypes_Sadick2022",
                 "AstAgeing_Krawczyk2022",
                 "AstAgeing_Palmer2021",
                 "AstMaturation_Zhang2016",
                 "AstMaturation_Krawczyk2022",
                 "AstActivation_iAstro_Leng2022",
                 "AstActivation_hiPSC_TCW_Leng2022",
                 "AstActivation_hiPSC_Li_Leng2022",
                 "AstActivation_hiPSC_Krenciko_Leng2022",
                 "AstDisease_AD_Sadick2022",
                 "AstDisease_MS_Jakel2019",
                 "AstDisease_ASD_Gandal2022",
                 "AstDisease_GBM_Zhang2016",
                 "AstDisease_Peritumour_Krawczyk2022",
                 "AstDisease_Epilepsy_Zhang2016")
  x <- x[,keep.cols]
  x[,-1] <- apply(x[,-1], 2, as.numeric)
  
  # clean column names
  colnames(x) <- gsub("^Ast", "", colnames(x)) %>%
    gsub("_Leng2022", "", .) %>%
    gsub("hiPSC_", "", .) %>%
    gsub("Disease_", "", .) %>%
    gsub("Krenciko", "Krencik", .) %>%
    gsub("Activation", "Activation", .) %>%
    gsub("_Zhong", "", .) %>%
    gsub("_", "\n", .) %>%
    splitter(., "2", 1) %>%
    gsub("iAstro_", "", .) 
    
  
  
## Plot
  # row side colours to label hits
  rsc <- as.factor(x$Hit)
  levels(rsc) <- c("grey90", "firebrick1")
  rsc <- as.character(rsc)
  
  # plot all genes
  pdf(file = "Genes/Heatmap of Logical - All Genes.pdf", height = 20, width = 5)
  gplots::heatmap.2(as.matrix(x[,-1]), RowSideColors = rsc, trace = "none", col = carto_pal(7, "ag_GrnYl"), margins = c(12, 5))
  dev.off()  
  
  # plot hit genes
  pdf(file = "Genes/Heatmap of Logical - Hit Genes.pdf", height = 20, width = 5)
  gplots::heatmap.2(as.matrix(x[which(x$Hit),-1]), trace = "none", col = carto_pal(7, "ag_GrnYl"), margins = c(10, 5))
  dev.off()  
  