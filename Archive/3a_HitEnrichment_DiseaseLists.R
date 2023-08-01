## This script annotates hit genes according to various public resources, including:

  # gene ontology
  # disease annotations from disgenet
  # (mouse) ko phenotype annotations
  # constraint within human populations
  # astrocyte-specific differential expression in a range of disorders and phenotypes
  # cell-type- and tissue-regulated genes



################################################################################################################################ #
## Setup ----


## Generic
# rm(list = ls()); gc()
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
  annotation <- data.frame(Gene = bg, Hit = bg %in% hit.genes)
  annotation <- annotation[order(!(annotation$Hit), annotation$Gene),]
  
  # a dataframe to store fisher test results for overenrichments
  
  
  
  
  
  # (and a function to fill this)
  run.fisher.vsHits <- function(data = annotation, data.column) {
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
                      Upper = f$conf.int[2])  
    
    return(out)
  }
    

  
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
  


################################################################################################################################ #
## Housekeeping genes ----
  
## From Lin 2019: Evaluating stably expressed genes in single cells
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748759
  # the Supp table in the paper is not pre-filtered, so download from here: http://www.maths.usyd.edu.au/u/pengyi/software/scHK/scHK_human.xlsx
  
## Read in and process
  hkg <- readxl::read_xlsx("../../../PublicData/HousekeepingGenes/scHK_human.xlsx", sheet = 1) # this is the human list
  hkg <- hkg$`Gene Symbol`
  
## Enrichment test
  hkg.test <- data.frame(Symbol = bg)
  hkg.test$Hit <- hkg.test$Symbol %in% hit.genes
  hkg.test$HKG <- hkg.test$Symbol %in% hkg
  
  table(hkg.test$Hit, hkg.test$HKG) %>% fisher.test()
  
################################################################################################################################ #
## Astrocyte-specific genes ----
  
## Setup dataframe
  AstMarkers <- data.frame(Symbol = bg,
                           EnsID = NA,
                           Hit = bg %in% hit.genes)
  
  m <- match(AstMarkers$Symbol, geneInfo$Symbol)
  AstMarkers$EnsID <- geneInfo$EnsID[m] # add EnsID because that's what the signature uses...  
  
## Load and process
  load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/Annotation/sigsBrain.rda") # brain signatures
  sigsBrain <- sigsBrain[c("IP", "DM", "VL", "NG", "CA", "LK")] # subset to human adult brain tissue samples
  markerCount <- lapply(sigsBrain, function(a) { 
    a <- a[intersect(rownames(a), spec$EnsID),]
    b <- apply(a, 1, which.max)
    names(b)[which(b == grep("Astrocytes", colnames(a)))]
  })
  
  markerCount <- do.call("c", markerCount)
  markerCount <- table(markerCount) # number of datasets in which the tf is ast-highest (zero excluded)
  
## Add number of times each gene is a marker
  m <- match(AstMarkers$EnsID, names(markerCount))
  AstMarkers$AstSignatureTop <- markerCount[m]
  AstMarkers$AstSignatureTop[which(is.na(AstMarkers$AstSignatureTop))] <- 0
  
## Define a cut-off, and stats
  table(AstMarkers$Hit, AstMarkers$AstSignatureTop >= 2) %>% fisher.test() # ms
  AstMarkers$Specific <- AstMarkers$AstSignatureTop >= 2
  
## Save
  write.csv(AstMarkers, file = "Genes/Ast Markers.csv")
  
## Plot
  p <- table(AstMarkers$Hit, AstMarkers$AstSignatureTop)
  p <- p / rowSums(p) # normalise to total number in category
  
  p <- as.data.frame(p)
  colnames(p) <- c("Hit", "MarkerCount", "Freq")
  levels(p$Hit) <- c("Non-hit Gene", "Hit Gene")
  
  pdf(file = "Genes/Ast Markers.pdf", height = 3, width = 4)
  ggplot(p, aes(x = MarkerCount, y = Freq, fill = Hit)) +
    geom_col(position = "dodge", width = 0.7) +
    theme_bw() +
    scale_fill_lancet() +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.7)) +
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = 2.5, linetype = 2) +
    labs(y = "Fraction of Gene Set", x = "Number of Datasets Where Gene\nis Highest Expressed in Ast")
  dev.off()
  
  
   
  
  
    
################################################################################################################################ #
## snRNA-seq DEGs ----
  
  
## ASD (Velmeshev et al, 2019)
  # load in
  asd <- readxl::read_xls("../../../../PublicData/snRNAseq/Velmeshev2019_ASD/NIHMS1053005-supplement-Data_S4.xls", sheet = 1)
  
  # filter to ast DEGs
  asd <- asd[grep("AST", asd$`Cell type`),]
  
  # match
  table(hCore %in% asd$`Gene name`) # nada....
  table(hPerm %in% asd$`Gene name`) # nada....
  

## ASD (Gandal et al, 2022)
  # read in
  asd2022 <- readxl::read_xlsx("../../../../PublicData/snRNAseq/Gandal2022_ASD/41586_2022_5377_MOESM10_ESM.xlsx", sheet = "DEA_ASDvCTL_sumstats", skip = 1)
  colnames(asd2022) <- c("Celltype", "Region", "Gene", "P", "logFC", "FDR")  
  
  # filter to ast degs
  asd2022 <- asd2022[grep("ASTRO", asd2022$Celltype),] # 8675?
  
  # filter to PFC
  asd2022 <- asd2022[asd2022$Region == "PFC",] # 8675?
  
  # overlap
  asd2022.overlap <- asd2022[which(asd2022$Gene %in% hPerm),]
  
## Multiple sclerosis (Jakel et al, 2019) 
  # read in
  ms <- list(Ast1 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes1")),
             Ast2 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes2")))
  
  ms <- do.call("rbind", ms)
  ms$Cluster <- splitter(rownames(ms), "\\.", 1)
  
  # overlap
  table(hPerm %in% ms$gene) # woo, seven!
  
  ms.overlap <- ms[which(ms$gene %in% hPerm),]
  

## MDD (Nagy 2020)
  # read in
  mdd  <- list(Ast2 = read_xlsx("../../../PublicData/snRNAseq/Nagy2020_MDD/41593_2020_621_MOESM3_ESM.xlsx", sheet = 5),
             Ast3 = read_xlsx("../../../PublicData/snRNAseq/Nagy2020_MDD/41593_2020_621_MOESM3_ESM.xlsx", sheet = 6))
  
  mdd <- lapply(mdd, function(x) {
    colnames(x) <- x[2,]
    x <- x[-c(1,2),]
  })
  
  mdd <- do.call("rbind", mdd)
  mdd$Cluster <- splitter(rownames(mdd), "\\.", 1)
  mdd <- mdd[which(mdd$p.adjust < 0.05),] # well, looks like there's nothing DE below FDR 0.1...
  
## Ageing (Palmer 2021)
  # read in
  ageing <- read_xlsx("../../../PublicData/snRNAseq/Palmer2021_Aging/pnas.2114326118.sd06.xlsx", sheet = "Ast")
  colnames(ageing)[1] <- "Gene"
  ageing <- ageing[which(ageing$p_val_adj < 0.05),]
  
  # note that a positive fold-change indicates higher expression in old vs. young brains.
  
  # overlap
  table(hPerm %in% ageing$Gene) # 6, woo!
  ageing.overlap <- ageing[which(ageing$Gene %in% hPerm),]
  

## AD (Sadick 2022)
  # read in
  ad <- list(Up = read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST6_AD_DEGs_Astro.xlsx", sheet = "Astro_upregulated_DEGs_dis"),
             Dn = read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST6_AD_DEGs_Astro.xlsx", sheet = "Astro_downregulated_DEGs_dis"))
  ad <- lapply(ad, function(x) {
    colnames(x)[1] <- "Gene"
    return(x)
  })
  ad <- do.call("rbind", ad)
  
  # overlap
  ad.overlap <- ad[which(ad$Gene %in% hPerm),]
  
## Astrocyte subtypes (Sadick 2022)
  # read in
  subtypeMarkers <- read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST3_ctMarkers.xlsx", sheet = "LEN_so_astro_r2_DEGs")
  subtypeMarkers <- subtypeMarkers[which(subtypeMarkers$p_val_adj < 0.05),] # removes 3? How curious...
  
  # overlap
  subtypeMarkers.overlap <- subtypeMarkers[which(subtypeMarkers$gene %in% hPerm),]
  
## Summarise
  sum <- data.frame(Gene = hPerm,
                    Subtype = hPerm %in% subtypeMarkers$gene,
                    AD = hPerm %in% ad$Gene,
                    Ageing = hPerm %in% ageing$Gene,
                    MS = hPerm %in% ms$gene,
                    ASD = hPerm %in% asd2022$Gene)
  write.csv(sum, file = "snRNAseq Overlaps Binary.csv")
  
################################################################################################################################ #
## Miscellaneous gene lists ----
  
  
## SFARI
  sfari <- read.csv("../../../../PublicData/SFARI-Gene_genes_10-12-2022release_10-13-2022export.csv")
  sfari <- sfari[which(sfari$gene.score == 1 | sfari$syndromic == 1),]
  
  table(hit.genes %in% sfari$gene.symbol) # just the one...
  
  sfari.hit <- sfari[which(sfari$gene.symbol %in% hPerm),]
  
  write.csv(sfari.hit, file = "SFARI.csv")
  
  
################################################################################################################################ #
## Brain markers ----
  
## Use GTEx V8, pre-transformed to median-per-tissue
  # read in
  gtex <- read.table("../../../PublicData/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", sep = "\t", skip = 2, header = TRUE)
  
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
  
  # genes
  gtex_genes <- gtex[,1:2]
  rownames(gtex) <- gtex$Name
  gtex <- gtex[,-c(1,2)]
  
  # save
  # save(gtex, file = "../../../PublicData/GTEx/GTEx_V8_AllSamps.rda")
  
## Model
  # setup expression data
  gtex <- log2(gtex + 0.5)
  gtex <- t(gtex)

  # remove cells rather than tissues
  remove <- grep("Cells", rownames(gtex))
  gtex <- gtex[-remove,]
  
  # metadata
  isBrain <- grepl("^Brain", rownames(gtex))
  
  # run
  brainSpecific <- lm(gtex~isBrain)
  brainSpecific <- summary(brainSpecific)
  
  brainSpecific <- lapply(brainSpecific, function(z) {
    data.frame(P = as.numeric(z$coefficients["isBrainTRUE", "Pr(>|t|)"]),
               log2fc = as.numeric(z$coefficients["isBrainTRUE", "Estimate"]))
  })
    
  brainSpecific <- do.call("rbind", brainSpecific)
  rownames(brainSpecific) <- splitter(rownames(brainSpecific), " ", 2)
  brainSpecific$FDR <- p.adjust(brainSpecific$P, method = "fdr")
  brainSpecific$BrainMarker <- brainSpecific$FDR < 0.05 & brainSpecific$log2fc > 1
  brainSpecific$Gene <- gtex_genes$Description
  brainSpecific <- brainSpecific[,c(5,2,1,3,4)]  
  
  write.csv(brainSpecific, file = "GTEx - Brain-specific (Median-derived).csv")
  
  brainSpecific <- brainSpecific[which(brainSpecific$Gene %in% res.final$Gene),]
  elfn2 <- geneInfo$EnsID[which(geneInfo$Symbol == "ELFN2")]
  brainSpecific <- brainSpecific[-which(brainSpecific$Gene == "ELFN2" & grepl(elfn2, rownames(brainSpecific))),] # this symbol has two rows for separate EnsIDs, so remove the EnsID which does not correspond to your records
  
## Now run stats
  brainSpecific$Hit <- brainSpecific$Gene %in% hit.genes
  
  table(brainSpecific$Hit, brainSpecific$BrainMarker) %>% fisher.test()
  
  
## Repeat the above, but use brain cortex and brain frontal cortex versus the rest, unaggregated
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
    a <- Sys.time()
    brainSpecific <- lm(e~isBrain)
    brainSpecific <- summary(brainSpecific)
    b  <- Sys.time()
    b-a
    
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
    
    write.csv(brainSpecific, "GTEx - Brain-specific (PFC vs non-brain).csv")
    

################################################################################################################################ #
## Astrocyte activation states ----
  
    
## Results from IL-1α+TNF+C1q stimulation of iAstrocytes (Leng 2022)
  ## Read in
    leng_data <- list()
    leng_path <- "../../../PublicData/Leng2022_iAstrocyteScreen/ST1_BulkActivationDEGs.xlsx"
    
    # now: this supplementary table contains activation results for 4 iPSC-derived astrocyte lines, as well as organoid data
    leng_data$iAstro <- read_xlsx("../../../PublicData/Leng2022_iAstrocyteScreen/ST1_BulkActivationDEGs.xlsx", sheet = "iAstrocytes")
    leng_data$hiPSC_TCW <- read_xlsx("../../../PublicData/Leng2022_iAstrocyteScreen/ST1_BulkActivationDEGs.xlsx", sheet = "TCW et al. astrocytes")
    leng_data$hiPSC_Li <- read_xlsx("../../../PublicData/Leng2022_iAstrocyteScreen/ST1_BulkActivationDEGs.xlsx", sheet = "Li et al. astrocytes")
    leng_data$hiPSC_Krencik <- read_xlsx("../../../PublicData/Leng2022_iAstrocyteScreen/ST1_BulkActivationDEGs.xlsx", sheet = "Krencik et al. astrocytes")
    
  ## Wrangle
    leng_data <- lapply(leng_data, function(x) {
      x <- x[,-c(7:8)] # removes two extraneous columns
      pCol <- grep("adj", colnames(x))
      x$Hit <- x[,pCol] < 0.05
      x <- x$gene[which(x$Hit)]
      return(x)
    })
    
    names(leng_data) <- paste0("Activation_", names(leng_data))
    
  ## Annotate our screened genes
    leng_overlap <- data.frame(Screen_Gene = bg)
    leng_overlap$Screen_Hit <- leng_overlap$Screen_Gene %in% hit.genes
    
    for (j in names(leng_data)) {
      leng_overlap[,j] <- leng_overlap$Screen_Gene %in% leng_data[[j]]
    }
  
    # how replicable are these calls across studies?
    leng_overlap$Activation_Sum <- rowSums(leng_overlap[,names(leng_data)])
    
    # sort
    leng_overlap <- leng_overlap[order(-leng_overlap[,"Screen_Hit"], -leng_overlap[,"Activation_Sum"]), ]
    
    write.csv(leng_overlap, file = "Activation - Leng2022.csv")  
    
  ## Overrepresentation within each dataset
    leng_summary <- apply(leng_overlap[,names(leng_data)], 2, function(x) {
      y <- leng_overlap$Screen_Hit
      
      fish <- table(x, y) %>% fisher.test()
      
      data.frame(nActivated = sum(x),
                 nActivated_Hit = sum(x & y),
                 p = fish$p.value,
                 OR = fish$estimate,
                 Lower = fish$conf.int[1],
                 Upper = fish$conf.int[2])
    })
    
    leng_summary <- do.call("rbind", leng_summary)
    write.csv(leng_summary, file = "Activation - Leng2022 Fisher.csv")
    
  ## Replicability of activation status as a function of screen hit
    # tab <- table(leng_overlap$Activation_Sum, leng_overlap$Screen_Hit)
    # tab <- as.data.frame(tab) %>% dcast(Var1~Var2)
    # colnames(tab) <- c("nDatasets_SigActivation", "Screen_ns", "Screen_Hit")
    
    leng_replicability <- list()
    for (j in 0:4) {
      
      if (j == 0) {
        x <- leng_overlap$Activation_Sum == j
        label <- "0"
      } else {
        x <- leng_overlap$Activation_Sum >= j  
        label <- paste0(j, "+")
      }
      
      
      y <- leng_overlap$Screen_Hit
      
      fish <- table(x, y) %>% fisher.test()
      
      leng_replicability[[as.character(j)]] <- data.frame(GeneActivatedIn = label,
                                                          nGenes = sum(x),
                                                          Also_Screen_Hit = sum(x & y),
                                                          p = fish$p.value,
                                                          OR = fish$estimate,
                                                          Lower = fish$conf.int[1],
                                                          Upper = fish$conf.int[2])
      
      
    }
    
    leng_replicability <- do.call("rbind", leng_replicability)
    write.csv(leng_replicability, file = "Activation - Leng2022 Replicability Fisher.csv") 
    
  ## Plots
    # for enrichments within each individual study
    p <- leng_summary
    p$Protocol <- gsub("hiPSC_", "", rownames(p)) %>% splitter(., "_", 2) %>% factor()
    p$Protocol <- factor(p$Protocol, levels = levels(p$Protocol)[c(1,3,4,2)])
    
    pdf(file = "Activation - Leng2022 Odds Ratio.pdf", height = 4, width = 3.5)
    ggplot(p, aes(x = Protocol, y = OR, ymin = Lower, ymax = Upper)) +
      geom_col(position = "dodge", width = 0.7, fill = pal_lancet(alpha = 0.7)(2)[2]) +
      geom_errorbar(width = 0.2) +
      theme_bw() +
      theme_gjs +
      scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
      geom_hline(yintercept = 1, linetype = 2) +
      theme(legend.position = "none") +
      labs(y = ("Odds Ratio for Hit Genes\nWithin IL-1alpha+TNF+C1q Treated Astrocytes"), x = "hiPSC Protocol")
    dev.off()
    
    # for enrichments within each individual study
    p <- leng_replicability
    p$GeneActivatedIn <- gsub("4\\+", "4", p$GeneActivatedIn)
    
    pdf(file = "Activation - Leng2022 Replicability Odds Ratio.pdf", height = 4, width = 3.5)
    ggplot(p, aes(x = GeneActivatedIn, y = OR, ymin = Lower, ymax = Upper)) +
      geom_col(position = "dodge", width = 0.7, fill = pal_lancet(alpha = 0.7)(2)[2]) +
      geom_errorbar(width = 0.2) +
      theme_bw() +
      theme_gjs +
      scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
      geom_hline(yintercept = 1, linetype = 2) +
      theme(legend.position = "none") +
      labs(y = ("Odds Ratio for Hit Genes\nWithin IL-1alpha+TNF+C1q Treated Astrocytes"), x = "Number of hiPSC Protocols in which Gene is Activated")
    dev.off()
    

## (Mouse) astrocyte activation by neuroactive compounts
  ## RNAseq results
    ## Read in
      sardar_rna_path <- "../../../PublicData/Sardar2021_NeuroactiveComponentsMouseAst/TableS2.xlsx"
      sardar_rna_sheets <- excel_sheets(sardar_rna_path)
      
      sardar_rna_data <- list()
      for (j in sardar_rna_sheets) sardar_rna_data[[j]] <- read_xlsx(sardar_rna_path, sheet = j)
      sardar_rna_data <- do.call("rbind", sardar_rna_data)
      sardar_rna_data$Compound <- splitter(rownames(sardar_rna_data), "_DEG", 1)
      sardar_rna_data <- as.data.frame(sardar_rna_data)
      
      
    ## Add human gene symbol
      sardar_rna_data$HumanSymbol <- convert_mouse2human(genes = sardar_rna_data$GeneSymbol, path.fix = TRUE, return.vector = TRUE)
      
    ## Save
      write.csv(sardar_rna_data, file = "Neuroactive Compound Induction - Sardar RNAseq Data.csv")
      
    ## Compare to our hit genes
      sardar_rna_comparison <- data.frame(Gene = bg)
      sardar_rna_comparison$Hit <- sardar_rna_comparison$Gene %in% hit.genes
      sardar_rna_comparison$NeuroactiveDE <- sardar_rna_comparison$Gene %in% sardar_rna_data$HumanSymbol

      write.csv(sardar_rna_comparison, file = "Neuroactive Compound Induction - Sardar RNAseq Overlaps.csv")
      
  ## ATAC results
    ## Read in
      sardar_atac_path <- "../../../PublicData/Sardar2021_NeuroactiveComponentsMouseAst/TableS4.xlsx"
      sardar_atac_sheets <- excel_sheets(sardar_atac_path)
      sardar_atac_sheets <- sardar_atac_sheets[-grep("GFP", sardar_atac_sheets)]
      # sardar_atac_sheets <- gsub(" only", "", sardar_atac_sheets)
      
      sardar_atac_data <- list()
      for (j in sardar_atac_sheets) sardar_atac_data[[j]] <- read_xlsx(sardar_atac_path, sheet = j)
      names(sardar_atac_data) <- gsub(" only", "", names(sardar_atac_data))
      sardar_atac_data <- do.call("rbind", sardar_atac_data)
      sardar_atac_data$Compound <- splitter(rownames(sardar_atac_data), "\\.", 1)
      sardar_atac_data <- as.data.frame(sardar_atac_data)
      
      
    ## Add human gene symbol
      sardar_atac_data$HumanSymbol <- convert_mouse2human(genes = sardar_atac_data$GeneSymbol, path.fix = TRUE, return.vector = TRUE)
      
    ## Save
      write.csv(sardar_atac_data, file = "Neuroactive Compound Induction - Sardar ATAC Data.csv")
      
    ## Compare to our hit genes
      sardar_atac_comparison <- data.frame(Gene = bg)
      sardar_atac_comparison$Hit <- sardar_atac_comparison$Gene %in% hit.genes
      sardar_atac_comparison$NeuroactiveDE <- sardar_atac_comparison$Gene %in% sardar_atac_data$HumanSymbol
  
      write.csv(sardar_atac_comparison, file = "Neuroactive Compound Induction - Sardar ATAC Overlaps.csv")
      
      table(sardar_atac_comparison$Hit, sardar_atac_comparison$NeuroactiveDE) %>% fisher.test() # ns, but suggestive OR of 1.55
      
      
################################################################################################################################ #
## Immunopanned astrocytes from the human brain ----
   
## Zhang 2016       
  ## A seminal work which pioneered the immunopanning of astrocytes from the human brain. Many comparisons within:
    # mouse vs human
    # foetal vs adult
    # glioblastoma vs healthy
    # epilepsy vs healthy
        
  ## Read in
    zhang_path <- "../../../PublicData/Zhang2016_ImmunpanningHumanAst/1-s2.0-S0896627315010193-mmc5.xlsx"
    zhang_sheets <- excel_sheets(zhang_path)
    
    zhang_data <- list()
    for (j in zhang_sheets) zhang_data[[j]] <- read_xlsx(zhang_path, sheet = j, skip = 1)
      
  ## Collect the genes
    zhang_data <- lapply(zhang_data, function(x) x$Gene) # column 1 always has gene symbol, but its first entry is always "Gene"
    
  ## Convert mouse symbol to human (where relevant)
    m <- grep("mouse", zhang_sheets) 
    zhang_data[m] <- lapply(zhang_data[m], convert_mouse2human, path.fix = TRUE)
    
  ## Rename levels
    names(zhang_data) <- gsub(" ", "_", names(zhang_data))
      
  ## Compare to hit genes
    zhang_overlap <- data.frame(Gene = bg)
    zhang_overlap$Hit <- zhang_overlap$Gene %in% hit.genes
    for (j in names(zhang_data)) zhang_overlap[,j] <- zhang_overlap$Gene %in% zhang_data[[j]]
    
    # combine up and down
    m1 <- grep("mouse", colnames(zhang_overlap))
    zhang_overlap$human_vs_mouse <- rowSums(zhang_overlap[,m1]) >= 1
    
    m2 <- grep("fetal", colnames(zhang_overlap))
    zhang_overlap$fetal_vs_adult <- rowSums(zhang_overlap[,m2]) >= 1
    
    m3 <- grep("GBM", colnames(zhang_overlap))
    zhang_overlap$GBM_vs_Healthy <- rowSums(zhang_overlap[,m3]) >= 1
    
    m4 <- grep("Epilepsy", colnames(zhang_overlap))
    zhang_overlap$Epilepsy_vs_Healthy <- rowSums(zhang_overlap[,m4]) >= 1
    
    # save
    write.csv(zhang_overlap, file = "Zhang 2016 - Overlaps.csv")
    
  ## Statistical comparisons
    zhang_stats <- apply(zhang_overlap[,-c(1,2)], 2, function(x) {
      total <- length(x)
      y <- zhang_overlap$Hit
      f <- table(x, y) %>% fisher.test()
      data.frame(nSig = sum(x),
                 nSigAndHit = sum(x & y),
                 p = f$p.value,
                 OR = f$estimate,
                 Lower = f$conf.int[1],
                 Upper = f$conf.int[2])
    })
    
    zhang_stats <- do.call("rbind", zhang_stats)
    write.csv(zhang_stats, file = "Zhang 2016 - Fisher.csv")

    
## Krawczyk 2022  
  # similar to the above study, but with a much large scope, with ~4x as many samples
    
  ## Read in
    krawczyk_data <- list()
    krawczyk_data$Peritumour <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-3_peritumour.xlsx", skip = 1)
    krawczyk_data$FocalCorticalDysplasia <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-9_FCD.xlsx", skip = 1)
    krawczyk_data$Maturation <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Up")
    krawczyk_data$Maturation <- rbind(krawczyk_data$Maturation, read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Down"))
    krawczyk_data$Ageing <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-12_Ageing.xlsx", skip = 1)
  
  ## Extract gene symbol
    krawczyk_data <- lapply(krawczyk_data, function(x) x$`Gene Name`)
    
  ## Overlap
    krawczyk_overlap <- data.frame(Gene = bg)
    krawczyk_overlap$Hit <- krawczyk_overlap$Gene %in% hit.genes
    for (j in names(krawczyk_data)) krawczyk_overlap[,j] <- krawczyk_overlap$Gene %in% krawczyk_data[[j]]
    write.csv(krawczyk_overlap, file = "Krawczyk 2022 - Overlaps.csv")
    
  ## Statistical comparisons
    krawczyk_stats <- apply(krawczyk_overlap[,-c(1,2)], 2, function(x) {
      total <- length(x)
      y <- krawczyk_overlap$Hit

        f <- table(x, y) %>% fisher.test()
        data.frame(nSig = sum(x),
                   nSigAndHit = sum(x & y),
                   p = f$p.value,
                   OR = f$estimate,
                   Lower = f$conf.int[1],
                   Upper = f$conf.int[2])  

    })
    
    krawczyk_stats <- do.call("rbind", krawczyk_stats)
    write.csv(krawczyk_stats, file = "Krawczyk 2022 - Fisher.csv")

  
    
      
################################################################################################################################ #
## New work from Nov 2022: a unified annotation ----
  
## Setup dataframe with basic annotation
  hits <- res.final[which(res.final$HitPermissive),]
  
  gene.annot <- data.frame(Symbol = unique(hits$Gene[which(hits$HitPermissive)]),
                                 EnsID = ".",
                           ManualSummary = "")
  
  # ensid
  m <- match(gene.annot$Symbol, geneInfo$Symbol)
  gene.annot$EnsID <- geneInfo$EnsID[m] # add EnsID because that's what the signature uses...  
  
  # gene type
  gene.annot$TranscriptType <- factor(geneInfo$Type[m])
  levels(gene.annot$TranscriptType) <- c("lncRNA", "PC")
  
  
  # hit type
  gene.annot$HitPermissive <- TRUE
  gene.annot$HitCore <- gene.annot$Symbol %in% hits$Gene[which(hits$HitCore)]
  
  # stats from differential expression
  gene.annot$LogFCvst <- gene.annot$Distance <- gene.annot$Enh <- "."
  for (j in 1:nrow(gene.annot)) {
    w <- which(hits$Gene == gene.annot$Symbol[j])
    gene.annot$Enh[j] <- paste(hits$Enh[w], collapse = " | ")
    gene.annot$Distance[j] <- paste(round(hits$Gene.Distance[w] / 1000, 1), collapse = " | ")
    gene.annot$LogFCvst[j] <- paste(round(hits$logfc.vst[w], 3), collapse = " | ")
  }
  
  # mean expression
  m <- match(gene.annot$Symbol, hits$Gene)
  gene.annot$Expression.Screen <- round(hits$Gene.Exp[m],3)
  
## Trend information
  load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/5.DevTimecourseAnalysis/ProcessedData/trends.rda") # dev trends
  m <- match(gene.annot$Symbol, trends$Symbol)
  gene.annot$AstDevTrend <- trends$Astro[m]
  gene.annot$AstDevTrend <- factor(gene.annot$AstDevTrend)
  levels(gene.annot$AstDevTrend) <- paste0(levels(gene.annot$AstDevTrend), " ", c("", rep("Up", 3), rep("Dn", 6)))

      
## Astrocyte specificity
  load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/Annotation/sigsBrain.rda") # brain signatures
  sigsBrain <- sigsBrain[c("IP", "DM", "VL", "NG", "CA", "LK")] # subset to human adult brain tissue samples
  spec <- lapply(sigsBrain, function(a) { # gives list of all tfs in the dataset which are most highly expressed in Ast
    a <- a[intersect(rownames(a), gene.annot$EnsID),]
    b <- apply(a, 1, which.max)
    names(b)[which(b == grep("Astrocytes", colnames(a)))]
  })
  
  spec <- do.call("c", spec)
  spec <- table(spec) # number of datasets in which the tf is ast-highest (zero excluded)
  
  m <- match(gene.annot$EnsID, names(spec))
  gene.annot$AstSignatureTop <- spec[m]
  gene.annot$AstSignatureTop[which(is.na(gene.annot$AstSignatureTop))] <- 0
  # gene.annot$AstSignatureTop <- paste0(gene.annot$AstSignatureTop, "/6")
  
  
## gnomAD
  gnomad <- read.delim("../../../PublicData/gnomAD_full_constraint_metrics.tsv")
  gnomad <- gnomad[which(gnomad$canonical == "true"),] # canonical transcript only
  m <- match(gene.annot$Symbol, gnomad$gene)
  
  # pLI: probability of loss of function intolerance, where >0.9 is the accepted cut-off for being intolerant
  gene.annot$pLI <- round(gnomad$pLI[m], 3)
  
  # the new gnomad metric is observed/expected loss of function mutations. they recommend the upper bound of 90% confidence interval, so called: LOEUF (0.35 is the accepted cut-off)
  gene.annot$LOEUF <- gnomad$oe_lof_upper[m]
  
## Mouse phenotype
  mgi <- read.delim("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/Annotation/MGIPhenotypes/MGIBatchReport_20220707_043045.txt")

  # homology
  hom <- read.delim("../../../PublicData/HOM_MouseHumanSequence.txt") # downloaded from: http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
  mgi <- mgi[which(mgi$Symbol %in% hom$Symbol),]
  m <- match(mgi$Symbol, hom$Symbol)
  m <- m + 1 # to reflect the fact that the human homologue of a mouse gene in row m is found in row m + 1
  mgi$HumanSymbol <- hom$Symbol[m]

  # match
  mgi <- mgi[which(mgi$HumanSymbol %in% gene.annot$Symbol),]
  
  gene.annot$MGI <- ""
    
    for (j in 1:nrow(gene.annot)) {
      g <- gene.annot$Symbol[j]
      if (g %in% mgi$HumanSymbol) {
        y <- mgi[which(mgi$HumanSymbol == g),]
        gene.annot$MGI[j] <- paste(y$Term, collapse = " | ")
      } else {
        next
      }
    }
  
  # gene.annot$MGI.Brain <- grepl("brain", gene.annot$MGI, ignore.case = TRUE)
  gene.annot$MGI.Lethal <- grepl("lethal", gene.annot$MGI, ignore.case = TRUE)
  
  
  
  
## DisGeNeT
  ## TEMPORARY ##
  disgenet.g2d <- read.csv("../Scratchspace/X.csv", row.names = 1)
  
  ## What are the disease classes herein?
    x <- strsplit(disgenet.g2d$disease_class_name, ";")
    x <- do.call("c", x)
    x <- unique(x) 
    x <- gsub(" ", "", x)
    x <- unique(x)
    
    # from these, use:
    # Nervous System Diseases, Mental Disorders, and Behavior and Behavior Mechanisms
    # Neoplasms
    
  ## Get neoplasm associations
    x <- disgenet.g2d[grep("Neoplasms", disgenet.g2d$disease_class_name),]
    
    gene.annot$DGNT.Neoplasm <- ""
    
    for (j in 1:nrow(gene.annot)) {
      g <- gene.annot$Symbol[j]
      if (g %in% x$gene_symbol) {
        y <- x[which(x$gene_symbol == g),]
        gene.annot$DGNT.Neoplasm[j] <- paste(y$disease_name, collapse = " | ")
      } else {
        next
      }
    }
  
  ## Get brain-associated phenotype associations
    x <- disgenet.g2d[grep("Nervous System Diseases|Mental Disorders|Behaviour", disgenet.g2d$disease_class_name),]
    
    gene.annot$DGNT.Brain <- ""
    
    for (j in 1:nrow(gene.annot)) {
      g <- gene.annot$Symbol[j]
      if (g %in% x$gene_symbol) {
        y <- x[which(x$gene_symbol == g),]
        gene.annot$DGNT.Brain[j] <- paste(y$disease_name, collapse = " | ")
      } else {
        next
      }
    }
    
## snRNAseq
  # read in results from an above section
  sum <- read.csv("snRNAseq Overlaps Binary.csv", row.names = 1)
  
  # for each gene, get the snRNA-seq datasets in which it is DE
  x <- apply(sum, 1, function(x) {
    x <- as.logical(x[-1])
    if (any (x)) {
      x <- paste(colnames(sum)[-1][x], collapse = " | ")
    
    } else {
      x <- ""
    }
    return(x)
  })
  
  m <- match(gene.annot$Symbol, sum$Gene)
  gene.annot$snRNAseq <- x[m]
  
  
## Save
  write.csv(gene.annot, "Genes/Gene Annotation V2.csv", row.names = FALSE)
  