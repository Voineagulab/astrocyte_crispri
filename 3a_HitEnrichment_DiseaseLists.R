## This script preprocesses all counts data output from CellRanger using R

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

## Load
  source("../../Scripts/Functions.R")
  load("../../Data/Preprocessed/NHA Pooled.rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist.csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]

## Data information
  pos <- guides$TargetID[which(guides$TargetCat == "Promoter")]
  pos <- unique(pos)
  enh <- guides$TargetID[which(guides$TargetCat == "Enh")]
  enh <- unique(enh)
  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  sig.colours <- c("black", "firebrick1")

## Results
  # s <- read.csv("Enhancers - All Characterisation.csv")
  load("../2_DE/Enhancers - All Results Summary.rda")

  
## Set up gene lists
  query <- unique(res$Gene[which(res$Hit)]) # hits
  bg <- unique(res$Gene) # highly-expressed genes as background
  
  
################################################################################################################################ #
## GO ----
  
library(gprofiler2)
  


## Run
  go <- gost(query = query, custom_bg = bg, significant = TRUE, evcodes = TRUE, organism = "hsapiens", user_threshold = 0.05, correction_method = "fdr")
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
  
  write.csv(g, file = "GO.csv", quote = FALSE)  
  
## Analyse and plot!
  # top hits
                                                           

  
  ## Curated hits
    curated <- c("animal organ development", 
                 "memory",
                 "positive regulation of neuron projection development",
                 "regulation of developmental process",
                 "calcium ion import")
    
    p <- go[which(go$term_name %in% curated),]
    
    p$term_name[3] <- gsub("neuron", "neuron\n", p$term_name[3])
    p$term_name[4] <- gsub("developmental", "developmental\n", p$term_name[4])
    p$FDR <- -log10(p$p_value)
    p$term_name <- factor(p$term_name, levels = p$term_name)
    
    pdf(file = "GO (Curated).pdf", height = 2.5, width = maxw)
    ggplot(p, aes(x = term_name, y = FDR)) +
      geom_col(width = 0.6, fill = "black", colour = "black") +
      theme_bw() +
      theme(panel.border = invis, axis.line.x = element_line(), axis.text.x = element_text()) +
      labs(y = "-log10(FDR)", x = "GO Term") +
      scale_y_continuous(expand = c(0,0)) +
      coord_flip()
  dev.off()
  
  

   
  
  
################################################################################################################################ #
## Disgenet ----
  
library(disgenet2r)
  
## Associated packages
  library(VennDiagram) 
  library(stringr)
  library(tidyr)
  library(SPARQL)
  library(RCurl)
  library(igraph)
  library(ggplot2)
  library(reshape2)
  
## Once-off authorisation
  disgenet_api_key <- get_disgenet_api_key(
    email = "i.voineagu@unsw.edu.au", 
    password = "UNSWlab123" )
  Sys.setenv(DISGENET_API_KEY = disgenet_api_key)
  
## Run
  # disease associations for each gene
  disgenet.g2d <- gene2disease(gene = query, 
                      database = "ALL", 
                      score = c(0.2, 1), 
                      verbose= TRUE)
  
  # disease set enrichments for the list of hits
  disgenet.enrich <- disease_enrichment(entities = query, 
                                        database = "ALL", 
                                        universe = "CUSTOM", 
                                        custom_universe = bg, 
                                        verbose = TRUE)
  
  
  # note: One or more of the genes in the list is not in DisGeNET ( 'ALL' ):
   # - COA7
   # - FAM107B
   # - PCGF5
   # - CTR9
   # - EHBP1L1
   # - RAB6A
   # - TMEM9B
   # - AC068790.7
   # - SYNE3
   # - COPRS
   # - CNN1
   # - AC020916.1
   # - FBXO17
   # - UQCR10
   # - NCKIPSD
   # - AMD1
   # - GTF3C6
   # - GLIPR2
   # - MAMDC2
   # - CAMK2N1
   # - RPLP1
  

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
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / 70
  p <- p[which(p$FDR < 0.05),]
  p <- p[1:10,]
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Disgenet Enrichment - Mental Disorders.pdf", height = 4, width = maxw)
  denrich(p)
  dev.off()
  
  # "Nervous System Diseases"
  p <- disgenet.enrich@qresult[grepl("Nervous System Diseases", disgenet.enrich@qresult$disease_class_name),]
  p <- p[!(grepl("Neoplasms|Mental|Cardiovascular|Musculo", p$disease_class_name)),]
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / 70
  p <- p[which(p$FDR < 0.05),]
  p <- p[c(1,2,4:11),]
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Disgenet Enrichment - Nervous System Diseases.pdf", height = 3.5, width = maxw)
  denrich(p)
  dev.off()
  
  # "Neoplasms"
  p <- disgenet.enrich@qresult[grepl("Neoplasms", disgenet.enrich@qresult$disease_class_name),]
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / 70
  p <- p[1:15,]
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Disgenet Enrichment - Neoplasms.pdf", height = 4, width = maxw)
  denrich(p)
  dev.off()
  
  
  x <- read.delim("../../../PublicData/DisGeNet_070422/all_gene_disease_associations.tsv")
  y <- x[which(x$geneSymbol %in% query),]
  y <- y[which(y$score >= 0.2),]
  
  z <- y[grep("Menta", y$diseaseSemanticType),]
  
  z <- y[grep("Neopl", y$diseaseSemanticType),]
  z <- z[grep("Glio|Neur", z$diseaseName),]
