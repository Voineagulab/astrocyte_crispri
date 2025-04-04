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

## Only tests for well-powered genes
filterForWellPowered <- TRUE
# if (filterForWellPowered) setwd("")

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
  library(clusterProfiler)
  library(enrichplot)
  

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
  
## Load power
  pow <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable3_DE/3F_Power.csv")
  wellpowered_egps <- pow$Pair[which(pow$WellPowered015  | pow$Hit)]
  wellpowered_enh <- pow$Enhancer[which(pow$WellPowered015 | pow$Hit)] %>% unique()
  wellpowered_genes <- pow$Gene[which(pow$WellPowered015 | pow$Hit)] %>% unique()
  
## Set up gene lists
  # the hit and background genes, useful for enrichment
  hit.genes <- unique(res.final$Gene[which(res.final$HitPermissive)]) # hit genes
  bg <- unique(res.final$Gene) # highly-expressed genes as background
  
  bg_wellpowered <- bg[which(bg %in% wellpowered_genes)]
  
  # a dataframe to store gene-level logical annotations
  annot.logi <- data.frame(Gene = bg, Hit = bg %in% hit.genes)
  annot.logi <- annot.logi[order(!(annot.logi$Hit), annot.logi$Gene),]
  annot.logi.up <- annot.logi.down <- annot.logi
  
  
  # a dataframe to store fisher test results for overenrichments
  # (actually, a function to fill initialise/fill this)
  run.fisher.vsHits <- function(data = annot.logi, 
                                data.column, 
                                rbind = FALSE, 
                                rbind.to = enrichments,
                                replace = FALSE,
                                filt.pow = filterForWellPowered,
                                alternative = "two.sided",
                                signed = "Unsigned") {
    
    # filtering to well-powered hits
    if (filt.pow) data <- data[which(data$Gene %in% wellpowered_genes),]
    
    # get information
    total <- nrow(data)
    x <- data[,data.column] %>% as.logical()
    y <- data$Hit
    
    # run stats
    f <- table(x, y) %>% fisher.test(alternative = alternative)
    
    # output stats
    out <- data.frame(Resource = data.column,
                      Total_TRUE = sum(x, na.rm = TRUE),
                      Fraction_Bg_TRUE = sum(x, na.rm = TRUE) / total,
                      Total_Hit_TRUE = sum(x & y, na.rm = TRUE),
                      Fraction_Hit_TRUE = sum(x & y, na.rm = TRUE) / sum(y, na.rm = TRUE),
                      p = f$p.value,
                      OR = f$estimate,
                      Lower = f$conf.int[1],
                      Upper = f$conf.int[2],
                      Signed = signed,
                      row.names = data.column)  
    
    if (rbind) {
      if (replace) {
        enrichments[data.column,] <- out
        return(enrichments  )
      } else {
        return(rbind(rbind.to, out))
      }
      
    } else {
      return(out)  
    }
    
  }
    
  
## For hit genes only, store detailed annotation information in this dataframe
  annot.detailed <- data.frame(Gene = sort(hit.genes))
  
## Quick load
  # load("Genes/Wellpowered/Final.rda")

################################################################################################################################ #
## GO ----
  
## Run gprofiler 2
  # run
  go <- gprofiler2::gost(query = hit.genes, custom_bg = bg_wellpowered, significant = TRUE, evcodes = TRUE, organism = "hsapiens", user_threshold = 0.05, correction_method = "fdr")
  go <- as.data.frame(go$result)
  go$parents <- sapply(go$parents, function(x) paste(x, collapse = "_"))

  # save
  g <- apply(go, 2, function(x) { # this code makes the dataframe suitable for csv format
    
    if (class(x) != "character") {
      return(x)
    } else {
      return(gsub(",", "_", x))
    }
    
  })
  
  write.csv(g, file = "Genes/Wellpowered/GO - gprofiler2.csv", quote = FALSE, row.names = FALSE)  
  
## Run Clusterprofiler
  cp <- enrichGO(gene = hit.genes,
                 universe = bg_wellpowered,
                 # universe = bg, 
                 OrgDb = "org.Hs.eg.db", 
                 keyType = "SYMBOL", 
                 ont = "all", 
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", 
                 readable = FALSE) # a new update
  
  save(cp, file = "Genes/Wellpowered/GO - Clusterprofiler.rda")
  x <- cp@result
  write.csv(x, "Genes/Wellpowered/GO - Clusterprofiler.csv")
  
  # test a few plots
  # goplot(cp) # not working
  dotplot(cp, split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY~., scale="free") +
    theme_bw()
  
  p <- cp@result[1:10,]
  
  
  
  cp_dotplot <- function(p) {
    p$Ratio <- as.numeric(splitter(p$GeneRatio, "/", 1)) / as.numeric(splitter(p$GeneRatio, "/", 2))
    p$Description <- factor(p$Description, levels = p$Description[order(-p$pvalue)])
    
    ggplot(p, aes(x = Description, y = -log10(p.adjust), size = Ratio)) +
      geom_segment(aes(y = 0, yend = -log10(p.adjust), xend = Description), size = 0.4, linetype = 2,
                   colour = "grey50", alpha = 0.5) +  
      geom_point( colour = "black", fill = "black") +
      coord_flip() +
      # scale_fill_viridis_c(option = "A") +
      theme_bw() +
      scale_y_continuous(limits = c(0, NA)) +
      guides(size = guide_legend(ncol = 1)) +
      geom_hline(yintercept = -log10(0.05)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 3.5)) +
      theme(axis.title.y = invis, panel.border = invis, axis.line.x = element_line(), panel.grid = invis,
            legend.position = "right") +
      labs(y = "FDR (-log10)")
  }
  
  
  
#   cnetplot(cp@, layout="dh", colorEdge = TRUE,
#          cex_label_gene=0.5,cex_label_category=0.7)
# 
#    cnetplot(cp, categorySize="pvalue", showCategory =20, font.size=12,label_format=30,cex_line =1,colorEdge = TRUE)
# 
#   
# #plot interesting terms
# pdf("GOenrichment plots.pdf", height=5, width=10)
# terms=c( "tissue development", "cell migration","cellular response to growth factor stimulus", "chemotaxis", "response to wounding", "regulation of cell differentiation")
# keep=which(cp@result$Description%in%terms)
# cp2=cp; cp2@result=cp2@result[keep,]
# 
# cnetplot(cp2,   layout="dh", colorEdge = TRUE,
#          cex_label_gene=0.5,cex_label_category=0.7)
# cnetplot(go2,   layout="dh", colorEdge = TRUE,
#          cex_label_gene=0.5,cex_label_category=0.7, circular=TRUE)
# terms=c( "cell migration","cellular response to growth factor stimulus", "chemotaxis", "response to wounding", "regulation of cell differentiation")
# keep=which(go@result$Description%in%terms)
# go2=go; go2@result=go@result[keep,]
# cnetplot(go2,   layout="dh", colorEdge = TRUE,
#          cex_label_gene=0.5,cex_label_category=0.7)
# cnetplot(go2,   layout="dh", colorEdge = TRUE,
#          cex_label_gene=0.5,cex_label_category=0.7, circular=TRUE)

  
  
################################################################################################################################ #
## Disgenet ----
  

## Once-off authorisation
  disgenet_api_key <- get_disgenet_api_key(
    email = "i.voineagu@unsw.edu.au", 
    password = "UNSWlab123" )
  Sys.setenv(DISGENET_API_KEY = disgenet_api_key)
  
  
## Run
  # disease associations for each gene
  # disgenet.g2d <- gene2disease(gene = hit.genes, 
  #                     database = "ALL", 
  #                     score = c(0, 1), 
  #                     verbose = TRUE)
  # 
  # x <- disgenet.g2d@qresult
  # x <- x[order(x$gene_symbol),]
  # x <- relocate(x, "gene_symbol", "disease_name")
  # write.csv(x, file = "Genes/Wellpowered/Disgenet - Gene2Disease.csv", row.names = FALSE)
  
  
  
  
  
  # disease set enrichments for the list of hits
  # disgenet.enrich <- disease_enrichment(entities = hit.genes,
  #                                      database = "ALL",
  #                                      universe = "DISGENET",
  #                                      # custom_universe = bg,
  #                                      vocabulary = "HGNC",
  #                                      verbose = TRUE)
  #
  # x_disgenet <- disgenet.enrich@qresult
  # write.csv(x, file = "Genes/Wellpowered/Disgenet - Enrichment (Curated).csv", row.names = FALSE)
  # GJS note: the function seems to be insensitive to the choice of universe, and defaults to all genes.
  
  ## Disease associations for each gene
    # read in  
    disg_all <- read.delim("../../../PublicData/DisGeNet_070422/all_gene_disease_associations.tsv")
    
    # filter to tested genes
    # g2d <- disg_all[which(disg_all$geneSymbol %in% res.final$Gene),]
    g2d <- disg_all[which(disg_all$geneSymbol %in% bg_wellpowered),]
    
    # clean table
    g2d <- data.frame(Gene = g2d$geneSymbol,
                      Hit = g2d$geneSymbol %in% hit.genes,
                      Disease = g2d$diseaseName,
                      Disease_Type = g2d$diseaseType,
                      Disease_SemanticType = g2d$diseaseSemanticType,
                      Score = g2d$score,
                      NStudies = g2d$NofPmids)
    
    # test enrichment
    g2d_annot <- annot.logi[,1:2]
    if (exists("g2d_enrich")) rm(g2d_enrich)
    
    for (j in unique(g2d$Disease)) {
      print(j)
      
      # filter
      w <- which(g2d$Disease == j)
      if (length(w) < 10) next 
      k <- g2d[w,]
      
      # enrich
      g2d_annot[,j] <- g2d_annot$Gene %in% k$Gene
      g2d_enrich <- run.fisher.vsHits(data.column = j, 
                                      rbind = exists("g2d_enrich"), 
                                      data = g2d_annot, rbind.to = 
                                        g2d_enrich, 
                                      signed = "Unsigned",
                                      alternative = "greater")  
      
      
    }
    
    g2d_enrich <- g2d_enrich[order(g2d_enrich$p),]
    g2d_enrich$FDR <- p.adjust(g2d_enrich$p, method = "fdr")
    
    # augment annotation
    colnames(g2d_enrich) <- c("Disease", "Total_count", "Total_ratio", "Hit_count", "Hit_ratio", "p", "OR", "Lower", "Upper", "Signed", "FDR")
    m <- match(g2d_enrich$Disease, g2d$Disease)
    g2d_enrich$Disease_Type <- g2d$Disease_Type[m]
    g2d_enrich$Disease_SemanticType <- g2d$Disease_SemanticType[m]
    
    g2d_enrich <- g2d_enrich[,c("Disease", "p", "FDR", "OR", "Disease_Type", "Disease_SemanticType", "Total_count", "Total_ratio", "Hit_count", "Hit_ratio", "Lower", "Upper")]
    rownames(g2d_enrich) <- 1:nrow(g2d_enrich)
    
    write.csv(g2d_enrich, "Genes/Wellpowered/Disgenet - Enrichment.csv")
    
  

## Save
  # save(disgenet.enrich, disgenet.g2d, file = "Genes/Wellpowered/Disgenet.rda")
  
  
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
  
  pdf(file = "Genes/Wellpowered/Disgenet Enrichment - Mental Disorders (Curated).pdf", height = 4, width = 8)
  denrich(p)
  dev.off()
  
  # "Nervous System Diseases"
  p <- disgenet.enrich@qresult[grepl("Nervous System Diseases", disgenet.enrich@qresult$disease_class_name),]
  p <- p[!(grepl("Neoplasms|Mental|Cardiovascular|Musculo", p$disease_class_name)),]
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / length(hit.genes)
  p <- p[which(p$FDR < 0.05),]
  p <- p[1:10,]
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Genes/Wellpowered/Disgenet Enrichment - Nervous System Diseases (Curated).pdf", height = 3.5, width = 8)
  denrich(p)
  dev.off()
  
  # "Neoplasms"
  p <- disgenet.enrich@qresult[grepl("Neoplasms", disgenet.enrich@qresult$disease_class_name),]
  p$Ratio <- as.numeric(splitter(p$Ratio, "/", 1)) / length(hit.genes)
  p <- p[1:15,]
  p$Description <- factor(p$Description, levels = rev(p$Description))
  
  pdf(file = "Genes/Wellpowered/Disgenet Enrichment - Neoplasms (Curated).pdf", height = 4, width = 8)
  denrich(p)
  dev.off()
  
## Add disgenet annotations to each gene
  annot.detailed$Disgenet <- "."
  for (j in annot.detailed$Gene) {
    x <- disgenet.g2d@qresult[which(disgenet.g2d@qresult$gene_symbol == j),]
    annot.detailed$Disgenet[which(annot.detailed$Gene == j)] <- paste(x$disease_name, collapse = "; ")
  }
  
  annot.detailed$Disgenet <- gsub(",", "", annot.detailed$Disgenet)
  
  
################################################################################################################################ #
## Miscellaneous gene lists ----
  
  
## SFARI
  sfari <- read.csv("../../../PublicData/SFARI-Gene_genes_10-12-2022release_10-13-2022export.csv")
  m <- match(annot.detailed$Gene, sfari$gene.symbol)
  
  annot.detailed$SFARI_Score <- sfari$gene.score[m]
  annot.detailed$SFARI_Syndromic <- sfari$syndromic[m]
  
  
## Phenotype in mouse knockouts
  # access from https://www.informatics.jax.org/batch
  write.table(bg, file = "../../../PublicData/MGIPhenotypes/BackgroundGenes.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  # read in
  mgi <- read.delim("../../../PublicData/MGIPhenotypes/MGIBatchReport_20230309_050302.txt")
  mgi <- mgi[which(mgi$Term != ""),]

  # extract annotation
  annot.detailed$MGI <- NA
  
    for (j in 1:nrow(annot.detailed)) {
      print(j)
      g <- annot.detailed$Gene[j]
      
      if (g %in% mgi$Input) {
        
        x <- mgi[which(mgi$Input == g),]
        annot.detailed$MGI[j] <- paste(x$Term, collapse = "; ")
        
      } else {
        
        next
        
      }
    }
  
  annot.detailed$MGI.Lethal <- grepl("lethal", annot.detailed$MGI, ignore.case = TRUE)
  
  
## gnomAD
  gnomad <- read.delim("../../../PublicData/gnomAD_full_constraint_metrics.tsv")
  gnomad <- gnomad[which(gnomad$canonical == "true"),] # canonical transcript only
  
  ## Add detailed annotation
    m <- match(annot.detailed$Gene, gnomad$gene)
    
    # pLI: probability of loss of function intolerance
    annot.detailed$pLI <- round(gnomad$pLI[m], 3)
    
    # the new gnomad metric is observed/expected loss of function mutations
    annot.detailed$LOEUF <- gnomad$oe_lof_upper[m]  
    
  ## Add logical annotation  
    m <- match(annot.logi$Gene, gnomad$gene)
    
    # pLI: where >0.9 is the accepted cut-off for being intolerant
    annot.logi$pLI <- gnomad$pLI[m] > 0.9
    enrichments <- run.fisher.vsHits(data.column = "pLI", rbind = FALSE)
    
    # observed/expected loss of function mutations: they recommend the upper bound of 90% confidence interval, so called: LOEUF (0.35 is the accepted cut-off)
    annot.logi$LOEUF <- gnomad$oe_lof_upper[m] < 0.35
    enrichments <- run.fisher.vsHits(data.column = "LOEUF", rbind = TRUE)  
    
    

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
      annot.logi$Housekeeping_Lin2019 <- annot.logi$Gene %in% hkg
      enrichments <- run.fisher.vsHits(data.column = "Housekeeping_Lin2019", rbind = TRUE)
      
  ## From Eisenberg 2013:
    # https://www.tau.ac.il/~elieis/HKG/
    
    ## Read in
      eisenberg <- read.table("../../../PublicData/HousekeepingGenes/Eisenberg2013_Genes.txt")
      eisenberg <- eisenberg$V1
        
    ## Add to universal annotation
      annot.logi$Housekeeping_Eisenberg2013 <- annot.logi$Gene %in% eisenberg
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
    
    
    
  ## Coarsely summarise the trend across tissues
    gtex_trends <- list()
    
    ## No tissue pooling
    gtex_trends$NoPool <- data.frame(Mean = rowMeans(gtex),
                                     sd = apply(gtex, 1, sd),
                                     cv = NA,
                                     PctExp0 = rowMeans(gtex > 0),
                                     PctExp1 = rowMeans(gtex > 1),
                                     PctExp2 = rowMeans(gtex > 2),
                                     PctExp5 = rowMeans(gtex > 5),
                                     PctExp10 = rowMeans(gtex > 10))
    
    gtex_trends$NoPool$cv <- gtex_trends$NoPool$sd / gtex_trends$NoPool$Mean
  
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
    
    write.csv(brainSpecific, "Genes/Wellpowered/Tissue-Specificity - GTEx PFC vs non-brain.csv")
    
  ## Add to annotation
    annot.logi$BrainSpecific_GTEx <- annot.logi$Gene %in% brainSpecific$Gene[which(brainSpecific$BrainMarker)] 
    enrichments <- run.fisher.vsHits(data.column = "BrainSpecific_GTEx", rbind = TRUE, signed = "Up")
    
    
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
    write.csv(herring_astSpecific, file = "Genes/Wellpowered/Astrocyte Markers - Herring 2022.csv")
    
  ## Add to annotation dataframes
    annot.logi$AstMarkers_AllAges_Herring2022 <- annot.logi$Gene %in% herring_astSpecific$Gene[which(herring_astSpecific$AstMarker)]
    enrichments <- run.fisher.vsHits(data.column = "AstMarkers_AllAges_Herring2022", rbind = TRUE, signed = "Up")
   
## Adult brain markers (Morabito 2021)
  # read in
  morabito_markers <- read_xlsx("../../../PublicData/Morabito2021_AD_snATAC_snRNA/41588_2021_894_MOESM5_ESM.xlsx", sheet = "Supplementary Data 1d", skip = 2)
  morabito_markers <- morabito_markers[grep("ASC", morabito_markers$cluster),]  
  morabito_markers <- morabito_markers[which(morabito_markers$avg_logFC > 0),]

  # add to annotation
  annot.logi$AstMarkers_Adult_Morabito2021 <- annot.logi$Gene %in% morabito_markers$gene
  
  enrichments <- run.fisher.vsHits(data.column = "AstMarkers_Adult_Morabito2021", rbind = TRUE, signed = "Up")
  

## Astrocyte markers (foetal)
  ## Read in
    zhong_data <- read_xlsx("../../../PublicData/snRNAseq/Zhong2018_FoetalSmartSeq/41586_2018_BFnature25980_MOESM2_ESM.xlsx", skip = 4, sheet = 1)
    zhong_data <- zhong_data[which(zhong_data$cluster == "Astrocytes"),] # 353
    zhong_data <- zhong_data[which(zhong_data$avg_diff > 0),] # 326
  
  ## Add to annotation dataframes
    annot.logi$AstMarkers_Foetal_Zhong2018 <- annot.logi$Gene %in% zhong_data$gene
    enrichments <- run.fisher.vsHits(data.column = "AstMarkers_Foetal_Zhong2018", rbind = TRUE, signed = "Up")
    
    
## Astrocyte subtypes (Sadick 2022)
  ## Read in
    subtypeMarkers <- read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST3_ctMarkers.xlsx", sheet = "LEN_so_astro_r2_DEGs")
    subtypeMarkers <- subtypeMarkers[which(subtypeMarkers$p_val_adj < 0.05),] # removes 3? How curious...
  
  ## Add to annotation dataframes
    annot.logi$AstMarkers_Subtypes_Sadick2022 <- annot.logi$Gene %in% subtypeMarkers$gene
    enrichments <- run.fisher.vsHits(data.column = "AstMarkers_Subtypes_Sadick2022", rbind = TRUE, signed = "Up")
  
  ## Detailed annotation too
    annot.detailed$AstSubtypes_Sadick2022 <- NA
  
    for (j in 1:nrow(annot.detailed)) {
      print(j)
      g <- annot.detailed$Gene[j]
      
      if (g %in% subtypeMarkers$gene) {
        
        x <- subtypeMarkers[which(subtypeMarkers$gene == g),]
        annot.detailed$AstSubtypes_Sadick2022[j] <- paste(x$LEN_so_astro_r2_cluster, collapse = "; ")
        
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
    annot.logi$AstMaturation_Herring2022 <- annot.logi$Gene %in% herring_trends$Symbol
    annot.logi.up$AstMaturation_Herring2022 <- annot.logi.up$Gene %in% herring_trends$Symbol[grep("up", herring_trends$Astro)]
    annot.logi.down$AstMaturation_Herring2022 <- annot.logi.down$Gene %in% herring_trends$Symbol[grep("down", herring_trends$Astro)]
    
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Herring2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Herring2022", rbind = TRUE, data = annot.logi.up, signed = "Up")
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Herring2022", rbind = TRUE, data = annot.logi.down, signed = "Down")
  
  ## Detailed annotation too
    annot.detailed$AstMaturation_Herring2022 <- NA
  
    for (j in 1:nrow(annot.detailed)) {
      print(j)
      g <- annot.detailed$Gene[j]
      
      if (g %in% herring_trends$Symbol) {
        m <- match(g, herring_trends$Symbol)
        annot.detailed$AstMaturation_Herring2022[j] <- herring_trends$Astro[m]
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
    zhang_data[m] <- lapply(zhang_data[m], convert_mouse2human, path.fix = TRUE)
    
  ## Rename levels
    names(zhang_data) <- gsub(" ", "_", names(zhang_data))
    
  # ## Combine up and down
  #   m1 <- grep("mouse", names(zhang_data))
  #   zhang_data$human_vs_mouse <- do.call("c", (zhang_data[m1]))
  #   
  #   m2 <- grep("fetal", names(zhang_data))
  #   zhang_data$fetal_vs_adult <- do.call("c", (zhang_data[m2]))
  #   
  #   m3 <- grep("GBM", names(zhang_data))
  #   zhang_data$GBM_vs_Healthy <- do.call("c", (zhang_data[m3]))
  #   
  #   m4 <- grep("Epilepsy", names(zhang_data))
  #   zhang_data$Epilepsy_vs_Healthy <- do.call("c", (zhang_data[m4]))
      
  ## Add to annotation (foetal vs adult only)
    annot.logi$AstMaturation_Zhang2016 <- annot.logi$Gene %in% c(zhang_data$up_fetal_vs_adult, zhang_data$down_fetal_vs_adult)
    annot.logi.up$AstMaturation_Zhang2016 <- annot.logi.up$Gene %in% zhang_data$up_fetal_vs_adult
    annot.logi.down$AstMaturation_Zhang2016 <- annot.logi.down$Gene %in% zhang_data$down_fetal_vs_adult
    
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Zhang2016", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Zhang2016", rbind = TRUE, data = annot.logi.up, signed = "Up")
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Zhang2016", rbind = TRUE, data = annot.logi.down, signed = "Down")
    
  
## Maturation: Krawczyk  
  # similar to the above study, but with a much large scope, with ~4x as many samples
    
  # again, read in all relevant data here, but only output to the annotation dataframes the relevant maturation work
    
  ## Read in all
    krawczyk_data <- list()
    krawczyk_data$Peritumour <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-3_peritumour.xlsx", skip = 1)
    krawczyk_data$FocalCorticalDysplasia <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-9_FCD.xlsx", skip = 1)
    krawczyk_data$Maturation_Up <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Up")
    krawczyk_data$Maturation_Down <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-10_Maturation.xlsx", skip = 1, sheet = "Maturation Down")
    krawczyk_data$Ageing <- read_xlsx("../../../PublicData/Krawczyk2022_ImmunopanningHumanAst/inline-supplementary-material-12_Ageing.xlsx", skip = 1)
  
  ## Extract gene symbol
    # krawczyk_data <- lapply(krawczyk_data, function(x) x$`Gene Name`)
    
  ## Add to annotation (maturation only)
    annot.logi$AstMaturation_Krawczyk2022 <- annot.logi$Gene %in% c(krawczyk_data$Maturation_Up$`Gene Name`, krawczyk_data$Maturation_Down$`Gene Name`)
    annot.logi.up$AstMaturation_Krawczyk2022 <- annot.logi.up$Gene %in% krawczyk_data$Maturation_Up$`Gene Name`
    annot.logi.down$AstMaturation_Krawczyk2022 <- annot.logi.down$Gene %in% krawczyk_data$Maturation_Down$`Gene Name`
    
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Krawczyk2022", rbind = TRUE, data = annot.logi.up, signed = "Up")
    enrichments <- run.fisher.vsHits(data.column = "AstMaturation_Krawczyk2022", rbind = TRUE, data = annot.logi.down, signed = "Down")
    
## Ageing: Krawczyk
  # already processed above
    
  ## Add to annotation (ageing only)
    annot.logi$AstAgeing_Krawczyk2022 <- annot.logi$Gene %in% krawczyk_data$Ageing$`Gene Name`
    annot.logi.up$AstAgeing_Krawczyk2022 <- annot.logi.up$Gene %in% krawczyk_data$Ageing$`Gene Name`[which(krawczyk_data$Ageing$`log2(Fold Change)` > 0)]
    annot.logi.down$AstAgeing_Krawczyk2022 <- annot.logi.down$Gene %in% krawczyk_data$Ageing$`Gene Name`[which(krawczyk_data$Ageing$`log2(Fold Change)` < 0)]
    
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Krawczyk2022", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Krawczyk2022", rbind = TRUE, data = annot.logi.up, signed = "Up")
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Krawczyk2022", rbind = TRUE, data = annot.logi.down, signed = "Down")
    
## Ageing: Palmer 2021
  # snRNA-seq of the ageing brain
    
  ## Read in
    ageing_palmer <- read_xlsx("../../../PublicData/snRNAseq/Palmer2021_Aging/pnas.2114326118.sd06.xlsx", sheet = "Ast")
    colnames(ageing_palmer)[1] <- "Gene"
    ageing_palmer <- ageing_palmer[which(ageing_palmer$p_val_adj < 0.05),]
    
  # note that a positive fold-change indicates higher expression in old vs. young brains.
  
  ## Add to annotation
    annot.logi$AstAgeing_Palmer2021 <- annot.logi$Gene %in% ageing_palmer$Gene
    annot.logi.up$AstAgeing_Palmer2021 <- annot.logi.up$Gene %in% ageing_palmer$Gene[which(ageing_palmer$avg_logFC > 0)]
    annot.logi.down$AstAgeing_Palmer2021 <- annot.logi.down$Gene %in% ageing_palmer$Gene[which(ageing_palmer$avg_logFC < 0)]
    
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Palmer2021", rbind = TRUE)
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Palmer2021", rbind = TRUE, data = annot.logi.up, signed = "Up")
    enrichments <- run.fisher.vsHits(data.column = "AstAgeing_Palmer2021", rbind = TRUE, data = annot.logi.down, signed = "Down")
    

################################################################################################################################ #
## Astrocyte biology: Disease DEGs ----
   
     
## ASD snRNA-seq (Velmeshev et al, 2019)
  # load in
  asd <- readxl::read_xls("../../../PublicData/snRNAseq/Velmeshev2019_ASD/NIHMS1053005-supplement-Data_S4.xls", sheet = 1)
  
  # filter to ast DEGs
  asd <- asd[grep("AST", asd$`Cell type`),]
  
  # add to annotation
  annot.logi$AstDisease_ASD_Velmeshev2019 <- annot.logi$Gene %in% asd$`Gene name`
  annot.logi.up$AstDisease_ASD_Velmeshev2019 <- annot.logi$Gene %in% asd$`Gene name`[which(asd$`Fold change` > 0)]
  annot.logi.down$AstDisease_ASD_Velmeshev2019 <- annot.logi.down$Gene %in% asd$`Gene name`[which(asd$`Fold change` < 0)]
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Velmeshev2019", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Velmeshev2019", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Velmeshev2019", rbind = TRUE, data = annot.logi.down, signed = "Down")
  

## ASD snRNA-seq (Gandal et al, 2022)
  # read in
  asd2022 <- readxl::read_xlsx("../../../PublicData/snRNAseq/Gandal2022_ASD/41586_2022_5377_MOESM10_ESM.xlsx", sheet = "DEA_ASDvCTL_sumstats", skip = 1)
  colnames(asd2022) <- c("Celltype", "Region", "Gene", "P", "logFC", "FDR")  
  
  # filter to ast degs
  asd2022 <- asd2022[grep("ASTRO", asd2022$Celltype),] 
  
  # filter to PFC
  asd2022 <- asd2022[asd2022$Region == "PFC",] 
  
  # add to annotation
  annot.logi$AstDisease_ASD_Gandal2022 <- annot.logi$Gene %in% asd2022$Gene
  annot.logi.up$AstDisease_ASD_Gandal2022 <- annot.logi$Gene %in% asd2022$Gene[which(asd2022$logFC > 0)]
  annot.logi.down$AstDisease_ASD_Gandal2022 <- annot.logi$Gene %in% asd2022$Gene[which(asd2022$logFC < 0)]
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Gandal2022", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Gandal2022", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_ASD_Gandal2022", rbind = TRUE, data = annot.logi.down, signed = "Down")
  
  
## Multiple sclerosis snRNA-seq (Jakel et al, 2019) 
  # read in
  ms <- list(Ast1 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes1")),
             Ast2 = read_xlsx("../../../PublicData/snRNAseq/Jakel2019_MS/41586_2019_903_MOESM5_ESM.xlsx", sheet = c("Astrocytes2")))
  
  ms <- do.call("rbind", ms)
  ms$Cluster <- splitter(rownames(ms), "\\.", 1)
  
  # add to annotation
  annot.logi$AstDisease_MS_Jakel2019 <- annot.logi$Gene %in% ms$gene
  annot.logi.up$AstDisease_MS_Jakel2019 <- annot.logi$Gene %in% ms$gene[which(ms$avg_logFC > 0)]
  annot.logi.down$AstDisease_MS_Jakel2019 <- NA # there are no genes; their DE pipeline excluded them
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_MS_Jakel2019", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_MS_Jakel2019", rbind = TRUE, data = annot.logi.up, signed = "Up")
  # enrichments <- run.fisher.vsHits(data.column = "AstDisease_MS_Jakel2019", rbind = TRUE, data = annot.logi.down, signed = "Down") # note: no 


## AD snRNA-seq (Morabito 2021)
  # read in
  ad2021 <- read_xlsx("../../../PublicData/Morabito2021_AD_snATAC_snRNA/41588_2021_894_MOESM5_ESM.xlsx", sheet = "Supplementary Data 1e", skip = 2)
  ad2021 <- ad2021[which(ad2021$celltype == "ASC"),]  
  
  # add to annotation
  annot.logi$AstDisease_AD_Morabito2021 <- annot.logi$Gene %in% ad2021$gene
  annot.logi.up$AstDisease_AD_Morabito2021 <- annot.logi$Gene %in% ad2021$gene[which(ad2021$avg_logFC > 0)]
  annot.logi.down$AstDisease_AD_Morabito2021 <- annot.logi$Gene %in% ad2021$gene[which(ad2021$avg_logFC < 0)]
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Morabito2021", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Morabito2021", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Morabito2021", rbind = TRUE, data = annot.logi.down, signed = "Down")
  
  

## AD snRNA-seq (Sadick 2022)
  # read in
  ad <- list(Up = read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST6_AD_DEGs_Astro.xlsx", sheet = "Astro_upregulated_DEGs_dis"),
             Dn = read_xlsx("../../../PublicData/snRNAseq/Sadick2022_AD_AstroEnriched/ST6_AD_DEGs_Astro.xlsx", sheet = "Astro_downregulated_DEGs_dis"))
  
  ad <- lapply(ad, function(x) {
    colnames(x)[1] <- "Gene"
    return(x)
  })
  
  # ad <- do.call("rbind", ad)
  
  # add to annotation
  annot.logi$AstDisease_AD_Sadick2022 <- annot.logi$Gene %in% c(ad$Up$Gene, ad$Dn$Gene)
  annot.logi.up$AstDisease_AD_Sadick2022 <- annot.logi$Gene %in% ad$Up$Gene
  annot.logi.down$AstDisease_AD_Sadick2022 <- annot.logi$Gene %in% ad$Dn$Gene
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Sadick2022", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Sadick2022", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_AD_Sadick2022", rbind = TRUE, data = annot.logi.down, signed = "Down")
  
## Glioblastoma (Zhang 2016)
  # data already processed above
  
  # add to annotation (GBM results only)
  annot.logi$AstDisease_GBM_Zhang2016 <- annot.logi$Gene %in% c(zhang_data$up_GBM_vs_Healthy, zhang_data$down_GBM_vs_Healthy)
  annot.logi.up$AstDisease_GBM_Zhang2016 <- annot.logi$Gene %in% zhang_data$up_GBM_vs_Healthy
  annot.logi.down$AstDisease_GBM_Zhang2016 <- annot.logi$Gene %in% zhang_data$down_GBM_vs_Healthy
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_GBM_Zhang2016", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_GBM_Zhang2016", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_GBM_Zhang2016", rbind = TRUE, data = annot.logi.down, signed = "Down")

  
## Peritumour (Krawczyk 2022)
  # data already processed above
  
  # add to annotation (GBM results only)
  annot.logi$AstDisease_Peritumour_Krawczyk2022 <- annot.logi$Gene %in% krawczyk_data$Peritumour$`Gene Name`
  annot.logi.up$AstDisease_Peritumour_Krawczyk2022 <- annot.logi$Gene %in% krawczyk_data$Peritumour$`Gene Name`[which(krawczyk_data$Peritumour$`log2(Fold Change)` > 0)]
  annot.logi.down$AstDisease_Peritumour_Krawczyk2022 <- annot.logi$Gene %in% krawczyk_data$Peritumour$`Gene Name`[which(krawczyk_data$Peritumour$`log2(Fold Change)` < 0)]
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_Peritumour_Krawczyk2022", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_Peritumour_Krawczyk2022", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_Peritumour_Krawczyk2022", rbind = TRUE, data = annot.logi.down, signed = "Down")

  
## Epilepsy (Zhang 2016)
  # data already processed above
  
  # add to annotation (epilepsy results only)
  annot.logi$AstDisease_Epilepsy_Zhang2016 <- annot.logi$Gene %in% c(zhang_data$up_Epilepsy_vs_Healthy, zhang_data$down_Epilepsy_vs_Healthy)
  annot.logi.up$AstDisease_Epilepsy_Zhang2016 <- annot.logi$Gene %in% zhang_data$up_Epilepsy_vs_Healthy
  annot.logi.down$AstDisease_Epilepsy_Zhang2016 <- annot.logi$Gene %in% zhang_data$down_Epilepsy_vs_Healthy
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_Epilepsy_Zhang2016", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_Epilepsy_Zhang2016", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_Epilepsy_Zhang2016", rbind = TRUE, data = annot.logi.down, signed = "Down")
  

## Focal Cortical Dysplasia (Krawczyk 2022)
  # data already processed above
  
  # add to annotation (FCD results only)
  # annot.logi$AstDisease_FocalCorticalDysplasia_Krawczyk2022 <- annot.logi$Gene %in% krawczyk_data$FocalCorticalDysplasia$`Gene Name`
  # annot.logi.up$AstDisease_FocalCorticalDysplasia_Krawczyk2022 <- annot.logi$Gene %in% krawczyk_data$FocalCorticalDysplasia$`Gene Name`[which(krawczyk_data$FocalCorticalDysplasia$`log2(Fold Change)` > 0)]
  # annot.logi.down$AstDisease_FocalCorticalDysplasia_Krawczyk2022 <- annot.logi$Gene %in% krawczyk_data$FocalCorticalDysplasia$`Gene Name`[which(krawczyk_data$FocalCorticalDysplasia$`log2(Fold Change)` < 0)]
  
  # enrichments <- run.fisher.vsHits(data.column = "AstDisease_FocalCorticalDysplasia_Krawczyk2022", rbind = TRUE)
  # enrichments <- run.fisher.vsHits(data.column = "AstDisease_FocalCorticalDysplasia_Krawczyk2022", rbind = TRUE, data = annot.logi.up, signed = "Up")
  # enrichments <- run.fisher.vsHits(data.column = "AstDisease_FocalCorticalDysplasia_Krawczyk2022", rbind = TRUE, data = annot.logi.down, signed = "Down")
  # the above line of code was commented out as there are no down DEGs in our screened gene set
  
## "Disease-associated astrocytes": cluster 4 in Habib et al Nat Neuro 2020
  # six astrocyte clusters were found, with pairwise DE performed.
  # cluster 6 was found to be a GFAP-high state
  # cluster 4 was termed DAA, because it was unique to their Alzheimer's population
  # we will focus on the the comparison of these to cluster 1, which they deem the GFAP-low homeostatic state
  
  # read in
  daa <- read_xlsx("../../../PublicData/snRNAseq/Habib2020_MouseDAA/41593_2020_624_MOESM2_ESM.xlsx", sheet = "Supplementary Table 2", skip = 5)
  
  # convert to human symbol
  daa$Human <- convert_mouse2human(daa$Gene)
  
  # filter
  daa <- daa[-which(is.na(daa$Human)),] # remove rows with no found human symbol
  
  # overlap
  annot.logi$AstDisease_DAA_Habib2020 <- annot.logi$Gene %in% daa$Human[which((daa$`Cluster ID #1` == 1) & (daa$`Cluster ID #2` == 4))]
  annot.logi.up$AstDisease_DAA_Habib2020 <- annot.logi$Gene %in% daa$Human[which((daa$`Cluster ID #1` == 1) & (daa$`Cluster ID #2` == 4) & daa$`Avg log fold-change #1/#2` < 0)] # note: this is < 0 for upregulation in Cluster 4
  annot.logi.down$AstDisease_DAA_Habib2020 <- annot.logi$Gene %in% daa$Human[which((daa$`Cluster ID #1` == 1) & (daa$`Cluster ID #2` == 4) & daa$`Avg log fold-change #1/#2` > 0)] # note: this is > 0 for upregulation in Cluster 4
  
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_DAA_Habib2020", rbind = TRUE)
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_DAA_Habib2020", rbind = TRUE, data = annot.logi.up, signed = "Up")
  enrichments <- run.fisher.vsHits(data.column = "AstDisease_DAA_Habib2020", rbind = TRUE, data = annot.logi.down, signed = "Down")
  
    
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

    
  ## Wrangle
    leng_data <- lapply(leng_data, function(x) {
      x <- x[,-c(7:8)] # removes two extraneous columns
      pCol <- grep("adj", colnames(x))
      fcCol <- grep("FC|FoldChange", colnames(x))
      x$Hit <- (x[,pCol] < 0.05) & (abs(x[,fcCol]) > log2(1.5))
      x <- x[which(x$Hit),]
      return(x)
    })
  
    names(leng_data) <- paste0("AstActivation_", names(leng_data))
  

  ## Add each to annotation dataframes
    for (j in names(leng_data)) {

      annot.logi[,j] <- annot.logi$Gene %in% leng_data[[j]]$gene
      annot.logi.up[,j] <- annot.logi$Gene %in% leng_data[[j]]$gene[which(leng_data[[j]]$log2FoldChange > 0)]
      annot.logi.down[,j] <- annot.logi$Gene %in% leng_data[[j]]$gene[which(leng_data[[j]]$log2FoldChange < 0)]
  
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE)
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE, data = annot.logi.up, signed = "Up")
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE, data = annot.logi.down, signed = "Down")
      
    }
    
    
  ## An aside: how replicable are these calls across studies?
    # get 
    leng_consistency <- data.frame(Gene = annot.logi$Gene,
                                   Hit = annot.logi$Hit,
                                   Unsigned = rowSums(annot.logi[,grep("AstActivation_TIC_hiPSC", colnames(annot.logi))]),
                                   Up = rowSums(annot.logi.up[,grep("AstActivation_TIC_hiPSC", colnames(annot.logi.up))]),
                                   Down = rowSums(annot.logi.down[,grep("AstActivation_TIC_hiPSC", colnames(annot.logi.down))]))
  
    write.csv(leng_consistency, file = "Genes/Wellpowered/Activation - Leng2022.csv")
    
    # test
    for (j in 1:4) {
      
      label <- paste0("AstActivation_TIC_hiPSC_Rep", j)
      
      annot.logi[,label] <- annot.logi$Gene %in% leng_consistency$Gene[which(leng_consistency$Unsigned >= j)]
      annot.logi.up[,label] <- annot.logi$Gene %in% leng_consistency$Gene[which(leng_consistency$Up >= j)]
      annot.logi.down[,label] <- annot.logi$Gene %in% leng_consistency$Gene[which(leng_consistency$Down >= j)]
  
      enrichments <- run.fisher.vsHits(data.column = label, rbind = TRUE)
      enrichments <- run.fisher.vsHits(data.column = label, rbind = TRUE, data = annot.logi.up, signed = "Up")
      enrichments <- run.fisher.vsHits(data.column = label, rbind = TRUE, data = annot.logi.down, signed = "Down")
      
      
    }
 


## LPS activation timecourse of mouse astrocytes in vivo
  ## From Hasel 2021
    
  ## Read in
    hasel_genes <- list()
    
    # bulk
    hasel_bulk <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST1_LPS_BulkTimecourse.xlsx", sheet = 2)
    
    get_hasel_bulkDEGs <- function(y, fix = FALSE) {
      
      # genes passing padjusted threshold
      padj <- hasel_bulk[[paste0("padj_", y, "_LPS")]]
      padj <- gsub("E", "e", padj) 
      padj <- as.numeric(padj)
      pass.padj <- padj < 0.05
      
      # genes passing log2fc threshold
      # fc <- hasel_bulk[[paste0("l2f_", y, "_LPS")]]
      # pass.fc <- fc > 1
      
      # get degs
      # z <- hasel_bulk$`Gene Name`[pass.padj & pass.fc]
      z <- hasel_bulk[which(pass.padj), c("Gene Name", paste0("l2f_", y, "_LPS"))]
      colnames(z) <- c("Gene", "logf2c")
      
      # convert to human symbol
      z$Gene <- convert_mouse2human(genes = z$Gene, path.fix = fix, return.vector = TRUE)
      z <- z[-which(is.na(z$Gene)),]
      z <- unique(z)
      
      # return
      return(z)
    }
    
    hasel_genes$MouseLPS_3h <- get_hasel_bulkDEGs("3h")
    hasel_genes$MouseLPS_24h <- rbind(get_hasel_bulkDEGs("24h_Male"), get_hasel_bulkDEGs("24h_Female"))
    hasel_genes$MouseLPS_72h <- get_hasel_bulkDEGs("72h")
    
    # rat
    hasel_rat <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST7_Rat.xlsx", sheet = 2)
    hasel_rat$Human <- convert_rat2human(hasel_rat$SYMBOL, path.fix = FALSE)
    hasel_rat <- hasel_rat[-which(is.na(hasel_rat$Human)),]
    
    get_hasel_ratDEGs <- function(y) {
    
      # genes passing padjusted threshold
      padj <- hasel_rat[[paste0(y, "_padj")]]
      padj <- gsub("E", "e", padj) 
      padj <- as.numeric(padj)
      pass.padj <- padj < 0.05
      
      # # genes passing log2fc threshold
      # fc <- hasel_rat[[paste0(y, "_l2f")]]
      # pass.fc <- fc > 1
      
      # get degs
      z <- hasel_rat[which(pass.padj), c("Human", paste0(y, "_l2f"))]
      colnames(z) <- c("Gene", "logf2c")
      z <- unique(z)
      
      # return
      return(z)
    }
    
    hasel_genes$Rat_IFNgamma <- get_hasel_ratDEGs("Ifng_vs_Cnt") # interferon gamma
    hasel_genes$Rat_IFNbeta <- get_hasel_ratDEGs("Ifnb_vs_Cnt") # interferon beta
    hasel_genes$Rat_TIC <- get_hasel_ratDEGs("TIC_vs_Cnt") # tnfa, ilf1 alpha, c1q
    
    # single cell
    hasel_sc <- read_xlsx("../../../PublicData/Hasel2021_MouseAstActivation/ST6_scRNAseq_DEGs.xlsx", skip = 1)
    hasel_sc$Human <- convert_mouse2human(hasel_sc$Gene_name, path.fix = FALSE)
    hasel_sc <- hasel_sc[-which(is.na(hasel_sc$Human)),]
    
    hasel_genes$MouseLPS_SC <- data.frame(Gene = hasel_sc$Human, logf2c = hasel_sc$log_2_fold)
    
    # rename
    names(hasel_genes) <- paste0("AstActivation_", names(hasel_genes), "_Hasel2021")
    
  ## Add each list to the dataframe
    for (j in names(hasel_genes)) {
      # annot.logi[,j] <- annot.logi$Gene %in% hasel_genes[[j]]
      # enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE)
      
      annot.logi[,j] <- annot.logi$Gene %in% hasel_genes[[j]]$Gene
      annot.logi.up[,j] <- annot.logi$Gene %in% hasel_genes[[j]]$Gene[which(hasel_genes[[j]]$logf2c > 0)]
      annot.logi.down[,j] <- annot.logi$Gene %in% hasel_genes[[j]]$Gene[which(hasel_genes[[j]]$logf2c < 0)]
  
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE)
      enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE, data = annot.logi.up, signed = "Up")
      
      if (j == "AstActivation_MouseLPS_72h_Hasel2021") next
      # enrichments <- run.fisher.vsHits(data.column = j, rbind = TRUE, data = annot.logi.down, signed = "Down")
      
    }
    
    
## Finally: from Leng 2022, looking at DEGs between IRAS1 and IRAS2
  # these are transcriptomically-distinct populations from scRNAseq of TIC-activated iAstrocytes
    
  ## Read in
    iras <- read_xlsx("../../../PublicData/Leng2022_iAstrocyteScreen/ST4_ClusterMarkers.xlsx", sheet = "DEGs between IRAS2 vs. IRAS1")
    
  ## Add to annotation matrices
    # annot.logi$AstActivation_TIC_hiPSC_IRAS <- annot.logi$Gene %in% iras$gene
    annot.logi$AstActivation_TIC_hiPSC_IRAS1 <- annot.logi$Gene %in% iras$gene[which(iras$avg_diff > 0)] # this is a < sign. it is because the paper implies that IRAS2 is group1:  The columns ‘pct.1’ and ‘pct.2’ contain the percent of cells expressing the gene of interest in the cluster of interest (‘pct.1’) or all other cells (‘pct.2’) or in IRAS2 (‘pct.1’) iAstrocytes or IRAS1 (‘pct.2’) iAstrocytes. 
    annot.logi$AstActivation_TIC_hiPSC_IRAS2 <- annot.logi$Gene %in% iras$gene[which(iras$avg_diff > 0)]
    
    # enrichments <- run.fisher.vsHits(data.column = "AstActivation_TIC_hiPSC_IRAS", rbind = TRUE, signed = "Unsigned")
    enrichments <- run.fisher.vsHits(data.column = "AstActivation_TIC_hiPSC_IRAS1", rbind = TRUE, signed = "Up")
    enrichments <- run.fisher.vsHits(data.column = "AstActivation_TIC_hiPSC_IRAS2", rbind = TRUE, signed = "Up")
  
  
    
  
################################################################################################################################ #
## Make a combined annotation file ----
    
## Make dataframe from base scaffold
  annot.combined <- annot.logi[,1:2]
    
## Combine information on marker calls
  # notes: exclude subtype markers; all calls are for positive markers
  isMarker <- c("AstMarkers_AllAges_Herring2022", "AstMarkers_Adult_Morabito2021", "AstMarkers_Foetal_Zhong2018")
  annot.combined$AstMarker <- rowSums(annot.logi[,isMarker])
    
## Combine information on maturation and ageing
  # notes: this scores upregulation as 1 and downregulation as 0; this assumes same column order in the up and down dataframes
  isAgeing <- grep("AstAgeing_", colnames(annot.logi.down))
  annot.combined$AstAgeing <- rowSums(annot.logi.up[,isAgeing]) - rowSums(annot.logi.down[,isAgeing])
  
  isMaturation <- grep("AstMaturation_", colnames(annot.logi.down))
  annot.combined$AstMaturation <- rowSums(annot.logi.up[,isMaturation]) - rowSums(annot.logi.down[,isMaturation])
  
## Combine information on astrocyte activation
  # notes: this scores upregulation as 1 and downregulation as 0; this assumes same column order in the up and down dataframes
  isActivationHuman <- which(grepl("AstActivation_TIC_hiPSC", colnames(annot.logi.down)) & !(grepl("hiPSC_Rep", colnames(annot.logi.down))))
  annot.combined$AstActivation_HumanOnly <- rowSums(annot.logi.up[,isActivationHuman]) - rowSums(annot.logi.down[,isActivationHuman])
  
  isActivationNonHuman <- grep("AstActivation_Mouse|AstActivation_Rat", colnames(annot.logi.down))
  annot.combined$AstActivation_NonHumanOnly <- rowSums(annot.logi.up[,isActivationNonHuman]) - rowSums(annot.logi.down[,isActivationNonHuman])
  
  annot.combined$AstActivation <- annot.combined$AstActivation_HumanOnly + annot.combined$AstActivation_NonHumanOnly
  
## And housekeeping
  isHkg <- c("Housekeeping_Lin2019", "Housekeeping_Eisenberg2013")
  annot.combined$Housekeeping <- rowSums(annot.logi[,isHkg])
  
## Add to enrichment dataframe
  run.fisher.vsHits.combined <- function(data = annot.combined, 
                                         data.column, 
                                         rbind = FALSE, 
                                         rbind.to = enrichments.combined,
                                         filt.pow = filterForWellPowered,
                                         replace = FALSE,
                                         signed = "Unsigned") {
    
    # filtering to well-powered hits
    if (filt.pow) data <- data[which(data$Gene %in% wellpowered_genes),]
    
    # get information
    total <- nrow(data)
    x <- data[,data.column]
    if (signed == "Unsigned") x <- x != 0
    if (signed == "Up") x <- x > 0
    if (signed == "Down") x <- x < 0
    
    y <- data$Hit
    
    # run stats
    f <- table(x, y) %>% fisher.test()
    
    # output stats
    out <- data.frame(Resource = data.column,
                      Total_TRUE = sum(x, na.rm = TRUE),
                      Fraction_Bg_TRUE = sum(x, na.rm = TRUE) / total,
                      Total_Hit_TRUE = sum(x & y, na.rm = TRUE),
                      Fraction_Hit_TRUE = sum(x & y, na.rm = TRUE) / sum(y, na.rm = TRUE),
                      p = f$p.value,
                      OR = f$estimate,
                      Lower = f$conf.int[1],
                      Upper = f$conf.int[2],
                      Signed = signed,
                      row.names = data.column)  
    
    if (rbind) {
      if (replace) {
        enrichments[data.column,] <- out
        return(enrichments  )
      } else {
        return(rbind(rbind.to, out))
      }
      
    } else {
      return(out)  
    }
    
  }
  
  
  
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstActivation_HumanOnly", rbind = FALSE, signed = "Unsigned")
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstActivation_HumanOnly", rbind = TRUE, signed = "Up")
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstActivation_HumanOnly", rbind = TRUE, signed = "Down")
  
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstMarker", rbind = TRUE, signed = "Up")
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "Housekeeping", rbind = TRUE, signed = "Up")
  
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstAgeing", rbind = TRUE, signed = "Unsigned")
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstAgeing", rbind = TRUE, signed = "Up")
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstAgeing", rbind = TRUE, signed = "Down")
  
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstMaturation", rbind = TRUE, signed = "Unsigned")
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstMaturation", rbind = TRUE, signed = "Up")
  enrichments.combined <- run.fisher.vsHits.combined(data.column = "AstMaturation", rbind = TRUE, signed = "Down")
  
  
  
    
################################################################################################################################ #
## Save ----
  
    
write.csv(enrichments, file = "Genes/Wellpowered/Final - Enrichments.csv")
write.csv(enrichments.combined, file = "Genes/Wellpowered/Final - Enrichments (Combined Annotations).csv", row.names = FALSE)  
write.csv(annot.logi, file = "Genes/Wellpowered/Final - Annotation Logical.csv", row.names = FALSE)
write.csv(annot.detailed, file = "Genes/Wellpowered/Final - Annotation Detailed.csv", row.names = FALSE)
write.csv(annot.logi.up, file = "Genes/Wellpowered/Final - Annotation Logical (Signed Up).csv", row.names = FALSE)
write.csv(annot.logi.down, file = "Genes/Wellpowered/Final - Annotation Logical (Signed Down).csv", row.names = FALSE)
write.csv(annot.combined, file = "Genes/Wellpowered/Final - Annotation Summed Per Phenotype Category.csv", row.names = FALSE)


save(enrichments, 
     enrichments.combined,
     annot.logi, 
     annot.logi.down,
     annot.logi.up,
     annot.combined,
     annot.detailed, file = "Genes/Wellpowered/Final.rda")

  
