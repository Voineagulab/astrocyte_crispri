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
  library(tidyr)
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
  s <- read.csv("../2_DE/Enhancers - All Characterisation.csv")

  
## Set up enhancer lists
  targets <- unique(s$Enh.Pos)
  hits <- unique(s$Enh.Pos[which(s$Hit)])
  
## Gene annotation
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  attr <- listAttributes(ensembl)
  genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"), mart = ensembl)
  

################################################################################################################################ #
## Enhancer overlap with dbSNP ----
  
## Paths
  db.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/dbSNP153_280422.bed"
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  db.out <- paste0(getwd(), "/dbSNP_Overlap.bed")
  
## Run
  # pe
  call <- paste("windowBed",
                "-a", nha_hg38,
                "-b", db.dir,
                "-w", "1000", # window of 1000bp
                ">", db.out)
  
    system(call, intern = FALSE, wait = TRUE) 
  
## Analyse
  db.res <- read.delim(db.out)
    
  # column names
  colnames(db.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "SNP.chr", "SNP.start", "SNP.end", "SNP.id")
    
  # add hit annotation
  db.res$Peak.id <- sub("_", ":", db.res$Peak.id)
  db.res$Peak.id <- sub("_", "-", db.res$Peak.id)
  db.res$Hit <- db.res$Peak.id %in% hits
  
################################################################################################################################ #
## Enhancer overlap to bulk eQTLs ----
 
## Paths
  gtex.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/GTEX_cis_eQTLs_050422_ucscDownload.bed"
  pe.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/PsychENCODE_eQTL/DER-08a_hg38_eQTL.significant.bed"
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  gtex.out <- paste0(getwd(), "/eQTL_GTEx_Overlap.bed")
  pe.out <- paste0(getwd(), "/eQTL_PE_Overlap.bed")
   
## Process eQTLs
  # GTEx 
  gtex <- read.table("../../../PublicData/GTEX_cis_eQTLs_050422_ucscDownload", sep = "\t", header = TRUE)
  gtex.bed <- gtex[,c("eqtlChrom", "eqtlStart", "eqtlEnd", "eqtlName", "geneName", "tissue")]
  write.table(gtex.bed, file = gtex.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
  
  # PSYCHEncode
  pe <- read.table("../../../PublicData/PsychENCODE_eQTL/DER-08a_hg38_eQTL.significant.txt", sep = "\t", header = TRUE)
  pe.bed <- pe[,c("SNP_chr", "SNP_start", "SNP_end", "SNP_id", "gene_id", "FDR")]
  pe.bed$gene_id <- splitter(pe.bed$gene_id, "\\.", 1)
  write.table(pe.bed, file = pe.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
 
## Run
  # pe
  call <- paste("windowBed",
                "-a", nha_hg38,
                "-b", pe.dir,
                "-w", "1000", # window of 1000bp
                ">", pe.out)
  
  system(call, intern = FALSE, wait = TRUE) 
  
  # gtex
    call <- paste("windowBed",
                "-a", nha_hg38,
                "-b", gtex.dir,
                "-w", "1000", # window of 1000bp
                ">", gtex.out)
  
    system(call, intern = FALSE, wait = TRUE) 
  
  
## Analyse PE 
  ## Annotate
    # load
    pe.res <- read.table(pe.out, sep = "\t", header = FALSE)
    
    # column names
    colnames(pe.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "eQTL.chr", "eQTL.start", "eQTL.end", "eQTL.id", "eQTL.gene", "eQTL.fdr")
    
    # add hit annotation
    pe.res$Peak.id <- sub("_", ":", pe.res$Peak.id)
    pe.res$Peak.id <- sub("_", "-", pe.res$Peak.id)
    pe.res$Hit <- pe.res$Peak.id %in% hits
    
    # add hit gene
    m <- match(pe.res$Peak.id, s$Enh.Pos[which(s$Hit)])
    pe.res$Peak.gene <- s[which(s$Hit),"Gene"][m]
  
    # convert gene symbol
    m <- match(pe.res$Peak.gene, genes$hgnc_symbol)
    pe.res$Peak.gene <- genes$ensembl_gene_id[m]
    
    # same target?
    pe.res$SameTarget <- pe.res$eQTL.gene == pe.res$Peak.gene
    
  ## Explore
    ## Rate of having eqtls
      p <- data.frame(Region = targets,
                      Hit = targets %in% hits,
                      Has.eQTL = targets %in% pe.res$Peak.id)
      
      pdf(file = "eQTL - PE - Peak Has eQTL.pdf", height = 3, width = maxw)
      ggplot(p, aes(x = Hit, fill = Has.eQTL)) +
        geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
        theme_bw() +
        scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 475)) +
        labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
        guides(fill = guide_legend(title = "Has PE eQTL")) +
        theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
              panel.grid = invis)
      dev.off()
    
    ## Same target
      p <- p[which(p$Hit),]
    
      x <- pe.res$Peak.id[which(pe.res$SameTarget)] # has same target
      p$SameTarget <- p$Region %in% x
    
      pdf(file = "eQTL - PE - Peak And eQTL Share Target.pdf", height = 3, width = maxw)
      ggplot(p, aes(x = Has.eQTL, fill = SameTarget)) +
        geom_bar(colour = "black", width = 0.7, position = position_stack()) +
        theme_bw() +
        scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
        labs(y = "Count of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
        guides(fill = guide_legend(title = "Same Target Gene")) +
        theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis)
      dev.off()
      
    
## Analyse GTEx 
  ## Annotate
    # load
    gtex.res <- read.table(gtex.out, sep = "\t", header = FALSE)
    
    # column names
    colnames(gtex.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "eQTL.chr", "eQTL.start", "eQTL.end", "eQTL.id", "eQTL.gene", "eQTL.tissue")
    
    # add hit annotation
    gtex.res$Peak.id <- sub("_", ":", gtex.res$Peak.id)
    gtex.res$Peak.id <- sub("_", "-", gtex.res$Peak.id)
    gtex.res$Hit <- gtex.res$Peak.id %in% hits
    
    # add hit gene
    m <- match(gtex.res$Peak.id, s$Enh.Pos[which(s$Hit)])
    gtex.res$Peak.gene <- s[which(s$Hit),"Gene"][m]
  
    # same target?
    gtex.res$SameTarget <- gtex.res$eQTL.gene == gtex.res$Peak.gene
    
  ## Explore
    ## Rate of having eqtls
      p <- data.frame(Region = targets,
                      Hit = targets %in% hits,
                      Has.eQTL = targets %in% gtex.res$Peak.id)
      
      pdf(file = "eQTL - GTEx - Peak Has eQTL.pdf", height = 3, width = maxw)
      ggplot(p, aes(x = Hit, fill = Has.eQTL)) +
        geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
        theme_bw() +
        scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 550)) +
        labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
        guides(fill = guide_legend(title = "Has GTEx eQTL")) +
        theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
              panel.grid = invis)
      dev.off()
    
    ## Same target
      p <- p[which(p$Hit),]
    
      x <- gtex.res$Peak.id[which(gtex.res$SameTarget)] # has same target
      p$SameTarget <- p$Region %in% x
    
      pdf(file = "eQTL - GTEx - Peak And eQTL Share Target.pdf", height = 3, width = maxw)
      ggplot(p, aes(x = Has.eQTL, fill = SameTarget)) +
        geom_bar(colour = "black", width = 0.7, position = position_stack()) +
        theme_bw() +
        scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
        labs(y = "Count of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
        guides(fill = guide_legend(title = "Same Target Gene")) +
        theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis)
      dev.off()

    ## When there is an overlap, what is the tissue?
      gtex.res$eQTL.tissue.type <- "."
      g <- grep("Brain", gtex.res$eQTL.tissue)
      gtex.res$eQTL.tissue.type[g] <- "Brain"
      gtex.res$eQTL.tissue.type[-g] <- "Non-brain"
      gtex.res$eQTL.tissue.type[grep("Nerve", gtex.res$eQTL.tissue)] <- "Nerve Tibial"
      
      pdf(file = "eQTL - GTEx - Tissue Types.pdf", height = 3, width = maxw)
      pA <- ggplot(gtex.res, aes(x = Hit, fill = eQTL.tissue.type)) +
        geom_bar(colour = "black", width = 0.7, position = position_stack()) +
        theme_bw() +
        scale_fill_manual(values = carto_pal(3, "Pastel")) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 1200)) +
        labs(y = "Fraction of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
        guides(fill = guide_legend(title = "Same Target Gene")) +
        theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
              legend.position = c(0.75, 0.8))
      
      pB <- ggplot(gtex.res, aes(x = Hit, fill = eQTL.tissue.type)) +
        geom_bar(colour = "black", width = 0.7, position = position_fill()) +
        theme_bw() +
        scale_fill_manual(values = carto_pal(3, "Pastel")) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
        labs(y = "Count of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
        guides(fill = guide_legend(title = "eQTL Tissue Type")) +
        theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
              legend.position = "none")
      
      plot_grid(pB, pA, rel_widths = c(1, 1.4))
      dev.off()
      
      ## An alternative view of the above: for a given peak, are any of its eQTLs in the brain?
      x <- unique(gtex.res$Peak.id[grep("Brain", gtex.res$eQTL.tissue)])
      p <- data.frame(Region = targets,
                      Hit = targets %in% hits,
                      Has.eQTL = targets %in% gtex.res$Peak.id,
                      Has.brain.eQTL = targets %in% x)
      
      p <- p[which(p$Has.eQTL),]
      
      pdf(file = "eQTL - GTEx - Peak Has eQTL In Brain.pdf", height = 3, width = maxw)
      pA <- ggplot(p, aes(x = Hit, fill = Has.brain.eQTL)) +
        geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
        theme_bw() +
        scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 300)) +
        labs(y = "Fraction of Peaks in Category", x = "Enhancer With eQTL Has Is Hit") +
        guides(fill = guide_legend(title = "At Least 1 GTEx\neQTL In Brain")) +
        theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
              panel.grid = invis)
      
      pB <- ggplot(p, aes(x = Hit, fill = Has.brain.eQTL)) +
        geom_bar(colour = "black", width = 0.7, position = position_fill()) +
        theme_bw() +
        scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
        labs(y = "Count of Peaks in Category", x = "Enhancer With eQTL Has Is Hit") +
        guides(fill = guide_legend(title = "At Least 1 GTEx eQTL In Brain")) +
        theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
              legend.position = "none")
      
      plot_grid(pB, pA)
      dev.off()
      
      ## Save
      write.csv(gtex.res, file = "eQTL - GTEx - Data.csv")
      write.csv(pe.res, file = "eQTL - PE - Data.csv")
      
  
################################################################################################################################ #
## Enhancer overlap to snRNA-seq eQTLs ----
 
## Paths
  bryois.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Bryois2021/"
  bryois.bed <- paste0(bryois.dir, "Processed.bed")
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  bryois.out <- paste0(getwd(), "/eQTL_snAst_Overlap.bed")
   
## Process eQTLs
  # read in
  l <- list.files(bryois.dir)
  l <- l[-c(23:24)] # meta
  
  bryois.raw <- lapply(l, function(x) {
    print(x)
    y <- read.table(paste0(bryois.dir, x), sep = " ")
    colnames(y) <- c("Gene", "SNP", "Distance", "P", "Beta")
    y$Gene <- splitter(y$Gene, "_", 1)
    y <- y[,c("Gene", "SNP", "P")]
    return(y)
  })
  
  names(bryois.raw) <- gsub("Astrocytes.quantile.txt.gz.", "", l) %>% splitter(":", 1) %>% paste0("Chr", .)
  
  bryois.raw <- do.call("rbind", bryois.raw)

  # filter to significant hits
  bryois.raw$FDR <- p.adjust(bryois.raw$P, method = "fdr")
  keep <- which(bryois.raw$FDR < 0.05)
  bryois.raw <- bryois.raw[keep,]
  
  # add hg38 position
  posit <- read.delim(paste0(bryois.dir, "snp_pos.txt"))
  m <- match(bryois.raw$SNP, posit$SNP)
  bryois.raw$Chr <- posit$chr[m]
  bryois.raw$Start <- posit$pos_hg38[m]
  bryois.raw$End <- posit$pos_hg38[m]
  
  # write to disk as bed
  bryois.raw <- bryois.raw[,c("Chr", "Start", "End", "SNP", "Gene", "FDR")]
  write.table(bryois.raw, bryois.bed, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  
## Overlap
  call <- paste("windowBed",
              "-a", nha_hg38,
              "-b", bryois.bed,
              "-w", "1000", # window of 1000bp
              ">", bryois.out)

   system(call, intern = FALSE, wait = TRUE) 
  
  
## Read in and annotate 
  # load
  bryois.res <- read.table(bryois.out, sep = "\t", header = FALSE)
  
  # column names
  colnames(bryois.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "eQTL.chr", "eQTL.start", "eQTL.end", "eQTL.id", "eQTL.gene", "eQTL.fdr")
  
  # add hit annotation
  bryois.res$Peak.id <- sub("_", ":", bryois.res$Peak.id)
  bryois.res$Peak.id <- sub("_", "-", bryois.res$Peak.id)
  bryois.res$Hit <- bryois.res$Peak.id %in% hits
  
  # add hit gene
  m <- match(bryois.res$Peak.id, s$Enh.Pos[which(s$Hit)])
  bryois.res$Peak.gene <- s[which(s$Hit),"Gene"][m]

  # same target?
  bryois.res$SameTarget <- bryois.res$eQTL.gene == bryois.res$Peak.gene
  
## Little to report, so just save csv
  write.csv(bryois.res, file = "eQTL - snAst - Data.csv")
    
  
        
################################################################################################################################ #
## Enhancer overlap to clinical variants from DisGeNet ----
  
  
## Paths
  disg.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/DisGeNet_070422/all_variant_disease_associations.bed"
  disg.out <- paste0(getwd(), "/Disgenet_Variant_Overlap.bed")
  
## Process table
  disg <- read.delim("../../../PublicData/DisGeNet_070422/all_variant_disease_associations.tsv", sep = "\t", header = TRUE)
  disg.bed <- disg[,c("chromosome", "position", "position", "snpId", "diseaseName","diseaseType")]
  disg.bed$chromosome <- paste0("chr", disg.bed$chromosome)
  write.table(disg.bed, file = disg.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
  
## Run
  call <- paste("windowBed",
                "-a", nha_hg38,
                "-b", disg.dir,
                "-w", "1000", # window of 1000bp
                ">", disg.out)
  
  system(call, intern = FALSE, wait = TRUE) 
  
## Read
  disg.res <- read.delim(disg.out, sep = "\t", header = FALSE)
  colnames(disg.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "Disg.chr", "Disg.start", "Disg.end", "Disg.id", "Disg.disease", "Disg.class")
  
## Output variants for hit peaks
  disg.res$Peak.id <- sub("_", ":", disg.res$Peak.id)
  disg.res$Peak.id <- sub("_", "-", disg.res$Peak.id)
  disg.res$Hit <- disg.res$Peak.id %in% enh.hits
  disg.hits <- disg.res[which(disg.res$Hit),]
  
## Variant enrichment within peaks
  p <- data.frame(Region = targets,
                  Hit = targets %in% hits,
                  Has.disg = targets %in% disg.res$Peak.id)
  
  pdf(file = "Disgenet - Peak Has Variant.pdf", height = 3, width = maxw)
  pA <- ggplot(p, aes(x = Hit, fill = Has.disg)) +
    geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
    theme_bw() +
    scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 820)) +
    labs(y = "Fraction of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
    guides(fill = guide_legend(title = "Has Disgenet Variant")) +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
          panel.grid = invis)
  
  pB <- ggplot(p, aes(x = Hit, fill = Has.disg)) +
    geom_bar(colour = "black", width = 0.7, position = position_fill()) +
    theme_bw() +
    scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
    labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
    guides(fill = guide_legend(title = "Has Disgenet Variant")) +
    theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
          panel.grid = invis)
  
  plot_grid(pB, pA)
  dev.off()
  
## Output
  write.csv(disg.res, file = "Disgenet - Peak Variant Overlap.csv")
  write.csv(disg.hits, file = "Disgenet - Peak Variant Overlap (Hits Only).csv")
  
  
################################################################################################################################ #
## Enhancer overlap to GWAS ----
  
  
## GWAS Catalogue
  
  ## Paths
    gwas.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/gwasCatalog_050422_ucscDownload.bed"
    gwas.out <- paste0(getwd(), "/GWAS_Catalogue_Variant_Overlap.bed")
  
  ## Process table
    gwas <- read.delim("../../../PublicData/gwasCatalog_050422_ucscDownload", sep = "\t", header = TRUE)
    gwas.bed <- gwas[,c("chrom", "chromStart", "chromEnd", "name", "trait","genes")]
    gwas.bed$genes[gwas.bed$genes == ""] <- "(Blank)"
    write.table(gwas.bed, file = gwas.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
    
  ## Run
    call <- paste("windowBed",
                  "-a", nha_hg38,
                  "-b", gwas.dir,
                  "-w", "1000", # window of 1000bp
                  ">", gwas.out)
    
    system(call, intern = FALSE, wait = TRUE) 
    
  ## Read
    gwas.res <- read.delim(gwas.out, sep = "\t", header = FALSE)
    colnames(gwas.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "GWAS.chr", "GWAS.start", "GWAS.end", "GWAS.id", "GWAS.disease", "GWAS.gene")
    
  ## Output variants for hit peaks
    gwas.res$Peak.id <- sub("_", ":", gwas.res$Peak.id)
    gwas.res$Peak.id <- sub("_", "-", gwas.res$Peak.id)
    gwas.res$Hit <- gwas.res$Peak.id %in% hits
    gwas.hits <- gwas.res[which(gwas.res$Hit),]
    
  ## Variant enrichment within peaks
    p <- data.frame(Region = targets,
                    Hit = targets %in% hits,
                    Has.gwas = targets %in% gwas.res$Peak.id)
    
    pdf(file = "GWAS Catalogue - Peak Has Variant.pdf", height = 3, width = maxw)
    pA <- ggplot(p, aes(x = Hit, fill = Has.gwas)) +
      geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
      theme_bw() +
      scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 820)) +
      labs(y = "Fraction of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
      guides(fill = guide_legend(title = "Has Disgenet Variant")) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
            panel.grid = invis)
    
    pB <- ggplot(p, aes(x = Hit, fill = Has.gwas)) +
      geom_bar(colour = "black", width = 0.7, position = position_fill()) +
      theme_bw() +
      scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
      labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
      guides(fill = guide_legend(title = "Has GWAS Variant")) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
            panel.grid = invis)
    
    plot_grid(pB, pA)
    dev.off()
    
  ## Output for hit peaks
    write.csv(gwas.res, file = "GWAS Catalogue - Peak Variant Overlap.csv")
    write.csv(gwas.hits, file = "GWAS Catalogue - Peak Variant Overlap (Hits Only).csv")
  
## MAGMA
  ## Load functions
    source("/mnt/Data0/PROJECTS/GWAS_Enrichment/GAVIN/Scripts/0_Setup.R")
    
  ## Write enh peaks as a MAGMA loc file
    # load
    magma.peaks <- read.table(nha_hg38, sep = "\t")
    magma.peaks$V4 <- sub("_", ":", magma.peaks$V4)
    magma.peaks$V4 <- sub("_", "-", magma.peaks$V4)
    magma.peaks$V5 <- magma.peaks$V4 %in% hits
    magma.peaks$V6 <- "."
    
    # sort in terminal
    unsorted.magma <- "MAGMA/MAGMA_peaks_unsorted.bed"
    sorted.magma <- "MAGMA/MAGMA_peaks_sorted.bed"
    write.table(magma.peaks, file = unsorted.magma, 
                sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    call <- paste("sort",
                  "-k1,1",
                  "-k2,2n",
                  unsorted.magma,
                  ">",
                  sorted.magma, 
                  sep = " ")
    
    system(call, intern = FALSE, wait = TRUE)  
    
    
    # read in and convert to a file MAGMA can parse
    magma.peaks.sorted <- read.table(sorted.magma, sep = "\t")
    magma.bg <- data.frame(Peak = magma.peaks$V4,
                           Chromosome = gsub("chr", "", magma.peaks$V1),
                           Start = magma.peaks$V2,
                           End = magma.peaks$V3,
                           Strand = ".",
                           Hit = magma.peaks$V5)
    
    write.table(magma.bg, file = "MAGMA/MAGMA_peaks_bg.loc",
                quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    
  ## Annotate peaks with annotation with SNPs
    magma.annotate(bg = "MAGMA/MAGMA_peaks_bg.loc",
                   output = "MAGMA/MAGMA_bg_annotation")
  
  ## GWAS background risk
    # function  
    run.magma.background <- function(annot, gwas, gwas.columns = c("SNP", "P"), N, output, snps = paste0(paths$snps, "g1000_eur")) {
      
      # convert arguments into a terminal call
      pval.statement <- paste0(gwas,
                               " ", "N=", N,
                               " ", "use=", gwas.columns[1], ",", gwas.columns[2])
      
      # output <- paste0("Background_Analysis/", output)
      # annot <- paste0("/mnt/Data0/PROJECTS/Lister/Results/Annotation/", annot)
      
      call <- paste(paths$magma,
                    "--bfile", snps,
                    "--gene-annot", annot,
                    "--pval", pval.statement,
                    "--out", output,
                    sep = " ")
      
      # run
      system(call, intern = FALSE, wait = TRUE)
      return("Done")
    }
    
    # run
    for (j in names(n)) {
      run.magma.background(annot = "MAGMA/MAGMA_bg_annotation.genes.annot", gwas = paths[[j]], gwas.columns = cols[[j]], N = n[[j]], output = paste0("MAGMA/GWAS_Backgrounds/", j, "_bg"))
    }
    
    
  ## Query for increased risk in hit set
    # write query set
    x <- paste0("Hit ", magma.bg$Peak[which(magma.bg$Hit)])
    write.table(x, file = "MAGMA/Hit_Enrichment/SetInput.sets", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # function
    run.magma.queryset <- function(gwas, queryset = "MAGMA/Hit_Enrichment/SetInput.sets") {
      
      # convert arguments into a terminal call
      # queryset <- paste0("/mnt/Data0/PROJECTS/Lister/Data/MAGMA_QuerySets/", queryset)
      output <- paste0("MAGMA/Hit_Enrichment/", gwas)
      gwas <- paste0("MAGMA/GWAS_Backgrounds/", gwas, "_bg.genes.raw")
      
      
      call <- paste(paths$magma,
                    "--gene-results", gwas,
                    "--set-annot", queryset,  "col=2,1", # the second statement tells magma that the queryset is in a column format
                    "--out", output,
                    sep = " ")
      
      # run
      system(call, intern = FALSE, wait = TRUE)
      return("Done")
    }
    
    # run
    for (j in names(n)) { # for GWAS j
      run.magma.queryset(gwas = j)
    }
    
  
  
################################################################################################################################ #
## A combined variant annotation file! ---- 

## Here, you aim to create a table which summarises the variants for each enhancer
  
## Create
  out <- data.frame(Region = targets,
                    Hit = targets %in% hits)
  
## Function
  add.annot <- function(x = out,  overlaps, name, combine = FALSE) {
    # is there an overlap
    x[,paste0(name, "_Overlap")] <- x$Region %in% overlaps[,4]
    
    # if overlap, what is its nature?
    x[,paste0(name, "_Info")] <- ""
    m <- match(overlaps[,4], x$Region)
    
    if (combine) {
      for (j in 1:length(m)) {
        x[m[j],paste0(name, "_Info")] <- paste0(x[m[j],paste0(name, "_Info")], " | ", paste0(overlaps[j,9]), " in ", overlaps[j,10])
      }
    } else {
      for (j in 1:length(m)) {
        x[m[j],paste0(name, "_Info")] <- paste0(x[m[j],paste0(name, "_Info")], " | ", overlaps[j,9])
      }
    }
    
    
    
    return(x)
  }
  
out <- add.annot(out, overlaps = gwas.res, name = "gwas")
out <- add.annot(out, overlaps = disg.res, name = "disgenet")
out <- add.annot(out, overlaps = pe.res, name = "eQTLpe")
out <- add.annot(out, overlaps = gtex.res, name = "eQTLgtex", combine = TRUE)
out <- add.annot(out, overlaps = bryois.res, name = "eQTLbryois")

write.csv(out, "Variant Annotation To Peaks.csv")
  

      
      
      
      