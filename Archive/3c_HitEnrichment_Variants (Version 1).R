## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----


## Generic
# rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Variants")
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
  source("../../../Scripts/Functions.R")
  load("../../../Data/Preprocessed/NHA Pooled.rda")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
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
  res.final <- read.csv("../../2_DE/Enh/Results Final.csv")

  
## Set up enhancer lists
  targets <- unique(res.final$Enh.Pos)
  hits <- unique(res.final$Enh.Pos[which(s$Sceptre.FDR < 0.1)])
  
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
  db.res <- read.delim(db.out, header = FALSE)
    
  # column names
  colnames(db.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "SNP.chr", "SNP.start", "SNP.end", "SNP.id")
    
  # add hit annotation
  db.res$Peak.id <- sub("_", ":", db.res$Peak.id)
  db.res$Peak.id <- sub("_", "-", db.res$Peak.id)
  # db.res$Hit <- db.res$Peak.id %in% hits
  
## Include LD variants
  # calculate LD for all 
  
  # requires ~24h
  LDproxy_batch(snp = unique(db.res$SNP.id), # unique is necessary as these peaks have a window appended
                  pop = "EUR", # european
                  r2d="r2", 
                  append = TRUE,
                  genome_build = "grch38",
                  token = "5988cfbf2737")
  
  # read in
  all.ld.snps <- read.table("combined_query_snp_list_grch38.txt", sep = "\t", header = TRUE, row.names = NULL)
  all.ld.snps <- all.ld.snps[,c(2,3,4,6,7,9)]
  colnames(all.ld.snps) <- c("Query_SNP", "Linked_SNP", "Position", "MAF", "Distance", "R2")
  
  # stats
  length(unique(all.ld.snps$Query_SNP)) / length(unique(db.res$SNP.id)) # 6982 / 9854 = 0.7
  
  # filter
  ld.snps <- all.ld.snps[which(all.ld.snps$R2 > 0.8 & all.ld.snps$Distance != 0),] # filters to r2 > 0.8, and not identical
  
  ## Annotate
    snp.annotation <- list()
  
    for (j in targets) {
      print(j)
      
      # skip if no snps nearby. this holds for one enhancer at a 1kb window
      if ( !( j %in% db.res$Peak.id)) {
        snp.annotation[[j]] <- ("None") 
        next
        }
      
      # collect these snps in the dataframe x
      x <- db.res[which(db.res$Peak.id == j),]
      peak.start <- x$Peak.start[1]; peak.end <- x$Peak.end[1]
      peak.centre <- mean(c(peak.start, peak.end))
      x <- x[,c(4,8,7)]
      colnames(x) <- c("Peak", "SNP", "Pos")
      
      x$R2 <- NA
      x$Category <- "."
      
      # annotate as overlapping or nearby
      w <- ((x$Pos > peak.start) & (x$Pos < peak.end))
      if (any(w)) x$Category[w] <- "Overlapping"
      if (any(!(w))) {
        x$Category[!(w)] <- "Nearby"
      }
      
      # annotate with distance rank
      x$DistanceFromCentre <- abs(x$Pos - peak.centre)
      x$DistanceRank <- rank(x$DistanceFromCentre) 
      
      ## Collect snps in ld
        # start by searching for ld for SNPs marked "Overlapping"
        if (any(w)) {
          if (any(x$SNP[w] %in% ld.snps$Query_SNP)) { # when these "Overlapping" SNPs have at least one SNP in LD
            
            y <- ld.snps[which(ld.snps$Query_SNP %in% x$SNP[w]),] # collect
            
            # reformat
            z <- data.frame(Peak = j,
                            SNP = y$Linked_SNP, 
                            Pos = splitter(y$Position, ":", 2) %>% as.numeric(),
                            R2 = y$R2,
                            Category = "LD_ToOverlapping",
                            DistanceFromCentre = NaN,
                            DistanceRank = NA)
            
            z$DistanceFromCentre <- (z$Pos - peak.centre) %>% abs()
            
            # remove LD SNPs which are also nearby! match based on pos rather than ID, as sometimes ID is "."
            if (any(z$Pos %in% x$Pos)) {
              dup <- which(z$Pos %in% x$Pos)
              if (length(dup) == length(w)) { # if removing as many rows as there are in the dataframe, may cause error. instead, finish loop
                snp.annotation[[j]] <- x
                next
              }
              z <- z[-which(z$Pos %in% x$Pos),]  
            }
            
            # combine
            x <- rbind(x, z)
            
          } 
        }
      
      # but, if either A) none Overlapping, or B) those Overlapping have no linked SNPs, 
      if (!(any(w)) | !(any(x$SNP[w] %in% ld.snps$Query_SNP))) { 
        
        # do any of nearby SNPs have linked SNPs?
        if (any(x$SNP[!(w)] %in% ld.snps$Query_SNP)) {
          y <- which(x$SNP %in% ld.snps$Query_SNP) # get the index of these guides
          y <- y[which.min(x$DistanceRank[y])] # which of these is the closest?
          y <- ld.snps[which(ld.snps$Query_SNP == x$SNP[y]),]
          
          # collect LD to this SNP noted in y, and reformat
          z <- data.frame(Peak = j,
                            SNP = y$Linked_SNP, 
                            Pos = splitter(y$Position, ":", 2) %>% as.numeric(),
                            R2 = y$R2,
                            Category = "LD_ToNearby",
                            DistanceFromCentre = NaN,
                            DistanceRank = NA)
            
            z$DistanceFromCentre <- (z$Pos - peak.centre) %>% abs()
            
            # remove LD SNPs which are also nearby! match based on pos rather than ID, as sometimes ID is "."
            if (any(z$Pos %in% x$Pos)) {
              dup <- which(z$Pos %in% x$Pos)
              if (length(dup) == length(w)) { # if removing as many rows as there are in the dataframe, may cause error. instead, finish loop
                snp.annotation[[j]] <- x
                next
              }
              z <- z[-which(z$Pos %in% x$Pos),]  
            }
            
            # combine
            x <- rbind(x, z)
          
          
        }
        
      }
        
      # return
      snp.annotation[[j]] <- x
      
    }
    
  snp.annotation <- do.call("rbind", snp.annotation)
    
## Other annotation 
  # remove the "none"?
  snp.annotation <- snp.annotation[-which(snp.annotation$Category == "None"),]
  
  # add the enhancer id
  m <- match(snp.annotation$Peak, guides$TargetCoord)
  snp.annotation$Enh <- guides$TargetID[m] 
  
  # is it a hit enhancer?
  snp.annotation$HitPermissive <- snp.annotation$Peak %in% res.final$Enh.Pos[which(res.final$HitPermissive)]
  snp.annotation$HitCore <- snp.annotation$Peak %in% res.final$Enh.Pos[which(res.final$HitCore)]
    
  # is the snp near the summit  of the peak? This is for the MPRA...
  summits <- read.delim("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/MPRA/Results/0_Design/Peak_Centres_259bp.bed", header = FALSE)  
  summits$V4 <- sub("_", ":", summits$V4) %>% sub("_", "-", .) 
  
  snp.annotation$WithinEnhSummit <- "."
  
  for (j in unique(snp.annotation$Peak)) {
    print(j)
    
    # summit coordinates for the given enh
    summit.start <- summits$V2[summits$V4 == j]
    summit.end <- summits$V3[summits$V4 == j]
    
    # row indices of snps for given enh
    w <- which(snp.annotation$Peak == j)
    
    # snp coordinate within summit window
    snp.annotation$WithinEnhSummit[w] <- (snp.annotation$Pos[w] >= summit.start & snp.annotation$Pos[w] <= summit.end )
  }
  
  # save
  write.table(snp.annotation, file = "SNP Annotation.txt", sep = "\t", quote = FALSE, col.names = TRUE)
  

################################################################################################################################ #
## SNPs that alter TF binding sites ----
  
## Use the SNP2TFBS webtool:
  # https://ccg.epfl.ch/snp2tfbs/snpselect.php
  # https://academic.oup.com/nar/article/45/D1/D139/2605727?login=true
  
## Prepare for upload
  # snp.annotation <- read.delim("SNP Annotation.txt")
  o <- snp.annotation[which(snp.annotation$Category %in% c("Overlapping") & snp.annotation$HitPermissive),]
  write.table(unique(o$SNP), file = "SNP2TFBS Webtool Input.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
## Download and wrangle
  snp2tfbs <- read.delim("SNP2TFBS Webtool Output.txt", header = FALSE)
  colnames(snp2tfbs) <- c("hg19chr", "hg19start", "hg19end", "Ref", "Alt", "Match", "rsID")
  snp2tfbs <- snp2tfbs[,c("rsID", "Ref", "Alt", "Match")]
  
  # map to enh
  # snp2tfbs <- snp2tfbs[which(snp2tfbs$rsID %in% x$SNP[which(x$Category == "Overlapping")]),]
  m <- match(snp2tfbs$rsID, o$SNP)
  snp2tfbs$Peak <- o$Peak[m]
  snp2tfbs$Enh <- o$Enh[m]
  snp2tfbs$Category <- o$Category[m]
  snp2tfbs$WithinEnhSummit <- o$WithinEnhSummit[m]
  
  # split the match column
  split.col <- strsplit(snp2tfbs$Match, ";")
  snp2tfbs$TF.number <- sapply(split.col, "[", 1) %>% gsub("MATCH=", "", .) %>% as.numeric() # number of TF-affecting variant matches for a given SNP
  snp2tfbs$TF.name <- sapply(split.col, "[", 2) %>% gsub("TF=", "", .)
  snp2tfbs$TF.scorediff <- sapply(split.col, "[", 3) %>% gsub("ScoreDiff=", "", .) # %>% as.numeric()
  
  # create multiple rows for snps with multiple affected tfs
  # mult.snps <- snp2tfbs[which(snp2tfbs$TF.number > 1),]  # get those with > 1 match, and split into many rows via the below loop
  # snp2tfbs <- snp2tfbs[-which(snp2tfbs$TF.number > 1),]
  # for (j in 1:nrow(mult.snps)) {
  #   x <- mult.snps[j,]
  #   mult.tf <- strsplit(x$TF.name, ",") %>% do.call("c", .)
  #   mult.score <- strsplit(x$TF.scorediff, ",") %>% do.call("c", .)
  # 
  #   y <- list()
  #   for (k in 1:length(mult.tf)) {
  #     y[[k]] <- x
  #     y[[k]]$TF.name <- mult.tf[k]
  #     y[[k]]$TF.scorediff <- mult.score[k]
  #   }
  #   y <- do.call("rbind", y)
  #   snp2tfbs <- rbind(snp2tfbs, y)
  # }
  # 
  
  snp2tfbs$HitCore <- snp2tfbs$HitPermissive <- "."
  
  # add target gene
  for (j in 1:nrow(snp2tfbs)) {
    x <- res.final[which(res.final$Enh == snp2tfbs$Enh[j]),]
    
    snp2tfbs$HitPermissive[j] <- paste(x$Gene[which(x$HitPermissive)], collapse = ",")
    snp2tfbs$HitCore[j] <- paste(x$Gene[which(x$HitCore)], collapse = ",")
    
  }
  
## Is SNP2TFBS supported by ReMap data?
  # since SNP2TFBS is based on _predicted_ binding sites, is there ChIP-seq evidence for binding?
  
  # load ReMap
  remap <- read.csv("../ReMap Processed.csv")
  
  # loop each row of snp2tfbs
  snp2tfbs$ReMap <- "."
  
  for (j in 1:nrow(snp2tfbs)) {
    # get row
    print(j)
    x <- snp2tfbs[j,]
    
    # check the TF(s)
    tf <- x$TF.name
    
    if (grepl(",", tf)) { # many TFs!
      # vector of tfs
      tf <- strsplit(tf, ",") %>% do.call("c", .)
      
      # get the SNP's position
      m <- match(x$rsID, snp.annotation$SNP)
      snp.chr <- snp.annotation$Peak[m] %>% splitter(":", 1) 
      snp.pos <- snp.annotation$Pos[m] %>% as.numeric()
      
      # loop check for each tf
      tf.checker <- list()
      for (k in tf) {
        if (!(k %in% remap$TF)) {
          tf.checker[[k]] <- FALSE
          next
        }
        
        # overlap with tf peaks
        y <- remap$TF.Coord[which(remap$TF == k)]
        y <- sub(":", "-", y)
        tf.peaks <- data.frame(chr = splitter(y, "-", 1),
                               start = splitter(y, "-", 2),
                               end = splitter(y, "-", 3))
        
        tf.checker[[k]] <- any(snp.chr == tf.peaks$chr & # same chr
                                   snp.pos >= tf.peaks$start & # greater or equal to the tf peak's start
                                   snp.pos <= tf.peaks$end) # lesser or equal to the tf peak's end
      }
      
      # combine output
      snp2tfbs$ReMap[j] <- do.call("c", tf.checker) %>% paste(., collapse = ",")
      
      
    } else { # 1 TF
      if (!(tf %in% remap$TF)) {
        snp2tfbs$ReMap[j] <- FALSE
        next
      }
      
      # get the SNP's position
      m <- match(x$rsID, snp.annotation$SNP)
      snp.chr <- snp.annotation$Peak[m] %>% splitter(":", 1) 
      snp.pos <- snp.annotation$Pos[m] %>% as.numeric()
      
      # overlap with tf peaks
      x <- remap$TF.Coord[which(remap$TF == tf)]
      x <- sub(":", "-", x)
      tf.peaks <- data.frame(chr = splitter(x, "-", 1),
                             start = splitter(x, "-", 2),
                             end = splitter(x, "-", 3))
      
      snp2tfbs$ReMap[j] <- any(snp.chr == tf.peaks$chr & # same chr
                                 snp.pos >= tf.peaks$start & # greater or equal to the tf peak's start
                                 snp.pos <= tf.peaks$end) # lesser or equal to the tf peak's end
    } 
    
    
    
  }
  
  
  
## Save
  write.csv(snp2tfbs, file = "SNP2TFBS Processed.csv")
  
## Add to the annotation dataframe
  snp.annotation$SNP2TFBS <- snp.annotation$SNP %in% snp2tfbs$rsID
  
# ################################################################################################################################ #
# ## Enhancer overlap to bulk eQTLs ----
#  
# ## Paths
#   gtex.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/GTEX_cis_eQTLs_050422_ucscDownload.bed"
#   pe.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/PsychENCODE_eQTL/DER-08a_hg38_eQTL.significant.bed"
#   nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
#   gtex.out <- paste0(getwd(), "/eQTL_GTEx_Overlap.bed")
#   pe.out <- paste0(getwd(), "/eQTL_PE_Overlap.bed")
#    
# ## Process eQTLs
#   # GTEx 
#   gtex <- read.table("../../../PublicData/GTEX_cis_eQTLs_050422_ucscDownload", sep = "\t", header = TRUE)
#   gtex.bed <- gtex[,c("eqtlChrom", "eqtlStart", "eqtlEnd", "eqtlName", "geneName", "tissue")]
#   write.table(gtex.bed, file = gtex.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
#   
#   # PSYCHEncode
#   pe <- read.table("../../../PublicData/PsychENCODE_eQTL/DER-08a_hg38_eQTL.significant.txt", sep = "\t", header = TRUE)
#   pe.bed <- pe[,c("SNP_chr", "SNP_start", "SNP_end", "SNP_id", "gene_id", "FDR")]
#   pe.bed$gene_id <- splitter(pe.bed$gene_id, "\\.", 1)
#   write.table(pe.bed, file = pe.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#  
# ## Run
#   # pe
#   call <- paste("windowBed",
#                 "-a", nha_hg38,
#                 "-b", pe.dir,
#                 "-w", "1000", # window of 1000bp
#                 ">", pe.out)
#   
#   system(call, intern = FALSE, wait = TRUE) 
#   
#   # gtex
#     call <- paste("windowBed",
#                 "-a", nha_hg38,
#                 "-b", gtex.dir,
#                 "-w", "1000", # window of 1000bp
#                 ">", gtex.out)
#   
#     system(call, intern = FALSE, wait = TRUE) 
#   
#   
# ## Analyse PE 
#   ## Annotate
#     # load
#     pe.res <- read.table(pe.out, sep = "\t", header = FALSE)
#     
#     # column names
#     colnames(pe.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "eQTL.chr", "eQTL.start", "eQTL.end", "eQTL.id", "eQTL.gene", "eQTL.fdr")
#     
#     # add hit annotation
#     pe.res$Peak.id <- sub("_", ":", pe.res$Peak.id)
#     pe.res$Peak.id <- sub("_", "-", pe.res$Peak.id)
#     pe.res$Hit <- pe.res$Peak.id %in% hits
#     
#     # add hit gene
#     m <- match(pe.res$Peak.id, s$Enh.Pos[which(s$Hit)])
#     pe.res$Peak.gene <- s[which(s$Hit),"Gene"][m]
#   
#     # convert gene symbol
#     m <- match(pe.res$Peak.gene, genes$hgnc_symbol)
#     pe.res$Peak.gene <- genes$ensembl_gene_id[m]
#     
#     # same target?
#     pe.res$SameTarget <- pe.res$eQTL.gene == pe.res$Peak.gene
#     
#   ## Explore
#     ## Rate of having eqtls
#       p <- data.frame(Region = targets,
#                       Hit = targets %in% hits,
#                       Has.eQTL = targets %in% pe.res$Peak.id)
#       
#       pdf(file = "eQTL - PE - Peak Has eQTL.pdf", height = 3, width = maxw)
#       ggplot(p, aes(x = Hit, fill = Has.eQTL)) +
#         geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
#         theme_bw() +
#         scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 475)) +
#         labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
#         guides(fill = guide_legend(title = "Has PE eQTL")) +
#         theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
#               panel.grid = invis)
#       dev.off()
#     
#     ## Same target
#       p <- p[which(p$Hit),]
#     
#       x <- pe.res$Peak.id[which(pe.res$SameTarget)] # has same target
#       p$SameTarget <- p$Region %in% x
#     
#       pdf(file = "eQTL - PE - Peak And eQTL Share Target.pdf", height = 3, width = maxw)
#       ggplot(p, aes(x = Has.eQTL, fill = SameTarget)) +
#         geom_bar(colour = "black", width = 0.7, position = position_stack()) +
#         theme_bw() +
#         scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
#         labs(y = "Count of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
#         guides(fill = guide_legend(title = "Same Target Gene")) +
#         theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis)
#       dev.off()
#       
#     
# ## Analyse GTEx 
#   ## Annotate
#     # load
#     gtex.res <- read.table(gtex.out, sep = "\t", header = FALSE)
#     
#     # column names
#     colnames(gtex.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "eQTL.chr", "eQTL.start", "eQTL.end", "eQTL.id", "eQTL.gene", "eQTL.tissue")
#     
#     # add hit annotation
#     gtex.res$Peak.id <- sub("_", ":", gtex.res$Peak.id)
#     gtex.res$Peak.id <- sub("_", "-", gtex.res$Peak.id)
#     gtex.res$Hit <- gtex.res$Peak.id %in% hits
#     
#     # add hit gene
#     m <- match(gtex.res$Peak.id, s$Enh.Pos[which(s$Hit)])
#     gtex.res$Peak.gene <- s[which(s$Hit),"Gene"][m]
#   
#     # same target?
#     gtex.res$SameTarget <- gtex.res$eQTL.gene == gtex.res$Peak.gene
#     
#   ## Explore
#     ## Rate of having eqtls
#       p <- data.frame(Region = targets,
#                       Hit = targets %in% hits,
#                       Has.eQTL = targets %in% gtex.res$Peak.id)
#       
#       pdf(file = "eQTL - GTEx - Peak Has eQTL.pdf", height = 3, width = maxw)
#       ggplot(p, aes(x = Hit, fill = Has.eQTL)) +
#         geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
#         theme_bw() +
#         scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 550)) +
#         labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
#         guides(fill = guide_legend(title = "Has GTEx eQTL")) +
#         theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
#               panel.grid = invis)
#       dev.off()
#     
#     ## Same target
#       p <- p[which(p$Hit),]
#     
#       x <- gtex.res$Peak.id[which(gtex.res$SameTarget)] # has same target
#       p$SameTarget <- p$Region %in% x
#     
#       pdf(file = "eQTL - GTEx - Peak And eQTL Share Target.pdf", height = 3, width = maxw)
#       ggplot(p, aes(x = Has.eQTL, fill = SameTarget)) +
#         geom_bar(colour = "black", width = 0.7, position = position_stack()) +
#         theme_bw() +
#         scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
#         labs(y = "Count of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
#         guides(fill = guide_legend(title = "Same Target Gene")) +
#         theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis)
#       dev.off()
# 
#     ## When there is an overlap, what is the tissue?
#       gtex.res$eQTL.tissue.type <- "."
#       g <- grep("Brain", gtex.res$eQTL.tissue)
#       gtex.res$eQTL.tissue.type[g] <- "Brain"
#       gtex.res$eQTL.tissue.type[-g] <- "Non-brain"
#       gtex.res$eQTL.tissue.type[grep("Nerve", gtex.res$eQTL.tissue)] <- "Nerve Tibial"
#       
#       pdf(file = "eQTL - GTEx - Tissue Types.pdf", height = 3, width = maxw)
#       pA <- ggplot(gtex.res, aes(x = Hit, fill = eQTL.tissue.type)) +
#         geom_bar(colour = "black", width = 0.7, position = position_stack()) +
#         theme_bw() +
#         scale_fill_manual(values = carto_pal(3, "Pastel")) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 1200)) +
#         labs(y = "Fraction of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
#         guides(fill = guide_legend(title = "Same Target Gene")) +
#         theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
#               legend.position = c(0.75, 0.8))
#       
#       pB <- ggplot(gtex.res, aes(x = Hit, fill = eQTL.tissue.type)) +
#         geom_bar(colour = "black", width = 0.7, position = position_fill()) +
#         theme_bw() +
#         scale_fill_manual(values = carto_pal(3, "Pastel")) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
#         labs(y = "Count of Peaks in Category", x = "Hit Enhancer Overlaps eQTL") +
#         guides(fill = guide_legend(title = "eQTL Tissue Type")) +
#         theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
#               legend.position = "none")
#       
#       plot_grid(pB, pA, rel_widths = c(1, 1.4))
#       dev.off()
#       
#       ## An alternative view of the above: for a given peak, are any of its eQTLs in the brain?
#       x <- unique(gtex.res$Peak.id[grep("Brain", gtex.res$eQTL.tissue)])
#       p <- data.frame(Region = targets,
#                       Hit = targets %in% hits,
#                       Has.eQTL = targets %in% gtex.res$Peak.id,
#                       Has.brain.eQTL = targets %in% x)
#       
#       p <- p[which(p$Has.eQTL),]
#       
#       pdf(file = "eQTL - GTEx - Peak Has eQTL In Brain.pdf", height = 3, width = maxw)
#       pA <- ggplot(p, aes(x = Hit, fill = Has.brain.eQTL)) +
#         geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
#         theme_bw() +
#         scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 300)) +
#         labs(y = "Fraction of Peaks in Category", x = "Enhancer With eQTL Has Is Hit") +
#         guides(fill = guide_legend(title = "At Least 1 GTEx\neQTL In Brain")) +
#         theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
#               panel.grid = invis)
#       
#       pB <- ggplot(p, aes(x = Hit, fill = Has.brain.eQTL)) +
#         geom_bar(colour = "black", width = 0.7, position = position_fill()) +
#         theme_bw() +
#         scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#         scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
#         labs(y = "Count of Peaks in Category", x = "Enhancer With eQTL Has Is Hit") +
#         guides(fill = guide_legend(title = "At Least 1 GTEx eQTL In Brain")) +
#         theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
#               legend.position = "none")
#       
#       plot_grid(pB, pA)
#       dev.off()
#       
#       ## Save
#       write.csv(gtex.res, file = "eQTL - GTEx - Data.csv")
#       write.csv(pe.res, file = "eQTL - PE - Data.csv")
  
################################################################################################################################ #
## Enhancer overlap to bulk eQTLs ----
 
## Versus GTEx
  ## Read in eQTLs
    # GTEx 
    gtex <- read.table("../../../../PublicData/GTEX_cis_eQTLs_050422_ucscDownload", sep = "\t", header = TRUE)
    gtex <- gtex[,c("eqtlChrom", "eqtlStart", "eqtlEnd", "eqtlName", "geneName", "tissue", "cpp")]
    colnames(gtex) <- c("Chr", "Start", "End", "rsID", "Gene", "Tissue", "CausalPosteriorProbability")
    
  ## Filter to those in our list of linked SNPs
    x <- gtex[which(gtex$rsID %in% snp.annotation$SNP),]
    
    x$Consistent <- x$EnhGene <- x$EnhHit <- x$LinkageCategory <- x$Enh <- "."
    
    h <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs
    
    for (j in 1:nrow(x)) {
      # get screen enhancer snp annotation of eqtls
      y <- snp.annotation[which(snp.annotation$SNP == x$rsID[j]),]
      y <- unique(y)
      
      # add the enhancer
      x$Enh[j] <- paste(y$Enh, collapse = "/")
      x$LinkageCategory[j] <- paste(y$Category, collapse = "/")
      
      # is the enhancer a hit?
      if (any(y$Enh %in% h$Enh)) {
        z <- h[which(h$Enh %in% y$Enh),]
        x$EnhHit[j] <- TRUE
        x$EnhGene[j] <- paste(z$Gene, collapse = "/")
        x$Consistent[j] <- any(z$Gene %in% x$Gene[j])
      }
    }
    
    gtex.eqtl.overlap <- data.frame(SNP = x$rsID,
                                    Enh = x$Enh,
                                    Category = x$LinkageCategory,
                                    Enh.Hit = x$EnhHit,
                                    Enh.Gene = x$EnhGene,
                                    eQTL.Gene = x$Gene,
                                    eQT.Tissue = x$Tissue,
                                    eQTL.CPP = x$CausalPosteriorProbability,
                                    Consistent = x$Consistent)
    
    
    
    write.csv(gtex.eqtl.overlap, file = "eQTL - GTEx.csv")
    

## Versus multi-ancestry brain metaanalysis from Zeng 2022
  ## Read in (fine-mapped)
    mmeqtl <- read.delim("../../../../PublicData/Zeng2022_eQTL/mmQTL_brain_meta_eqtl_finemap.tsv")
    
  ## Overlap
    x <- mmeqtl[which(mmeqtl$Variant %in% snp.annotation$SNP),]
    
    m <- match(x$Gene, geneInfo$EnsID)
    x$Gene <- geneInfo$Symbol[m]
    x <- x[-which(is.na(x$Gene)),]
    
    x$Consistent <- x$EnhGene <- x$EnhHit <- x$LinkageCategory <- x$Enh <- "."
    
    h <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs
    
    for (j in 1:nrow(x)) {
      # get screen enhancer snp annotation of eqtls
      y <- snp.annotation[which(snp.annotation$SNP == x$Variant[j]),]
      y <- unique(y)
      
      # add the enhancer
      x$Enh[j] <- paste(y$Enh, collapse = "/")
      x$LinkageCategory[j] <- paste(y$Category, collapse = "/")
      
      # is the enhancer a hit?
      if (any(y$Enh %in% h$Enh)) {
        z <- h[which(h$Enh %in% y$Enh),]
        x$EnhHit[j] <- TRUE
        x$EnhGene[j] <- paste(z$Gene, collapse = "/")
        x$Consistent[j] <- any(z$Gene %in% x$Gene[j])
      }
    }
    
    mmeqtl.overlap <- data.frame(SNP = x$Variant,
                                    Enh = x$Enh,
                                    Category = x$LinkageCategory,
                                    Enh.Hit = x$EnhHit,
                                    Enh.Gene = x$EnhGene,
                                    eQTL.Gene = x$Gene,
                                    eQL.Tissue = "Brain",
                                    eQTL.CPP = x$PP,
                                 Consistent = x$Consistent)
    
  ## Save
    write.csv(mmeqtl.overlap, file = "eQTL - mmeQTL.csv")
    
## Versus Metabrain, from de Klein 2023
  ## Read in
    
    # this resource has a similar scope to mmeQTL
    metabrain <- read_xlsx("../../../../PublicData/deKlein2023_eQTL/Supplementary Table 2 - Cis-eQTLs.xlsx", sheet = "Cortex-EUR", skip = 1)  
      
    # I note that the authors report 16619 hits a q < 0.05, so filter to this 
    metabrain$qval <- gsub("E", "e", metabrain$qval) %>% as.numeric() 
    metabrain <- metabrain[which(metabrain$qval < 0.05),]# all already are?
    
    # extract rs id
    metabrain$rsID <- splitter(metabrain$SNP, ":", 3) # I have confirmed that all have an associated rsID
  
   ## Overlap
    x$Consistent <- x$EnhGene <- x$EnhHit <- x$LinkageCategory <- x$Enh <- "."
    
    h <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs
    
    for (j in 1:nrow(x)) {
      # get screen enhancer snp annotation of eqtls
      y <- snp.annotation[which(snp.annotation$SNP == x$rsID[j]),]
      y <- unique(y)
      
      # add the enhancer
      x$Enh[j] <- paste(y$Enh, collapse = "/")
      x$LinkageCategory[j] <- paste(y$Category, collapse = "/")
      
      # is the enhancer a hit?
      if (any(y$Enh %in% h$Enh)) {
        z <- h[which(h$Enh %in% y$Enh),]
        x$EnhHit[j] <- TRUE
        x$EnhGene[j] <- paste(z$Gene, collapse = "/")
        x$Consistent[j] <- any(z$Gene %in% x$GeneSymbol[j])
      }
    }
    
    metabrain.overlap <- data.frame(SNP = x$rsID,
                                    Enh = x$Enh,
                                    Category = x$LinkageCategory,
                                    Enh.Hit = x$EnhHit,
                                    Enh.Gene = x$EnhGene,
                                    eQTL.Gene = x$GeneSymbol,
                                    eQTL.Tissue = "Brain",
                                    Consistent = x$Consistent)
    
  ## Save
    write.csv(metabrain.overlap, file = "eQTL - Metabrain.csv")
  
## Annotate snp.annotation object
  snp.annotation$eQTL.GTEx <- snp.annotation$SNP %in% gtex.eqtl.overlap$SNP[which(gtex.eqtl.overlap$Consistent == "TRUE")]
  snp.annotation$eQTL.mmeQTL <- snp.annotation$SNP %in% mmeqtl.overlap$SNP[which(mmeqtl.overlap$Consistent == "TRUE" & mmeqtl.overlap$eQTL.CPP > 0.1)]

    
  
################################################################################################################################ #
## Enhancer overlap to snRNA-seq eQTLs ----
 
# ## Paths
#   bryois.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Bryois2021_eQTL_snRNAseq/"
#   # bryois.bed <- paste0(bryois.dir, "Processed.bed")
#   # nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
#   # bryois.out <- paste0(getwd(), "/eQTL_snAst_Overlap.bed")
#    
# ## Process eQTLs
#   # read in
#   l <- list.files(bryois.dir)
#   l <- l[-c(23:24)] # that's metadata
#   
#   bryois.raw <- lapply(l, function(x) {
#     print(x)
#     y <- read.table(paste0(bryois.dir, x), sep = " ")
#     colnames(y) <- c("Gene", "SNP", "Distance", "P", "Beta")
#     y$Gene <- splitter(y$Gene, "_", 1)
#     y <- y[,c("Gene", "SNP", "P", "Beta")]
#     return(y)
#   })
#   
#   # names(bryois.raw) <- gsub("Astrocytes.quantile.txt.gz.", "", l) %>% splitter(":", 1) %>% paste0("Chr", .)
#   
#   bryois.raw <- do.call("rbind", bryois.raw)
# 
#   # filter to significant hits
#   bryois.raw$FDR <- p.adjust(bryois.raw$P, method = "fdr")
#   keep <- which(bryois.raw$FDR < 0.05)
#   bryois.raw <- bryois.raw[keep,]
#   
#   # add hg38 position
#   posit <- read.delim(paste0(bryois.dir, "snp_pos.txt"))
#   m <- match(bryois.raw$SNP, posit$SNP)
#   bryois.raw$Chr <- posit$chr[m]
#   bryois.raw$Start <- posit$pos_hg38[m]
#   bryois.raw$End <- posit$pos_hg38[m]
# 
#   # write to disk as bed
#   bryois.raw <- bryois.raw[,c("Chr", "Start", "End", "SNP", "Gene", "FDR")]
#   write.table(bryois.raw, bryois.bed, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#   
#   
# ## Overlap
#   call <- paste("windowBed",
#               "-a", nha_hg38,
#               "-b", bryois.bed,
#               "-w", "1000", # window of 1000bp
#               ">", bryois.out)
# 
#    system(call, intern = FALSE, wait = TRUE)
#   
#   
# ## Read in and annotate 
#   load
#   bryois.res <- read.table(bryois.out, sep = "\t", header = FALSE)
#   
#   # column names
#   colnames(bryois.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "eQTL.chr", "eQTL.start", "eQTL.end", "eQTL.id", "eQTL.gene", "eQTL.fdr")
# 
#   # add hit annotation
#   bryois.res$Peak.id <- sub("_", ":", bryois.res$Peak.id)
#   bryois.res$Peak.id <- sub("_", "-", bryois.res$Peak.id)
#   bryois.res$Hit <- bryois.res$Peak.id %in% hits
# 
#   # add hit gene
#   m <- match(bryois.res$Peak.id, s$Enh.Pos[which(s$Hit)])
#   bryois.res$Peak.gene <- s[which(s$Hit),"Gene"][m]
# 
#   # same target?
#   bryois.res$SameTarget <- bryois.res$eQTL.gene == bryois.res$Peak.gene
# 
# ## Little to report, so just save csv
#   write.csv(bryois.res, file = "eQTL - snAst - Data.csv")
    
  
## New overlapping method
  bryois <- readxl::read_xlsx("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/Bryois2021_eQTL_snRNAseq/41593_2022_1128_MOESM3_ESM.xlsx", sheet = "Table S3", skip = 3) # fine mapped eQTLs
  # bryois <- bryois[which(bryois$cell_type == "Astrocytes"),]
  
  ## Overlap
    x <- bryois[which(bryois$SNP %in% snp.annotation$SNP),] %>% as.data.frame()
    
    x$EnhGene <- x$EnhHit <- x$LinkageCategory <- x$Enh <- "."
    
    h <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs
    
    for (j in 1:nrow(x)) {
      # get screen enhancer snp annotation of eqtls
      y <- snp.annotation[which(snp.annotation$SNP == x$SNP[j]),]
      y <- unique(y)
      
      # add the enhancer
      x$Enh[j] <- paste(y$Enh, collapse = "/")
      x$LinkageCategory[j] <- paste(y$Category, collapse = "/")
      
      # is the enhancer a hit?
      if (any(y$Enh %in% h$Enh)) {
        z <- h[which(h$Enh %in% y$Enh),]
        x$EnhHit[j] <- TRUE
        x$EnhGene[j] <- paste(z$Gene, collapse = "/")
      }
    }
    
    bryois.overlap <- data.frame(SNP = x$SNP,
                                    Enh = x$Enh,
                                    Category = x$LinkageCategory,
                                    Enh.Hit = x$EnhHit,
                                    Enh.Gene = x$EnhGene,
                                    eQTL.Gene = x$symbol,
                                    eQTL.Tissue = x$cell_type,
                                    eQTL.CPP = x$Probability)
    
  ## Save
    write.csv(bryois.overlap, file = "eQTL - Bryois snRNAseq.csv")
  
        
################################################################################################################################ #
## Enhancer overlap to clinical variants from DisGeNet ----

    
## Old version below  
  
# ## Paths
#   disg.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/DisGeNet_070422/all_variant_disease_associations.bed"
#   disg.out <- paste0(getwd(), "/Disgenet_Variant_Overlap.bed")
#   
# ## Process table
#   disg <- read.delim("../../../PublicData/DisGeNet_070422/all_variant_disease_associations.tsv", sep = "\t", header = TRUE)
#   disg.bed <- disg[,c("chromosome", "position", "position", "snpId", "diseaseName","diseaseType")]
#   disg.bed$chromosome <- paste0("chr", disg.bed$chromosome)
#   write.table(disg.bed, file = disg.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
#   
# ## Run
#   call <- paste("windowBed",
#                 "-a", nha_hg38,
#                 "-b", disg.dir,
#                 "-w", "1000", # window of 1000bp
#                 ">", disg.out)
#   
#   system(call, intern = FALSE, wait = TRUE) 
#   
# ## Read
#   disg.res <- read.delim(disg.out, sep = "\t", header = FALSE)
#   colnames(disg.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "Disg.chr", "Disg.start", "Disg.end", "Disg.id", "Disg.disease", "Disg.class")
#   
# ## Output variants for hit peaks
#   disg.res$Peak.id <- sub("_", ":", disg.res$Peak.id)
#   disg.res$Peak.id <- sub("_", "-", disg.res$Peak.id)
#   disg.res$Hit <- disg.res$Peak.id %in% enh.hits
#   disg.hits <- disg.res[which(disg.res$Hit),]
#   
# ## Variant enrichment within peaks
#   p <- data.frame(Region = targets,
#                   Hit = targets %in% hits,
#                   Has.disg = targets %in% disg.res$Peak.id)
#   
#   pdf(file = "Disgenet - Peak Has Variant.pdf", height = 3, width = maxw)
#   pA <- ggplot(p, aes(x = Hit, fill = Has.disg)) +
#     geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
#     theme_bw() +
#     scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#     scale_y_continuous(expand = c(0,0), limits = c(0, 820)) +
#     labs(y = "Fraction of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
#     guides(fill = guide_legend(title = "Has Disgenet Variant")) +
#     theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
#           panel.grid = invis)
#   
#   pB <- ggplot(p, aes(x = Hit, fill = Has.disg)) +
#     geom_bar(colour = "black", width = 0.7, position = position_fill()) +
#     theme_bw() +
#     scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
#     scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
#     labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
#     guides(fill = guide_legend(title = "Has Disgenet Variant")) +
#     theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
#           panel.grid = invis)
#   
#   plot_grid(pB, pA)
#   dev.off()
#   
# ## Output
#   write.csv(disg.res, file = "Disgenet - Peak Variant Overlap.csv")
#   write.csv(disg.hits, file = "Disgenet - Peak Variant Overlap (Hits Only).csv")
  
## Read in
  disg <- read.delim("../../../../PublicData/DisGeNet_070422/all_variant_disease_associations.tsv", sep = "\t", header = TRUE)
    
## Overlap to our enhancers' SNPs
  x <- disg[which(disg$snpId %in% snp.annotation$SNP),]
  
  x$EnhGene <- x$EnhHit <- x$LinkageCategory <- x$Enh <- "."
    
  h <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs
    
  for (j in 1:nrow(x)) {
    # get screen enhancer snp annotation of eqtls
    y <- snp.annotation[which(snp.annotation$SNP == x$snpId[j]),]
    y <- unique(y)
    
    # add the enhancer
    x$Enh[j] <- paste(y$Enh, collapse = "/")
    x$LinkageCategory[j] <- paste(y$Category, collapse = "/")
    
    # is the enhancer a hit?
    if (any(y$Enh %in% h$Enh)) {
      z <- h[which(h$Enh %in% y$Enh),]
      x$EnhHit[j] <- TRUE
      x$EnhGene[j] <- paste(z$Gene, collapse = "/")
    }
  }
  
## Extract salient information
  disg.overlap <- data.frame(SNP = x$snpId,
                             x[,17:20],
                             Disease = x$diseaseName,
                             DiseaseType = x$diseaseSemanticType,
                             DiseaseScore = x$score,
                             DSI = x$DSI,
                             DPI = x$DPI)
  
  write.csv(disg.overlap, "Disgenet.csv")
  
## Add to annotation dataframe
  snp.annotation$Disgenet <- snp.annotation$SNP %in% disg.overlap$SNP
  
################################################################################################################################ #
## Enhancer overlap to GWAS ----
  
  
## GWAS Catalogue
  
    # ## Paths
    #   gwas.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/gwasCatalog_050422_ucscDownload.bed"
    #   gwas.out <- paste0(getwd(), "/GWAS_Catalogue_Variant_Overlap.bed")
    # 
    # ## Process table
    #   gwas <- read.delim("../../../PublicData/gwasCatalog_050422_ucscDownload", sep = "\t", header = TRUE)
    #   gwas.bed <- gwas[,c("chrom", "chromStart", "chromEnd", "name", "trait","genes")]
    #   gwas.bed$genes[gwas.bed$genes == ""] <- "(Blank)"
    #   write.table(gwas.bed, file = gwas.dir, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
    #   
    # ## Run
    #   call <- paste("windowBed",
    #                 "-a", nha_hg38,
    #                 "-b", gwas.dir,
    #                 "-w", "1000", # window of 1000bp
    #                 ">", gwas.out)
    #   
    #   system(call, intern = FALSE, wait = TRUE) 
    #   
    # ## Read
    #   gwas.res <- read.delim(gwas.out, sep = "\t", header = FALSE)
    #   colnames(gwas.res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.id", "GWAS.chr", "GWAS.start", "GWAS.end", "GWAS.id", "GWAS.disease", "GWAS.gene")
    #   
    # ## Output variants for hit peaks
    #   gwas.res$Peak.id <- sub("_", ":", gwas.res$Peak.id)
    #   gwas.res$Peak.id <- sub("_", "-", gwas.res$Peak.id)
    #   gwas.res$Hit <- gwas.res$Peak.id %in% hits
    #   gwas.hits <- gwas.res[which(gwas.res$Hit),]
    #   
    # ## Variant enrichment within peaks
    #   p <- data.frame(Region = targets,
    #                   Hit = targets %in% hits,
    #                   Has.gwas = targets %in% gwas.res$Peak.id)
    #   
    #   pdf(file = "GWAS Catalogue - Peak Has Variant.pdf", height = 3, width = maxw)
    #   pA <- ggplot(p, aes(x = Hit, fill = Has.gwas)) +
    #     geom_bar(colour = "black", width = 0.7, position = position_dodge()) +
    #     theme_bw() +
    #     scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
    #     scale_y_continuous(expand = c(0,0), limits = c(0, 820)) +
    #     labs(y = "Fraction of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
    #     guides(fill = guide_legend(title = "Has Disgenet Variant")) +
    #     theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.7, 0.8),
    #           panel.grid = invis)
    #   
    #   pB <- ggplot(p, aes(x = Hit, fill = Has.gwas)) +
    #     geom_bar(colour = "black", width = 0.7, position = position_fill()) +
    #     theme_bw() +
    #     scale_fill_manual(values = rev(carto_pal(2, "Geyser"))) +
    #     scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
    #     labs(y = "Count of Peaks in Category", x = "Enhancer Has Significantly Associated Gene") +
    #     guides(fill = guide_legend(title = "Has GWAS Variant")) +
    #     theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
    #           panel.grid = invis)
    #   
    #   plot_grid(pB, pA)
    #   dev.off()
    #   
    # ## Output for hit peaks
    #   write.csv(gwas.res, file = "GWAS Catalogue - Peak Variant Overlap.csv")
    #   write.csv(gwas.hits, file = "GWAS Catalogue - Peak Variant Overlap (Hits Only).csv")
  
  ## Redo the above
    gwas.cat <- read.delim("../../../../PublicData/gwasCatalog_050422_ucscDownload", sep = "\t", header = TRUE)  
  
    x <- gwas.cat[which(gwas.cat$name %in% snp.annotation$SNP),]
    
    x$EnhGene <- x$EnhHit <- x$LinkageCategory <- x$Enh <- "."
      
    h <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs
      
    for (j in 1:nrow(x)) {
      # get screen enhancer snp annotation of eqtls
      y <- snp.annotation[which(snp.annotation$SNP == x$name[j]),]
      y <- unique(y)
      
      # add the enhancer
      x$Enh[j] <- paste(y$Enh, collapse = "/")
      x$LinkageCategory[j] <- paste(y$Category, collapse = "/")
      
      # is the enhancer a hit?
      if (any(y$Enh %in% h$Enh)) {
        z <- h[which(h$Enh %in% y$Enh),]
        x$EnhHit[j] <- TRUE
        x$EnhGene[j] <- paste(z$Gene, collapse = "/")
      }
    }
    
    gwas.cat.overlap <- data.frame(SNP = x$name,
                             Trait = x$trait,
                             x[,(ncol(x)-4):ncol(x)])
    
    write.csv(gwas.cat.overlap, "GWAS Catalogue.csv")
    
  ## Add to annotation dataframe
    snp.annotation$GWASCat <- snp.annotation$SNP %in% gwas.cat.overlap$SNP
  
# ## MAGMA
#   ## Load functions
#     source("/mnt/Data0/PROJECTS/GWAS_Enrichment/GAVIN/Scripts/0_Setup.R")
#     
#   ## Write enh peaks as a MAGMA loc file
#     # load
#     magma.peaks <- read.table(nha_hg38, sep = "\t")
#     magma.peaks$V4 <- sub("_", ":", magma.peaks$V4)
#     magma.peaks$V4 <- sub("_", "-", magma.peaks$V4)
#     magma.peaks$V5 <- magma.peaks$V4 %in% hits
#     magma.peaks$V6 <- "."
#     
#     # sort in terminal
#     unsorted.magma <- "MAGMA/MAGMA_peaks_unsorted.bed"
#     sorted.magma <- "MAGMA/MAGMA_peaks_sorted.bed"
#     write.table(magma.peaks, file = unsorted.magma, 
#                 sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#     
#     call <- paste("sort",
#                   "-k1,1",
#                   "-k2,2n",
#                   unsorted.magma,
#                   ">",
#                   sorted.magma, 
#                   sep = " ")
#     
#     system(call, intern = FALSE, wait = TRUE)  
#     
#     
#     # read in and convert to a file MAGMA can parse
#     magma.peaks.sorted <- read.table(sorted.magma, sep = "\t")
#     magma.bg <- data.frame(Peak = magma.peaks$V4,
#                            Chromosome = gsub("chr", "", magma.peaks$V1),
#                            Start = magma.peaks$V2,
#                            End = magma.peaks$V3,
#                            Strand = ".",
#                            Hit = magma.peaks$V5)
#     
#     write.table(magma.bg, file = "MAGMA/MAGMA_peaks_bg.loc",
#                 quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
#     
#   ## Annotate peaks with annotation with SNPs
#     magma.annotate(bg = "MAGMA/MAGMA_peaks_bg.loc",
#                    output = "MAGMA/MAGMA_bg_annotation")
#   
#   ## GWAS background risk
#     # function  
#     run.magma.background <- function(annot, gwas, gwas.columns = c("SNP", "P"), N, output, snps = paste0(paths$snps, "g1000_eur")) {
#       
#       # convert arguments into a terminal call
#       pval.statement <- paste0(gwas,
#                                " ", "N=", N,
#                                " ", "use=", gwas.columns[1], ",", gwas.columns[2])
#       
#       # output <- paste0("Background_Analysis/", output)
#       # annot <- paste0("/mnt/Data0/PROJECTS/Lister/Results/Annotation/", annot)
#       
#       call <- paste(paths$magma,
#                     "--bfile", snps,
#                     "--gene-annot", annot,
#                     "--pval", pval.statement,
#                     "--out", output,
#                     sep = " ")
#       
#       # run
#       system(call, intern = FALSE, wait = TRUE)
#       return("Done")
#     }
#     
#     # run
#     for (j in names(n)) {
#       run.magma.background(annot = "MAGMA/MAGMA_bg_annotation.genes.annot", gwas = paths[[j]], gwas.columns = cols[[j]], N = n[[j]], output = paste0("MAGMA/GWAS_Backgrounds/", j, "_bg"))
#     }
#     
#     
#   ## Query for increased risk in hit set
#     # write query set
#     x <- paste0("Hit ", magma.bg$Peak[which(magma.bg$Hit)])
#     write.table(x, file = "MAGMA/Hit_Enrichment/SetInput.sets", col.names = FALSE, row.names = FALSE, quote = FALSE)
#     
#     # function
#     run.magma.queryset <- function(gwas, queryset = "MAGMA/Hit_Enrichment/SetInput.sets") {
#       
#       # convert arguments into a terminal call
#       # queryset <- paste0("/mnt/Data0/PROJECTS/Lister/Data/MAGMA_QuerySets/", queryset)
#       output <- paste0("MAGMA/Hit_Enrichment/", gwas)
#       gwas <- paste0("MAGMA/GWAS_Backgrounds/", gwas, "_bg.genes.raw")
#       
#       
#       call <- paste(paths$magma,
#                     "--gene-results", gwas,
#                     "--set-annot", queryset,  "col=2,1", # the second statement tells magma that the queryset is in a column format
#                     "--out", output,
#                     sep = " ")
#       
#       # run
#       system(call, intern = FALSE, wait = TRUE)
#       return("Done")
#     }
#     
#     # run
#     for (j in names(n)) { # for GWAS j
#       run.magma.queryset(gwas = j)
#     }
    
    
## Top SNPs from various GWAS
  ## Read in
    risk.variants <- list()
    
    
    ## Summary statistics
    source("/Volumes/shar/mnt/Data0/PROJECTS/GWAS_Enrichment/GAVIN/Scripts/0_Setup.R")  
    cols$INS[2] <- "P"
    cols$SZ[2] <- "PVAL"
    snp.thresh <- 5e-8
    
    
    for (ss in names(n)[1:10]) {
        print(ss)
        corrected.path <- gsub("~/PROGRAMS/magma_v1.09b_static/GWAS_SS/", "~/REFERENCE/GWAS/", paths[[ss]])
        if (corrected.path == "~/REFERENCE/GWAS/Insomnia_sumstats_Jansenetal.txt") corrected.path <- gsub("GWAS", "GWAS/Unmodified", corrected.path)
        if (corrected.path == "~/REFERENCE/GWAS/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix") corrected.path <- "~/REFERENCE/GWAS/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv"
        
        
        if (ss %in% c("MDD")) {
          x <- read.delim(corrected.path, sep = " ")
        } else {
          x <- read.delim(corrected.path)
        }
        
        x <- x[which(x[,cols[[ss]][2]] < snp.thresh),]
        
        risk.variants[[ss]] <- x
    }
      

    risk.variants <- lapply(risk.variants, function(x) {
      g <- grep("^rs", x[1,])
      x <- as.character(x[,g])
      return(x)
    })
      
    
    
    risk.variants.filt <- lapply(risk.variants, function(x) {
      snp.annotation[which(snp.annotation$SNP %in% x),1:11]
    })
    
    save(risk.variants, risk.variants.filt, file = "Summary Statistic Sig SNPs.rda")
    
    risk.variants.filt.sum <- sapply(risk.variants.filt, function(x) table(x$Category)) %>% melt() %>% dcast(., L1~Var1)
    colnames(risk.variants.filt.sum) <- c("GWAS", "LD_ToNearest", "LD_ToOverlapping", "Nearby", "Overlapping")
    write.csv(risk.variants.filt.sum, file = "Summary Statistic Sig SNPs Overlapping Count.csv", row.names = FALSE)
      
  ## Append to SNP annotation object
    for (j in names(risk.variants.filt)) {
      snp.annotation[,j] <- snp.annotation$SNP %in% risk.variants[[j]]
    }
    
    ## Glioma is done separately as it is an xlsx of signicant hits
      # read in
      x <- read_xlsx("../../../../../PublicData/GWAS_TopHits/Glioma_Melin2017.xlsx", sheet = "1")
      
      # set colnames
      colnames(x) <- x[2,]
      x <- x[,1:16]
      for (j in 2:4) {
        i <- (((4*j)-3):(4*j))
        colnames(x)[i] <- paste0(x[1,i[1]], "_", colnames(x)[i]) 
      }
      colnames(x) <- gsub(" ", "_", colnames(x))
      colnames(x) <- gsub("-", "", colnames(x))
      
      # filter columns to this study only
      x <- x[-c(1:3, 29),]
      
      # significance?
      x$Sig_AllGlioma <- gsub("x10", "e", x$All_glioma_P) %>% as.numeric(.) < 5e-8
      x$Sig_GBM <- gsub("x10", "e", x$GBM_glioma_P) %>% as.numeric(.) < 5e-8
      x$Sig_NonGBM <- gsub("x10", "e", x$NonGBM_glioma_P) %>% as.numeric(.) < 5e-8
      
      # some are false for all of these, but have been found in previous studies
      risk.variants$Glioma_2020 <- x
      
## Overlap
  # table(risk.variants$Glioma_2020$SNP %in% snp.annotation$SNP)
  
################################################################################################################################ #
## A combined variant annotation file! ---- 

    # ## Here, you aim to create a table which summarises the variants for each enhancer
    #   
    # ## Create
    #   out <- data.frame(Region = targets,
    #                     Hit = targets %in% hits)
    #   
    # ## Function
    #   add.annot <- function(x = out,  overlaps, name, combine = FALSE) {
    #     # is there an overlap
    #     x[,paste0(name, "_Overlap")] <- x$Region %in% overlaps[,4]
    #     
    #     # if overlap, what is its nature?
    #     x[,paste0(name, "_Info")] <- ""
    #     m <- match(overlaps[,4], x$Region)
    #     
    #     if (combine) {
    #       for (j in 1:length(m)) {
    #         x[m[j],paste0(name, "_Info")] <- paste0(x[m[j],paste0(name, "_Info")], " | ", paste0(overlaps[j,9]), " in ", overlaps[j,10])
    #       }
    #     } else {
    #       for (j in 1:length(m)) {
    #         x[m[j],paste0(name, "_Info")] <- paste0(x[m[j],paste0(name, "_Info")], " | ", overlaps[j,9])
    #       }
    #     }
    #     
    #     
    #     
    #     return(x)
    #   }
    #   
    # out <- add.annot(out, overlaps = gwas.res, name = "gwas")
    # out <- add.annot(out, overlaps = disg.res, name = "disgenet")
    # out <- add.annot(out, overlaps = pe.res, name = "eQTLpe")
    # out <- add.annot(out, overlaps = gtex.res, name = "eQTLgtex", combine = TRUE)
    # out <- add.annot(out, overlaps = bryois.res, name = "eQTLbryois")
    # 
    # write.csv(out, "Variant Annotation To Peaks.csv")
  
## Calculate the number of biological signals for each SNP
snp.annotation$nBiologicalSignals <- rowSums(snp.annotation[,which(colnames(snp.annotation) == "SNP2TFBS") :which(colnames(snp.annotation) == "INS") ])
      
## Remove duplicated rows
snp.annotation <- snp.annotation[-which(duplicated(paste0(snp.annotation$SNP, snp.annotation$Pos, snp.annotation$Enh))),]      
# y <- snp.annotation[-which(duplicated(snp.annotation)),]      

write.csv(snp.annotation, file = "SNP Annotation With Biological Signals.csv", row.names = FALSE)
write.csv(snp.annotation[which(snp.annotation$HitPermissive & as.logical(snp.annotation$WithinEnhSummit)),], file = "SNP Annotation With Biological Signals (Filtered To Summits Of Interest).csv", row.names = FALSE)
      
      
      
      