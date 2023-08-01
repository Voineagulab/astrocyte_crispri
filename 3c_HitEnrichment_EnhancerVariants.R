## This script explores how variants associate with our peaks

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Variants")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  # no additional libraries versus the function script
  

## Load
  source("../../../Scripts/Functions.R")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]

  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  sig.colours <- c("black", "firebrick1")

## Results dataframe
  res.final <- read.csv("../../2_DE/Enh/Results Final.csv")
  hit.enh <- res.final$Enh[which(res.final$HitPermissive)]
  hit.pairs <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs

  
## Main paths
  nha_dir_38 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  nha_dir_19 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed" 
  

################################################################################################################################ #
## First: annotate enhancers with SNPs from dbSNP ----
  
## Paths
  db_dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/dbSNP/dbSNP153_280422.bed"
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  db_out <- "dbSNP_Overlap.bed"
  
## Intersect
  call <- paste("windowBed",
                "-a", nha_dir_38,
                "-b", db_dir,
                "-w", "1000", # window of 1000bp
                ">", db_out)
  
    system(call, intern = FALSE, wait = TRUE) 
  
## Analyse
  db_res <- read.delim(db_out, header = FALSE)
    
  # column names
  colnames(db_res) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.coord", "SNP.chr", "SNP.start", "SNP.end", "SNP.id")
    
  # add hit annotation
  db_res$Peak.coord <- sub("_", ":", db_res$Peak.coord) %>% sub("_", "-", .)
  
## Collect variants in LD to these variants
  ## Calculate LD
    # requires ~24h
  
    run.LDcode <- FALSE
    
    if (run.LDcode) {
      LDproxy_batch(snp = unique(db_res$SNP.id), # unique is necessary as these peaks have a window appended
                    pop = "EUR", # european
                    r2d = "r2",
                    append = TRUE,
                    genome_build = "grch38",
                    token = "5988cfbf2737")  
    }
    
  
  ## Read in
    all.ld.snps <- read.table("combined_query_snp_list_grch38.txt", sep = "\t", header = TRUE, row.names = NULL)
    all.ld.snps <- all.ld.snps[,c("query_snp", "RS_Number", "Coord", "MAF", "Distance", "R2")]
    colnames(all.ld.snps)[1:3] <- c("Query_SNP", "Linked_SNP", "Position")
  
  ## Filter
    ld.snps <- all.ld.snps[which(all.ld.snps$R2 > 0.8 & all.ld.snps$Distance != 0),] # subset to r2 > 0.8, and not same-snp
  

## Combine peak annotations to snps
  snp2peak <- list()
  tested.peaks <- res.final$Enh.Pos %>% unique()

  for (j in tested.peaks) { # for each tested peak (957), rather than all peaks (979)
    
    print(j)
    
    ## First, look at nearby (within 1kb) SNPs
      
      # skip if no snps within 1kb window
      if ( ! ( j %in% db_res$Peak.coord)) {
       
         snp2peak[[j]] <- ("None") 
        
        next
         
      }
    
      # collect these snps in the dataframe x
      x <- db_res[which(db_res$Peak.coord == j),]
      peak.start <- x$Peak.start[1]; peak.end <- x$Peak.end[1]
      peak.centre <- mean(c(peak.start, peak.end)) %>% round()
      x <- x[,c("Peak.coord", "SNP.id", "SNP.start")]
      colnames(x) <- c("Peak", "SNP", "SNP.Pos")
      
      x$R2 <- NA
      x$Category <- "."
    
      # annotate as overlapping or nearby
      w <- ((x$SNP.Pos > peak.start) & (x$SNP.Pos < peak.end)) # this is technically based on SNP.start
      if (any(w)) x$Category[w] <- "Overlapping"
      if (any(!(w))) x$Category[!(w)] <- "Nearby"
      
      # annotate with distance rank
      x$DistanceFromCentre <- (x$SNP.Pos - peak.centre)
      x$DistanceRank <- rank(abs(x$DistanceFromCentre))
      
      
    ## Collect SNPs in ld 
      
      if (any(x$SNP %in% ld.snps$Query_SNP)) { 
        
        # get
        y <- ld.snps[which(ld.snps$Query_SNP %in% x$SNP),] # collect
        
        # annotate by category-of-linked-SNP
        y$Linked_SNP_Category <- "Nearby"
        y$Linked_SNP_Category[which(y$Query_SNP %in% x$SNP[w])] <- "Overlapping"
        
        # reformat
        z <- data.frame(Peak = j,
                        SNP = y$Linked_SNP, 
                        SNP.Pos = splitter(y$Position, ":", 2) %>% as.numeric(),
                        R2 = y$R2,
                        Category = paste0("LD_to_", y$Linked_SNP_Category),
                        DistanceFromCentre = NaN,
                        DistanceRank = NA)
        
        z$DistanceFromCentre <- (z$SNP.Pos - peak.centre)
        
        # remove LD SNPs which are also nearby! match based on pos rather than ID, as sometimes ID is "."
        if (any(z$SNP.Pos %in% x$SNP.Pos)) {
          
          dup <- which(z$SNP.Pos %in% x$SNP.Pos)
          
          if (length(dup) == length(w)) { # if removing as many rows as there are in the dataframe, may cause error. instead, finish loop
            snp2peak[[j]] <- x
            next
          }
          
          z <- z[-which(z$Peak.coord %in% x$SNP.Pos),]  
          
        }
        
  
        
        # combine
        x <- rbind(x, z)
        
      } 
      
  
    ## Return
      snp2peak[[j]] <- x
    
  }
    
  
## Process
  snp2peak <- do.call("rbind", snp2peak)
    
  # remove the "none", i.e. when there is a peak with no nearby SNPs
  snp2peak <- snp2peak[-which(snp2peak$Category == "None"),]
  
  # add the enhancer id
  m <- match(snp2peak$Peak, guides$TargetCoord)
  snp2peak$Enh <- guides$TargetID[m] 
  
  # is it a hit enhancer?
  snp2peak$Hit <- snp2peak$Enh %in% hit.enh
  
  ## Add snp synonyms 
    snp2peak$Synonym <- NA
  
    # read from: # http://genehopper.ifis.cs.tu-bs.de/downloads, Ensembl VariationDB 84
    syn <- read.table("../../../../PublicData/dbSNP/variation_synonyms.txt", sep = "\t", header = TRUE) 
    syn <- syn[-grep("N", syn$synonyms),] # removes all those without a synonym
    
    # match your id to main id, extract the archive
    m <- match(snp2peak$SNP, syn$main)
    snp2peak$Synonym <- syn$synonyms[m]
    
    # match your id to archive id, extract the main
    syn.split <- strsplit(syn$synonyms, ",")
    names(syn.split) <- paste0(syn$main, "_") # the paste helps with matching downstream
    syn.split <- do.call("c", syn.split)
    syn.split <- data.frame(Main = splitter(names(syn.split), "_", 1),
                    Archive = syn.split)
    
    m <- match(snp2peak$SNP, syn.split$Archive)
    snp2peak$Synonym <- paste0(snp2peak$Synonym, ",", syn.split$Main[m], ",")
    
    # clean
    snp2peak$Synonym <- gsub("NA,", "", snp2peak$Synonym) %>%
      # gsub(",NA", "|", .) %>% # commas don't suit csv
      gsub(",", "|", .) # separator
    
    
## Final clean
  snp2peak <- relocate(snp2peak, "Enh", "Hit", "SNP", "Synonym")
  rownames(snp2peak) <- 1:nrow(snp2peak)
    
  
## Save
  write.csv(snp2peak, file = "Final - SNP-to-Peak List.csv", row.names = FALSE)

   
  
################################################################################################################################ #
## Setup annotation process ----
   
## Now start annotating
  snp.annot <- read.csv("Final - SNP-to-Peak List.csv")
  # to skip ahead: snp.annot <- read.csv("Final - SNP Annotation.csv", row.names = 1)
  
  
## List of all SNPs for easier matching
  x <- gsub("^\\|", "", snp.annot$Synonym) %>%
    strsplit("|", fixed = TRUE) %>%
    do.call("c", .) %>% 
    unique()

  all.snps <- c(x, snp.annot$SNP) %>% unique()  
  
## Function to annotate our SNPs by a given SNP list
  check.snps <- function(snps, colname, data = snp.annot) {
    
    # match to main
    in_main <- data$SNP %in% snps
    
    # match to synonyms
    in_syn <- list()
    for (j in snps) { in_syn[[j]] <- grepl(paste0(j, "|"), data$Synonym, fixed = TRUE) }
    in_syn <- do.call("cbind", in_syn)
    in_syn <- rowSums(in_syn) > 0
    
    # combine
    data[,colname] <- in_main | in_syn
    
    # return
    return(data)
  
  }
  

## Function to annotate eQTL (and other similarly-formatted data) with more detailed information
  
  # requires columns for "Tissue", "Gene", and "rsID". a user-specified score column is optional.
  detailed.overlap <- function(snps, score.col = NA, score.makecolname = "SNP.CPP", check.consistency = TRUE, gwas.format = FALSE) {
    
    # filter to overlapping snps
    x <- snps[which(snps$rsID %in% all.snps),] %>% as.data.frame()
    
    # add columns
    x$Consistent <- x$EnhGene <- x$EnhHit <- x$LinkageCategory <- x$Enh <- "."
    
    # annotate
    for (j in 1:nrow(x)) {
      
      # match snp
      w <- (snp.annot$SNP == x$rsID[j]) | grepl(paste0(x$rsID[j], "|"), snp.annot$Synonym, fixed = TRUE)
      
      # get screen enhancer snp annotation of eqtls
      y <- snp.annot[which(w),]
      y <- unique(y)
      
      # add the enhancer
      x$Enh[j] <- paste(y$Enh, collapse = "/")
      x$LinkageCategory[j] <- paste(y$Category, collapse = "/")
      
      # is the enhancer a hit?
      if (any(y$Enh %in% hit.pairs$Enh)) {
        z <- hit.pairs[which(hit.pairs$Enh %in% y$Enh),]
        x$EnhHit[j] <- TRUE
        x$EnhGene[j] <- paste(z$Gene, collapse = "/")
        x$Consistent[j] <- any(z$Gene %in% x$Gene[j])
      }
    }
    
    # output
    if (gwas.format) { 
      return(x) 
    } else {
      overlap <- data.frame(SNP = x$rsID,
                            LinkedEnh = x$Enh,
                            Category = x$LinkageCategory,
                            SNP.Association = x$Gene,
                            SNP.Tissue = x$Tissue,
                            LinkedEnh.Hit = x$EnhHit,
                            LinkedEnh.Gene = x$EnhGene)
      
      if (!(anyNA(score.col))) {
        overlap[,score.makecolname] <- x[,score.col]
        overlap <- relocate(overlap, "SNP", "LinkedEnh", "Category", "SNP.Association", "SNP.Tissue", all_of(score.makecolname))
      }
      if (check.consistency) overlap$Consistent <- x$Consistent
      
      return(overlap)
    }
    
  }
  
 
################################################################################################################################ #
## SNPs that alter TF binding sites ----
  
## Use the SNP2TFBS webtool:
  # https://ccg.epfl.ch/snp2tfbs/snpselect.php
  # https://academic.oup.com/nar/article/45/D1/D139/2605727?login=true
  
## Prepare for upload
  x <- snp.annot[which(snp.annot$Category == "Overlapping"),]
  y <- gsub("^\\|", "", x$Synonym) %>%
    strsplit("|", fixed = TRUE) %>%
    do.call("c", .) %>% 
    unique()
  z <- c(x$SNP, y) %>% unique
  write.table(z, file = "SNP2TFBS Webtool Input.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
## Download and wrangle
  snp2tfbs <- read.delim("SNP2TFBS Webtool Output - Variant Annotation.txt", header = FALSE)
  colnames(snp2tfbs) <- c("hg19chr", "hg19pos", "Ref", "Alt", "Annotation", "Genes", "Score", "Unused", "rsID", "TF")
  snp2tfbs <- snp2tfbs[,c("rsID", "TF", "Score", "Ref", "Alt")]
  
  # map to enh
  m <- match(snp2tfbs$rsID, snp.annot$SNP) # can use match as it's overlapping SNPs and thus unique, noting that all map to $SNP
  snp2tfbs$Enh <- snp.annot$Enh[m]
  snp2tfbs$Hit <- snp.annot$Hit[m]
  snp2tfbs$Category <- snp.annot$Category[m]

  
## Categorise by overlap to TOBIAS motifs
  ## Read in TOBIAS data
    tob <- read.delim("../../../../EnhancerPredictionModels/Results/Tobias/Footprint/BINDetect/bound_overlaps.bed", header = FALSE)
    tob$TF <- splitter(tob$V8, "_", 1) %>% splitter("\\.", 3)
    tob <- tob[,c("TF", "V5", "V6", "V7")]
    colnames(tob) <- c("TF", "Chr", "Start", "End")
    
  # Overlap
    snp2tfbs$TOBIAS <- "Not Evaluated"
    
    for (j in unique(snp2tfbs$TF)) { # for each TF
      
      print(j)
      
      # skip if not tested
      if (!(j %in% tob$TF)) next
      
      # get tobias bound coordinates for the given TF
      x <- tob[which(tob$TF == j),]
      write.bed(x[,-1], "../../Scratchspace/Temp_tobias.bed")
      
      # get coordinates for the given snp2tfbs
      i <- snp2tfbs$TF == j
      y <- snp2tfbs[i,] 
      m <- match(y$rsID, snp.annot$SNP) # match works as I have checked that all snp2tfbs$rsID are in snp.annot$SNP rather than $Synonym
      y <- data.frame(Chr = splitter(snp.annot$Peak[m], ":", 1),
                      Start = snp.annot$SNP.Pos[m],
                      End = snp.annot$SNP.Pos[m] + 1,
                      ID = y$rsID)
      write.bed(y, "../../Scratchspace/Temp_snp2tfbs.bed")
      
      # intersect
      call <- paste("intersectBed",
                    "-a", "../../Scratchspace/Temp_snp2tfbs.bed",
                    "-b", "../../Scratchspace/Temp_tobias.bed",
                    "-c", 
                    ">", "../../Scratchspace/Temp_intersect.bed")
      
      system(call, intern = FALSE, wait = TRUE)   
      
      # annotate
      z <- read.delim("../../Scratchspace/Temp_intersect.bed", header = FALSE)
      z$TOBIAS <- "Not Bound"
      z$TOBIAS[which(z$V5 > 0)] <- "Bound"
      
      # final output
      snp2tfbs$TOBIAS[i] <- z$TOBIAS
 
    }


## Save
  write.csv(snp2tfbs, file = "SNP2TFBS Processed.csv", row.names = FALSE)

  
## Add to the annotation dataframe
  snp.annot <- check.snps(snp2tfbs$rsID, "snp2tfbs")
  snp.annot <- check.snps(snp2tfbs$rsID[which(snp2tfbs$TOBIAS == "Bound")], "snp2tfbs_TOBIAS")
  

  
################################################################################################################################ #
## Enhancer overlap to eQTLs ----
 
  
## Last check of GTEx
  # get the values from dapg_eqtl: https://zenodo.org/record/3517189#.ZGXGfHVByWR
  # found this from the UCSC track: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&c=chrX&g=gtexEqtlHighConf
  # and try bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/hg38/gtex/eQtl/gtexCaviar.bb -chrom=chr16 -start=34990190 -end=36727467 stdout

  # * `{tissue}.variants_pip.txt.gz` contains the variants' posterior inclusion probabilities at being causal for every gene.
  #   * gene: gene id (or intron id)
  #   * rank: ranking of the variant according to its PIP (see below)
  #   * variant_id: gtex variant id
  #   * pip: posterior inclusion probability of the variant in the causal models
  #   * log10_abf: approximate Bayes factor (-log10)
  #   * cluster_id: id of cluster to which the variant belongs 
  
  
  
  
  
## Versus GTEx
  ## Read in eQTLs
    # GTEx 
    gtex <- read.table("../../../../PublicData/GTEx/GTEX_cis_eQTLs_050422_ucscDownload", sep = "\t", header = TRUE)
    gtex <- gtex[,c("eqtlChrom", "eqtlStart", "eqtlEnd", "eqtlName", "geneName", "tissue", "cpp")]
    colnames(gtex) <- c("Chr", "Start", "End", "rsID", "Gene", "Tissue", "CPP")
    
  ## Overlap
    gtex.overlap <- detailed.overlap(gtex, score.col = "CPP")

  ## Save
    write.csv(gtex.overlap, file = "eQTL - GTEx.csv", row.names = FALSE)
    

## Versus multi-ancestry brain metaanalysis from Zeng 2022
  ## Read in (fine-mapped)
    mmqtl <- read.delim("../../../../PublicData/Zeng2022_eQTL/mmQTL_brain_meta_eqtl_finemap.tsv")
    colnames(mmqtl)[4] <- "rsID"
    mmqtl$Tissue <- "Brain"
    
    m <- match(mmqtl$Gene, geneInfo$EnsID)
    mmqtl$Gene <- geneInfo$Symbol[m]
    
    
  ## Overlap
    mmqtl.overlap <- detailed.overlap(mmqtl, score.col = "PP")
    mmqtl.overlap.all <-mmqtl.overlap # retain unfiltered version for downstream plotting purposes
    mmqtl.overlap <-  mmqtl.overlap[which(mmqtl.overlap$SNP.CPP > 0.1),] # removes those with CPP < 0.1
    
  ## Save
    write.csv(mmqtl.overlap, file = "eQTL - mmQTL.csv", row.names = FALSE)
    
 
## Versus Metabrain, from de Klein 2023
  ## Read in
    
    # this resource has a similar scope to mmqtl
    metabrain <- read_xlsx("../../../../PublicData/deKlein2023_eQTL/Supplementary Table 2 - Cis-eQTLs.xlsx", sheet = "Cortex-EUR", skip = 1)  
      
    # I note that the authors report 16619 hits a q < 0.05, so filter to this 
    metabrain$qval <- gsub("E", "e", metabrain$qval) %>% as.numeric() 
    metabrain <- metabrain[which(metabrain$qval < 0.05),]# all already are?
    
    # extract rs id
    metabrain$rsID <- splitter(metabrain$SNP, ":", 3) # I have confirmed that all have an associated rsID
    
    # colnames
    metabrain$Gene <- metabrain$GeneSymbol
    metabrain$Tissue <- "Brain"
    
  ## Overlap
    metabrain.overlap <- detailed.overlap(metabrain, score.col = c("MetaP", "MetaBeta", "qval"), score.makecolname = c("P", "Beta", "Qval"))
    
  ## Save
    write.csv(metabrain.overlap, file = "eQTL - Metabrain.csv", row.names = FALSE)
  

## snRNA-seq from Bryois et al 2018. 
    bryois <- read_xlsx("../../../../PublicData/Bryois2021_eQTL_snRNAseq/41593_2022_1128_MOESM3_ESM.xlsx", sheet = "Table S3", skip = 3) # fine mapped eQTLs
    colnames(bryois) <- c("Tissue", "Gene", "EnsID", "rsID", "Chr", "Pos", "CaVEMaN", "CPP") 
    
  ## Overlap
    bryois.overlap <- detailed.overlap(bryois, score.col = "CPP")
    
  ## Save
    write.csv(bryois.overlap, file = "eQTL - Bryois snRNAseq.csv")

## Annotate snp.annot object
  snp.annot <- check.snps(gtex.overlap$SNP, "eqtl_gtex")
  snp.annot <- check.snps(gtex.overlap$SNP[which(gtex.overlap$Consistent == "TRUE")], "eqtl_gtex_consistent")
  
  snp.annot <- check.snps(mmqtl.overlap$SNP, "eqtl_mmqtl")
  snp.annot <- check.snps(mmqtl.overlap$SNP[which(mmqtl.overlap$Consistent == "TRUE")], "eqtl_mmqtl_consistent")
  
  snp.annot <- check.snps(metabrain.overlap$SNP, "eqtl_metabrain")
  snp.annot <- check.snps(metabrain.overlap$SNP[which(metabrain.overlap$Consistent == "TRUE")], "eqtl_metabrain_consistent")
  
  snp.annot <- check.snps(bryois.overlap$SNP, "eqtl_bryois")
  snp.annot <- check.snps(bryois.overlap$SNP[which(bryois.overlap$Consistent == "TRUE")], "eqtl_bryois_consistent")
  
  
## Some basic explorations of trends
  # ## Create combined dataframe
  #   eqtl <- list(GTEx = gtex.overlap,
  #                GTEx_Brain = gtex.overlap[grep("Brain", gtex.overlap$SNP.Tissue),],
  #                Metabrain = metabrain.overlap,
  #                mmQTL = mmqtl.overlap,
  #                snRNAseq = bryois.overlap,
  #                snRNAseq_Astrocytes = bryois.overlap[which(bryois.overlap$SNP.Tissue == "Astrocytes"),])
  #   
  #   # base dataframe
  #   eqtl.sum <- data.frame(Enh = unique(res.final$Enh))
  #   eqtl.sum$Hit <- eqtl.sum$Enh %in% hit.enh
  #   
  #   # binary annotation
  #   list.all.enh <- function(overlaps) {
  #     # all enh with a snp in the database
  #     x <- strsplit(overlaps$LinkedEnh, "/") %>% do.call("c", .)
  #     
  #     # check which are nearby and overlapping
  #     y <- strsplit(overlaps$Category, "/") %>% do.call("c", .)
  #     keep <- y %in% c("Overlapping", "Nearby")
  #     
  #     # filter
  #     x <- x[which(keep)]
  #     x <- unique(x)
  #     
  #     return(x)
  #   }
  #   
  # ## Are hit enhancers likelier to be near eQTLs?
  #   eqtl.sum$GTEx_Brain <- eqtl.sum$Enh %in% list.all.enh(gtex.overlap[grep("Brain", gtex.overlap$SNP.Tissue),])
  #   eqtl.sum$GTEx_Nonbrain <- eqtl.sum$Enh %in% list.all.enh(gtex.overlap[-grep("Brain", gtex.overlap$SNP.Tissue),])
  #   
  #   eqtl.sum$mmQTL <- eqtl.sum$Enh %in% list.all.enh(mmqtl.overlap)
  #   eqtl.sum$Metabrain <- eqtl.sum$Enh %in% list.all.enh(metabrain.overlap)
  #   
  #   eqtl.sum$snRNAseq_Astrocytes <- eqtl.sum$Enh %in% list.all.enh(bryois.overlap[grep("Astro", bryois.overlap$SNP.Tissue),])
  #   eqtl.sum$snRNAseq_Nonastro <- eqtl.sum$Enh %in% list.all.enh(bryois.overlap[-grep("Astro", bryois.overlap$SNP.Tissue),])
  
  ## For each screen pair, what does the eQTL evidence say?
    scaffold <- res.final[which(res.final$HitPermissive), c("Pair", "Enh", "Gene")]
    
    # function
    check.pair <- function(overlaps, data_in = scaffold) {
      w <- overlaps[which(overlaps$LinkedEnh.Hit == "TRUE"),]
      
      # all enh with a snp in the database
      x <- strsplit(w$LinkedEnh, "/") 
      names(x) <- paste0(w$SNP.Association, "_")
      x <- do.call("c", x)
      
      # check which are nearby and overlapping
      y <- strsplit(w$Category, "/") %>% do.call("c", .)
      y[grep("LD", y)] <- "LD"
      
      # combine
      z <- data.frame(Enh = x,
                      eQTL = splitter(names(x), "_", 1),
                      Linkage = y)
      
      z <- unique(z)
      rm(x, y)
      
      # return categorisation
      data_in$eQTL_Nearby_Consistent <- data_in$eQTL_LD_Consistent <- data_in$eQTL_Nearby <- data_in$eQTL_LD <- data_in$eQTL_Consistent <- data_in$eQTL_Present <- NA
      
      for (j in 1:nrow(data_in)) {
        
        k <- data_in$Enh[j]
        
        # check if there's an eqtl for this enhancer
        if (k %in% z$Enh) {
          
          # note that there is an eqtl
          data_in$eQTL_Present[j] <- TRUE
          
          # is there a match in gene target?
          targets <- z[which(z$Enh == k),]
          data_in$eQTL_Consistent[j] <- data_in$Gene[j] %in% targets$eQTL
          
          # and when looking at nearby SNPs only
          data_in$eQTL_Nearby[j] <- any(targets$Linkage != "LD")
          if (data_in$eQTL_Nearby[j]) data_in$eQTL_Nearby_Consistent[j] <- (data_in$Gene[j] %in% targets$eQTL[which(targets$Linkage != "LD")])
          
          # and when looking at LD SNPs only
          data_in$eQTL_LD[j] <- any(targets$Linkage == "LD")
          if (data_in$eQTL_LD[j]) data_in$eQTL_LD_Consistent[j] <- (data_in$Gene[j] %in% targets$eQTL[which(targets$Linkage == "LD")]) 
          
        } else {
          
          data_in$eQTL_Present[j] <- FALSE
          next 
          
        }
        
        
      }
      
      return(data_in)
    }
    
    # apply function
    checked.pairs <- list()
    
    g <- grep("Brain", gtex.overlap$SNP.Tissue)
    checked.pairs$GTEx_Brain <- check.pair(gtex.overlap[g,])
    checked.pairs$GTEx_Nonbrain <- check.pair(gtex.overlap[-g,])
    
    checked.pairs$mmQTL <- check.pair(mmqtl.overlap)
    checked.pairs$Metabrain <- check.pair(metabrain.overlap)
    
    g <- grep("Astro", bryois.overlap$SNP.Tissue)
    checked.pairs$snRNAseq_Astro <- check.pair(bryois.overlap[g,])
    checked.pairs$snRNAseq_Nonastro <- check.pair(bryois.overlap[-g,])
    
    # stats within each eqtl resource
    stats.within <- sapply(checked.pairs, function(x) colSums(x[,-c(1:3)], na.rm = TRUE) )
    
    stats.within <- t(stats.within) %>% as.data.frame()
    stats.within <- stats.within[,1:2]
    stats.within$Study <- rownames(stats.within)
    stats.within$eQTL_Inconsistent <- stats.within$eQTL_Present - stats.within$eQTL_Consistent
    stats.within$eQTL_None <- length(hit.enh) - stats.within$eQTL_Present
    stats.within <- melt(stats.within[,-1])
    
    stats.within$variable <- gsub("eQTL_", "", stats.within$variable)
    stats.within$Study <- gsub("_", " ", stats.within$Study)
    stats.within$variable <- factor(stats.within$variable, levels = (c("None", "Inconsistent", "Consistent")))
    levels(stats.within$variable) <- (c("No eQTL", "Inconsistent eQTL(s)", "1+ Supporting eQTL"))
    
    # cols <- carto_pal(3, "Earth")[c(1,3)] %>% c("grey80", .)
    cols <- c("grey80", carto_pal(3, "ArmyRose")[c(1,3)])
    pdf(file = "eQTL - Support Rate.pdf", height = 3, width = 5)
    ggplot(stats.within, aes(x = Study, y = value, fill = variable, alpha = variable)) +
      geom_col(position = "stack", width = 0.8, colour = "black", size = 0.5) +
      theme_bw() +
      # scale_y_continuous(expand = c(0,0)) +
      labs(y = "Screen Enhancer-Gene Pairs", x = "eQTL Resource") +
      scale_fill_manual(values = cols) +
      scale_alpha_manual(values = c(0.85, 0.85, 1)) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "right", 
            legend.title = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.title.x = invis, panel.grid = invis, axis.ticks.x = invis)
    dev.off()
    
    # stats across eqtl resources
    stats.across <- sapply(checked.pairs, function(x) {
      x <- x[,c("eQTL_Present", "eQTL_Consistent")]
      x$Category <- "No eQTL"
      x$Category[which(x$eQTL_Present)] <- as.character(x$eQTL_Consistent[which(x$eQTL_Present)])
      x$Category <- factor(x$Category, levels = c("No eQTL", "FALSE", "TRUE"))
      levels(x$Category) <- c("No eQTL", "Inconsistent eQTL(s)", "1+ Supporting eQTL")
      return(x$Category)
    })
    
    ord <- scaffold$Pair[order(rowSums(stats.across == "1+ Supporting eQTL"),
                               rowSums(stats.across == "Inconsistent eQTL(s)"))]
    stats.across <- as.data.frame(stats.across)
    stats.across$Pair <- scaffold$Pair
    stats.across <- melt(stats.across, id.vars = "Pair")
    stats.across$Pair <- factor(stats.across$Pair, levels = ord)
    stats.across$variable <- factor(stats.across$variable, levels = levels(stats.across$variable)[c(2,1,4,3,6,5)])
    
    pdf(file = "eQTL - Support Heatmap.pdf", height = 11, width = 6)
    ggplot(stats.across, aes(x = variable, y = Pair, fill = value, alpha = value)) +
      geom_tile() +
      scale_fill_manual(values = rev(cols)) +
      scale_alpha_manual(values = c(1, 0.6, 0.5)) +
      theme_bw() +
      theme(panel.border = invis, legend.position = "right", 
            legend.title = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8),
            # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
            axis.title.x = invis, panel.grid = invis, axis.ticks = invis,
            axis.text.y = element_text(size = 6))
    dev.off()
          
  ## And the effect of a stat like Causal Prior Probability?
    # use mmqtl, as that has the range from 0-1.
    p <- mmqtl.overlap.all[which(mmqtl.overlap.all$LinkedEnh.Hit == "TRUE"),]
    p <- data.frame(Consistent = p$Consistent, CPP = p$SNP.CPP)
    
    pdf(file = "eQTL - Prior Probability in mmQTL.pdf", height = 3.5, width = 4)
    ggplot(p, aes(x = Consistent, y = CPP)) +
      geom_quasirandom(size = 1) +
      scale_y_continuous(trans = "log10", limits = c(NA, 1)) +
      stat_summary(fun = median, geom = "point", shape = "+", colour = "red", size = 10) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) +
      labs(x = "eQTL and Enhancer Sharing Target Gene", y = "eQTL Prior Probability")
    dev.off()
    
    table(p$CPP > 0.1, p$Consistent) %>% fisher.test() # p = 0.002, OR = 3.4
    wilcox.test(p$CPP ~ p$Consistent) # p = 1e-12
    
    # and, in metabrain, pvalue and beta
    p <- metabrain.overlap[which(metabrain.overlap$LinkedEnh.Hit == "TRUE"),]
    p$P <- gsub("E", "e", p$P) %>% as.numeric()
    
    pdf(file = "eQTL - Qvalue and Beta in Metabrain.pdf", height = 3.5, width = 4)
    ggplot(p, aes(x = Consistent, y = -log10(Qval))) +
      geom_quasirandom(size = 1) +
      # scale_y_continuous(trans = "log10", limits = c(NA, 1)) +
      stat_summary(fun = median, geom = "point", shape = "+", colour = "red", size = 10) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) +
      labs(x = "eQTL and Enhancer Sharing Target Gene", y = "eQTL -log10(Q-value)")
    
    ggplot(p, aes(x = Consistent, y = abs(Beta))) +
      geom_quasirandom(size = 1) +
      # scale_y_continuous(trans = "log10", limits = c(NA, 1)) +
      stat_summary(fun = median, geom = "point", shape = "+", colour = "red", size = 10) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) +
      labs(x = "eQTL and Enhancer Sharing Target Gene", y = "Beta (Absolute)")
    dev.off()
    
    table(-log10(p$Qval) > 5, p$Consistent) %>% fisher.test() # OR = 2.57, p = 0.38
    wilcox.test(abs(p$Qval) ~ p$Consistent) # p = 0.88
    
    # gtex cpp
    p <- gtex.overlap[which(gtex.overlap$LinkedEnh.Hit == "TRUE"),]
    p <- data.frame(Consistent = p$Consistent, CPP = p$SNP.CPP, Brain = grepl("Brain", p$SNP.Tissue))
    
    pdf(file = "eQTL - Prior Probability in GTEx.pdf", height = 3.5, width = 4)
    ggplot(p[which(p$Brain),], aes(x = Consistent, y = CPP)) +
      geom_quasirandom(size = 1) +
      scale_y_continuous(trans = "log10", limits = c(NA, 1)) +
      stat_summary(fun = median, geom = "point", shape = "+", colour = "red", size = 10) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) +
      labs(x = "eQTL and Enhancer Sharing Target Gene", y = "GTEx Brain\neQTL Prior Probability")
    
    ggplot(p[-which(p$Brain),], aes(x = Consistent, y = CPP)) +
      geom_quasirandom(size = 1) +
      scale_y_continuous(trans = "log10", limits = c(NA, 1)) +
      stat_summary(fun = median, geom = "point", shape = "+", colour = "red", size = 10) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis) +
      labs(x = "eQTL and Enhancer Sharing Target Gene", y = "GTEx NonBrain\neQTL Prior Probability")
    dev.off()
    
    q <- p[which(p$Brain),]
    table(q$CPP > 0.2, q$Consistent) %>% fisher.test() # p = 0.03, OR = 2.4
    wilcox.test(q$CPP ~ q$Consistent) # p = 0.007
    
################################################################################################################################ #
## Enhancer overlap to clinical variants from DisGeNet ----


## Read in
  disg <- read.delim("../../../../PublicData/DisGeNet_070422/all_variant_disease_associations.tsv", sep = "\t", header = TRUE)

## Wrangle
  disg$rsID <- disg$snpId
  disg$Gene <- disg$diseaseName
  disg$Tissue <- disg$diseaseSemanticType
  
## Overlap
  disg.overlap <- detailed.overlap(disg, score.col = "score", score.makecolname = "Disg.Score", check.consistency = FALSE)
      
## Save
  write.csv(disg.overlap, "Disgenet.csv")
  
## Add to annotation dataframe
  snp.annot <- check.snps(disg.overlap$SNP, "disgenet")
  
  
################################################################################################################################ #
## Enhancer overlap to GWAS ----
  

## Read in
  gwascat <- read.delim("../../../../PublicData/gwasCatalog_050422_ucscDownload", sep = "\t", header = TRUE)  

## Wrangle
  gwascat$rsID <- gwascat$name
  gwascat$Gene <- gwascat$trait
  gwascat$Tissue <- NA
    
## Overlap
  gwascat.overlap <- detailed.overlap(gwascat, score.col = "pValue", score.makecolname = "P", check.consistency = FALSE)
  gwascat.overlap <- gwascat.overlap[,-grep("Tissue", colnames(gwascat.overlap))]
      
## Save
  write.csv(gwascat.overlap, "GWAS Catalogue.csv")
  

## Add to annotation dataframe
  snp.annot <- check.snps(gwascat.overlap$SNP, "gwasCatalogue")
  
  
################################################################################################################################ #
## Enhancer overlap to AD GWAS ----
  

## As our hit genes are strongly associated with neurodegeneration and Alzheimer's, we shall annotate our SNPs with AD GWAS effects
  
  
## Dataset
  ## The latest is Bellenguez 2022 Nat Genet: https://www.nature.com/articles/s41588-022-01024-z
    # stage 1 analyses: based on 39,106 clinically diagnosed AD cases, 46,828 proxy-ADD cases, and 401,577 control. 
    # keep all variants at p < 1e-5, collect all non-overlapping regions around, and run Stage II analyses (25,392 AD cases and 276,086 controls)
    # significant if nominal < 1e-5 in both stages, same direction in both stages, and <5e-8 in meta-analysis of both stages
  
    bell <- list()
    bell$Hits <- read.delim("../../../../PublicData/Bellenguez2022_AD_GWAS_Latest/gwas-association-downloaded_2023-07-31-pubmedId_35379992.tsv")
    bell$SumStat <- read.delim("../../../../PublicData/Bellenguez2022_AD_GWAS_Latest/GCST90027158_buildGRCh38.tsv")
    colnames(bell$SumStat)[1] <- "rsID"
  
    
  
  ## Schwarzentruber presented a fine-mapped 2021 analysis of the UKBB: https://www.nature.com/articles/s41588-020-00776-w
  
  
  
## Raw p-value for all SNPs
  ad_allSNP <- detailed.overlap(bell$SumStat, score.col = "p_value", check.consistency = FALSE, score.makecolname = "P", gwas.format = TRUE) # requires ~5min

## Lead variants
  ad_leadSNP <- read.csv("Final - SNP-to-Peak List.csv")
  ad_leadSNP <- check.snps(bell$Hits$SNPS, "Bellenguez2023_Top")
  
## Save
  save(ad_allSNP, ad_leadSNP, file = "AD GWAS Bellenguez 2022.rda")
  
## Process
  x <- ad_allSNP[which(ad_allSNP$p_value < 5e-8),]
  x <- x[,c("rsID", "Enh", "LinkageCategory", "p_value", "odds_ratio", "beta", "n_cases", "n_controls")]
  
## Write out
  
################################################################################################################################ #
## Output ---- 

  
## Output
  write.csv(snp.annot, file = "Final - SNP Annotation.csv")
      
      