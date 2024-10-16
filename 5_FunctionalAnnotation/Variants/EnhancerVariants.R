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
  source("../../../../Manuscript/Figs/FinalFigureFunctions.R")
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

 
## Function
  # get synonyms for rsIDs
  get_synonyms <- function(dat) {
    dat$Synonym <- NA
    
    # match your id to main id, extract the archive
    m <- match(dat$SNP, syn$main)
    dat$Synonym <- syn$synonyms[m]
    
    # match your id to archive id, extract the main
    m <- match(dat$SNP, syn.split$Archive)
    dat$Synonym <- paste0(dat$Synonym, ",", syn.split$Main[m], ",")
    
    # clean
    dat$Synonym <- gsub("NA,", "", dat$Synonym) %>%
      gsub(",", "|", .) # separator
    
    # return
    return(dat)
  }
  
  # calculate distance to edge
  get_distanceToEdge <- function(dat, l = peak.start, r = peak.end) {
    apply(dat, 1, function(x) {
      g <- grep("SNP.pos", colnames(dat))
      left <- as.numeric(x[g]) - l
      right <- as.numeric(x[g]) - r
      dist <- c(left, right)
      w <- which.min(abs(dist))  # what is the index of the shortest absolute distance?
      return(dist[w]) # get shortest distance, including its sign
    })
  }
  

################################################################################################################################ #
## First: annotate enhancers with SNPs from dbSNP ----
  
## This section will use dbSNP153 to tag SNPs to enhancers, based on proximity (overlapping or within 1kb)
## as well as linkage disequilibrium
## This allows us to annotate the SNPs by possible functionality, as well as explore phenotypic effects

## Paths
  db_dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/dbSNP/dbSNP153_280422.bed"
  nha_hg38 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  db_out <- "dbSNP_Overlap.bed"
  
  
## Read in large SNP annotation files  
   saved_snp_files <- "dbSNP and synonyms, formatted.rda"
  
  if (! file.exists(saved_snp_files)) {
    
    message("Creating and saving")
    
    # dbsnp 153
    db153 <- read.delim(db_dir, header = FALSE)

    # snp synonyms, from: # http://genehopper.ifis.cs.tu-bs.de/downloads, Ensembl VariationDB 84
    syn <- read.table("../../../../PublicData/dbSNP/variation_synonyms.txt", sep = "\t", header = TRUE) 
    syn <- syn[-grep("N", syn$synonyms),] # removes all those without a synonym
    
    # snp synonyms, alternative formatting
    syn.split <- strsplit(syn$synonyms, ",")
    names(syn.split) <- paste0(syn$main, "_") # the paste helps with matching downstream
    syn.split <- do.call("c", syn.split)
    syn.split <- data.frame(Main = splitter(names(syn.split), "_", 1),
                    Archive = syn.split)
  
    # save
    save(db153, syn, syn.split, file = saved_snp_files)
      
  } else {
    
    message("Loading")
    load(saved_snp_files, verbose = TRUE)
    
  }
  
  
## Intersect our enhancers with dbSBP153
    call <- paste("windowBed",
                "-a", nha_dir_38,
                "-b", db_dir,
                "-w", "1000", # window of 1000bp
                ">", db_out)
  
    system(call, intern = FALSE, wait = TRUE) 
  
## Analyse
  db_overlap <- read.delim(db_out, header = FALSE)
    
  # column names
  colnames(db_overlap) <- c("Peak.chr", "Peak.start", "Peak.end", "Peak.coord", "SNP.chr", "SNP.start", "SNP.end", "SNP.id")
    
  # add hit annotation
  db_overlap$Peak.coord <- sub("_", ":", db_overlap$Peak.coord) %>% sub("_", "-", .)
  
## Collect variants in LD to these variants
  ## Calculate LD
    # requires ~24h
  
    run.LDcode <- FALSE
    
    if (run.LDcode) {
      
      library(LDlinkR)
      
      # get all snps nearby
      useSNP <- unique(db_overlap$SNP.id)
      useSNP_syn <- syn$synonyms[which(syn$main %in% useSNP)]
      useSNP_syn <- strsplit(useSNP_syn, ",")
      useSNP_syn <- do.call("c", useSNP_syn) %>% unique()
      use <- union(useSNP, useSNP_syn)
      
      LDproxy_batch(snp = use, # unique is necessary as these peaks have a window appended
                    pop = "EUR", # european
                    r2d = "r2",
                    append = TRUE,
                    genome_build = "grch38",
                    token = "5988cfbf2737")  
    }
    
  
  ## Read in
      if (! file.exists("combined_query_snp_list_grch38_filtered.rda")) {
        message("Loading and processing LD SNPs")
        
        # all.ld.snps <- read.table("combined_query_snp_list_grch38.txt", sep = "\t", header = TRUE, row.names = NULL)
        all.ld.snps <- read.table("LDLinkR/combined_query_snp_list_grch38.txt", sep = "\t", header = TRUE, row.names = NULL)
        all.ld.snps <- all.ld.snps[,c("query_snp", "RS_Number", "Coord", "MAF", "Distance", "R2")]
        colnames(all.ld.snps)[1:3] <- c("Query_SNP", "Linked_SNP", "Position")
      
        # filter to high ld snps
        ld.snps <- all.ld.snps[which(all.ld.snps$R2 > 0.8 & all.ld.snps$Distance != 0),] # subset to r2 > 0.8, and not same-snp
        save(ld.snps, file = "combined_query_snp_list_grch38_filtered.rda")
        
      } else {
        
        message("Loading LD SNPs")
        load("combined_query_snp_list_grch38_filtered.rda", verbose = TRUE)
        
      }
    
    
    

## Combine peak annotations to snps
  snp2peak <- list()
  tested.peaks <- res.final$Enh.Pos %>% unique()

  a <- Sys.time()
  for (j in tested.peaks) { # for each tested peak (957), rather than all peaks (979)
    
    print(j)
    
    ## First, look at nearby (within 1kb of edge) SNPs
      # skip if no snps within 1kb window
      if ( ! ( j %in% db_overlap$Peak.coord)) {
       
         snp2peak[[j]] <- ("None") 
        
        next
         
      }
    
      # collect snps within 1kb window
      j_near <- db_overlap[which(db_overlap$Peak.coord == j),]
      peak.start <- j_near$Peak.start[1]
      peak.end <- j_near$Peak.end[1]
      peak.centre <- mean(c(peak.start, peak.end)) %>% round()
      j_near <- j_near[,c("Peak.coord", "SNP.id", "SNP.start")]
      colnames(j_near) <- c("Peak", "SNP", "SNP.pos")
      
      j_near$R2 <- NA
      j_near$Category <- "."
    
      # annotate as overlapping or nearby
      j_overlaps <- ((j_near$SNP.pos >= peak.start) & (j_near$SNP.pos <= peak.end)) 
      
      if (any(j_overlaps)) j_near$Category[j_overlaps] <- "Overlapping"
      if (any(!(j_overlaps))) j_near$Category[!(j_overlaps)] <- "Nearby"
      
      # annotate with distance 
      j_near$DistanceFromCentre <- (j_near$SNP.pos - peak.centre) # uses start as a simplification
      j_near$DistanceFromEdge <- get_distanceToEdge(j_near)
      
    ## Get synonyms from other versions of dbSNP
      # apply function
      j_near <- get_synonyms(j_near)
      
      # correspondence between main and synonyms
      j_hasSyns <- j_near$Synonym != ""

      if (any(j_hasSyns)) {
        j_syn.split <- j_near[which(j_hasSyns),]
        j_syn.split <- strsplit(j_syn.split$Synonym, "|", fixed = TRUE)
        names(j_syn.split) <- paste0(j_near$SNP[which(j_hasSyns)], "_")
        j_syn.split <- do.call("c", j_syn.split)
        j_syn.split <- data.frame(Main = splitter(names(j_syn.split), "_", 1),
                        Archive = j_syn.split)
        colnames(j_syn.split)[2] <- "Archive"
      }
      
    ## Collect SNPs in ld 
      # vector of rsID j_near
      j_rsMain <- j_near$SNP
      j_rsSyns <- strsplit(j_near$Synonym, "|", fixed = TRUE)
      j_rsAll <- unique(j_rsMain, do.call("c", j_rsSyns))
      
      # check if there are any with ld data
      if (any(j_rsAll %in% ld.snps$Query_SNP)) { 
        
        # get
        j_ld <- ld.snps[which(ld.snps$Query_SNP %in% j_rsAll),] # collect
        
        # reformat
        j_ld <- data.frame(Peak = j,
                        SNP = j_ld$Linked_SNP, 
                        SNP.pos = splitter(j_ld$Position, ":", 2) %>% as.numeric(),
                        R2 = j_ld$R2,
                        Category = "LD",
                        DistanceFromCentre = NaN,
                        DistanceFromEdge = NaN)
        
        j_ld <- get_synonyms(j_ld)
        j_ld$DistanceFromCentre <- (j_ld$SNP.pos - peak.centre)
        j_ld$DistanceFromEdge <- get_distanceToEdge(j_ld)
        
        # reannotate ld snps that are nearby/overlapping as such
        actually_overlapping <- ((j_ld$SNP.pos >= peak.start) & (j_ld$SNP.pos <= peak.end))   
        if (any (actually_overlapping)) j_ld$Category[which(actually_overlapping)] <- "Overlapping"
        
        actually_near <- abs(j_ld$DistanceFromEdge) <= 1000 & !(actually_overlapping) 
        if (any (actually_near)) {
          j_ld$Category[which(actually_near)] <- "Nearby"
          j_ld$R2[which(actually_near)] <- NA
        } 
        
        # filter 1: remove LD SNPs that are duplicated
        j_ld <- unique(j_ld) # both steps are necessary
        j_ld <- j_ld[order(j_ld$R2, decreasing = TRUE),] # this step ensures that the highest LD is kept first in case of duplicates
        j_dup1 <- (duplicated(j_ld$SNP)) & (j_ld$SNP != ".")
        if (any(j_dup1)) j_ld <- j_ld[-which(j_dup1),]
        
        # filter 2: remove LD SNPs which are also nearby! 
        i <- grep("Synonym", colnames(j_ld))

        j_dup2 <- (j_ld$SNP.pos %in% j_near$SNP.pos) | # same position
          (j_ld$SNP %in% j_rsAll) | # same main rsid
          apply(j_ld, 1, function(x) { any(strsplit(x[i], "|", fixed = TRUE) %in% j_rsAll) }) # any of the synonyms are the same
        
        if (any(j_dup2)) {
          
          j_dup2 <- which(j_dup2)
          
          if (length(j_dup2) == length(j_overlaps)) { # if removing as many rows as there are in the dataframe, may cause error. instead, finish loop
            snp2peak[[j]] <- j_near
            next
          }
          
          j_ld <- j_ld[-j_dup2,]  
          
        }
        

        # combine
        j_out <- rbind(j_near, j_ld)
        
        
      } else { # if no LD snps
        j_out <- j_near
      }
      
  
    ## Return
      snp2peak[[j]] <- j_out
    
  }
  b <- Sys.time()
    
  
## Process
  snp2peak <- do.call("rbind", snp2peak)
    
  # remove the "none", i.e. when there is a peak with no nearby SNPs
  snp2peak <- snp2peak[-which(snp2peak$Category == "None"),]
  
  # remove SNPs with no rsID, as rsID is used for matching 
  snp2peak <- snp2peak[-which(snp2peak$SNP == "."),] # remove snps without an rsid
  
  # remove SNPs that are not in dbSNP153
  # snp2peak <- snp2peak[which(snp2peak$SNP %in% db153[,4]),] 
  # the above will be done later, after extracting synonyms
  
  # add the enhancer id
  m <- match(snp2peak$Peak, guides$TargetCoord)
  snp2peak$Enh <- guides$TargetID[m] 
  
  # is it a hit enhancer?
  snp2peak$Hit <- snp2peak$Enh %in% hit.enh
  
  
    
## Remove any SNP without an entry in dbSNP153
  # here, your original calls of nearby were based on dbSNP153, but LD was calculated from matches to dbSNP151
  # remove any that are exclusive to 151 for which you couldn't find the 153 synonym

  # concatenate main and synonym for each row
  keepcheck1 <- apply(snp2peak, 1, function(x) {

    y <- paste(x[which(colnames(snp2peak) == "SNP")], x[which(colnames(snp2peak) == "Synonym")], sep = "|")
    y <- strsplit(y, "\\|")[[1]]

  })

  # pre-filter dbSNP153
  keepcheck2 <- do.call("c", keepcheck1)
  keepcheck3 <- db153$V4
  keepcheck3 <- keepcheck3[which(keepcheck3 %in% keepcheck2)]

  # check if each row has at least one main/synonym (keepcheck1) in db153 (keepcheck3)
  in153 <- sapply(keepcheck1, function(x) {
    any(x %in% keepcheck3)
  })

  # snp2peak <- snp2peak[which(in153),]
  snp2peak$in153 <- in153


  
## Final clean
  snp2peak <- relocate(snp2peak, "Enh", "Hit", "SNP", "Synonym")
  rownames(snp2peak) <- 1:nrow(snp2peak)
  
## Save
  write.csv(snp2peak, file = "Final - SNP-to-Peak List.csv", row.names = FALSE)

   
  
################################################################################################################################ #
## Setup annotation process ----
   
## Now start annotating
  snp.annot <- read.csv("Final - SNP-to-Peak List.csv")

  
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
## General characterisation ----
    
## Number of SNPs per enhancer
  p <- snp.annot
  p <- p[which(p$Hit),]
  
  tab <- table(p$Enh, p$Category)
  tab <- as.data.frame.matrix(tab)
  tab$Total <- rowSums(tab)
  tab$Enh <- rownames(tab)
  tab$Enh <- factor(tab$Enh, levels = tab$Enh[order(tab$Total)])
  tab <- melt(tab, id.vars = c("Enh", "Total"))
  tab$variable <- factor(tab$variable, levels = rev(levels(tab$variable)))
  
  offset <- 1
  
  # as violin
  pdf(file = "Count of SNPs per enhancer.pdf", height = 3, width = 3)
  ggplot(tab, aes(x = variable, y = value+offset, fill = variable, colour = variable)) +
    geom_violin(draw_quantiles = 0.5, scale = "width") +
    labs(x = "SNP category", y = "SNPs per\nfunctional enhancer") +
    theme_bw() +
    scale_fill_manual(values = pals$Primary[c(4,5,7)]) +
    scale_colour_manual(values = pals$Primary_Darker[c(4,5,7)]) +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(),
          legend.position = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(1, 1025),
                       trans = "log2", labels = function(x) { comma(x - offset) }, 
                       breaks = c(0, 1, 4, 16, 64, 256, 1024) + offset) 
  
  ggplot(tab, aes(x = ".", y = Total+offset)) +
      geom_violin(draw_quantiles = 0.5, scale = "width", width = 0.3) +
      labs(x = "SNP category", y = "SNPs per\nfunctional enhancer") +
      theme_bw() +
      theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(),
            legend.position = "none") +
      scale_y_continuous(expand = c(0,0), limits = c(1, 1025),
                         trans = "log2", labels = function(x) { comma(x - offset) }, 
                         breaks = c(0, 1, 4, 16, 64, 256, 1024) + offset) 
    dev.off()
    
   # as histogram
  tab$Bin <- cut(tab$Total, breaks = c(0,1,5,10,20,50,100,500,1000), include.lowest = TRUE,
                 labels = c("0", "1-4", "5-9", "10-19", "20-49", "50-99", "100-499", "500-999"))
  
  p <- table(tab$Bin) %>% as.data.frame.table()
  
  pdf(file = "Count of SNPs per enhancer (Hist).pdf", height = 2.5, width = 4)
  ggplot(p, aes(x = Var1, y = Freq)) +
    geom_col(fill = pals$One) +
    labs(y = "Count of\nfunctional enhancers", x = "Assigned SNPs") +
    theme_bw() +
    scale_y_continuous(expand  = c(0,0), limits = c(0, 160), breaks = c(0, 50, 100, 150)) +
    theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line(),
          legend.position = "none", axis.text.y = text90) 
  dev.off()
 
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
  snp2tfbs <- read.delim("SNP2TFBS Webtool Output.txt", header = FALSE)
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
                      Start = snp.annot$SNP.pos[m],
                      End = snp.annot$SNP.pos[m] + 1,
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
## Enhancer overlap to eQTLs ----
 

## Versus GTEx, high confidence fine-mapped
  ## Read in eQTLs
    # gtex <- read.delim("../../../../PublicData/GTEx/CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.txt")
    gtex_cav <- read.delim("../../../../PublicData/GTEx/CAVIAR/CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF.txt")
    colnames(gtex_cav) <- c("Tissue", "Gene", "eQTL", "Chr", "Coord", "CPP")
  
  ## Process
    # add gene symbol
    gtex_cav$Gene <- splitter(gtex_cav$Gene, "\\.", 1)
    m <- match(gtex_cav$Gene, geneInfo$EnsID)
    gtex_cav$Gene <- geneInfo$Symbol[m]
    gtex_cav <- gtex_cav[-which(is.na(gtex_cav$Gene)),]
    
    # rename tissue
    gtex_cav$Tissue <- gsub("_", " ", gtex_cav$Tissue)
    
    # add rsid (requires a few minutes to read in the gtexAnnot variable)
    # gtexAnnot <- read.delim("../../../../PublicData/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
    # gtexAnnot <- data.frame(Variant = paste0(gtexAnnot$chr, "_", gtexAnnot$variant_pos),
    #                 rsID = gtexAnnot$rs_id_dbSNP151_GRCh38p7)
    # save(gtexAnnot, file = "../../../../PublicData/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.rda")
    
    load("../../../../PublicData/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.rda", verbose = TRUE)
    gtex_cav$eQTL <- paste0("chr", gtex_cav$eQTL)
    m <- match(gtex_cav$eQTL, gtexAnnot$Variant)
    gtex_cav$rsID <- gtexAnnot$rsID[m]
    gtex_cav <- gtex_cav[grep("rs", gtex_cav$rsID),] # removes the non-matches, or matches where the variant has no rsID

    # save
    save(gtex_cav, file = "eQTLs - GTEx data, processed.rda")
    
  ## Filter 
    # to high confidence (CPP > 0.1)
    gtex_cav <- gtex_cav[which(gtex_cav$CPP > 0.1),]
    
    # to brain
    gtex_cav <- gtex_cav[grep("Brain", gtex_cav$Tissue),]
    
  ## Overlap
    gtex_cav.overlap <- detailed.overlap(gtex_cav, score.col = c("CPP"), score.makecolname = c("CPP"))
    write.csv(gtex_cav.overlap, file = "eQTL - GTEx (Caviar High Confidence).csv", row.names = FALSE)
    
    
## Versus GTEx, significant eQTLs
  ## Read in
    path <- "../../../../PublicData/GTEx/GTEX_Analysis_v8_eQTL/"
    l <- list.files(path, pattern = "Brain")
    l <- l[grep("signif_pairs", l)]

    gtex_sig <- lapply(l, function(x) {
      x <- paste0(path, x)
      print(x)
      read.delim(x)
    })

  ## Process
    names(gtex_sig) <- splitter(l, "\\.", 1)
    gtex_sig <- do.call("rbind", gtex_sig)
    gtex_sig <- gtex_sig[,c("phenotype_id", "variant_id", "pval_nominal")]
    gtex_sig$variant_id <- splitter(gtex_sig$variant_id, "_A|_C|_T|_G", 1)
    gtex_sig$Tissue <- splitter(rownames(gtex_sig), "\\.", 1) %>% gsub("_", " ", .)
    gtex_sig$phenotype_id <- splitter(gtex_sig$phenotype_id, "\\.", 1)
    
    # convert gene symbol
    m <- match(gtex_sig$phenotype_id, geneInfo$EnsID)
    gtex_sig$Symbol <- geneInfo$Symbol[m]
    gtex_sig <- gtex_sig[-which(is.na(gtex_sig$Symbol)),]
    
    # add rsID
    load("../../../../PublicData/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.rda", verbose = TRUE)
    m <- match(gtex_sig$variant_id, gtexAnnot$Variant)
    gtex_sig$rsID <- gtexAnnot$rsID[m]
    gtex_sig <- gtex_sig[grep("rs", gtex_sig$rsID),] # removes the non-matches, or matches where the variant has no rsID
    
    # filter columns
    rownames(gtex_sig) <- 1:nrow(gtex_sig)
    gtex_sig <- gtex_sig[,c("Symbol", "rsID", "Tissue", "pval_nominal")]
    
    # save
    save(gtex_sig, file = "eQTLs - GTEx data, processed (Sig).rda")
    
    # further filtering post-save: to SNPs in the screen, for memory reasons
    gtex_sig <- gtex_sig[which(gtex_sig$rsID %in% all.snps),]
    
    # rename columns
    colnames(gtex_sig)[1] <- "Gene"
    
  ## Overlap to our screen
    gtex_sig.overlap <- detailed.overlap(gtex_sig, score.col = c("pval_nominal"), score.makecolname = c("pval_nominal"))
    write.csv(gtex_sig.overlap, file = "eQTL - GTEx (Sig).csv", row.names = FALSE)
    
    
## Versus multi-ancestry brain metaanalysis from Zeng 2022
  ## Read in all p-value output
    mmqtl_fdr <- read.delim("../../../../PublicData/Zeng2022_eQTL/mmQTL_brain_meta_eqtl_all.tsv")
    mmqtl_fdr <- mmqtl_fdr[,c("Gene", "Variant", "eQTL_order", "z_score_random", "p_random")] # filter columns
    m <- match(mmqtl_fdr$Gene, geneInfo$EnsID)
    mmqtl_fdr$Gene <- geneInfo$Symbol[m]
    mmqtl_fdr <- mmqtl_fdr[-which(is.na(mmqtl_fdr$Gene)),]
    
    # calculate FDR
    mmqtl_fdr$FDR_random <- p.adjust(mmqtl_fdr$p_random)
    mmqtl_fdr <- mmqtl_fdr[which(mmqtl_fdr$FDR_random < 0.05),]
    
    # save
    save(mmqtl_fdr, file = "eQTL - mmQTL Processed Data (FDR sig).rda")
    
    # further filtering post-save: to SNPs in the screen, for memory reasons
    mmqtl_fdr <- mmqtl_fdr[which(mmqtl_fdr$Variant %in% all.snps),]
    
    # further processing
    colnames(mmqtl_fdr)[2] <- "rsID"
    mmqtl_fdr$Tissue <- "Brain"
    
    # mmqtl_fdr$HasFinemap <- mmqtl_fdr$Id %in% mmqtl$Id
    # mmqtl_fdr$Id <- paste0(mmqtl_fdr$Gene, mmqtl_fdr$Variant) 

  ## Read in fine-mapped output (passing significance and CAVIAR thresholds)
    mmqtl_cav <- read.delim("../../../../PublicData/Zeng2022_eQTL/mmQTL_brain_meta_eqtl_finemap.tsv")
    colnames(mmqtl_cav)[4] <- "rsID"
    mmqtl_cav$Tissue <- "Brain"
    
    m <- match(mmqtl_cav$Gene, geneInfo$EnsID)
    mmqtl_cav$Gene <- geneInfo$Symbol[m]
    mmqtl_cav$Id <- paste0(mmqtl_cav$Gene, mmqtl_cav$rsID)
    
    
  ## Overlap
    mmqtl_cav.overlap <- detailed.overlap(mmqtl_cav, score.col = c("PP", "eQTL_order"), score.makecolname = c("SNP.CPP", "eQTL_order"))
    mmqtl_fdr.overlap <- detailed.overlap(mmqtl_fdr, score.col = c("FDR_random", "eQTL_order"), score.makecolname = c("SNP.FDR", "eQTL_order"))

  ## Save
    write.csv(mmqtl_cav.overlap, file = "eQTL - mmQTL (Fine-mapped).csv", row.names = FALSE)
    write.csv(mmqtl_fdr.overlap, file = "eQTL - mmQTL (FDR).csv", row.names = FALSE)

 
## Versus Metabrain, from de Klein 2023
  ## Read in
    # this resource has a similar scope to mmqtl
    metabrain <- read_xlsx("../../../../PublicData/deKlein2023_eQTL/Supplementary Table 2 - Cis-eQTLs.xlsx", sheet = "Cortex-EUR", skip = 1)  
      
    # I note that the authors report 16619 hits a q < 0.05, so filter to this 
    metabrain$qval <- gsub("E", "e", metabrain$qval) %>% as.numeric() 
    metabrain <- metabrain[which(metabrain$qval < 0.05),] # all already are?
    
    # extract rs id
    metabrain$rsID <- splitter(metabrain$SNP, ":", 3) # I have confirmed that all have an associated rsID
    
    # colnames
    metabrain$Gene <- metabrain$GeneSymbol
    metabrain$Tissue <- "Brain"
    
  ## Overlap
    metabrain.overlap <- detailed.overlap(metabrain, score.col = c("MetaP", "MetaBeta", "qval", "QTLRank"), score.makecolname = c("P", "Beta", "Qval", "Rank"))
    
  ## Save
    write.csv(metabrain.overlap, file = "eQTL - Metabrain.csv", row.names = FALSE)
    
    
## Versus Metabrain, from de Klein 2023, but using their celltype-interaction eQTLs
  ## Read in
    # this resource has a similar scope to mmqtl
    deconQTL <- read_xlsx("../../../../PublicData/deKlein2023_eQTL/Supplementary Table 8 - Deconvoluted ieQTLs cis.xlsx", sheet = "Decon-eQTL Cortex EUR")  
      
    # annotate significance
    sigCols <- grep("BH-FDR", colnames(deconQTL))
    names(sigCols) <- colnames(deconQTL)[sigCols]
    names(sigCols) <- splitter(names(sigCols), " ", 1)
    deconQTL$Tissue <- apply(deconQTL[,sigCols], 1, function(x) {
      w <- which(as.numeric(x) < 0.05)
      paste0(names(sigCols)[w], collapse = "/")
    })
    
    
    # filter to those significant in Ast
    # deconQTL <- deconQTL[which(deconQTL$`Astrocyte BH-FDR` < 0.05),] # although both BH and permutation-based FDRs are present, the manuscript reports BH
    
    # extract rs id
    deconQTL$rsID <- splitter(deconQTL$SNP, ":", 3) # I have confirmed that all have an associated rsID
    
    # colnames
    deconQTL$Gene <- deconQTL$`Gene symbol`
    # deconQTL$Tissue <- "Brain"
    
  ## Overlap
    deconQTL.overlap <- detailed.overlap(deconQTL, score.col = c("Astrocyte BH-FDR"), score.makecolname = c("FDR_Ast"))
    
  ## Save
    write.csv(deconQTL, file = "eQTL - Metabrain (Deconvolved, Ast).csv", row.names = FALSE)
  

## snRNA-seq from Bryois et al 2021 
    bryois <- read_xlsx("../../../../PublicData/Bryois2021_eQTL_snRNAseq/41593_2022_1128_MOESM3_ESM.xlsx", sheet = "Table S3", skip = 3) # fine mapped eQTLs
    colnames(bryois) <- c("Tissue", "Gene", "EnsID", "rsID", "Chr", "Pos", "CaVEMaN", "CPP") 
    
  ## Overlap
    bryois.overlap <- detailed.overlap(bryois, score.col = "CPP")
    
  ## Save
    write.csv(bryois.overlap, file = "eQTL - Bryois snRNAseq.csv")
    
## Save all
  save(bryois.overlap,
       gtex_cav.overlap, gtex_sig.overlap,
       mmqtl_cav.overlap, mmqtl_fdr.overlap,
       metabrain.overlap, deconQTL.overlap,
       file = "eQTL - All Overlaps.rda")

## Annotate snp.annot object
  snp.annot <- check.snps(gtex.overlap$SNP, "eqtl_gtex")
  snp.annot <- check.snps(gtex.overlap$SNP[which(gtex.overlap$Consistent == "TRUE")], "eqtl_gtex_consistent")
  
  snp.annot <- check.snps(mmqtl.overlap$SNP, "eqtl_mmqtl")
  snp.annot <- check.snps(mmqtl.overlap$SNP[which(mmqtl.overlap$Consistent == "TRUE")], "eqtl_mmqtl_consistent")
  
  snp.annot <- check.snps(metabrain.overlap$SNP, "eqtl_metabrain")
  snp.annot <- check.snps(metabrain.overlap$SNP[which(metabrain.overlap$Consistent == "TRUE")], "eqtl_metabrain_consistent")
  
  snp.annot <- check.snps(bryois.overlap$SNP, "eqtl_bryois")
  snp.annot <- check.snps(bryois.overlap$SNP[which(bryois.overlap$Consistent == "TRUE")], "eqtl_bryois_consistent")
  
################################################################################################################################ #
## Analysis of eQTLs ----
 
    
   
## Plot 1: for each EGP, is it supported by a replicated eQTL?
  # this operates at a pooled level across datasets

  # get list of all rs-gene pairs
    fun1 <- function(w) {
      # convert rsID-gene pairs to enh-gene pairs
      x <- strsplit(w$LinkedEnh, "/")
      names(x) <- paste0(w$SNP.Association, "_")
      x <- do.call("c", x)
      
      x <- data.frame(Enh = x,
                      eQTL = splitter(names(x), "_", 1))
      x$Pair <- paste0(x$Enh, "_", x$eQTL)
      x$Expressed <- x$eQTL %in% names(meanExp_NHA)[which(meanExp_NHA > meanExp_NHA_thresh)]
      
      x <- unique(x)
    }
    
    rsPairs <- list(GTEx = fun1(gtex_sig.overlap),
                    mmQTL = fun1(mmqtl_fdr.overlap),
                    Metabrain = fun1(metabrain.overlap),
                    Bryois = fun1(bryois.overlap))
    rsPairs <- do.call("rbind", rsPairs)
    rsPairs$Dataset <- splitter(rownames(rsPairs), "\\.", 1)
  
  # categorise
    eQTL_pair_pooled <- res.final[which(res.final$HitPermissive), c("Pair", "Enh", "Gene")]
    eQTL_pair_pooled$OtherGenes_NotExpressed <- eQTL_pair_pooled$OtherGenes_Expressed <- eQTL_pair_pooled$eQTL_Category <- eQTL_pair_pooled$nRep <- eQTL_pair_pooled$Bryois <- eQTL_pair_pooled$GTEx <- eQTL_pair_pooled$mmQTL <- eQTL_pair_pooled$Metabrain <- NA

    for (j in 1:nrow(eQTL_pair_pooled)) {

      e <- eQTL_pair_pooled$Enh[j]
      g <- eQTL_pair_pooled$Gene[j]


      # check if there's an eqtl for this enhancer
      if (e %in% rsPairs$Enh) {

        # is the gene present as an eQTL?
        eqtl <- rsPairs[which(rsPairs$Enh == e),]

        # check in each study
        eQTL_pair_pooled$Metabrain[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "Metabrain")]
        eQTL_pair_pooled$mmQTL[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "mmQTL")]
        eQTL_pair_pooled$GTEx[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "GTEx")]
        eQTL_pair_pooled$Bryois[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "Bryois")]

        eQTL_pair_pooled$nRep[j] <- sum(eQTL_pair_pooled[j,c("Metabrain", "GTEx", "mmQTL", "Bryois")])

        if (eQTL_pair_pooled$nRep[j] >= 1) { # if the pair is supported by at least one study

          # eQTL_pair_pooled$eQTL_Category[j] <- paste0("Same gene (", eQTL_pair_pooled$nRep[j], ")")
          eQTL_pair_pooled$eQTL_Category[j] <- "Same gene"

          if (any(eqtl$eQTL != g)) { # other eQTLs for this enhancer, but targeting different genes
            
            sameButDiff <- eqtl[-which(eqtl$eQTL == g),]
            
            if (any(sameButDiff$Expressed)) {
              tab <- table(sameButDiff$eQTL[which(sameButDiff$Expressed)])
              exp <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
              eQTL_pair_pooled$OtherGenes_Expressed[j] <- exp
              
            }
            
            if (any(!sameButDiff$Expressed)) {
              tab <- table(sameButDiff$eQTL[which(!sameButDiff$Expressed)])
              not <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
              eQTL_pair_pooled$OtherGenes_NotExpressed[j] <- not
  
            }
            
          } # no else statement, they are left as NA

        } else { # if the pair is not found in any study...

          if (any(eqtl$Expressed)) {
            tab <- table(eqtl$eQTL[which(eqtl$Expressed)])
            exp <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
            eQTL_pair_pooled$OtherGenes_Expressed[j] <- exp

            m <- max(tab)

            if (m == 1) { # cannot be zero
              eQTL_pair_pooled$eQTL_Category[j] <- "Different gene (non-replicably)"
            } else {
              eQTL_pair_pooled$eQTL_Category[j] <- "Different gene (replicably)"
            }
          }

          if (any(!eqtl$Expressed)) {
            tab <- table(eqtl$eQTL[which(!eqtl$Expressed)])
            not <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
            eQTL_pair_pooled$OtherGenes_NotExpressed[j] <- not

            if(!(any(eqtl$Expressed))) {
              m <- max(tab)

              if (m == 1) { # cannot be zero
                eQTL_pair_pooled$eQTL_Category[j] <- "Different non-expressed gene (non-replicably)"
              } else {
                eQTL_pair_pooled$eQTL_Category[j] <- "Different non-expressed gene (replicably)"
              }
            }
          }

        }

      } else {

        eQTL_pair_pooled$eQTL_Category[j] <- "No eQTL"
        next

      }


    }
    
    write.csv(eQTL_pair_pooled, "eQTL Final - EGPs pooled.csv")

  ## Plot
    p <- eQTL_pair_pooled
    p$eQTL_Category <- gsub(" non-expressed", "", p$eQTL_Category) %>%
      gsub("1", "non-replicably", .) %>%
      gsub("2|3", "replicably", .)
    # pal <- c("grey95", pals$grn2orng[c(6:9, 3:1)])
    pal <- c("grey95", pals$grn2orng[c(7,8,2)])
  
  
    # p$eQTL_Category <- factor(p$eQTL_Category, levels = unique(p$eQTL_Category)[c(2,3,8,1,7,4,5,6)])
    p$eQTL_Category <- factor(p$eQTL_Category, levels = unique(p$eQTL_Category)[c(2,3,1,4)])
    pdf(file = "eQTL Final - EGPs Pooled.pdf", height = 2.2, width = 4)
    ggplot(p, aes(x = ".", fill = eQTL_Category)) +
      geom_bar(position = "stack", width = 0.8, colour = "black", linewidth = 0.5) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 158)) +
      labs(y = "Experimentally-derived EGPs") +
      scale_fill_manual(values = pal) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "right",
            legend.title = invis, axis.text.x = invis,
            axis.title.x = invis, panel.grid = invis, axis.ticks.x = invis)
    dev.off()
    
    
## Repeat the above in the fine-mapped data
    fmPairs <- list(GTEx = fun1(gtex_cav.overlap),
                    mmQTL = fun1(mmqtl_cav.overlap))
    fmPairs <- do.call("rbind", fmPairs)
    fmPairs$Dataset <- splitter(rownames(fmPairs), "\\.", 1)
  
  # categorise
    fm_res <- res.final[which(res.final$HitPermissive), c("Pair", "Enh", "Gene")]
    fm_res$OtherGenes_NotExpressed <- fm_res$OtherGenes_Expressed <- fm_res$eQTL_Category <- fm_res$nRep <-  fm_res$GTEx <- fm_res$mmQTL <- NA

    # redo the above to collect different gene in cases of same gene
    for (j in 1:nrow(fm_res)) {

      e <- fm_res$Enh[j]
      g <- fm_res$Gene[j]


      # check if there's an eqtl for this enhancer
      if (e %in% fmPairs$Enh) {

        # is the gene present as an eQTL?
        eqtl <- fmPairs[which(fmPairs$Enh == e),]

        # check in each study
        fm_res$mmQTL[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "mmQTL")]
        fm_res$GTEx[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "GTEx")]

        fm_res$nRep[j] <- sum(fm_res[j,c("GTEx", "mmQTL")])

        if (fm_res$nRep[j] >= 1) { # if the pair is supported by at least one study

          # fm_res$eQTL_Category[j] <- paste0("Same gene (", fm_res$nRep[j], ")")
          fm_res$eQTL_Category[j] <- "Same gene"

          if (any(eqtl$eQTL != g)) { # other eQTLs for this enhancer, but targeting different genes
            
            sameButDiff <- eqtl[-which(eqtl$eQTL == g),]
            
            if (any(sameButDiff$Expressed)) {
              tab <- table(sameButDiff$eQTL[which(sameButDiff$Expressed)])
              exp <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
              fm_res$OtherGenes_Expressed[j] <- exp
              
            }
            
            if (any(!sameButDiff$Expressed)) {
              tab <- table(sameButDiff$eQTL[which(!sameButDiff$Expressed)])
              not <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
              fm_res$OtherGenes_NotExpressed[j] <- not
  
            }
            
          } # no else statement, they are left as NA

        } else { # if the pair is not found in any study...

          if (any(eqtl$Expressed)) {
            tab <- table(eqtl$eQTL[which(eqtl$Expressed)])
            exp <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
            fm_res$OtherGenes_Expressed[j] <- exp

            m <- max(tab)

            if (m == 1) { # cannot be zero
              fm_res$eQTL_Category[j] <- "Different gene (non-replicably)"
            } else {
              fm_res$eQTL_Category[j] <- "Different gene (replicably)"
            }
          }

          if (any(!eqtl$Expressed)) {
            tab <- table(eqtl$eQTL[which(!eqtl$Expressed)])
            not <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
            fm_res$OtherGenes_NotExpressed[j] <- not

            if(!(any(eqtl$Expressed))) {
              m <- max(tab)

              if (m == 1) { # cannot be zero
                fm_res$eQTL_Category[j] <- "Different non-expressed gene (non-replicably)"
              } else {
                fm_res$eQTL_Category[j] <- "Different non-expressed gene (replicably)"
              }
            }
          }

        }

      } else {

        fm_res$eQTL_Category[j] <- "No eQTL"
        next

      }


    }
    
    
    write.csv(fm_res, "eQTL Final - EGPs pooled, Fine-mapped.csv")

  ## Plot
    p <- fm_res
    p$eQTL_Category <- gsub(" non-expressed", "", p$eQTL_Category) %>%
      gsub("1", "non-replicably", .) %>%
      gsub("2|3", "replicably", .)
    # pal <- c("grey95", pals$grn2orng[c(6:9, 3:1)])
    pal <- c("grey95", pals$grn2orng[c(7,8,2)])
  
  
    # p$eQTL_Category <- factor(p$eQTL_Category, levels = unique(p$eQTL_Category)[c(2,3,8,1,7,4,5,6)])
    p$eQTL_Category <- factor(p$eQTL_Category, levels = unique(p$eQTL_Category)[c(1,2,4,3)])
    pdf(file = "eQTL Final - EGPs Pooled.pdf", height = 2.2, width = 4)
    ggplot(p, aes(x = ".", fill = eQTL_Category)) +
      geom_bar(position = "stack", width = 0.8, colour = "black", linewidth = 0.5) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 158)) +
      labs(y = "Experimentally-derived EGPs") +
      scale_fill_manual(values = pal) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "right",
            legend.title = invis, axis.text.x = invis,
            axis.title.x = invis, panel.grid = invis, axis.ticks.x = invis)
    dev.off()
    
## Repeat the above in the two astrocyte-specific datasets
    astPairs <- list(Metabrain = fun1(deconQTL.overlap[grep("Astrocyte", deconQTL.overlap$SNP.Tissue),]),
                     Bryois = fun1(bryois.overlap[which(bryois.overlap$SNP.Tissue == "Astrocytes"),]))
    
    astPairs <- do.call("rbind", astPairs)
    astPairs$Dataset <- splitter(rownames(astPairs), "\\.", 1)
  
  # categorise
    ast_res <- res.final[which(res.final$HitPermissive), c("Pair", "Enh", "Gene")]
    ast_res$OtherGenes_NotExpressed <- ast_res$OtherGenes_Expressed <- ast_res$eQTL_Category <- ast_res$nRep <-  ast_res$Bryois <- ast_res$Metabrain <- NA

    # redo the above to collect different gene in cases of same gene
    for (j in 1:nrow(ast_res)) {

      e <- ast_res$Enh[j]
      g <- ast_res$Gene[j]


      # check if there's an eqtl for this enhancer
      if (e %in% astPairs$Enh) {

        # is the gene present as an eQTL?
        eqtl <- astPairs[which(astPairs$Enh == e),]

        # check in each study
        ast_res$Bryois[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "Bryois")]
        ast_res$Metabrain[j] <- g %in% eqtl$eQTL[which(eqtl$Dataset == "Metabrain")]

        ast_res$nRep[j] <- sum(ast_res[j,c("Metabrain", "Bryois")])

        if (ast_res$nRep[j] >= 1) { # if the pair is supported by at least one study

          # ast_res$eQTL_Category[j] <- paste0("Same gene (", ast_res$nRep[j], ")")
          ast_res$eQTL_Category[j] <- "Same gene"

          if (any(eqtl$eQTL != g)) { # other eQTLs for this enhancer, but targeting different genes
            
            sameButDiff <- eqtl[-which(eqtl$eQTL == g),]
            
            if (any(sameButDiff$Expressed)) {
              tab <- table(sameButDiff$eQTL[which(sameButDiff$Expressed)])
              exp <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
              ast_res$OtherGenes_Expressed[j] <- exp
              
            }
            
            if (any(!sameButDiff$Expressed)) {
              tab <- table(sameButDiff$eQTL[which(!sameButDiff$Expressed)])
              not <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
              ast_res$OtherGenes_NotExpressed[j] <- not
  
            }
            
          } # no else statement, they are left as NA

        } else { # if the pair is not found in any study...

          if (any(eqtl$Expressed)) {
            tab <- table(eqtl$eQTL[which(eqtl$Expressed)])
            exp <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
            ast_res$OtherGenes_Expressed[j] <- exp

            m <- max(tab)

            if (m == 1) { # cannot be zero
              ast_res$eQTL_Category[j] <- "Different gene (non-replicably)"
            } else {
              ast_res$eQTL_Category[j] <- "Different gene (replicably)"
            }
          }

          if (any(!eqtl$Expressed)) {
            tab <- table(eqtl$eQTL[which(!eqtl$Expressed)])
            not <- paste(names(tab), "(", tab, ")", sep = "", collapse = "|")
            ast_res$OtherGenes_NotExpressed[j] <- not

            if(!(any(eqtl$Expressed))) {
              m <- max(tab)

              if (m == 1) { # cannot be zero
                ast_res$eQTL_Category[j] <- "Different non-expressed gene (non-replicably)"
              } else {
                ast_res$eQTL_Category[j] <- "Different non-expressed gene (replicably)"
              }
            }
          }

        }

      } else {

        ast_res$eQTL_Category[j] <- "No eQTL"
        next

      }


    }
    
    
    write.csv(ast_res, "eQTL Final - EGPs pooled, Ast-specific.csv")

  ## Plot
    p <- ast_res
    p$eQTL_Category <- gsub(" non-expressed", "", p$eQTL_Category) %>%
      gsub("1", "non-replicably", .) %>%
      gsub("2|3", "replicably", .)
    # pal <- c("grey95", pals$grn2orng[c(6:9, 3:1)])
    pal <- c("grey95", pals$grn2orng[c(7,8,2)])
  
  
    # p$eQTL_Category <- factor(p$eQTL_Category, levels = unique(p$eQTL_Category)[c(2,3,8,1,7,4,5,6)])
    p$eQTL_Category <- factor(p$eQTL_Category, levels = unique(p$eQTL_Category)[c(1,2,4,3)])
    pdf(file = "eQTL Final - EGPs Pooled.pdf", height = 2.2, width = 4)
    ggplot(p, aes(x = ".", fill = eQTL_Category)) +
      geom_bar(position = "stack", width = 0.8, colour = "black", linewidth = 0.5) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 158)) +
      labs(y = "Experimentally-derived EGPs") +
      scale_fill_manual(values = pal) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "right",
            legend.title = invis, axis.text.x = invis,
            axis.title.x = invis, panel.grid = invis, axis.ticks.x = invis)
    dev.off()
  
## Plot 2: examine how  eQTL pairs within each dataset are supportive of the screen.
     rsPairs2 <- list(GTEx = fun1(gtex_sig.overlap),
                    mmQTL = fun1(mmqtl_fdr.overlap),
                    Metabrain = fun1(metabrain.overlap),
                    Bryois = fun1(bryois.overlap),
                    AstBryois = fun1(bryois.overlap[which(bryois.overlap$SNP.Tissue == "Astrocytes"),]),
                    AstMetabrain = fun1(deconQTL.overlap[grep("Astrocyte", deconQTL.overlap$SNP.Tissue),]),
                    GTEx_FM = fun1(gtex_cav.overlap),
                    mmQTL_FM = fun1(mmqtl_cav.overlap))
    rsPairs2 <- do.call("rbind", rsPairs2)
    rsPairs2$Dataset <- splitter(rownames(rsPairs2), "\\.", 1)
  
    
  # categorise
    u <- unique(rsPairs2$Dataset)
    eQTL_pair_across <- res.final[which(res.final$HitPermissive), c("Pair", "Enh", "Gene")]
    eQTL_pair_across[,u] <- NA

    for (j in 1:nrow(eQTL_pair_across)) {
      
      e <- eQTL_pair_across$Enh[j]
      g <- eQTL_pair_across$Gene[j]
      
      # check if there's an eqtl for this enhancer
      if (e %in% rsPairs2$Enh) {
        
        # get eqtls for this gene
        eqtl <- rsPairs2[which(rsPairs2$Enh == e),]
        
        # are the eqtls replicable? that is, the gene is found in at least two datasets
        tab <- table(eqtl$eQTL)
        rep <- names(tab)[which(tab > 1)]
        eqtl$Rep <- eqtl$eQTL %in% rep
        
        # for each study
        for (k in u) {
          
          # check if the same gene is identified
          same <- g %in% eqtl$eQTL[which(eqtl$Dataset == k)]
          
          if (same) {
          
            eQTL_pair_across[j,k] <- "Same gene" 
            
          } else { # if not the same gene, are there other genes and are they replicably identified?
            
            if (k %in% eqtl$Dataset) { # if there is at least one eqtl for the dataset for this enhancer
              
              l <- eqtl[which(eqtl$Dataset == k),]
              
              if (any(l$Rep)) {
                
                eQTL_pair_across[j,k] <- "Different gene (replicable)"  
                
              } else {
                
                eQTL_pair_across[j,k] <- "Different gene (non-replicable)"  
                
              }
              
            } else { # no eqtls for this study for this enhancer
              
              eQTL_pair_across[j,k] <- "No eQTL"  
              
            }
          }
        }
        
        
      } else {
        
        eQTL_pair_across[j,u] <- "No eQTL" # if the enhancer has no eqtls in any study
        
        
      }
      
      
      
      
    }
    
  ## Save
    write.csv(eQTL_pair_across, file = "eQTL Final - EGPs across studies.csv")
    
  ## Plot
      p <- eQTL_pair_across[,c("Pair", u)]
      p <- melt(p, id.vars = "Pair")
      colnames(p) <- c("Pair", "Study", "Cat")
      # p$Study <- factor(p$Study, levels = unique(p$Study)[c(1,2,3,4,6,5)])
      p$Astro <- grepl("Ast", p$Study) %>% factor()
      p$Finemapped <- grepl("FM", p$Study) %>% factor()
      levels(p$Astro) <- c("Brain eQTL", "eQTL within astrocytes")
      levels(p$Finemapped) <- c("Significant set", "Fine-mapped set")
      levels(p$Study) <- gsub("Ast|_FM", "", levels(p$Study))
      
      p$Cat <- splitter(p$Cat, " \\(", 1)
      p$Cat <- factor(p$Cat, levels = unique(p$Cat)[1:3])
      pal <- c("grey95", pals$grn2orng[c(7,2)])
      p$Group <- paste0(p$Finemapped, "\n", p$Astro)
      
      pdf(file = "eQTL Final - EGPs across studies.pdf", height = 3, width = 7)
      ggplot(p, aes(x = Study, fill = Cat)) +
        geom_bar(position = "stack", width = 0.8, linewidth = 0.5) +
        theme_bw() +
        scale_y_continuous(expand = c(0,0)) +
        facet_grid(.~Group, scales = "free", space = "free") +
        labs(y = "Experimentally-validated EGPs", x = "eQTL Resource") +
        scale_fill_manual(values = pal) +
        scale_alpha_manual(values = c(0.85, 0.85, 1)) +
        theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none",
              legend.title = invis, axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              axis.title.x = invis, panel.grid = invis, axis.ticks.x = invis)
      dev.off()
      
      
  
################################################################################################################################ #
## Output ---- 

  
## Output as is
  write.csv(snp.annot, file = "Final - SNP Annotation.csv", row.names = FALSE)
  
  
## Wrangle to enhancer-level scores
  # function
  read.detailedOverlap <- function(x) {
    
    score.col <- which(grepl("Score|^P$", colnames(x)))
    
    y <- apply(x, 1, function(y) {
      
      if (grepl("/", y[2])) { # if the SNP is associated with many enhanceres
        z <- data.frame(Enh = do.call("c", strsplit(y[2], "/")),
                        Hit = NA,
                        SNP = y[1],
                        Resource = NA,
                        Category = do.call("c", strsplit(y[3], "/")),
                        SNP.Association = y[4],
                        Score = y[score.col])
      } else {
        z <- data.frame(Enh = y[2],
                        Hit = NA,
                        SNP = y[1],
                        Resource = NA,
                        Category = y[3],
                        SNP.Association = y[4],
                        Score = y[score.col])
      }
      
     
      
      return(z)
    })
    
    y <- do.call("rbind", y)

    
    # output
    return(y)
  }
  

  # apply
  final_Enh2SNP <- list(GWAS_Catalogue = read.detailedOverlap(gwascat.overlap),
                        Disgenet = read.detailedOverlap(disg.overlap))
  
  final_Enh2SNP <- do.call("rbind", final_Enh2SNP)
  
  # add GWAS SNP information
  cols <- colnames(snp.annot)[grep("GWAS_", colnames(snp.annot))]
  for (j in cols) {
    if (any(snp.annot[,j])) {
      w <- which(snp.annot[,j])
      x <- snp.annot[w,]
      y <- data.frame(Enh = x$Enh,
                      Hit = NA,
                      Resource = NA,
                      SNP = x$SNP,
                      Category = x$Category,
                      SNP.Association = splitter(j, "_", 2),
                       Score = NA,
                      row.names = paste0("GWAS.", 1:nrow(x)))
      final_Enh2SNP <- rbind(final_Enh2SNP, y)
    }
    
  }
  
  # clean
  ord <- sub("Enh", "", final_Enh2SNP$Enh) %>% as.numeric() %>% rank() %>% order()
  final_Enh2SNP <- final_Enh2SNP[ord,]
  final_Enh2SNP$Hit <- final_Enh2SNP$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
  final_Enh2SNP$Resource <- splitter(rownames(final_Enh2SNP), "\\.", 1)
  final_Enh2SNP <- unique(final_Enh2SNP)
  # final_Enh2SNP <- relocate(final_Enh2SNP, c("Enh", "Hit", "Resource", "SNP"))
  
  # save
  write.csv(final_Enh2SNP, file = "Final - SNPs of Interest.csv", row.names = FALSE)

  
    
      