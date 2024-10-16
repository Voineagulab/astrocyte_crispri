## This script explores enhancer transcription using in house TTseq data and public FANTOM5 CAGE datasets

################################################################################################################################ #
## Setup ----

rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/4_EnhancerTranscription/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Rsubread)


## Load
  source("../../Scripts/Functions.R")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  atac <- read.csv("../../../FullLibrary_Selection/Results/Peaks_Annotated.csv")

  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  sig.colours <- c("black", "firebrick1")

## Results dataframe
  res.final <- read.csv("../2_DE/Enh/Results Final.csv")
  hit.enh <- res.final$Enh[which(res.final$HitPermissive)]
  hit.pairs <- res.final[which(res.final$HitPermissive),] # significant enhancer-gene pairs

  
## Main paths
  nha_dir_38 <- "../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  nha_dir_19 <- "../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed" 
  
## Functions
  # read in a bed file
  read.bed <- function(dir) read.delim(dir, header = FALSE)
  
  # convert coordinate id to enhancer id
  convert.id <- function(x) {
    m <- match(x, guides$TargetCoord)
    y <- guides$TargetID[m]
    return(y)
  }
  
  # bedtools intersect, with maximal information retained
  wawb <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-wa", "-wb", 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
  
  # intersect with a “left outer join”. that is, for each feature in A report each overlap with B. if no overlaps are found, report a NULL feature for B
  # useful for binary calls
  loj <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-loj", 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
  
  # enrichment
  run.eRNAFisher <- function(data, data.column, rbind = FALSE, rbind.to = eRNA.enrich) {
    # filter to enhancers tested against at least one gene (957 of 979)
    if ("Tested" %in% colnames(data)) {
      data <- data[which(data$Tested),]  
    }
    
    
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
      return(rbind(rbind.to, out))
    } else {
      return(out)
    }

  }

  
################################################################################################################################ #
## Start by exploring enhancer expression in FANTOM5 ----
  
## A CAGE resource of enhancer transcription across the human body
  
## Read in 
  ## Sample list
    fantom.samps <- read.delim("../../../PublicData/FANTOM5/Human.sample_name2library_id.txt", header = FALSE) 
    colnames(fantom.samps) <- c("Name", "ID")
    fantom.sampsAst <- fantom.samps$ID[grep("Astrocyte", fantom.samps$Name)] # astrocytes
    fantom.sampsGBM <- fantom.samps$ID[grep("astrocytoma|glioblastoma|glioma", fantom.samps$Name)] # astrocytes
    
  ## Enhancer binary calls of usage across celltypes
    fantom.usage <- read.delim("../../../PublicData/FANTOM5/F5.hg38.enhancers.expression.usage.matrix") # one row per enhancer id, one column per sample 
    fantom.usage <- fantom.usage[,c(fantom.sampsAst, fantom.sampsGBM)]
    # fantom.usage <- rowSums(fantom.usage) # number of astrocyte samples in which it is expressed (out of 3)
  
    
  ## Enhancer coordinates
    fantom.coord <- read.bed("../../../PublicData/FANTOM5/F5.hg38.enhancers.bed") # bed12 format
    fantom.coord <- fantom.coord[,1:4] # bed 4 format
    fantom.coord$UsedAst <- rowSums(fantom.usage[,fantom.sampsAst])
    fantom.coord$UsedGBM <- rowSums(fantom.usage[,fantom.sampsGBM])
    
    fantom.dir <- "../../../PublicData/FANTOM5/F5.hg38.enhancers.ast.bed"
    write.bed(fantom.coord, fantom.dir)
    
## Compare to our peaks
  ## Overlap
    loj(nha_dir_38, fantom.dir, "FANTOM5/ScreenIntersect.bed")
    fantom.overlap <- read.bed("FANTOM5/ScreenIntersect.bed")
    fantom.overlap <- data.frame(Enh = convert.id(fantom.overlap$V4),
                                 Hit = NA,
                                 FANTOM5 = fantom.overlap$V6 != -1,
                                 UsedAst = fantom.overlap$V9,
                                 UsedGBM = fantom.overlap$V10)
    fantom.overlap$Hit <- fantom.overlap$Enh %in% hit.enh 
    
    g <- grep("Used", colnames(fantom.overlap))
    fantom.overlap[,g] <- apply(fantom.overlap[,g], 2, function(x) {
      x[which(x == ".")] <- 0
      x <- as.numeric(x)
      return(x)
    })
   
  ## There are 20 screen peaks that overlap >1 (i.e., 2) F5 enhancers
    # collect the max for these
    x <- list()  
    
    for (j in fantom.overlap$Enh) {
      
      if (length(which(fantom.overlap$Enh == j)) == 1) { # only one overlap
        
        x[[j]] <- fantom.overlap[which(fantom.overlap$Enh == j),] 
        
      } else { # multiple overlaps, collect max
        
        y <- fantom.overlap[which(fantom.overlap$Enh == j),]
        
        for (k in g) y[,k] <- max(y[,k]) # g is defined before the loop, as columns containing "Used"
          
        y <- unique(y)
        x[[j]] <- y
        
      }
      
    }
    
    fantom.overlap <- do.call("rbind", x)
      
  ## Save
    write.csv(fantom.overlap, "FANTOM5/Peak Annotation.csv", row.names = FALSE)
  
    
## Statistics
  ## Test
    table(fantom.overlap$Hit, fantom.overlap$FANTOM5) %>% fisher.test() # 1.7, p=2e-3
    table(fantom.overlap$Hit, fantom.overlap$UsedAst > 0) %>% fisher.test() # 2.11, p=5e-4
    table(fantom.overlap$Hit, fantom.overlap$UsedGBM > 0) %>% fisher.test() # 2.3, p=1e-3
    
    fantom.overlap$FANTOM5_Any <- fantom.overlap$FANTOM5
    fantom.overlap$FANTOM5_Ast <- fantom.overlap$UsedAst > 0
    fantom.overlap$FANTOM5_GBM <- fantom.overlap$UsedGBM > 0
    
    eRNA.enrich <- run.eRNAFisher(data = fantom.overlap, data.column = "FANTOM5_Any", rbind = FALSE)
    eRNA.enrich <- run.eRNAFisher(data = fantom.overlap, data.column = "FANTOM5_Ast", rbind = TRUE)
    eRNA.enrich <- run.eRNAFisher(data = fantom.overlap, data.column = "FANTOM5_GBM", rbind = TRUE)
    
 
  
################################################################################################################################ #
## Explore TTseq data: quantify expression at enhancer using feature counts ----
  
## Prepare input files
  # for feature counts, the peak bed file needs to be converted to a .saf
  create_safs <- function (directory, convertid = T) {
    bed <- read.bed(directory)
    colnames(bed) <- c("chr", "start", "end", "id")
    if (convertid == T) {
      bed$id <- convert.id(bed$id)
    }
    bed <- data.frame(GeneID = bed$id,
                        Chr = bed$chr, 
                        Start = bed$start,
                        End = bed$end,
                        Strand = ".")
    
    # apparently remove chr so it works with the reference genome used for generating the bam files, peaks$Chr=gsub("chr", "", peaks$Chr)
    
    # add flanking window
    window <- 1000 # total window
    midpoint <- round(((bed$End + bed$Start) / 2))
    bed$Start <- midpoint  - (window / 2) # take midpoint, and centre the window around it
    bed$End <- midpoint + (window / 2)
    
    saf <- list()
    saf$Unstranded <- bed
    saf$Pos <- bed; saf$Pos$Strand = "+"
    saf$Neg <- bed; saf$Neg$Strand = "-"
    return(saf)
  }
  
  saf  <- create_safs(nha_dir_38)
  saf_intergenic  <- create_safs("/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Intergenic_Peaks_chr.bed", convertid = F) #/Volumes/share
  
## Run feature counts
  ## Setup functions and directories
    bam_tt <- "/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/TTseq/NFcore_RNAseq_Basic/star_salmon/GOK10844A1.markdup.sorted.bam"
    bam_rnaseq <- "/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/TTseq/NFcore_RNAseq_Basic/star_salmon/GOK10856A1.markdup.sorted.bam"
      
    run_fc <- function(saf, saf_dir, bam_tt, bam_rnaseq, stranded = 2) {
      options(scipen = 10) #This stops numbers from being printed as scientific notation
      write.table(saf, saf_dir, sep = "\t", row.names = FALSE, quote = FALSE, )
      options(scipen=0)
      featureCounts(c(bam_tt, bam_rnaseq), 
                    annot.ext = saf_dir, isGTFAnnotationFile = FALSE, # a SAF of filtered bins
                    useMetaFeatures = FALSE, # feature-level, rather than pooling features-to-genes
                    allowMultiOverlap = TRUE, # a change versus Irina
                    isPairedEnd = TRUE, # self-explanatory
                    nthreads = 8, # self-explanatory 
                    strandSpecific = stranded, 
                    minMQS = 10, # minimum quality score at either end of the pair
                    checkFragLength = FALSE, 
                    countMultiMappingReads = FALSE, 
                    requireBothEndsMapped = FALSE,                              
                    countChimericFragments = FALSE)
      
      # the change in this versus Irina is that allowMultOverlap = TRUE. Because some of our features are overlapping, this is useful
      # further, because countMultiMappingReads = FALSE, we know they aren't multimappers
    }
    
  ## Run
    counts <- list()
    counts$Unstranded <- run_fc(saf$Unstranded, "TTseq/SAF_Unstranded_w1000.saf", 
                                bam_tt = bam_tt, bam_rnaseq = bam_rnaseq, stranded = 0)
    counts$Pos <- run_fc(saf$Pos, "TTseq/SAF_Pos_w1000.saf",
                         bam_tt = bam_tt, bam_rnaseq = bam_rnaseq,stranded = 2)
    counts$Neg <- run_fc(saf$Neg, "TTseq/SAF_Neg_w1000.saf",
                         bam_tt = bam_tt, bam_rnaseq = bam_rnaseq,stranded = 2)
    save(counts, file = "TTseq/FeatureCounts.rda")
    
    counts_intergenic <- list()
    counts_intergenic$Unstranded <- run_fc(saf_intergenic$Unstranded, 
                                           "TTseq/SAF_Unstranded_w1000_intergenic.saf",
                                           bam_tt = bam_tt, bam_rnaseq = bam_rnaseq,
                                           stranded = 0)
    counts_intergenic$Pos <- run_fc(saf_intergenic$Pos, "TTseq/SAF_Pos_w1000_intergenic.saf",
                                    bam_tt = bam_tt, bam_rnaseq = bam_rnaseq, stranded = 2)
    counts_intergenic$Neg <- run_fc(saf_intergenic$Neg, "TTseq/SAF_Neg_w1000_intergenic.saf",
                                    bam_tt = bam_tt, bam_rnaseq = bam_rnaseq,stranded = 2)
    save(counts_intergenic, file = "TTseq/FeatureCounts_intergenic.rda")


  ## Explore stats
    seqStats <- counts$Unstranded$stat # same for all levels, so use unstranded
    colnames(seqStats) <- c("Status", "TTseq", "RNAseq")
    libSize <- colSums(seqStats[,-1]) / 10 ^ 6
    
  
## Process
process_TTseq <- function (counts, libSize, min.exp.thresh = 3) {
    x <- lapply(counts, function(y) {
      z <- as.data.frame(y$counts)
      colnames(z) <- c("TTseq", "RNAseq")
      return(z)
    })
      
    x <- do.call("cbind", x)
  ## Calculations
    # note: because the total library size in RNA-seq and TT-seq is similar (66M and 63M, respectively)
    # we shall use counts rather than adjust for this factor
    
    ## Formatting
      x$Hit <- rownames(x) %in% hit.enh
      colnames(x) <- gsub("Unstranded", "Total", colnames(x))
      colnames(x) <- paste0(splitter(colnames(x), "\\.", 2), "_", splitter(colnames(x), "\\.", 1))
      colnames(x) <- gsub("NA_", "", colnames(x))
      x <- x[,c(7,3,5,1,4,6,2)]
      
    ## Libsize normalisation
      x$TPM_TTseq_Total <- x$TTseq_Total / libSize["TTseq"] # note that no need to correct for length, as the size is 1kb
      x$TPM_TTseq_Pos <- x$TTseq_Pos / libSize["TTseq"]
      x$TPM_TTseq_Neg <- x$TTseq_Neg / libSize["TTseq"]
      x$TPM_RNAseq_Total <- x$RNAseq_Total / libSize["RNAseq"] 
      x$TPM_RNAseq_Pos <- x$RNAseq_Pos / libSize["RNAseq"] 
      x$TPM_RNAseq_Neg <- x$RNAseq_Neg / libSize["RNAseq"] 
   
    ## Comparisons across libraries
      # strand bias
      x$FractionPos_TTseq <- x$TTseq_Pos / (x$TTseq_Total)
      x$FractionPos_RNAseq <- x$RNAseq_Pos / (x$RNAseq_Total)
      x$FractionPos_TTseq[which(is.na(x$FractionPos_TTseq))] <- 0 # this and line below: as NA was returned when both strands 0, set to 0
      x$FractionPos_RNAseq[which(is.na(x$FractionPos_RNAseq))] <- 0
      
      # TTseq to RNAseq ratio
      x$Ratio_TTversusRNA <- (x$TTseq_Total + 0.01) / (x$RNAseq_Total + 0.01)
      x$PassThresh_TT <- x$TTseq_Total >= min.exp.thresh
      x$PassThresh_RNA <- x$RNAseq_Total >= min.exp.thresh
      x$TT_Enrich <- (x$Ratio_TTversusRNA > 1) & (x$PassThresh_TT)
  
    ## Save
      return(x)
    }
    
    featCounts <- process_TTseq(counts, libSize)
    write.csv(featCounts, file = "TTseq/Results Table.csv")
    
    seqStats_intergenic <- counts_intergenic$Unstranded$stat
    colnames(seqStats_intergenic) <- c("Status", "TTseq", "RNAseq")
    libSize_intergenic <- colSums(seqStats_intergenic[,-1]) / 10 ^ 6
    featCounts_intergenic <- process_TTseq(counts_intergenic, libSize_intergenic)
    write.csv(featCounts_intergenic, file = "TTseq/Results_Table_intergenic.csv")
    
    
    #Repeat these analyses for K562 ENCODE_rE2G peaks 
    #@Sam bagot
    k562_bam_tt <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562TTseq/TTseq_R1.markdup.sorted.bam"
    k562_bam_rnaseq <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562TTseq/RNAseq_R1.markdup.sorted.bam"
    saf_K562_ENCODE <- create_safs("/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/ENCODE_Candidates_chr.bed", convertid = F)
    counts_K562_ENCODE <- list()
    counts_K562_ENCODE$Unstranded <- run_fc(saf_K562_ENCODE$Unstranded, 
                                           "K562TTseq/SAF_Unstranded_w1000_K562_ENCODE.saf",
                                           bam_tt = k562_bam_tt, bam_rnaseq = k562_bam_rnaseq,
                                           stranded = 0)
    counts_K562_ENCODE$Pos <- run_fc(saf_K562_ENCODE$Pos, "K562TTseq/SAF_Pos_w1000_K562_ENCODE.saf",
                                    bam_tt = k562_bam_tt, bam_rnaseq = k562_bam_rnaseq, stranded = 2)
    counts_K562_ENCODE$Neg <- run_fc(saf_K562_ENCODE$Neg, "K562TTseq/SAF_Neg_w1000_K562_ENCODE.saf",
                                    bam_tt = k562_bam_tt, bam_rnaseq = k562_bam_rnaseq,stranded = 2)
    save(counts_K562_ENCODE, file = "K562TTseq/FeatureCounts_K562_ENCODE.rda")
    seqStats_K562_ENCODE <- counts_K562_ENCODE$Unstranded$stat
    colnames(seqStats_K562_ENCODE) <- c("Status", "TTseq", "RNAseq")
    libSize_K562_ENCODE <- colSums(seqStats_K562_ENCODE[,-1]) / 10 ^ 6
    featCounts_K562_ENCODE <- process_TTseq(counts_K562_ENCODE, libSize_K562_ENCODE)
    write.csv(featCounts_K562_ENCODE, file = "K562TTseq/Results_Table_K562_ENCODE.csv")
    
  
################################################################################################################################ #
## Binary classification of peaks to be expressing eRNA, based on TTseq and RNAseq data (not CAGE) ----  
    
## Classify
  transcribed <- featCounts[,c("Hit", "TTseq_Pos", "TTseq_Neg" )]
  colnames(transcribed) <- c("Hit", "Pos", "Neg")
  transcribed$NascentEnriched <- featCounts$TTseq_Total > featCounts$RNAseq_Total
  transcribed$Enh <- rownames(transcribed)
  transcribed <- relocate(transcribed, "Enh")
  
  
  threshes <- c(2, 3, 5, 10)
  
  for (use.thresh in threshes) {
    u <- paste0("Unidirectional_", use.thresh) 
    b <- paste0("Bidirectional_", use.thresh)
    c <- paste0("Category_", use.thresh)
    
    transcribed[,u] <- (rowSums(transcribed[,c("Pos", "Neg")] >= use.thresh) >= 1) & transcribed$NascentEnriched
    transcribed[,b] <- (rowSums(transcribed[,c("Pos", "Neg")] >= use.thresh) == 2) & transcribed$NascentEnriched
    
    transcribed[,c] <- "Not transcribed"
    transcribed[which(transcribed[,u]),c] <- "Unidirectional"
    transcribed[which(transcribed[,b]),c] <- "Bidirectional"
    transcribed[,c] <- factor(transcribed[,c], levels = c("Not transcribed", "Unidirectional", "Bidirectional"))
  }
  
## Save
  write.csv(transcribed, file = "TTseq/Transcriptional classification.csv", row.names = FALSE)
    
## Enrichments
  x <- transcribed
  x$TTseq_Unidirectional <- x$Unidirectional_3
  x$TTseq_Bidirectional <- x$Bidirectional_3
  x$TTseq_Transcribed <- x$TTseq_Unidirectional | x$TTseq_Bidirectional
  x <- x[which(x$Enh %in% res.final$Enh),]
  
  eRNA.enrich <- run.eRNAFisher(data = x, data.column = "TTseq_Unidirectional", rbind = TRUE)
  eRNA.enrich <- run.eRNAFisher(data = x, data.column = "TTseq_Bidirectional", rbind = TRUE)
  eRNA.enrich <- run.eRNAFisher(data = x, data.column = "TTseq_Transcribed", rbind = TRUE)
  
  write.csv(eRNA.enrich, file = "Hit Enrichments.csv")
  
  
################################################################################################################################ #
## Comparison of CAGE to TTseq ----
 
     
## Setup
  x <- data.frame(Enh = rownames(featCounts), # same order confirmed with: match(rownames(featCounts), rownames(fantom.overlap))
                  TPM = featCounts$TPM_TTseq_Total,
                  Ratio = featCounts$Ratio_TTversusRNA,
                  Enrich = featCounts$TT_Enrich,
                  FANTOM5 = fantom.overlap$FANTOM5,
                  FANTOM5_Ast = fantom.overlap$UsedAst)
  
## Stats
  # TTSeq enriched peaks versus FANTOM5
  table(x$Enrich, x$FANTOM5) %>% fisher.test()  # OR = 2.32, p = 3e-10
  table(x$Enrich, x$FANTOM5_Ast > 0) %>% fisher.test() # OR = 4.01, p = 9e-15
  
  # TTseq above threshold versus FANTOM5
  table(x$TPM > 0.05, x$FANTOM5) %>% fisher.test() # OR = 2.18, p = 5e-9  
  table(x$TPM > 0.05, x$FANTOM5_Ast > 0) %>% fisher.test()  # OR = 4.29, p = 5e-16
  
  # FANTOM5 expression 
  wilcox.test(x$TPM ~ x$FANTOM5) # p=1e-10

  
  
  
  