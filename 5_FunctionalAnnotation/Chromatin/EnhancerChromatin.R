## This script annotates hit enhancers by overlap to "interesting" coordinates from various resources, including:

  # enhancer annotations across cell-types

  # regions of different accessibility in various brain disorders

 
################################################################################################################################ #
## Setup ----


## Generic
  rm(list = ls()); gc()
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/")
  options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Rsamtools)
  library(rtracklayer)
  
  source("../../../Scripts/Functions.R")

## Load
  # load("../../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  res.final <- read.csv("../../2_DE/Enh/Results Final.csv")  
  atac <- read.csv("../../../../FullLibrary_Selection/Results/Peaks_Annotated.csv")

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  sig.colours <- c("black", "firebrick1")

  
## Set up enhancer lists
  # a dataframe to store enhancer-level logical annotations for those in our screen
  hit.enh <- unique(res.final$Enh[which(res.final$HitPermissive)]) # hit enhancers
  bg <- unique(guides$TargetID[which(guides$TargetCat == "Enh")]) # all enhancers
  bg.tested <- unique(res.final$Enh) # some enhancers were not tested in the screen due to filtering of genes
  
  candidate.annot <- data.frame(Enh = bg, 
                                Tested = bg %in% bg.tested, 
                                Hit = bg %in% hit.enh, 
                                row.names = bg)
  
  s <- sub("Enh", "", bg) %>% as.numeric() %>% rank() %>% order()
  candidate.annot <- candidate.annot[s,]
  m <- match(candidate.annot$Enh, guides$TargetID)
  candidate.annot$Coord <- guides$TargetCoord[m]
 
  
  # a dataframe to store fisher test results for overenrichments
  # (actually, a function to fill initialise/fill this)
  run.EnhancerFisher <- function(data = candidate.annot, data.column, rbind = FALSE, rbind.to = candidate.enrich) {
    # filter to enhancers tested against at least one gene (957 of 979)
    data <- data[which(data$Tested),]
    
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
    
  
## Functions
  # intersect with a “left outer join”. that is, for each feature in A report each overlap with B. if no overlaps are found, report a NULL feature for B
  # useful for binary calls
  loj <- function(a, b, out, fun = "-loj") {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                fun, 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
  

  
  # intersect with full report
  wawb <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-wa", "-wb", 
                ">", out)

    system(call, intern = FALSE, wait = TRUE)    
    print("Complete!")
  }
 
  # read in a bed file
  read.bed <- function(dir) read.delim(dir, header = FALSE)
  
  # convert coordinate id to enhancer id
  convert.id <- function(x) {
    m <- match(x, guides$TargetCoord)
    y <- guides$TargetID[m]
    return(y)
  }
  
  # calculate coverage from bigwig
  read_bw <- function(bw, importedBed, ids = "standard") {

    # get coverage
    y <- import(bw, selection = importedBed, as = "NumericList")

    # get ids
    if (ids == "get") {
      i <- as.data.frame(y@metadata$ranges)
      i <- paste0(i$seqnames,
                    ":",
                    i$start-1, # the -1 is to account for the indexing differences
                    "-",
                    i$end)
    }

    if (ids == "standard") i <- importedBed$name

    # return summary
    z <- data.frame(id = i,
                    Mean = sapply(y, mean),
                    Max = sapply(y, max))
    # y <- sapply(y, mean) # gets the mean coverage
    # names(y) <- ids # $name is the GRCh38 coordinate!
    return(z)
  }
  
  
## Commonly used paths
  # candidate coordinates
  nha_dir_38 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed" 
  nha_dir_19 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed" 
  nha_dir_19_w1000 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19_window1000.bed"
  nha_dir_38_w1000 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg38_window1000.bed"
  
## Create a windowed version of the NHA peaks
  ## Function
    create_window <- function(bed_in, bed_out, window = 1000) {
      x <- read.delim(bed_in, header = FALSE)
      colnames(x)[1:4] <- c("chr", "start", "end", "id")
      
      # extend from midpoint
      midpoints <- round(((x$end + x$start) / 2))
      x$start <- midpoints  - (window / 2)
      x$end <- midpoints + (window / 2)
      
      ## Save
      write.bed(x, bed_out)
    }
    
  ## Apply
    create_window(nha_dir_19, nha_dir_19_w1000, 1000)
    create_window(nha_dir_38, nha_dir_38_w1000, 1000)
  
    
## Import beds
  impBed <- list()
  impBed$hg38 <- import.bed(nha_dir_38)
  impBed$hg38_w1000 <- import.bed(nha_dir_38_w1000)
  impBed$hg19 <- import.bed(nha_dir_19)
  impBed$hg19_w1000 <- import.bed(nha_dir_19_w1000)
  
  # and create a correspondence dataframe
  grch <- read.bed(nha_dir_19)
  grch$hg19 <- paste0(grch$V1, ":", grch$V2, "-", grch$V3)
  grch <- data.frame(hg19 = grch$hg19,
                     hg38 = grch$V4,
                     Enh = convert.id(grch$V4))
  grch <- grch[order(as.numeric(splitter(grch$Enh, "h", 2))),]
  
    
## And for speed, read in existing annotations
  # load("Final.rda", verbose = TRUE)

 
    
################################################################################################################################ #
## Analyse in-house ATAC-seq data ----
  
## Compare peak pileup
  # filter to candidates only
  p <- atac[which(atac$id %in% res.final$Enh.Pos),]
  
  # annotate as hit
  p$Hit <- p$id %in% res.final$Enh.Pos[which(res.final$HitPermissive)]
  
  # reformat
  p$Hit <- factor(p$Hit)
  levels(p$Hit) <- c("ns Peak", "Hit Enhancer")
  
  # plot
  pdf(file = "Inhouse ATAC Pileup.pdf", height = 2, width = 2)
  ggplot(p, aes(y = cpm, x = Hit, fill = Hit, colour = Hit)) +
    geom_violin(scale = "width", alpha = 0.75) +
    geom_boxplot(width = 0.2, outlier.shape = NA, colour = "black") +
    scale_fill_lancet() +
    scale_colour_lancet() +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    theme(axis.title.x = invis, panel.grid = invis, panel.border = invis,
          axis.line.y = element_line(), legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)) +
    labs(y = "ATAC pileup (CPM)")
  dev.off()
  


################################################################################################################################ #
## Enhancer classification using data from ENCODE V3 (Meuleman 2020), which we've called the Altius dataset ----
  
# This Altius work was exploratory in nature and does not feature heavily in the final manuscript


## Overlap peak coordinates to our screen
  # at FDR = 0.1%
  # typically, 100k per sample in 733 samples, and 76M total
  # however, they called consensus peaks
  # there are 3.5M consensus peaks (average width of ~200bp)
  # based on using the summit as centroids of all 76M peaks and aligning with local clustering
  # the final width is based on he "Full-Width at Half Maximum" on the histogram of peak signals across cell-types
  # hg38
  
## Read in
 # dat: expression matrix (DHS) for each peak
  load("../../../../PublicData/Meuleman_AltiusDHS_2020/dat_FDR01_hg38.RData", verbose = TRUE)
  
  # peak coordinates
  altius_coord <- read.delim("../../../../PublicData/Meuleman_AltiusDHS_2020/DHS_Index_and_Vocabulary_WM20190703.txt")

  # biosample metadata
  altius_meta <- read_xlsx("../../../../PublicData/Meuleman_AltiusDHS_2020/DHS_Index_and_Vocabulary_metadata.xlsx")
  altius_meta <- altius_meta[-which(is.na(altius_meta$`library order`)),] # this is the last row, and is entirely NA
  altius_meta <- as.data.frame(altius_meta)
  rownames(altius_meta) <- paste0(altius_meta$`Biosample name`, ".", altius_meta$`Altius Biosample ID`) # matches with colnames


## Intersect
  # write out
  altius_coord_filt <- altius_coord[,c(1:4, 7, 10)]
  colnames(altius_coord_filt) <- c("chr", "start", "end", "id", "summit", "NMF")
  write.bed(altius_coord_filt, "Altius_Coord.bed")
  
  # intersect
  call <- paste("intersectBed",
                "-a", nha_dir_38,
                "-b", "Altius_Coord.bed",
                "-wo", # number of bp overlap
                ">", "Altius_Intersect.bed")
  system(call, intern = FALSE, wait = TRUE)    
  
## Process
  # read in 
  altius_overlap <- read.delim("Altius_Intersect.bed", header = FALSE)
  
  # rename existing columns and add new ones related to the screen information
  colnames(altius_overlap) <- c(paste0("Screen_", c("Chr", "Start", "End", "ID")),
                      paste0("Altius_", colnames(altius_coord_filt)),
                      "Overlap")
  m <- match(altius_overlap$Screen_ID, guides$TargetCoord)
  altius_overlap$Screen_Enh <- guides$TargetID[m]
  # altius_overlap <- altius_overlap[which(altius_overlap$Screen_Enh %in% res.final$Enh),] # filter to only those enhancers we tested in the screen
  
  # where there are multiple encode peaks for a screen peak, pick the closest based on summit
  m <- match(altius_overlap$Screen_ID, atac$id)
  altius_overlap$Screen_Summit <- atac$abs_summit[m]
  
  for (j in unique(altius_overlap$Screen_Enh)) {
    if (length(which(altius_overlap$Screen_Enh == j)) == 1) next
    
    i <- which(altius_overlap$Screen_Enh == j)
    k <- which.min(abs(altius_overlap$Screen_Summit[i] - altius_overlap$Altius_summit[i])) # minimum distance between summits
    i <- i[-k] 
    altius_overlap <- altius_overlap[-i,] 
    
  }
  
  altius_overlap$Screen_Hit <- altius_overlap$Screen_Enh %in% hit.enh
  
  # finalise order
  altius_overlap$Altius_id <- paste0(altius_overlap$Altius_id, "_", altius_overlap$Altius_chr, ":", altius_overlap$Altius_start, "-", altius_overlap$Altius_end)
  altius_overlap <- altius_overlap[,c("Screen_Enh", "Screen_ID", "Screen_Hit", "Altius_id", "Altius_NMF", "Overlap", "Screen_Summit", "Altius_summit")]
  
  write.csv(altius_overlap, file = "Altius - Peak Annotation.csv")


## Filter expression matrix 
  # filter peaks to those that intersect our screen's
  altius_exp <- dat[splitter(altius_overlap$Altius_id, "_", 2),]
  rownames(altius_exp) <- paste0(altius_overlap$Screen_Enh, "_", rownames(altius_exp)) # no need to match order as the above row does that implicitly
  write.csv(altius_exp, file = "Altius - Summit-matched Peak Expression.csv")  
  # rm(dat); gc()
    
  # focus on primary biosamples
  prim <- which(altius_meta$`Biosample type` == "Primary")
  altius_meta_prim <- altius_meta[prim,]
  altius_exp_prim <- altius_exp[,rownames(altius_meta)]
  
  tab <- table(altius_meta_prim$Organ, altius_meta_prim$`Growth stage`) %>% as.data.frame.matrix()
  tab <- tab[,c(3,4,2,1)]
  write.csv(tab, file = "Altius - Primary Biosample Count.csv")
  
## Pool expression based on tissue/stage
  # a grouping variable
  altius_meta_prim$Group <- paste0(altius_meta_prim$Organ, "_", altius_meta_prim$`Growth stage`)
  
  # average across biosamples per group
  altius_pool <- list(Mean = list(),
                      FracNonzero = list(),
                      FracAbove0.5 = list())
  
  for (j in unique(altius_meta_prim$Group)) {
    w <- which(altius_meta_prim$Group == j)

    if (length(w) == 1) { # if only one sample in group
      
      altius_pool$Mean[[j]] <- altius_exp_prim[,w]
      altius_pool$FracNonzero[[j]] <- altius_exp_prim[,w] > 0
      altius_pool$FracAbove0.5[[j]] <- altius_exp_prim[,w] > 0.5
      
      
    } else {
      
      altius_pool$Mean[[j]] <- rowMeans(altius_exp_prim[,w])
      altius_pool$FracNonzero[[j]] <- rowMeans(altius_exp_prim[,w] > 0)
      altius_pool$FracAbove0.5[[j]] <- rowMeans(altius_exp_prim[,w] > 0.5)
      
    }

  }
 
  altius_pool <- lapply(altius_pool, function(x) {
    x <- do.call("cbind", x)
    x <- x[,sort(colnames(x))]
    x <- as.data.frame(x)
    rownames(x) <- rownames(altius_exp_prim)
    return(x)
  })

  # save
  save(altius_pool, file = "Altius - Pooled Expression Matrices.rda")
  write.csv(altius_pool$Mean, file = "Altius - Pooled Expression Matrix (Mean).csv")
  write.csv(altius_pool$FracNonzero, file = "Altius - Pooled Expression Matrix (Fraction Nonzero).csv")
  write.csv(altius_pool$FracAbove0.5, file = "Altius - Pooled Expression Matrix (Fraction Above 0.5).csv")
  
  
## Classify enhancers based on their expression patterns
  # calculate coefficient of variation
  sd <- apply(altius_pool$Mean, 1, sd)
  mean <- rowMeans(altius_pool$Mean)
  cv <- sd / mean
  
  # calculate the fraction of tissues in which the peak is open
  thresh <- 0.5
  onRate <- apply(altius_pool$Mean, 1, function(x) length(which(x > thresh))) / ncol(altius_pool$Mean)
  # y <- apply(altius_pool$Mean, 2, function(x) length(which(x > thresh))) / nrow(altius_pool$Mean)

  altius_stats <- data.frame(SD = sd,
                  Mean = mean,
                  CV = cv,
                  OnRate = onRate)
  
  pdf(file = "Altius - CV vs Onrate.pdf", height = 3, width = 3)
  ggplot(altius_stats, aes(x = OnRate, y = CV)) +
    geom_point() +
    scale_y_continuous(trans = "log2") +
    geom_hline(yintercept = 1, colour = "red") +
    geom_vline(xintercept = 0.5, colour = "red")
  dev.off()
  
  altius_stats$Classification <- "."
  altius_stats$Classification[which(altius_stats$OnRate <= 0.5)] <- "Tissue-specific"
  altius_stats$Classification[which(altius_stats$OnRate > 0.5 & altius_stats$CV > 1)] <- "Ubiquitous/Variable"
  altius_stats$Classification[which(altius_stats$OnRate > 0.5 & altius_stats$CV < 1)] <- "Ubiquitous/Constant"
  altius_stats$Classification[which(altius_stats$OnRate == 0)] <- "Not detected"
  
  altius_stats$Meuleman_NMF <- altius_overlap$Altius_NMF
  
  # save
  altius_stats$Enh <- splitter(rownames(altius_stats), "_", 1)
  altius_stats$Hit <- altius_stats$Enh %in% hit.enh
  altius_stats <- relocate(altius_stats, "Enh")
  write.csv(altius_stats, file = "Altius - CV vs Onrate.csv")
  
  # quick plots
  wilcox.test(altius_stats$OnRate ~ altius_stats$Hit) # p = 0.0018
  pdf(file = "Altius - Onrate in Hit Enh.pdf", height = 3, width = 3)
  ggplot(altius_stats, aes(x = Hit, y = OnRate)) +
    geom_violin(scale = "width", draw_quantiles = 0.5) +
    geom_quasirandom()
  dev.off()
  
 
## Cellular patterning of DNA accessibility
  # well, the authors acknowledge that cellular patterning is complex
  # and not cell-type specific
  # nnmf used to create what they call a "vocabulary" (k=16)
  # individual enhacners can be biologically annotaed by linear combinations of its NNMF components
  # very specific enhancers can be described using a single component, whilst constitutive enhancers use all
  # you can also summarise enhanceres using the single strongest component - this is the data we have.
  
  # tabulate
  tab_altius <- table(altius_coord$component) / nrow(altius_coord)
  tab_cand <- table(altius_overlap$Altius_NMF) / nrow(altius_overlap)
    
  ## Plot
    pdf(file = "Altius - NMF Classification.pdf", height = 4, width = 4.5)
    
    # plot fraction in ENCODE background versus candidates
    p <- data.frame(Altius = as.numeric(tab_altius),
                       Candidates = as.numeric(tab_cand),
                       Factor = names(tab_cand))
    
    ggplot(melt(p), aes(x = Factor, y = value, fill = variable)) +
      geom_col(position = "dodge") +
      theme_bw() +
      coord_flip() +
      guides(fill = guide_legend(title = "Dataset")) +
      labs(y = "Fraction of Peaks in Factor")
      
    # plot overrepresentation of hits
    p <- table(altius_overlap$Screen_Hit, altius_overlap$Altius_NMF) 
    p <- p / rowSums(p)
    p <- as.data.frame(p)
    colnames(p) <- c("Hit", "Factor", "Freq")

    ggplot(p, aes(x = Factor, y = Freq, fill = Hit)) +
      geom_col(position = "dodge") +
      theme_bw() +
      coord_flip() +
      guides(fill = guide_legend(title = "Peak is Hit")) +
      labs(y = "Fraction of Peaks in Factor")
    
    # above, but as odds ratio
    p <- list()
    for (j in unique(altius_overlap$Altius_NMF)) {
      fish <- table(altius_overlap$Altius_NMF == j, altius_overlap$Screen_Hit) %>% fisher.test()
      p[[j]] <- data.frame(OR = fish$estimate, p = fish$p.value, Sig = fish$p.value < 0.05)
    }

    p <- do.call("rbind", p)    
    p$Factor <- rownames(p)
    
    ggplot(p, aes(x = Factor, y = OR, fill = Sig, colour = Sig)) +
      geom_col() +
      coord_flip() +
      geom_hline(yintercept = 1) +
      theme_bw() +
      scale_y_continuous(trans = "log2") +
      labs(y = "Odds Ratio for Hit Enhancers")
    dev.off()
    
## Compare Altius annotation to our annotation
  tab <- table(altius_overlap$Altius_NMF, altius_stats$Classification)
  tab <- as.data.frame.matrix(tab)
  write.csv(tab, file = "Altius - Classification Comparison.csv")
 
  
  
## Add to final enhancer annotation dataframe
  candidate.annot$Altius_UbiquitousConstant <- candidate.annot$Enh %in% altius_stats$Enh[which(altius_stats$Classification == "Ubiquitous/Constant")]
  candidate.annot$Altius_UbiquitousVariable <- candidate.annot$Enh %in% altius_stats$Enh[which(altius_stats$Classification == "Ubiquitous/Variable")]
  candidate.annot$Altius_TissueSpecific <- candidate.annot$Enh %in% altius_stats$Enh[which(altius_stats$Classification == "Tissue-specific")]
  
  # enrichment table
  candidate.enrich <- run.EnhancerFisher(data.column = "Altius_UbiquitousConstant", rbind = FALSE)
  candidate.enrich <- run.EnhancerFisher(data.column = "Altius_UbiquitousVariable", rbind = TRUE)
  candidate.enrich <- run.EnhancerFisher(data.column = "Altius_TissueSpecific", rbind = TRUE)
  
################################################################################################################################ #
## Analyses of the Nott 2019 sorted ATAC and H3K27ac data ----
  
## These are antibody-sorted nuclei from the human brain, derived from samples from childhood and adolescence (4-18) 
    
## Paths
  # bigwig
  bw_nott <- list.files("../../../../PublicData/Nott2019/BigWigs", pattern = "bigWig", full.names = TRUE) # note hg19
  
  # bed
  bed_nott <- list.files("../../../../PublicData/Nott2019/hg38/", pattern = "bed", full.names = TRUE) # note hg38
  
  # samples
  nott_ct <- c("Astro", "Oligo", "Neuro", "Micro")
  
## Intersect
  ints_nott <- list()
  
  for (j in bed_nott) {
    out <- "../../Scratchspace/Temp.bed"
    
    loj(nha_dir_38, j, out, fun = "-c")
    
    name <- splitter(j, "/", 9) %>%
      gsub("optimal_peak", "", .) %>%
      gsub(".bed", "", .) %>% 
      gsub("IDR_ENCODE", "", .) %>% 
      gsub("__", "_", .) %>% 
      gsub("\\.", "", .)
    
    ints_nott[[name]] <- read.bed(out)
  }
  
    
## Call coverage
  # call
  covNott <- lapply(bw_nott, function(x) {
    y <- read_bw(x, impBed$hg19, ids = "standard")
    rownames(y) <- convert.id(y$id)
    return(y)
  })
  
  
  # name
  names(covNott) <- splitter(bw_nott, "/", 8) %>%
      gsub("human_", "", .) %>%
      gsub("_epilepsy_pooled_hg19.ucsc.bigWig", "", .) %>% 
      gsub("nuclei", "", .)
  

## Combine
  nott <- list()
  
  # intersects
  nott$Intersects <- sapply(ints_nott, function(x) x$V5)
  rownames(nott$Intersects) <- convert.id(ints_nott$LHX2_ATAC$V4)
  
  # coverage max
  nott$CoverageMax <- sapply(covNott, function(x) x$Max)
  rownames(nott$CoverageMax) <- rownames(covNott$LHX2_atac)
  
  # coverage mean 
  nott$CoverageMean <- sapply(covNott, function(x) x$Mean)
  rownames(nott$CoverageMean) <- rownames(covNott$LHX2_atac)
  
  # same formatting on colnames 
  colnames(nott$Intersects) <- gsub("LHX2", "Astro", colnames(nott$Intersects)) %>%
    gsub("NeuN", "Neuro", .) %>% 
    gsub("Olig2", "Oligo", .) %>%
    gsub("H3K27", "H3K27ac", .) %>%
    gsub("PU1", "Micro", .) 
  
  colnames(nott$CoverageMax) <- colnames(nott$CoverageMean) <- gsub("LHX2", "Astro", colnames(nott$CoverageMax)) %>%
    gsub("NEUN", "Neuro", .) %>% 
    gsub("OLIG2", "Oligo", .) %>%
    gsub("atac", "ATAC", .) %>%
    gsub("PU1", "Micro", .) 
  
  # filtered coverage
  nott$CoverageMax_PeakCalledEither <- nott$CoverageMax_PeakCalled <- nott$CoverageMax
  
  for (j in colnames(nott$CoverageMax_PeakCalled)) {
    # remove expression where there is no peak call
    nott$CoverageMax_PeakCalled[,j] <- nott$CoverageMax_PeakCalled[,j] * (nott$Intersects[,j] > 0) 
    
    # remove expression where there is no peak call in either ATAC or H3K27ac for the given cell-type
    k <- splitter(j, "_", 1)
    nott$CoverageMax_PeakCalledEither[,j] <- nott$CoverageMax_PeakCalledEither[,j] * (apply(nott$Intersects[,grep(k, colnames(nott$Intersects))], 1, any))
  }
  
  # average coverage across H3K27ac and ATAC
  nott$CoverageMax_PeakCalledEither_GeomMean <- list()
  geom_mean <- function(x) exp(mean(log(x)))
  
  for (j in nott_ct) {
    x <- nott$CoverageMax_PeakCalledEither # use the filtered version for being called in either cell-type
    x <- x[,grep(j, colnames(x))]
    y <- apply(x, 1, geom_mean) # geometric mean! Limits the skew driven by extreme values
    
    nott$CoverageMax_PeakCalledEither_GeomMean[[j]] <- y
  }
  
  nott$CoverageMax_PeakCalledEither_GeomMean <- do.call("cbind", nott$CoverageMax_PeakCalledEither_GeomMean)
  

  
## Save
  save(nott, file = "Coverage/Nott - Intersections and Coverage.rda")

  
## How many peaks intersect per dataset?
  p <- apply(nott$Intersects, 2, function(x) length(which(x > 0)))
  p <- as.data.frame(p)
  write.csv(p, file = "Coverage/Nott - Intersection Count to 979 Screen Enh.csv")
  
## Compare coverage within peak and non-peak calls
  p <- cbind(melt(nott$Intersects[,colnames(nott$CoverageMax)]),
             melt(nott$CoverageMax))
  
  p <- p[,c(1,2,3,6)]
  colnames(p) <- c("Enh", "Sample", "Intersect", "Coverage")  
  p$Intersect <- p$Intersect > 0

  pdf(file = "Coverage/Nott - Coverage in Intersecting Peaks.pdf", height = 4, width = 5)
  ggplot(p, aes(x = Sample, y = Coverage, colour = Intersect)) +
    geom_violin(draw_quantiles = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(trans = "log2") +
    labs(y = "Max Coverage")
  dev.off()
  
## Heatmap
  ## Scaled heatmap
    # function
    heatmap.scaled <- function(cov, do.scale = TRUE, ttl = "Blarg") {
      
      # setup
      x <- as.matrix(cov)
      if (do.scale) x <- apply(x, 1, function(y) y / max(y)) %>% t()
      colnames(x) <- colnames(cov)

      if (anyNA(x)) x <- x[-which(is.na(x[1,])),]

      # label rows by hit status
      rsc <- rownames(x) %in% res.final$Enh[which(res.final$HitPermissive)] %>% as.factor()
      levels(rsc) <- c("grey90", "firebrick1")
      rsc <- as.character(rsc)

      # label columns by cell-type
      csc <- splitter(colnames(x), "_", 1) %>% as.factor()
      # levels(csc)[2:6] <- "Neuronal"
      levels(csc) <- pal_lancet()(4)
      csc <- as.character(csc)


      gplots::heatmap.2(x, trace = "none", dendrogram = "row", Colv = FALSE, margins = c(10, 5),
                        RowSideColors = rsc, ColSideColors = csc, col = carto_pal(7, "Geyser"), symbreaks = FALSE, main = ttl)
      
      # gplots::heatmap.2(x, trace = "none", dendrogram = "row",, margins = c(10, 5),
      #                   col = carto_pal(7, "Geyser"), symbreaks = FALSE, main = ttl)

    }

    # get data
    p <- as.data.frame(nott$CoverageMax)
    p <- p[candidate.annot$Enh,] # reorder
    
    p <- list(ATAC = p[,grep("ATAC", colnames(p))],
              H3K27ac = p[,grep("H3K27ac", colnames(p))])
    
    
    # plot
    pdf(file = "Coverage/Nott - Heatmap.pdf", height = 10, width = 6)
    heatmap.scaled(p$ATAC, do.scale = TRUE, ttl = "ATAC (Scaled)")
    heatmap.scaled(log2(p$ATAC + 0.1), do.scale = FALSE, ttl = "ATAC (log2+0.1)")
    
    heatmap.scaled(p$H3K27ac, do.scale = TRUE, ttl = "H3K27ac (Scaled)")
    heatmap.scaled(log2(p$H3K27ac + 0.1), do.scale = FALSE, ttl = "H3K27ac (log2+0.1)")

    dev.off()

  
# ## Cell-type-specificity based on binary intersection
#   p <- nott$CoverageMax_Filt
  
## Cell-type-specificity based on fold-change
  # can use all peaks, regardless of peak call, due to the sanity check passing
  
  # get data
  p <- as.data.frame(nott$CoverageMax_Filt)
  p <- p[candidate.annot$Enh,] # reorder
  
  p <- list(ATAC = p[,grep("ATAC", colnames(p))],
            H3K27ac = p[,grep("H3K27ac", colnames(p))])
  
  # function
  calc.ctSpec <- function(covMx) {

    out <- apply(covMx, 1, function(a) {

      # highest coverage stage
      w <- which.max(a)

      # output
      out <- data.frame(MaxID = NaN,
                        Coverage = a[w],
                        RatioVsMedian = a[w] / median(a[-w]),
                        RatioVsRank2 = a[w] / max(a[-w]),
                        RatioVsLast = a[w] / min(a[-w]),
                        Specific = NA)

      out <- round(out, 2)
      out$MaxID <- colnames(covMx)[w]
      out$Specific <- out$RatioVsRank2 > 2
      return(out)

    })

    out <- do.call("rbind", out)

  }
  
  # apply function to call specificity as two-fold higher than all other ct
  nott_specific <- lapply(p, calc.ctSpec)
  nott_specific <- lapply(nott_specific, function(x) {
    # assay <- splitter(x$MaxID, "_", 2) %>% unique()
    x$MaxCt <- splitter(x$MaxID, "_", 1)
    x$SpecificCt <- "ns"
    w <- which(x$Specific)
    x$SpecificCt[w] <- x$MaxCt[w]
    
    x <- x[,c("MaxCt", "Coverage", "RatioVsMedian", "RatioVsRank2", "Specific", "SpecificCt")]
    # colnames(x) <- paste0(assay, "_", colnames(x))
    
    return(x)
  })
  nott_specific <- do.call("cbind", nott_specific)
  
  # save
  write.csv(nott_specific, "Coverage/Nott - Cell-type Specificity.csv")
  
  
## Output
  # annotation table 
  candidate.annot$AstSpecific_NottAtac <- candidate.annot$Enh %in% rownames(nott_specific)[which(nott_specific$ATAC.SpecificCt == "Astro")]
  candidate.annot$AstSpecific_NottH3K27ac <- candidate.annot$Enh %in% rownames(nott_specific)[which(nott_specific$H3K27ac.SpecificCt == "Astro")]

  # enrichment table
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_NottAtac", rbind = TRUE)
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_NottH3K27ac", rbind = TRUE)

  
  
################################################################################################################################ #
## Analyses of the Herring 2022 snATAC-seq data ----

## This is snATAC-seq from the human brain across maturation (foetal to early adulthood)

## Paths
  # bigwigs
  bw_herring <- list.files("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig", pattern = "bigwig", full.names = TRUE)

  # peaks
  bed_herring <- list.files("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/Bed/", pattern = "hg38", full.names = TRUE)


## Intersect.
  # this simply looks at whether the peaks called in each cell-type overlap out candidates
  ints_herring <- list()

  for (j in bed_herring) {
    out <- "../../Scratchspace/Temp.bed"

    loj(nha_dir_38, j, out, fun = "-c")

    name <- splitter(j, "bed_files_", 2) %>%
      splitter("_peaks", 1)

    ints_herring[[name]] <- read.bed(out)
  }


## Call coverage
  # this considers the mean coverage
  covHerring <- lapply(bw_herring, function(x) {
    print(x)

    y <- read_bw(x, impBed$hg19, ids = "get")
    
    m <- match(y$id, grch$hg19)
    y$id <- paste0("hg19_", y$id)
    rownames(y) <- grch$Enh[m]
    
    return(y)
  })
  

  # name
  names(covHerring) <- splitter(bw_herring, "BigWig/", 2) %>%
      splitter("_atac", 1)


## Combine
  herring <- list()

  # intersects
  herring$Intersects <- sapply(ints_herring, function(x) x$V5)
  rownames(herring$Intersects) <- convert.id(ints_herring$Astro$V4)

  # coverage max
  herring$CoverageMax <- sapply(covHerring, function(x) x$Max)
  rownames(herring$CoverageMax) <- rownames(covHerring$Astro_Adolescence)

  # coverage mean
  herring$CoverageMean <- sapply(covHerring, function(x) x$Mean)
  rownames(herring$CoverageMean) <- rownames(covHerring$Astro_Adolescence)

  # remove ChrX from intersects, and also set the same order
  herring$Intersects <- herring$Intersects[rownames(herring$CoverageMax),] # removes 5 enh on ChrX, as these were pre-filtered when getting coverage but not in the bed intersect

  # filtered coverage
  herring$CoverageMean_Filt <- herring$CoverageMean
  for (j in colnames(herring$CoverageMean_Filt)) {
    k <- gsub("_Adolescence|_Adult|_Childhood|_Fetal|_Infancy|_Neonatal", "", j)
    herring$CoverageMean_Filt[,j] <- herring$CoverageMean_Filt[,j] * (herring$Intersects[,k] > 0)
  }
  
  ## Pool coverage (averaged for Exc and Inh)
    # get coverage for glial celltypes, as these don't require further averaging
    x <- herring$CoverageMean # temp shorthand
    herring$CoverageMean_Pooled <- x[,grep("^Astro|^Oligo|^Micro", colnames(x))] %>% as.data.frame()
    
    # pool exc and inh using average coverage per stage
    pools <- list(Exc = c("L2_3", "L4", "L5_6"),
                  Inh = c("CGE_der", "MGE_der"))
    pools <- lapply(pools, paste, collapse = "|")
    stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

    for (j in names(pools)) {
      
      for (k in stages) {
        
        id <- paste0(j, "_", k)
        
        which_j <- grep(pools[[j]], colnames(x))
        which_k <- grep(paste0("_", k), colnames(x))
        w <- intersect(which_j, which_k)
        
        herring$CoverageMean_Pooled[,id] <- rowMeans(x[,w]) # the mean of depth-normalised coverage across the subtypes
        
      }
    }
    


## Save
  save(herring, file = "Coverage/Herring - Intersections and Coverage.rda")


## Number of intersecting peaks
  p <- apply(herring$Intersects, 2, function(x) length(which(x > 0)))
  p <- as.data.frame(p)
  write.csv(p, file = "Coverage/Herring - Intersection Count to 979 Screen Enh.csv")

## Compare coverage within peak and non-peak calls
  p <- list()

  for (j in colnames(herring$CoverageMean)) {
    k <- gsub("_Adolescence|_Adult|_Childhood|_Fetal|_Infancy|_Neonatal", "", j)
    p[[j]] <- data.frame(Intersect = (herring$Intersects[,k] > 0),
                         Coverage = herring$CoverageMean[,j],
                         Celltype = k,
                         Stage = gsub(paste0(k, "_"), "", j))
  }

  p <- do.call("rbind", p)

  pdf(file = "Coverage/Herring - Coverage in Intersecting Peaks.pdf", height = 8, width = 8)
  ggplot(p, aes(x = Stage, y = Coverage, colour = Intersect)) +
    geom_violin(draw_quantiles = 0.5, scale = "width") +
    theme_bw() +
    facet_wrap(~Celltype) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # scale_y_continuous(trans = "log2") +
    labs(y = "Mean Coverage")
  
  
  ggplot(p, aes(x = Stage, y = Coverage, colour = Intersect)) +
    geom_violin(draw_quantiles = 0.5, scale = "width") +
    theme_bw() +
    facet_wrap(~Celltype, scales= "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # scale_y_continuous(trans = "log2") +
    labs(y = "Mean Coverage")
  dev.off()

  
## Perform DE for cell-type-specificity (using repooled data)
  ## Get data
    offset <- 0.1
    hExp <- log2(herring$CoverageMean_Pooled + offset) %>% t() 
    
    hMeta <- data.frame(Sample = rownames(hExp))
    hMeta$Ct <- gsub("_Adolescence|_Adult|_Childhood|_Fetal|_Infancy|_Neonatal", "", hMeta$Sample) 
    hMeta$Stage <- gsub(paste(hMeta$Ct, collapse = "|"), "", hMeta$Sample) %>%
      sub("_", "", .)
    hMeta$Astro <- hMeta$Ct == "Astro"
    hMeta$Glia <- hMeta$Ct %in% c("Astro", "Micro", "Oligo")
    
  ## Run model
    hMods <- list()
    
    sumMod <- function(mod, var) {
      
      mod <- summary(mod)
      mod <- lapply(mod, function(x) {
        x <- x$coefficients
        y <- data.frame(P = x[paste0(var, "TRUE"), "Pr(>|t|)"],
                        EffectSize = x[paste0(var, "TRUE"), "Estimate"],
                        FDR = NaN)
        return(y)
      })
      mod <- do.call("rbind", mod)
      rownames(mod) <- sub("Response ", "", rownames(mod))
      
      # call hits
      mod$FDR <- p.adjust(mod$P, method = "fdr")
      mod[,paste0("Enriched")] <- mod$FDR < 0.05 & mod$EffectSize > 0
      mod[,paste0("Depleted")] <- mod$FDR < 0.05 & mod$EffectSize < 0
      
      # return
      return(mod)
    }
    
    hMods$Astro <- lm(hExp~Stage+Astro, data = hMeta) %>% sumMod(var = "Astro")
    hMods$Glia <- lm(hExp~Stage+Glia, data = hMeta) %>% sumMod(var = "Glia")
    
  # combine
    hMods <- do.call("cbind", hMods)
    write.csv(hMods, file = "Coverage/Herring - Cell-type Specificity Models.csv")
    

## What is a threshold to separate peak from non-peak?
  # calculate the FDR per cell-type
  # thresh_coverage <- list()
  test_threshes <- c(seq(0, 1, 0.1), 0.25)
  
  x <- list()
  for (j in colnames(herring$CoverageMean)) {

    # get cell-type calls
    k <- gsub("_Adolescence|_Adult|_Childhood|_Fetal|_Infancy|_Neonatal", "", j)
    i <- (herring$Intersects[,k]) > 0

    x[[j]] <- data.frame(Coverage = herring$CoverageMean[,j],
                    Peak = i,
                    Ct = j)
    
  }
  
  x <- do.call("rbind", x)
  x_peak <- x[which(x$Peak),]
  x_nonpeak <- x[-which(x$Peak),]
  
  aggregate(x$Coverage ~ x$Peak + x$Ct, FUN = summary)
  aggregate(x$Coverage ~ x$Peak, FUN = function(x) quantile(x, probs = c(0, 0.5, 0.9, 0.95, 0.99, 1)))
  
  y <- sapply(test_threshes, function(z) {
    data.frame(PeaksAboveThresh = length(which(x_peak$Coverage > z)) / nrow(x_peak),
               NonPeaksAboveThresh = length(which(x_nonpeak$Coverage > z)) / nrow(x_nonpeak),
               TPR = length(which(x_peak$Coverage > z)) / length(which(x$Coverage > z)))
  })
  
  colnames(y) <- paste0("CovThresh_", test_threshes)
  
    
## Add to annotation dataframes
  candidate.annot$AstSpecific_Herring <- candidate.annot$Enh %in% rownames(hMods_original)[which(hMods_original$Astro.Enriched)]
  candidate.annot$AstSpecific_Herring_Pooled <- candidate.annot$Enh %in% rownames(hMods)[which(hMods$Astro.Enriched)]
  candidate.annot$AstSpecific_Glia_Pooled <- candidate.annot$Enh %in% rownames(hMods)[which(hMods$Glia.Enriched)]

  # enrichment table
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_Herring", rbind = TRUE)
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_Herring_Pooled", rbind = TRUE)
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_Glia_Pooled", rbind = TRUE)

    
################################################################################################################################ #
## Cell-type-specificity in smaller resources with explicit calls  ----
    
## Using a different resource of calls from Morabito 2021
  
  # read in and process
  morabito <- read_xlsx("../../../../PublicData/Morabito2021_AD_snATAC_snRNA/41588_2021_894_MOESM7_ESM.xlsx", skip = 2)
  morabito <- morabito[which(morabito$cell_type == "ASC"),]
  morabito <- morabito[which(morabito$avg_logFC > 0),]
  morabito <- data.frame(chr = splitter(morabito$Peak, ":", 1),
                         start = splitter(splitter(morabito$Peak, ":", 2), "-", 1),
                         end = splitter(splitter(morabito$Peak, ":", 2), "-", 2),
                         id = morabito$Peak)
  write.bed(morabito, "../../../../PublicData/Morabito2021_AD_snATAC_snRNA/AstSpecific.bed")
  
  # overlap
  wawb(nha_dir_38, "../../../../PublicData/Morabito2021_AD_snATAC_snRNA/AstSpecific.bed", "AstSpecific_Morabito.bed")
  
  # read in
  morabitoAst <- read.bed("AstSpecific_Morabito.bed")
  
  # enrichments
  candidate.annot$AstSpecific_Morabito <- candidate.annot$Coord %in% morabitoAst$V4
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_Morabito", rbind = TRUE)

  

################################################################################################################################ #
## Superenhancers ----
  
  
## Use data from Hnisz et al 2013
  ## Directories
    hnisz_dir_in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/PROCESSED/liftOver_hg38/Astrocytes.hg38.super.bed" # note: this has super
    hnisz_dir_out <- "Superenhancers.bed"
    
    
  ## Read in and process
    hnisz_data <- read.table("/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/PROCESSED/liftOver_hg38/Astrocytes.hg38.bed") # note: this is not super
    hnisz_meta <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/1-s2.0-S0092867413012270-mmc7/Astrocytes.csv")
    
    # columns are: enhancer ID, chromosome, start, end, associated gene, enhancer rank, is enhancer a super-enhancer (1:yes, 0:no), H3K27ac ChIP-seq density (rpm/bp), read density in corresponding input sample (rpm/bp)
    colnames(hnisz_meta) <- c("ID", "chr", "start", "end", "gene", "rank", "super", "H3K27ac_density", "DensityInCorrespondingInputSample")
    
    hnisz_data <- hnisz_data[which(hnisz_data$V4 %in% hnisz_meta$ID[which(hnisz_meta$super == 1)]),] # filter to annotated superenhancers
    
    # save
    write.bed(hnisz_data, dir = hnisz_dir_in)
    
  
  ## Overlap
    loj(nha_dir_38, hnisz_dir_in, hnisz_dir_out)
  
  ## Analyse
    # read in
    hnisz_overlap <- read.bed(hnisz_dir_out)
    hnisz_overlap$V4 <- convert.id(hnisz_overlap$V4)
    
    # filter to peaks with overlap
    hnisz_overlap <- hnisz_overlap[which(hnisz_overlap$V9 != -1),]
    
    # save
    write.csv(hnisz_overlap, "Superenhancer - Annotation.csv")
  
  ## Output
    # annotation table 
    candidate.annot$Superenhancer <- candidate.annot$Enh %in% hnisz_overlap$V4
  
    # enrichment table
    candidate.enrich <- run.EnhancerFisher(data.column = "Superenhancer", rbind = TRUE)
    
  
################################################################################################################################ #
## Regions with human-specific activity ----
  
  
## Read in
  ## Regions with differential-acetylation in the human brain versus macaque and chimpanzee (Vermunt et al 2016)
    # samples are generally female, age 50+ in humans and but much younger for other species...
    vermunt_data <- list()
    vermunt_data$IncreaseVsMacaque_PFC <- read_xlsx("../../../../PublicData/Vermunt2016_HumanGainedH3K27ac/ST7_HumanAboveMacaque.xlsx", range = "A3:D2856")
    vermunt_data$IncreaseVsMacaque_WM <- read_xlsx("../../../../PublicData/Vermunt2016_HumanGainedH3K27ac/ST7_HumanAboveMacaque.xlsx", range = "M3:P5614")
    
    vermunt_data$IncreaseVsChimp_PFC <- read_xlsx("../../../../PublicData/Vermunt2016_HumanGainedH3K27ac/ST12_HumanAboveChimp.xlsx", range = "A3:D438")
    vermunt_data$IncreaseVsChimp_WM <- read_xlsx("../../../../PublicData/Vermunt2016_HumanGainedH3K27ac/ST12_HumanAboveChimp.xlsx", range = "M3:P1400")
    
    vermunt_data$GainVsMacaque <- read_xlsx("../../../../PublicData/Vermunt2016_HumanGainedH3K27ac/ST17_NewlyGainedEnhancers.xlsx", range = "A3:D1402")
    vermunt_data$GainVsChimp <- read_xlsx("../../../../PublicData/Vermunt2016_HumanGainedH3K27ac/ST17_NewlyGainedEnhancers.xlsx", range = "E3:H142")
    
    vermunt_data <- do.call("rbind", vermunt_data)
    vermunt_data$Category <- splitter(rownames(vermunt_data), "\\.", 1)
    
    vermunt_dir_in <- "../../../../PublicData/Vermunt2016_HumanGainedH3K27ac/Bed_combined.bed"
    vermunt_dir_out <- "Vermunt.bed"
    
    write.bed(vermunt_data, vermunt_dir_in)
  
## Intersect
  # full intersect information for Vermunt
  wawb(nha_dir_38, vermunt_dir_in, vermunt_dir_out)
  
## Read in and wrangle
 
    vermunt_intersect <- read.bed(vermunt_dir_out)
    
    # add to annotation
    for (j in unique(vermunt_intersect$V9)) {
      x <- vermunt_intersect[which(vermunt_intersect$V9 == j),]
      y <- paste0("HumanH3K27ac_", j, "_Vermunt2016")
      
      candidate.annot[,y] <- candidate.annot$Coord %in% x$V4
      candidate.enrich <- run.EnhancerFisher(data.column = y, rbind = TRUE) 
    }
 
    
################################################################################################################################ #
## Comparison to experimentally-validated enhancers ----

## Yao 2022 compiled this list from a range of CRISPR and CRISPRi studies in K562 cells
  # this includes Gasperini 2019

## Read in
  known_dir_in <- "../../../../PublicData/Yao_NatBiotech2022/ST2_KnownEnhancers.bed"
  known_dir_out <- "KnownEnhancerOverlap.bed"
  known_data <- read_xlsx("../../../../PublicData/Yao_NatBiotech2022/41587_2022_1211_MOESM3_ESM.xlsx",
                          sheet = "SupplementaryTable2",
                          skip = 1)
  write.bed(known_data, known_dir_in)

## Intersect
  # note that the above has no window to either side; this consistent with other analyses in this study
  wawb(nha_dir_38, known_dir_in, known_dir_out)

## Analyse
  known <- read.bed(known_dir_out)
  candidate.annot$ValidatedEnh_K562_Yao2022 <- candidate.annot$Coord %in% known$V4
  candidate.enrich <- run.EnhancerFisher(data.column = "ValidatedEnh_K562_Yao2022", rbind = TRUE)
    

    
################################################################################################################################ #
## Save ----

write.csv(candidate.enrich, file = "Final - Enrichments.csv")
write.csv(candidate.annot, file = "Final - Annotation Logical.csv", row.names = FALSE)
save(candidate.enrich, candidate.annot, file = "Final.rda")
# write.csv(annotation.detailed, file = "Final - Annotation Detailed.csv", row.names = FALSE)
  


