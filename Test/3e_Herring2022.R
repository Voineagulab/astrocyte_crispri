## This script looks at the relationship between our enhancers and a developmental brain maturation timecourse (Herring et al. 2022 Cell)

################################################################################################################################ #
## Setup ----


## Generic
# rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Herring2022")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(readxl)
  library(tidyr)
  library(liftOver)
  library(rtracklayer)
  library(rcartocolor)

  source("../../../Scripts/Functions.R")

## Load
  load("../../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  sig.colours <- c("black", "firebrick1")

## Results
  res.final <- read.csv("../../2_DE/Enh/Results Final.csv")
  

  
################################################################################################################################ #
## Some stats about the Herring data ----
  
  
## Get filenames and read in
  # beds
  bed_herring <- list.files("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/Bed/")
  names(bed_herring) <- splitter(bed_herring, "_peaks_final", 1) %>% splitter(., "files_", 2)
  bed_herring <- paste0("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/Bed/", bed_herring)
  names(bed_herring) <- splitter(bed_herring, "_peaks_final", 1) %>% splitter(., "files_", 2)
  
  peaks_herring <- lapply(bed_herring, function(x) read.delim(x, header = FALSE))
  names(peaks_herring) <- names(bed_herring)
  
  # meta data
  meta_herring <- read.csv("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/Processed_data_ATAC_BCs-meta-data.csv")
  
## Collect stats
  stats_herring <- data.frame(Ct = names(peaks_herring),
                              nCells = NA,
                              nPeaks = NA,
                              MeanWidth = NA,
                              MedianWidth = NA,
                              row.names = names(peaks_herring))
  
  # number of cells per cell-type
  tab <- table(meta_herring$predictedGroup) # these are the cell-types at a fine resolution; we need to pool these to the 9 main cell0types
  fun <- function(x) { grep(x, names(tab)) }
  stats_herring["Astro", "nCells"] <- sum(tab[fun("Astro")])
  stats_herring["CGE_der", "nCells"] <- sum(tab[fun("CGE")])
  stats_herring["L2_3", "nCells"] <- sum(tab[fun("L2")])
  stats_herring["L4", "nCells"] <- sum(tab[fun("L4")])
  stats_herring["L5_6", "nCells"] <- sum(tab[fun("L5")])
  stats_herring["MGE_der", "nCells"] <- sum(tab[fun("MGE")])
  stats_herring["Micro", "nCells"] <- sum(tab[fun("Micro")])
  stats_herring["Oligo", "nCells"] <- sum(tab[fun("Oligo")])
  stats_herring["Vas", "nCells"] <- sum(tab[fun("Vas")])
  
  # number of peaks
  stats_herring$nPeaks <- sapply(peaks_herring, nrow)
  
  # average size
  stats_herring$MeanWidth <- sapply(peaks_herring, function(x) mean(x$V3 - x$V2))
  stats_herring$MedianWidth <- sapply(peaks_herring, function(x) median(x$V3 - x$V2))
  

## Plot
  p <- melt(stats_herring)
  pdf(file = "Herring Peak Stats.pdf", height = 6, width = 8)
  ggplot(p, aes(x = Ct, y = value)) +
    geom_col() +
    facet_wrap(~variable, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = invis)
  dev.off()

  
################################################################################################################################ #
## Convert Herring data from hg19 to hg38 ----
  
  library()
  for (j in names(bed_herring)) {
    print(j)
    lift.bed(bed.in.19 = bed_herring[j], 
             bed.out.38 = gsub("\\.bed", "_hg38.bed", bed_herring[j]),
             isList = FALSE)
  }
  
    
################################################################################################################################ #
## Bed intersect ----
  
## Here, I use bedtools to find correspondence between our peaks and those defined in each of 9 cell-types in Herring 2022
  
## Intersect
  intersects <- list()
  
  # get filenames, where bed_herring was already generated above
  bed_nha <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  # bed_herring <- list.files("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/Bed/")
  # names(bed_herring) <- splitter(bed_herring, "_peaks_final", 1) %>% splitter(., "files_", 2)
  # bed_herring <- paste0("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/Bed/", bed_herring)
  # names(bed_herring) <- splitter(bed_herring, "_peaks_final", 1) %>% splitter(., "files_", 2)
  bed_herring_hg38 <- gsub("\\.bed", "_hg38.bed", bed_herring)
  
  # loop intersection
  for (j in 1:length(bed_herring_hg38)) {
    ct <- names(bed_herring_hg38)[j]
    print(ct)
    out <- paste0("intersectBed_", ct, ".bed")
    call <- paste("intersectBed",
                "-a", bed_nha,
                "-b", bed_herring_hg38[j],
                "-wa -wb", 
                ">", out)
  
    system(call, intern = FALSE, wait = TRUE)
    
    intersects[[ct]] <- read.delim(out, header = FALSE)
  } 
  
  # add Enh id
  intersects <- lapply(intersects, function(x) {
    m <- match(x$V4, guides$TargetCoord)
    x$V4 <- guides$TargetID[m]
    return(x)
  })
  
## Of our tested peaks, which have an overlap?
  # setup dataframe
  i <- data.frame(Enh = unique(guides$TargetID[which(guides$TargetCat == "Enh")]))
  i$Tested <- i$Enh %in% res.final$Enh
  i$HitCore <- i$Enh %in% res.final$Enh[which(res.final$HitCore)]
  i$HitPermissive <- i$Enh %in% res.final$Enh[which(res.final$HitPermissive)]

  # do they intersect?
  for (j in names(intersects)) {
    i[,j] <- i$Enh %in% intersects[[j]]$V4
  }

  # stats
  apply(i[,-c(1:4)], 2, function(x) table(x))
  
  # remove vasculate
  i <- i[,-grep("Vas", colnames(i))]
  
  # write
  write.csv(i, "Intersection Enhancer-level.csv")
  
  
## Report cell-type-specific peaks
  i$Glial <- rowSums(i[,c("Astro", "Micro", "Oligo")]) > 0
  i$OtherGlial <- rowSums(i[,c("Micro", "Oligo")]) > 0
  i$Neuronal <- rowSums(i[,c("CGE_der", "L2_3", "L4",  "L5_6", "MGE_der")]) > 0
  i$Broad <- i$Glial & i$Neuronal
  i$Absent <- rowSums(i[,5:12]) == 0
  
  ## Look at astrocyte peaks

    
    x <- data.frame(Ast = factor(i$Astro, levels = c("TRUE", "FALSE")),
                    Hit = factor(i$HitPermissive, levels = c("FALSE", "TRUE")),
                    Category = ".")
    
    x$Category[which(i$Astro & rowSums(i[,6:12]) == 0)] <- " Ast-Specific"
    x$Category[which(i$Astro & i$OtherGlial & !(i$Neuronal))] <- "Ast & Glia"
    x$Category[which(i$Astro & i$Neuronal & !(i$OtherGlial))] <- "Ast & Neu"
    x$Category[which(i$Astro & i$Neuronal & i$OtherGlial)] <- "Ast & Neu & Glia"
    x$Category[which(!(i$Astro) & i$OtherGlial & !(i$Neuronal))] <- "Non-Ast\nBut Glial"
    x$Category[which(!(i$Astro) & !(i$OtherGlial) & i$Neuronal)] <- "Non-Ast\nBut Neu"
    x$Category[which(!(i$Astro) & (i$OtherGlial) & i$Neuronal)] <- "Non-Ast\nBut Neu & Glial"
    # x$Category[which(!(i$Astro) & !(i$Absent))] <- "Other Ct Only"
    x$Category[which(i$Absent)] <- "Zero Ct"
    
    levels(x$Hit) <- c("Non-hit", "Permissive hit")
    levels(x$Ast) <- c("Peak Intersects Herring Ast", "Peak Not In Herring Ast")
  
    pdf(file = "Binary Intersections - Peak Categories.pdf", height = 4, width = 9)
    ggplot(x, aes(x = Category, fill = Hit)) +
      geom_bar() +
      theme_bw() +
      facet_wrap(~Ast, scales = "free_x") +
      theme() +
      scale_fill_lancet() +
      labs(y = "Count of Peaks", x = "Peak Category Based on Binary Intersection") +
      scale_y_continuous(expand = c(0,0))
    dev.off()
    
    
    ##  THEN STATS IRINA WANTS. THEN RNA
  
  
################################################################################################################################ #
## Calculate coverage directly for our peaks ----
  
  
## This section has now been transferred to Script 3b  
    
# ## Load required libraries
#   library(Rsamtools)
#   library(rtracklayer)
# 
# 
# ## Directories
#   bed_nha_hg19 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed" # it is easier to convert these than the bigwigs
#   bed_nha_hg19_file <- import.bed(bed_nha_hg19)
#   bed_nha_hg19_file <- bed_nha_hg19_file
#   bw_herring <- list.files("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig/")
#   bw_herring <- bw_herring[grep("bigwig", bw_herring)]
#   bw_herring <- paste0("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig/", bw_herring)
#   
# ## Aside: create hg38 bigwigs
#   # these are currently not used in the analysis (instead, peaks are converted to hg19 when relevant), but are useful for IGV
#   
#     # chain <- rtracklayer::import.chain(chainFile_19to38)
#     # x <- rtracklayer::import(x, as = "GRanges")
#     # y <- rtracklayer::liftOver(x, chain)
#     # export(y@unlistData, "../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig/hg38/Astro_Adolescence_atac_insertions_hg38.bigwig")
#   
# 
#   
#   
# 
# ## First: for each our our enhancers, quantify the coverage in each of the snATAC-seq bigwigs
#   coverage_nhaPeaks <- lapply(bw_herring, function(x) {
#     print(x)  
#     y <- import(x, selection = bed_nha_hg19_file, as = "NumericList")
#     z <- sapply(y, mean)
#     names(z) <- bed_nha_hg19_file$name[-grep("chrX", bed_nha_hg19_file$name)] # $name is the GRCh38 coordinate!
#     return(z)
#   })
#   
#   # requires ~20s
#   
#   # this is the mean of coverages across the bases, and is thus normalised to peak width
#   
#   # seems to not do chrX
#   
#   names(coverage_nhaPeaks) <- splitter(bw_herring, "BigWig\\/", 2) %>% sub("_atac_insertions.bigwig", "", .)
#   coverage_nhaPeaks <- do.call("cbind", coverage_nhaPeaks)
#   
#   # get the enhancer id
#   m <- match(rownames(coverage_nhaPeaks), guides$TargetCoord)
#   rownames(coverage_nhaPeaks) <- guides$TargetID[m]
#   
#   # add call annotation
#   m <- match(rownames(coverage_nhaPeaks), i$Enh)
#   coverage_nhaPeaks <- cbind(i[m,], coverage_nhaPeaks)
#   
#   # rename columns
#   colnames(coverage_nhaPeaks) <- gsub("_der", "", colnames(coverage_nhaPeaks))
#   colnames(coverage_nhaPeaks) <- gsub("L5_6_", "L56_", colnames(coverage_nhaPeaks))
#   colnames(coverage_nhaPeaks) <- gsub("L2_3_", "L23_", colnames(coverage_nhaPeaks))
#   
#   # save
#   write.csv(coverage_nhaPeaks, "NHA Peak Quantification in Herring (No ChrX).csv")
# 
# ## Sanity check plot
#   # when a peak is called, is there higher coverage?
#   
#   for (j in colnames(peaks_herring)) {
#     w <- which(colnames(coverage_nhaPeaks) == j)
#     p <- melt(coverage_nhaPeaks[,c(w, 14:58)])
#     
#     p$Ct <- splitter(p$variable, "_", 1)
#     p$Stage <- splitter(p$variable, "_", 2) 
#     p$Stage <- factor(p$Stage, levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
#     
#     ggplot(p, aes_string(x = j, y = "value", colour = j)) +
#       geom_quasirandom() +
#       facet_grid(Stage~Ct, scales = "free")
#     
#     ggplot(p, aes_string(x = "Ct", y = "value", colour = "Stage")) +
#       geom_violin() +
#       facet_wrap(as.formula(paste0("~", j)))
#   }
# 
#     
# ## Plots: heatmap
#   p <- coverage_nhaPeaks[,-c(1:13)]
#   
#   pdf(file = "NHA Peak Quantification in Herring - Heatmaps.pdf", height = 10, width = 10)
#   
#   # of all peaks
#   pp <- t(p)
#   pp <- pp[,-which(colSums(pp) == 0)]
#   pp <- cor(pp)
#   heatmap.2(pp, trace = "none")
#   
#   # of hit peaks
#   pp <- t(p[which(coverage_nhaPeaks$HitPermissive),])
#   pp <- pp[,-which(colSums(pp) == 0)]
#   pp <- cor(pp)
#   heatmap.2(pp, trace = "none")
#   
#   dev.off()
  

################################################################################################################################ #
## Calculate coverage for linked peaks ----

# this differs from the above section in that coverage is calculated for the peaks that intersect with our NHA peaks
 
## Loop
  
coverage_herring <- list()
  
  
  for (j in names(bed_herring)) { # this list was generated two sections up. it contains the correspondence between our peaks and the Herring peaks of cell-type j
    
    print(j)
    
    # skip vas
    if (j == "Vas") next
    
    # which bigwigs to import
    x <- bw_herring[grep(j, bw_herring)]
    
    # and which bed to use
    ct.bed <- import.bed(bed_herring[j])
    
    # get coverage
    y <- lapply(x, function(y) {
      print(y)
      y <- import(y, selection = ct.bed, as = "NumericList")
      z <- sapply(y, mean)
      names(z) <- ct.bed$name # $name is the GRCh38 coordinate!
      return(z)
    })
    
    
    # single matrix
    names(y) <- splitter(x, "BigWig\\/", 2) %>% sub("_atac_insertions.bigwig", "", .)
    y <- do.call("cbind", y)
    y <- as.data.frame(y)
    
    # and revert the organisation to the NHA overlap
    z <- data.frame(Enh = unique(guides$TargetID[which(guides$TargetCat == "Enh")]))
    z$Tested <- z$Enh %in% res.final$Enh
    z$HitCore <- z$Enh %in% res.final$Enh[which(res.final$HitCore)]
    z$HitPermissive <- z$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
    
    z$HerringID <- NA
    z[,colnames(y)] <- NA
    
    for (k in z$Enh) {
       if (k %in% intersects[[j]]$V4) {
        w <- which(intersects[[j]]$V4 == k)
        w <- intersects[[j]]$V8[w] # the id(s) intersecting the NHA peak k
        
        # note that it is rare for a NHA peak to intersect many peaks in the Herring data
        
        z[which(z$Enh == k), "HerringID"] <- paste(paste0(j, "_", w), collapse = " / ")
        z[which(z$Enh == k), colnames(y)] <- colMeans(y[w,]) # mean coverage at intersecting peak(s)
        
        
      } else {
        
        next
        
      }
    }
      
    # annotate Herring enhancers by their NHA overlap
    y$NHA_Enh <- NA
    for (k in rownames(y)) {
      if (k %in% intersects[[j]]$V8) {
        w <- which(intersects[[j]]$V8 == k)
        w <- paste(intersects[[j]]$V4[w], collapse = "/")
        y[k,"NHA_Enh"] <- w
      } else {
        
        next
        
      }
    }
    
    
    # output
    coverage_herring[[j]] <- list(By_HerringPeak = y,
                                  By_NHAPeak = z)
    
  }

## Save
  # rda
  save(coverage_herring, file = "Herring Peak Quantification.rda")
  
  # as a single dataframe
  cov <- lapply(coverage_herring, function(x) {
    x <- x$By_NHAPeak
    x <- x[,-c(1:4)]
  })

  cov <- do.call("cbind", cov)  
  cov <- cbind(coverage_herring$Astro$By_NHAPeak[,1:4], cov)

  write.csv(cov, file = "Herring Peak Quantification.csv")  
  
## Plot
  ## Setup
    p <- cov
  
    colnames(p) <- gsub("_der", "", colnames(p))
    p <- p[which(p$HitPermissive), grep("_",colnames(p))]
    p <- p[,-grep("Herring", colnames(p))]
    rownames(p) <- paste0("Enh", rownames(p))
    
    # replace NA with 0
    p[is.na(p)] <- 0
    
  ## Heatmap of correlations
    pdf(file = "Herring-NHA Peak Heatmap.pdf", height = 10, width = 10)
    cor <- cor(t(p))
    rem <- which(is.na(cor[1,]))
    cor <- cor[-rem,-rem]
    
    heatmap.2(cor, trace = "none", main = "Correlation Heatmap")
    
  ## Heatmap of raw values
    
    
    scaled <- apply(p, 1, scale)
    rem <- which(is.na(scaled[1,]))
    scaled <- scaled[,-rem]
    rownames(scaled) <- colnames(p)[-rem]
    heatmap.2(scaled, trace = "none", main = "Coverage Scaled-Per-CelltypeStage")
    dev.off()
    
  
## Run statistics
  ## Prepare coverage data
    p <- cov
  
    colnames(p) <- gsub("_der", "", colnames(p))
    p <- p[which(p$HitPermissive), grep("_",colnames(p))]
    p <- p[,-grep("Herring", colnames(p))]
    rownames(p) <- paste0("Enh", rownames(p))
    
    # replace NA with 0
    p[is.na(p)] <- 0
    
    # subset columns
    colnames(p) <- gsub("_der", "", colnames(p))
    
    # filter rows with more than 50% zeroes (when there is no corresponding peak)
    remove <- apply(p, 1, function(x) length(which(x == 0)) / length(x))
    p <- p[remove < 0.5,]
    
    # transpose
    p <- t(p)
    
    
  ## Prepare metadata
    y <- data.frame(Ct = splitter(rownames(p), "\\.", 1))
    y$IsAst <- y$Ct == "Astro"
    y$IsGlial <- y$Ct %in% c("Astro", "Micro", "Oligo")
    y$Stage <- gsub(paste(unique(y$Ct), collapse = "|"), "", rownames(p)) %>% splitter(., "_", 2) %>% factor(levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
    y$StageNumeric <- as.numeric(y$Stage)
    
  
  ## First: effect of cell-type in ANOVA
    anova <- list()
    
    for (j in colnames(p)) {
      print(j)
      a <- lm(p[,j] ~ Ct + StageNumeric, data = y) 
      anova[[j]] <- tryCatch(expr = {anova(a)},
                        error = NA)
    }
    
    anova <- sapply(anova, function(x) x["Ct", "Pr(>F)"])
   

  ## Second: effect of ast vs other, and glial vs other
    lm.AstEffect <- lm(p ~ IsAst + StageNumeric, data = y) %>% summary()
    lm.AstEffect <- sapply(lm.AstEffect, function(a) a$coefficients["IsAstTRUE", "Pr(>|t|)"]) # p-value for effect of the peak being in Ast vs. all other cell-types
    
    lm.GlialEffect <- lm(p ~ IsGlial + StageNumeric, data = y) %>% summary()
    lm.GlialEffect <- sapply(lm.GlialEffect, function(a) a$coefficients["IsGlialTRUE", "Pr(>|t|)"]) # p-value for effect of the peak being in Ast vs. all other cell-ty;es
    
  ## Third: effect of stage within astrocytes
    # setup
    pp <- cov
    
    pp <- pp[which(pp$HitPermissive),grep("Astro", colnames(pp))]
    pp <- pp[,-grep("Herring", colnames(pp))]
    pp <- pp[-which(is.na(pp$Astro.Astro_Adolescence)),]
    rownames(pp) <- paste0("Enh", rownames(pp))
    pp <- t(pp)
    
    # linear model
    use.samps <- which(y$Ct == "Astro")
    lm.StageEffect <- lm(pp ~ StageNumeric, data = y[use.samps,]) %>% summary()
    lm.StageEffect <- sapply(lm.StageEffect, function(a) a$coefficients["StageNumeric", c("Estimate", "Pr(>|t|)")]) %>% t() %>% as.data.frame()
    
    
    # an alternative: fold-changes.
    fc.StageEffect <- apply(pp, 2, function(a) {
      # highest coverage stage
      w <- which.max(a)
      
      # output
      out <- data.frame(MaxStage = NA,
                        Coverage = a[w],
                        RatioVsMedian = a[w] / median(a[-w]),
                        RatioVsRank2 = a[w] / max(a[-w]),
                        RatioVsLast = a[w] / min(a[-w]))
      
      out <- round(out, 2)
      out$MaxStage <- splitter(rownames(pp)[w], "_", 2)
      return(out)
      
    })
    
    fc.StageEffect <- do.call("rbind", fc.StageEffect)
      
    
    
  ## For the first and second and points, a plot!
    # side-by-side heatmaps. Both have enh on the y. Left is raw scores, right is 3x p-values
    
    ## Heatmap 1: raw scores
      p1 <- as.data.frame(t(p))
      p1$Enh <- rownames(p1)
      p1$Anova <- (anova)
      p1$AstEffect <- (lm.AstEffect)
      p1$GlialEffect <- (lm.GlialEffect)
    
      p2 <- melt(p1[,-c(47:49)])
      p2$Ct <- splitter(p2$variable, "\\.", 1) 
      p2$Ct <- factor(p2$Ct, levels = c("Astro", "Oligo", "Micro", "L2_3", "L4", "L5_6", "CGE", "MGE"))
      p2$variable <- splitter(p2$variable, "\\.", 2) 
      p2$Enh <- sub("Enh", "", p2$Enh)
      p2$Stage <- gsub(paste(unique(p2$Ct), collapse = "|"), "", p2$variable)
      p2$Stage <- sub("_", "", p2$Stage)
      p2$Stage <- factor(p2$Stage, levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
      p2$Coverage <- cut(p2$value, c(0, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 2), include.lowest = TRUE) # why is 0.15 the lowest bin? Because all but 2 enh have at least one cell greater than this
      levels(p2$Coverage) <- c("<15", "15-20", "20-25", "25-30", "30-40", "40-50", "50+")
      
      # and clustering!
      clust <- hclust(dist(t(p)))
      clust <- colnames(p)[clust$order]
      clust <- sub("Enh", "", clust)
      
      p2 <- ggplot(p2, aes(y = Enh, x = Stage, fill = Coverage, label = round(value*100,0))) +
        geom_tile(colour = "black") +
        geom_text(size = 2) +
        theme_bw() +
        theme(panel.grid = invis, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.title.x = invis, axis.ticks = invis) +
        labs() +
        scale_y_discrete(expand = c(0,0), limits = clust) +
        scale_x_discrete(expand = c(0,0)) +
        facet_wrap(~Ct, nrow = 1, strip.position = "bottom") +
        scale_fill_carto_d(palette = "Temps")
      
    ## Heatmap 2: P-values
      p3 <- (p1[,c(46:49)])
      p3$Enh <- sub("Enh", "", p3$Enh)
      p3$Enh <- factor(p3$Enh, levels = clust)
      p3 <- melt(p3, id.vars = "Enh") 
      p3$Bin <- cut(-log10(p3$value), c(0,1,2,3,4,5,7,10,20))
      levels(p3$Bin) <- c("ns", "< e-1", "< e-2", "< e-3", "< e-4", "< e-5", "< e-7", "< e-10")
      levels(p3$variable) <- c("All Ct Anova", "Stage-In-Ast", "Glia vs Neu")
      cols <- c((carto_pal(7, "Burg")), "black")
      
      p3 <- ggplot(p3, aes(y = Enh, x = "Adolescence", fill = Bin)) +
        geom_tile(colour = "black") +
        theme_bw() +
        theme(panel.grid = invis, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "white"),
              axis.title.x = invis, axis.ticks = invis) +
        labs() +
        scale_y_discrete(expand = c(0,0), limits = clust) +
        scale_x_discrete(expand = c(0,0)) +
        facet_wrap(~variable, nrow = 1, strip.position = "bottom") +
        scale_fill_manual(values = cols)
      
    ## Phew, output that...
      pdf(file = "Heatmap of Raw Values and Statistics.pdf", height = 10, width = 13)
      plot_grid(p2, p3, rel_widths = c(1, 0.45))
      dev.off()
  
      
      
        
    
## Calculate statistics for all enhancers. It would be a
  ## Prepare coverage data
    q <- cov
  
    colnames(q) <- gsub("_der", "", colnames(q))
    q <- q[, grep("_",colnames(q))]
    q <- q[,-grep("Herring", colnames(q))]
    rownames(q) <- paste0("Enh", rownames(q))
    
    # replace NA with 0
    q[is.na(q)] <- 0
    
    # subset columns
    colnames(q) <- gsub("_der", "", colnames(q))
    
    # filter rows with more than 50% zeroes (when there is no corresponding peak)
    remove <- apply(q, 1, function(x) length(which(x == 0)) / length(x))
    q <- q[remove < 0.5,]
    
    # transpose
    q <- t(q)
    
  ## Prepare metadata
    # y <- data.frame(Ct = splitter(rownames(q), "\\.", 1))
    # y$IsAst <- y$Ct == "Astro"
    # y$IsGlial <- y$Ct %in% c("Astro", "Micro", "Oligo")
    # y$Stage <- gsub(paste(unique(y$Ct), collapse = "|"), "", rownames(q)) %>% splitter(., "_", 2) %>% factor(levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
    # y$StageNumeric <- as.numeric(y$Stage)
    
    # can reuse y  
  
  ## List
    all.enh.stats <- list()
      
  ## First: effect of cell-type in ANOVA
    all.enh.stats$anova <- list()
    
    for (j in colnames(q)) {
      print(j)
      a <- lm(log2(q[,j] + 0.01) ~ Ct*StageNumeric, data = y) 
      all.enh.stats$anova[[j]] <- tryCatch(expr = {anova(a)},
                        error = NA)
    }
    
    all.enh.stats$anova <- sapply(all.enh.stats$anova, function(x) x[c("Ct", "StageNumeric"), "Pr(>F)"]) %>% t()
   
    # for the anova, Irina would like all ct pvalues (that's in the lm, NOT anova)
    
    
  ## Second: effect of ast vs other, and glial vs other
    all.enh.stats$AstVsOther <- all.enh.stats$GliaVsNeu <- list()
    
    all.enh.stats$AstVsOther <- lm(log2(q + 0.01) ~ IsAst + StageNumeric, data = y) %>% summary()
    all.enh.stats$AstVsOther <- sapply(all.enh.stats$AstVsOther, function(a) a$coefficients["IsAstTRUE", c("Estimate", "Pr(>|t|)")]) %>% t() # p-value for effect of the peak being in Ast vs. all other cell-types
    
    all.enh.stats$GliaVsNeu <- lm(log2(q + 0.01) ~ IsGlial + StageNumeric, data = y) %>% summary()
    all.enh.stats$GliaVsNeu <- sapply(all.enh.stats$GliaVsNeu, function(a) a$coefficients["IsGlialTRUE", c("Estimate", "Pr(>|t|)")]) %>% t()  # p-value for effect of the peak being in Ast vs. all other cell-types
    
  ## Third: effect of stage within astrocytes
    # setup
    q <- cov
    
    q <- q[,grep("Astro", colnames(q))]
    q <- q[,-grep("Herring", colnames(q))]
    q <- q[-which(is.na(q$Astro.Astro_Adolescence)),]
    rownames(q) <- paste0("Enh", rownames(q))
    q <- t(q)
    
    # linear model
    use.samps <- which(y$Ct == "Astro")
    
    all.enh.stats$StageLm <- lm(log2(q + 0.01) ~ StageNumeric, data = y[use.samps,]) %>% summary()
    all.enh.stats$StageLm <- sapply(all.enh.stats$StageLm, function(a) a$coefficients["StageNumeric", c("Estimate", "Pr(>|t|)")]) %>% t() %>% as.data.frame()
    
    
    # an alternative: fold-changes.
    all.enh.stats$StageFC <- apply(q, 2, function(a) {
      # highest coverage stage
      w <- which.max(a)
      
      # output
      out <- data.frame(MaxStage = NA,
                        Coverage = a[w],
                        RatioVsMedian = a[w] / median(a[-w]),
                        RatioVsRank2 = a[w] / max(a[-w]),
                        RatioVsLast = a[w] / min(a[-w]))
      
      out <- round(out, 2)
      out$MaxStage <- splitter(rownames(q)[w], "_", 2)
      return(out)
      
    })
    
    all.enh.stats$StageFC <- do.call("rbind", all.enh.stats$StageFC)
      
    
  ## Analyse: effect of cell-type
    x <- do.call("cbind", all.enh.stats[1:3])
    x <- as.data.frame(x)
    colnames(x) <- c("ANOVA_Ct", "ANOVA_Stage", "Glia.Coef", "Glia.P", "Ast.Coef", "Ast.P")
    x$Hit <- rownames(x) %in% res.final$Enh[which(res.final$HitPermissive)]
    
    x$ANOVA_Ct.FDR <- p.adjust(x$ANOVA_Ct, method = "fdr")
    x$ANOVA_Stage.FDR <- p.adjust(x$ANOVA_Stage, method = "fdr")
    x$Glia.FDR <- p.adjust(x$Glia.P, method = "fdr")
    x$Ast.FDR <- p.adjust(x$Ast.P, method = "fdr")
    
    table(x$Glia.FDR < 0.05, x$Hit) %>% fisher.test() # 72 hits are Glia+. ns versus non-hits.
    table(x$Ast.FDR < 0.05, x$Hit) %>% fisher.test() # 47 hits are Ast-different. ns versus non-hits.
    table(x$ANOVA_Ct.FDR < 0.05, x$Hit) %>% fisher.test() # 112 of 114 hits are ANOVA+. ns versus non-hits.
    table(x$ANOVA_Stage.FDR < 0.05, x$Hit) %>% fisher.test() # 51 of 114 hits are ANOVA+. ns versus non-hits.
    
    write.csv(x, file = "Herring Peak Quantification Stats - Ct Effect.csv")
    
  ## Analyse: stage
    x <- do.call("cbind", all.enh.stats[4:5])
    rownames(x) <- sub("Response ", "", rownames(x))
    x$Hit <- rownames(x) %in% res.final$Enh[which(res.final$HitPermissive)]
    x$StageLM.FDR <- p.adjust(x$`StageLm.Pr(>|t|)`, method = "fdr")
  
    table(x$StageLM.FDR < 0.05, x$Hit) # 1 is sig, but not a hit
    
    write.csv(x, file = "Herring Peak Quantification Stats - Stage Effect.csv")
    
    
  ## Save
    save(all.enh.stats, file = "Herring Peak Quantification Stats.rda")


## Plot Ast coverage
  ## Just the distribution within each stage
    x <- data.frame(Enh = unique(guides$TargetID[which(guides$TargetCat == "Enh")]))
    x$Tested <- x$Enh %in% res.final$Enh
    x$HitCore <- x$Enh %in% res.final$Enh[which(res.final$HitCore)]
    x$HitPermissive <- x$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
    
    m <- match(x$Enh, cov$Enh)
    x <- cbind(x, cov[m,6:11])
    x$Max <- apply(x[,5:10], 1, max)
  
    colnames(x) <- gsub("Astro\\.Astro_", "", colnames(x))
    x <- melt(x[,-c(1:3)])
    
    x <- x[-which(is.na(x$value)),]
    x$variable <- factor(x$variable, levels = levels(x$variable)[c(4:6, 3, 1, 2, 7)])
    
    pdf(file = "Hit vs Non-hit Coverage Comparison Across Stages.pdf", height = 5, width = 8)
    ggplot(x, aes(x = variable, y = value, colour = HitPermissive, alpha = HitPermissive)) +
      geom_quasirandom(dodge.width = 0.8, alpha = 0.5) +
      stat_summary(fun = "median", position = position_dodge(width = 0.8), shape = "-", colour = "black", size = 2) +
      scale_alpha_manual(values = c(0.99,1)) +
      labs(y = "Coverage At Herring Peak") +
      scale_colour_manual(values = pal_lancet()(2)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line()) 
    dev.off()
  
  ## And on stage specificity
    # setup dataframe
    x <- all.enh.stats$StageFC
    x$Hit <- rownames(x) %in% res.final$Enh[which(res.final$HitPermissive)]
    x$Specific <- factor(x$RatioVsMedian > 2)
    levels(x$Specific) <- c("No", "Yes")
    x$MaxStage <- factor(x$MaxStage, levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
    
    # plot 1: highest stage
    pdf(file = "Stage Specificity.pdf", height = 3, width = 8)
    ggplot(x, aes(x = MaxStage, fill = Specific)) +
      geom_bar(colour = "black", position = "dodge") +
      # scale_fill_manual() +
      scale_alpha_manual(values = c(1, 0.3)) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      guides(fill = guide_legend(title = "Specific\n(2-fold > median)")) +
      # facet_wrap(~Specific) +
      labs(y = "Enhancer Has Highest Coverage\nin the Given Stage") +
      theme(panel.border = invis, axis.line.y = element_line())
    
    # plot 2: stage specificity of hit enhancers
    ggplot(x, aes(x = Hit, y = log2(RatioVsMedian), colour = MaxStage)) +
      geom_quasirandom() +
      scale_colour_lancet() +
      theme_bw() +
      scale_y_continuous(limits = c(0, 4), expand = c(0,0)) +
      theme(panel.border = invis, axis.line.y = element_line()) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(y = "Coverage in Highest Stage vs Median Other Stage\n(log2fc)", x = "Enhancer Is Hit")
  dev.off()      
  
  table(x$Specific, x$Hit) %>% fisher.test() # ns
    
  
################################################################################################################################ #
## Trends in RNA ----
    
## Read in
  load("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/DE_Trends.rda")
  load("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/Pseudobulk.rda")
  
## Heatmap of expression for hit genes.
  ## Collect data
    samps <- paste(herring_pb$Meta$Ct, herring_pb$Meta$Stage) %>% unique()
    x <- list()
    
    for (j in samps) {
      print(j)
      w <- which(paste(herring_pb$Meta$Ct, herring_pb$Meta$Stage) == j)
      if (length(w) == 1) {
        x[[j]] <- herring_pb$CPM[,w]
      } else {
        x[[j]] <- rowMeans(herring_pb$CPM[,w])  
      }
    }
    
    x <- do.call("cbind", x)
  
  ## Filter
    # to hit genes
    x <- x[rownames(x) %in% (res.final$Gene[which(res.final$HitPermissive)]),]
    
    # to key ct
    x <- x[,-grep("Vas", colnames(x))]
    
  ## Plot 1 - all stages
    clust <- hclust(dist(log2(x+0.5)))
    clust <- rownames(x)[clust$order]
    cols <- c(carto_pal(7, "Temps"), carto_pal(2, "Fall")[2])
    
    p <- as.data.frame(x)
    p$Symbol <- rownames(p)
    p <- melt(p)
    p$Ct <- factor(splitter(p$variable, " ", 1))
    p$Ct <- factor(p$Ct, levels = levels(p$Ct)[c(1,8:10,11:14,2:7)])
    p$Stage <- factor(splitter(p$variable, " ", 2), levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
    p$Bin <- cut(p$value, c(0, 1, 10, 25, 50, 100, 500, 1000, 10000), include.lowest = TRUE)
    levels(p$Bin) <- c("<1", "1-10", "10-25", "25-50", "50-100", "100-500", "500-1000", "1000+")
    
    pdf(file = "RNA - All Ct Heatmap.pdf", height = 12, width = 15)
    ggplot(p, aes(y = Symbol, x = Stage, fill = Bin, label = round(value,0))) +
      geom_tile(colour = "black") +
      # geom_text(size = 2) +
      theme_bw() +
      theme(panel.grid = invis, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = invis, axis.ticks = invis) +
      labs() +
      scale_y_discrete(expand = c(0,0), limits = clust) +
      scale_x_discrete(expand = c(0,0)) +
      guides(fill = guide_legend(title = "CPM")) +
      facet_wrap(~Ct, nrow = 1, strip.position = "bottom") +
      scale_fill_manual(values = cols)
    dev.off()
        
    
    
  ## Next plot: just ast, and facet_wrap by trend!!
    p <- x[,grep("Ast", colnames(x))] %>% as.data.frame()
    clust <- hclust(dist(log2(p+0.5)))
    clust <- rownames(p)[clust$order]
    p$Symbol <- factor(rownames(p), levels = clust)
    m <- match(p$Symbol, herring_trends$Symbol)
    p$Trend <- herring_trends$Astro[m]
    p$Trend[is.na(p$Trend)] <- " No Trend"
    p$Trend <- gsub("trans_", "T", p$Trend)
    p$Trend <- splitter(p$Trend, "_", 1)
    # p$Trend <- gsub("_", "\n", p$Trend)
    # p$Trend <- gsub("nG", "n", p$Trend)
    
    
    p <- melt(p)
    p$Stage <- factor(splitter(p$variable, " ", 2), levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
    p$Bin <- cut(p$value, c(0, 1, 10, 25, 50, 100, 500, 1000, 10000), include.lowest = TRUE)
    levels(p$Bin) <- c("<1", "1-10", "10-25", "25-50", "50-100", "100-500", "500-1000", "1000+")
    
    pdf(file = "RNA - Astro Stage Comparison.pdf", height = 12, width = 6)
    ggplot(p, aes(y = Symbol, x = Stage, fill = Bin, label = round(value,0))) +
      geom_tile(colour = "black") +
      geom_text(size = 2) +
      theme_bw() +
      theme(panel.grid = invis, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = invis, axis.ticks = invis) +
      guides(fill = guide_legend(title = "CPM")) +
      labs() +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      # facet_wrap(~Trend, ncol = 1, strip.position = "right", scales = "free_y") +
      facet_grid(Trend~., space = "free_y", scales = "free_y", switch = "y") +
      scale_fill_manual(values = cols)
    dev.off()
   
    
## The visualisation raises a few interesting ideas to test statistically
  

## 1) ANOVA for stage within ast. And 1.5) lm for stage within Ast. 2) ANOVA for cell-type and stage.  

## Run expression statistics
  # get data
  e <- herring_pb$log2CPM
  e <- e[which(rownames(e) %in% res.final$Gene),] 
  e <- t(e)
    
  # note that metadata are in herring_pb$Meta
  m <- herring_pb$Meta
  m$Astro <- m$Ct == "Astro"
  m$StageNumeric <- as.numeric(m$Stage) # works because the latter is an already ordered factor

  # remove celltypes
  rem <- grep("Vas", m$Ct) # just vasculature
  e <- e[-rem,]
  m <- m[-rem,]
  
  ## Run stats!
  all.gene.stats <- data.frame(Gene = colnames(e),
                               Hit = colnames(e) %in% res.final$Gene[which(res.final$HitPermissive)],
                               ANOVA_CtBroad = NA,
                               ANOVA_Stage = NA,
                               ANOVA_StageNumeric = NA,
                               ANOVA_AstRegulation = NA,
                               LM_AstRegulation = NA,
                               LM_AstRegulation_Coef = NA,
                               LM_AstStageNumeric = NA,
                               ANOVA_AstStage = NA,
                               row.names = colnames(e))
  
  for (j in colnames(e)) {
    w <- which(all.gene.stats$Gene == j)
    print(j)
    
    # ANOVA for effect of numeric stage
    a <- lm(e[,j] ~ Ct_Broad*StageNumeric, data = m) %>% anova() # control for ct broad
    all.gene.stats$ANOVA_StageNumeric[w] <- a["StageNumeric", "Pr(>F)"]
    
    # ANOVA for effect of broad cell-type, and stage
    a <- lm(e[,j] ~ Ct_Broad*Stage, data = m) %>% anova() # control for ct broad
    all.gene.stats$ANOVA_CtBroad[w] <- a["Ct_Broad", "Pr(>F)"]
    all.gene.stats$ANOVA_Stage[w] <- a["Stage", "Pr(>F)"]
    
    # ANOVA for effect of ast vs other cell-types. this measures evidence for whether the mean in ast differs from the mean of other ct, controlling for stage
    a <- lm(e[,j] ~ Astro*Stage, data = m) %>% anova() # control for ct broad
    all.gene.stats$ANOVA_AstRegulation[w] <- a["Astro", "Pr(>F)"]
    
    # lm for effect of ast vs other cell-types. this measures evidence for whether the mean in ast differs from the mean of other ct, controlling for stage
    a <- lm(e[,j] ~ Astro*Stage, data = m) %>% summary() # control for ct broad
    all.gene.stats$LM_AstRegulation[w] <- a$coefficients["AstroTRUE", "Pr(>|t|)"]
    # all.gene.stats$LM_AstRegulation_Coef[w] <- a$coefficients["AstroTRUE", "Estimate"]
    
    # lm and anova for stage effect within Ast
    g <- grep("Astro", m$Ct_Broad)
    a <- lm(e[g,j] ~ StageNumeric, data = m[g,]) %>% summary()
    all.gene.stats$LM_AstStageNumeric[w] <- a$coefficients["StageNumeric", "Pr(>|t|)"]
    a <- lm(e[g,j] ~ Stage, data = m[g,]) %>% anova()
    all.gene.stats$ANOVA_AstStage[w] <- a["Stage", "Pr(>F)"]
  }
  
  # fdr adjustment of ANOVAs and the lm
  fdr.clmns <- colnames(all.gene.stats)[c(3:7, 9, 10)]
  for (k in fdr.clmns) all.gene.stats[,k] <- p.adjust(all.gene.stats[,k], method = "fdr")

  # numbers, and enrichments
  all.gene.stats.sum <- data.frame(Stat = fdr.clmns,
                                   TotalGenes = nrow(all.gene.stats),
                                   nSig = NA,
                                   nSigAndHit = NA,
                                   Fisher.OR = NA,
                                   Fisher.p = NA)
  
  for (j in fdr.clmns) {
    g <- grep(j, all.gene.stats.sum$Stat)
    k <- all.gene.stats[,j]
    thresh.fdr <- 0.05
    k <- k < thresh.fdr
    
    # number of significant on that test
    all.gene.stats.sum$nSig[g] <- length(which(k))
    all.gene.stats.sum$nSigAndHit[g] <- length(which(k & all.gene.stats$Hit))
    
    # fisher test
    f <- table(k, all.gene.stats$Hit) %>% fisher.test()
    all.gene.stats.sum$Fisher.OR[g] <- f$estimate
    all.gene.stats.sum$Fisher.p[g] <- f$p.value
  }
  
  write.csv(all.gene.stats, file = "RNA - All Gene-wise Stats.csv")
  write.csv(all.gene.stats.sum, file = "RNA - All Gene-wise Stats Summary.csv")
  
## For our ast stage analysis, concordance wth Lister DE
  ## Create dataframe
    within.stage.stats <- all.gene.stats[,c(1,2,9,10)]
    within.stage.stats$LM_AstStageNumeric <- within.stage.stats$LM_AstStageNumeric < 0.05
    within.stage.stats$ANOVA_AstStage <- within.stage.stats$ANOVA_AstStage < 0.05
    
    within.stage.stats$Trend <- herring_trends$Astro[match(within.stage.stats$Gene, herring_trends$Symbol)]
    within.stage.stats$HasTrend <- !(is.na(within.stage.stats$Trend))
    
  ## Fisher tests of concordance between the different methods of calling
    table(within.stage.stats$ANOVA_AstStage, within.stage.stats$HasTrend) %>% fisher.test()
    #         FALSE TRUE
    # FALSE  2322   38
    # TRUE    220  186
    # data:  .
    # p-value < 2.2e-16
    # alternative hypothesis: true odds ratio is not equal to 1
    # 95 percent confidence interval:
    #   22.42471 50.32286
    # sample estimates:
    #   odds ratio 
    # 33.16437
    
    table(within.stage.stats$LM_AstStageNumeric, within.stage.stats$HasTrend) %>% fisher.test() 
    #           FALSE TRUE
    # FALSE  2081  279
    # TRUE    117  289
    # data:  .
    # p-value < 2.2e-16
    # alternative hypothesis: true odds ratio is not equal to 1
    # 95 percent confidence interval:
    #   14.26586 23.82123
    # sample estimates:
    #   odds ratio 
    # 18.39798 
    
    table(within.stage.stats$LM_AstStageNumeric, within.stage.stats$ANOVA_AstStage) %>% fisher.test() 
    #           FALSE TRUE
    # FALSE  2322   38
    # TRUE    220  186
    # data:  .
    # p-value < 2.2e-16
    # alternative hypothesis: true odds ratio is not equal to 1
    # 95 percent confidence interval:
    #   35.11701 77.02796
    # sample estimates:
    #   odds ratio 
    # 51.49502 
    
    table(within.stage.stats$LM_AstStageNumeric, within.stage.stats$Hit) %>% fisher.test() 
    table(within.stage.stats$ANOVA_AstStage, within.stage.stats$Hit) %>% fisher.test() 
    table(within.stage.stats$HasTrend, within.stage.stats$Hit) %>% fisher.test() 
    
    
  ## Are hits likelier to have stage-specific regulation?
    stage.effects.hits <- within.stage.stats[which(within.stage.stats$Hit),]
    
    # number of significant stage-stage comparison methods
    stage.effects.hits$SignificantCalls <- rowSums(stage.effects.hits[,c(3,4,6)])
    
  ## Write
    write.csv(stage.effects.hits, file = "RNA - Stage-Specificty Method Concordance.csv")
    
  
  ## Err, can you reverse-engineer the FDR cut-off?
    # the issue is that not all genes without a significant pairwise comparison are omitted.
    # so I cannot determine the FDR cut-off for significant,
    
    # I will instead reverse-engineer it as follows:
      # assume that every gene must have at least one pairwise comparison below the FDR threshold. 
      # picked the lowest such value where every gene has at least one comparison below that point
    
    y <- herring_de$Astro$pval
    thresh <- max(rowMin(as.matrix(y))) # 0.00945... Not very stringent?
    
    # thus, I assume that any p-value below 0.00945 is significant
    y <- y <= thresh
    
    # look at the genes you've tested in the screen
    y <- y[which(rownames(y) %in% res.final$Gene),]
    table(rowSums(y), within.stage.stats$LM_AstStageNumeric[match(rownames(y), within.stage.stats$Gene)])
    table(rowSums(y), within.stage.stats$ANOVA_AstStage[match(rownames(y), within.stage.stats$Gene)])
  
## How about this for stage specificity: get the median of each stage, then compare medians?
  ##  Subset data to Ast
    e <- herring_pb$log2CPM
    e <- e[which(rownames(e) %in% res.final$Gene),] 
    e <- as.matrix(e)
      
    # note that metadata are in herring_pb$Meta
    m <- herring_pb$Meta
    m$Astro <- m$Ct == "Astro"
    m$StageNumeric <- as.numeric(m$Stage) # works because the latter is an already ordered factor
  
    # remove celltypes
    keep <- grep("Astro", m$Ct) # just vasculature
    e <- e[,keep]
    m <- m[keep,]
    
  ## Collect medians
    stage.medians <- list()
    for (j in levels(meta$Stage)) {
      stage.medians[[j]] <- rowMedians(e[,which(m$Stage == j)])
    }
    
    stage.medians <- do.call("cbind", stage.medians)
    rownames(stage.medians) <- rownames(e)
    stage.medians <- as.data.frame(stage.medians)
    
    
  ## Pick the max
    # stage.medians$MaxStage <- apply(stage.medians, 1, function(x) {
    #   w <- which.max(x)
    #   return(colnames(stage.medians)[w])
    # })
    
  ## Ratios
    x <- apply(stage.medians, 1, function(x) {
      w <- which.max(x)
      
      out <- data.frame(MaxStage = colnames(stage.medians)[w],
                        MaxLog2Exp = x[w],
                        Log2FCVsMedian = x[w] - median(x[-w]),
                        Log2FCVsRank2 = x[w] - max(x[-w]),
                        Log2FCVsLast = x[w] - min(x[-w]))
      
      return(out)
    })
  
    x <- do.call("rbind", x)
    stage.medians <- cbind(stage.medians, x)
    
    
  ## Call stage specificity
    stage.medians$SpecificVsMedian <- stage.medians$Log2FCVsMedian > 1
    stage.medians$SpecificVsRank2 <- stage.medians$Log2FCVsRank2 > 1
    
  ## Save
    write.csv(stage.medians, file = "RNA - Stage-specificity Median Calculation.csv")
    

################################################################################################################################ #
## Comparison of RNA to ATAC trends ----
      
    
## Combine data
  # get screen hit pairs
  x <- res.final[which(res.final$HitPermissive), 1:3]  
    
  # annotate with the timecourse atac astrocyte data
  m <- match(x$Enh, rownames(all.enh.stats$StageFC))
  x$ATAC_HighestStage <- all.enh.stats$StageFC$MaxStage[m]
  x$ATAC_HighestStageFCvsRank2 <- all.enh.stats$StageFC$RatioVsRank2[m]
    
    
  # annotate with the timecourse rna astrocyte data
    
    
  
## Quick check of stage-specificity
  # filter data
  x <- res.final[which(res.final$HitPermissive),1:3]
  
  # add ATAC stage
  m <- match(x$Enh, rownames(all.enh.stats$StageFC))
  
  

################################################################################################################################ #
## CREs ----
    
## Briefly, CREs were determined by:
  # pick 400 ATAC nuclei at random
  # pseudobulk them with their nearest 50 ATAC neighbours.
  # Use Seurat to integrate with RNA data and find 50 nearest RNA neighbours to the ATAC seed neighbour. Pseudobulk these.
  # Correlate peak and RNA, for all genes within 250kb of peak. 
  
    
## Read in  
  cre <- readxl::read_xlsx("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/SuppTable5_CREs.xlsx", sheet = 1, skip = 9) %>% as.data.frame()
    
  # convert ensid to symbol (using their annotation)
  m <- match(cre$Gene, pb$GeneInfo$gene_ids)  
  cre$Gene <- pb$GeneInfo$X[m]
    
  # convert to bed
  cre <- cre[,c(5:7,9,1,2)]
  write.bed(cre, "../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/CREs.bed")

    
## Intersect with our candidate peaks (in hg19)
  call <- paste("intersectBed",
                "-a", "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed",
                "-b", "../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/CREs.bed",
                "-wa", "-wb", 
                ">", "CRE_Overlap.bed")
  
    system(call, intern = FALSE, wait = TRUE)  
    
  call <- paste("windowBed",
                "-a", "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed",
                "-b", "../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/CREs.bed",
                "-w", "1000", 
                ">", "CRE_Overlap_w1000.bed")
  
  system(call, intern = FALSE, wait = TRUE)  

  
## Read in
  x <- read.delim("CRE_Overlap.bed", header = FALSE)
  
## Annotate
  # just Ast
  # x <- x[which(x$V10 == "Astro"),]
  
  # screen enhancer id
  m <- match(x$V4, guides$TargetCoord)
  x$ScreenEnh <- guides$TargetID[m]
  
  # screen target gene
  x$ScreenHit <- x$ScreenEnh %in% res.final$Enh[which(res.final$HitPermissive)]
  x$ScreenGene <- NA
  
  for (j in 1:nrow(x)) {
    if (x$ScreenHit[j]) {
      y <- res.final[which(res.final$Enh == x$ScreenEnh[j] & res.final$HitPermissive),]
      x$ScreenGene[j] <- paste(y$Gene, collapse = "/")
    }
  }
  
## Output
  x <- x[,c("ScreenEnh", "ScreenHit", "ScreenGene","V7","V8","V9","V10","V11")]
  colnames(x)[4:8] <- c("HerringChr", "HerringStart", "HerringEnd", "HerringCt", "HerringGene")
  write.csv(x, file = "CRE Annotated.csv")
  
  
 x <- read.delim("CRE_Overlap.bed", header = FALSE)
  
## Now look at the window
  # just Ast
  # x <- x[which(x$V10 == "Astro"),]
  
  # screen enhancer id
  m <- match(x$V4, guides$TargetCoord)
  x$ScreenEnh <- guides$TargetID[m]
  
  # screen target gene
  x$ScreenHit <- x$ScreenEnh %in% res.final$Enh[which(res.final$HitPermissive)]
  x$ScreenGene <- NA
  
  for (j in 1:nrow(x)) {
    if (x$ScreenHit[j]) {
      y <- res.final[which(res.final$Enh == x$ScreenEnh[j] & res.final$HitPermissive),]
      x$ScreenGene[j] <- paste(y$Gene, collapse = "/")
    }
  }
  
## Output
  x <- x[,c("ScreenEnh", "ScreenHit", "ScreenGene","V7","V8","V9","V10","V11")]
  colnames(x)[4:8] <- c("HerringChr", "HerringStart", "HerringEnd", "HerringCt", "HerringGene")
  write.csv(x, file = "CRE Annotated.csv")
  
  
################################################################################################################################ #
## Recalculate astrocyte markers using new pseudobulk ----
  
## Load
  # pseudobulk, called pb
  load("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/Pseudobulk_byGJS.rda")
  
  # Seurat object, called brainMat
  load("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/scRNAseq_UMIs_Seurat.rda")
  
  
# ## Process the Seurat data
#   ## Adult
#     # subset
#     adult <- subset(x = brainMat, cells = which(brainMat$stage_id %in% c("Adolescence", "Adult")))
#     adult <- subset(x = adult, cells = which(adult$cell_type != "Poor-Quality"))
#     
#     # cell-type identity
#     adult$Ct <- adult$cell_type
#     w <- which(adult$Ct == "Non-Neu")
#     adult$Ct[w] <- adult$major_clust[w]  
#     Idents(adult) <- "Ct"
#     
#     # for plotting
#     adult <- NormalizeData(adult)
#     adult <- FindVariableFeatures(object = adult, selection.method = "vst", nfeatures = 2000)
#     adult <- ScaleData(object = adult, vars.to.regress = c("batch", "log1p_total_counts", "percent_mito", "percent_ribo")) 
#     adult <- RunPCA(adult)
#     ElbowPlot(adult)
#     adult <- clust.bySeurat(adult, method = "UMAP", dim = 15, resolution = 1)
#     Idents(adult) <- "Ct"
#     DimPlot(adult)
#     
#   ## Foetal
#     # subset
#     foetal <- subset(x = brainMat, cells = which(brainMat$stage_id == "Fetal"))
#     foetal <- subset(x = foetal, cells = which(foetal$cell_type != "Poor-Quality"))
#     
#     # cell-type identity
#     foetal$Ct <- foetal$cell_type
#     w <- which(foetal$Ct == "Non-Neu")
#     foetal$Ct[w] <- foetal$major_clust[w]  
#     Idents(foetal) <- "Ct"
#     
#     # for plotting
#     foetal <- NormalizeData(foetal)
#     foetal <- FindVariableFeatures(object = foetal, selection.method = "vst", nfeatures = 2000)
#     foetal <- ScaleData(object = foetal, vars.to.regress = c("batch", "log1p_total_counts", "percent_mito", "percent_ribo")) 
#     foetal <- RunPCA(foetal)
#     ElbowPlot(foetal)
#     foetal <- clust.bySeurat(foetal, method = "UMAP", dim = 15, resolution = 1)
#     Idents(foetal) <- "Ct"
#     DimPlot(foetal)
#   
#   ## Save
#     save(foetal, adult, file = "Seurat Objects.rda")
    
## Approach 1: linear model versus all other samples in pseudobulk
  # function
    ## Get samples
      use.samps <- which(pb$Meta$nCells >= 10 & pb$Meta$Stage %in% c("Adolescence", "Adult"))
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
      y$Age <- splitter(y$OriginalID, "_", 2) %>% sub("yr", "", .) %>% as.numeric()
      y$log10UMI <- log10(y$nUMI)
      
    ## Linear model
      # run
      z <- lm(x~IsAstrocyte+Individual+Age+log10UMI, data = y)
      z <- summary(z)
      
      # extract effect of astrocytes
      z <- lapply(z, function(z) {
        data.frame(P = as.numeric(z$coefficients["IsAstrocyteTRUE", "Pr(>|t|)"]),
                   log2fc = as.numeric(z$coefficients["IsAstrocyteTRUE", "Estimate"]))
      })
    
      # return
      z <- do.call("rbind", z)
      rownames(z) <- splitter(rownames(z), " ", 2)
      z$FDR <- p.adjust(z$P, method = "fdr")
      z$AstMarker <- z$FDR < 0.05 & z$log2fc > 1
      z$Gene <- rownames(z)
      z <- z[,c(5,2,1,3,4)]  
      write.csv(z, file = "Astrocyte Markers - Pseudobulk, Adult and Adolescent.csv")
      
    ## Hit gene enrichment
      astMarkers <- data.frame(Gene = unique(res.final$Gene))
      astMarkers$InLister <- astMarkers$Gene %in% rownames(pb$Exp)
      astMarkers <- astMarkers[which(astMarkers$InLister),]
      astMarkers$AstMarker <- astMarkers$Gene %in% z$Gene[which(z$AstMarker)]
      astMarkers$Hit <- astMarkers$Gene %in% res.final$Gene[which(res.final$HitPermissive)]
      
      table(astMarkers$AstMarker, astMarkers$Hit)  %>% fisher.test()
      
      # p-value = 0.000109
      # alternative hypothesis: true odds ratio is not equal to 1
      # 95 percent confidence interval:
      #   1.575038 4.007459
      # sample estimates:
      #   odds ratio 
      # 2.541878 
    

#     
# ## Approach 2: Seurat's FindMarkers on single cells
#   gene.subset <- rownames(adult)[rownames(adult) %in% res.final$Gene]
#     
#   ## Adult samples
#     # run pairwise
#     fm.adult <- list()
#     
#     for (j in unique(adult$Ct)) {
#       print(j)
#       if (j == "Astro") next
#       
#       fm.adult[[paste0("Ast_vs_", j)]] <- FindMarkers(adult, ident.1 = "Astro", ident.2 = j,
#                                                       logfc.threshold = 0,
#                                                       features = gene.subset)
#     }
#     
#     save(fm.adult, file = "Astrocyte Markers - Seurat, Adult Pairwise, Screened Genes.rda")
#     
#     
#     
#     # combine ranked lists
#     
#   ## 
# 
#     
#     
#     # simplest process, comparing astrocytes to all other cells
#     test.fm <- FindMarkers(adult, ident.1 = "Astro")
#     write.csv(test.fm, "Astrocyte Markers - Seurat, Adult Ast vs Adult Other.csv")
#     test.fm.filt <- test.fm[which(test.fm$avg_log2FC > 0),]
#     test.fm.filt <- test.fm[which(test.fm$avg_log2FC > 1),]
#     
#   ## Foetal samples
#     foetal <- subset(x = brainMat, cells = which(brainMat$stage_id == "Fetal"))
#     foetal$Ct <- foetal$cell_type
#     w <- which(foetal$Ct == "Non-Neu")
#     foetal$Ct[w] <- foetal$major_clust[w]  
#     Idents(foetal) <- "Ct"
#     
#     foetal <- NormalizeData(foetal)
#     
#     
#     # simplest process, comparing astrocytes to all other cells
#     test.fm.foetal <- FindMarkers(foetal, ident.1 = "Astro")
#     write.csv(test.fm.foetal, "Astrocyte Markers - Seurat, Foetal Ast vs Foetal Other.csv")
#     test.fm.foetal.filt <- test.fm.foetal[which(test.fm.foetal$avg_log2FC > 1),]
#     
  
################################################################################################################################ #
## Visualisation of single-cell data ----
    
    
## For genes with developmental trajectories in astrocytes
    
## Get Seurat data
  astro_dev <- subset(x = brainMat, cells = which(brainMat$major_clust == "Astro"))
  Idents(astro_dev) <- "batch"
    
  # # for plotting
  # astro_dev <- NormalizeData(astro_dev)
  # astro_dev <- FindVariableFeatures(object = astro_dev, selection.method = "vst", nfeatures = 2000)
  # astro_dev <- ScaleData(object = astro_dev, vars.to.regress = c("batch", "log1p_total_counts", "percent_mito", "percent_ribo")) 
  # astro <- RunPCA(astro_dev)
  # ElbowPlot(astro, ndims = 50)
  # astro <- clust.bySeurat(astro, method = "UMAP", dim = 20, resolution = 1)
  # DimPlot(astro)
  # DimPlot(astro, group.by = "batch")
  
## Get pb data
  astro_dev_pb <- pb$Exp
  astro_dev_pb <- apply(astro_dev_pb, 2, function(x) x / (sum(x) / 10^6))
  astro_dev_pb <- log2(astro_dev_pb + 1)
  
## Move on to plotting
  library(ggbeeswarm)
  library(cowplot)
  
  # function
  plot.gene.devTrend <- function(gene, check.trend = TRUE) {
    # a trend label
    if (check.trend) {
      m <- match(gene, ast_trends$gene_name)
      gene.trend <- ast_trends$trend_class[m]
      ylab <- paste0(gene, " (", gene.trend, ")")
    } else {
      ylab <- gene
    }
    
    ## Single cell
      a <- data.frame(Exp = as.numeric(astro_dev@assays$RNA@data[gene,]),
                      Stage = astro_dev$stage_id,
                      Individual = astro_dev$RL.,
                      Age = astro_dev$numerical_age)
      a$Stage <- factor(a$Stage,levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
      
      # to order individuals...
      x <- unique(a[,c("Individual", "Age")])
      x <- x[order(x$Age),]
      a$Individual <- factor(a$Individual, levels = x$Individual)
      
      # plot
      pA <- ggplot(a, aes(x = Individual, y = Exp, colour = Stage, fill = Stage)) +
        # geom_violin(scale = "width", alpha = 0.2, colour = "white") +
        geom_quasirandom(size = 0.5, alpha = 0.5) +
        stat_summary(fun = "median", shape = "-", colour = "black", size = 2) +
        theme_bw() +
        scale_color_lancet() +
        scale_fill_lancet() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.border = invis, axis.line.y = element_line(),
              panel.grid.major.x = invis, legend.position = "none") +
        scale_y_continuous(expand = c(0,0), limits = c(0, (max(a$Exp)*1.05))) +
        labs(y = paste0(ylab, "\nSingle-cell Expression"))
      
    ## Pseudobulk
      w <- which(pb$Meta$Celltype == "Astro")
      b <- data.frame(Exp = astro_dev_pb[gene,w],
                      Stage = pb$Meta$Stage[w],
                      Individual = pb$Meta$Individual[w])
      b$Stage <- factor(b$Stage,levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
      b$Individual <- factor(b$Individual, levels = levels(a$Individual))
      
      pB <- ggplot(b, aes(x = Stage, y = Exp, colour = Stage, fill = Stage)) +
        # geom_violin(scale = "width", alpha = 0.5) +
        geom_quasirandom(size = 2, alpha = 1) +
        stat_summary(fun = "median", shape = "-", colour = "black", size = 2, show.legend = FALSE) +
        theme_bw() +
        scale_color_lancet() +
        scale_fill_lancet() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.border = invis, axis.line.y = element_line(),
              panel.grid.major.x = invis, axis.title.x = invis) +
        scale_y_continuous(expand = c(0,0), limits = c(0, (max(b$Exp)*1.1))) +
        labs(y = paste0(ylab, "\nPseudobulk Expression"))
      
      plot_grid(pA, pB, ncol = 1, rel_widths = c(1,1))
      
    
  }

  # apply for genes with trends in astrocytes
  ast_trends <- readxl::read_xlsx("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/SuppTable3_DE.xlsx", sheet  = 1, skip = 2)
  
  use.genes <- intersect(res.final$Gene[which(res.final$HitPermissive)],
                         ast_trends$gene_name)
  
  pdf(file = paste0("Gene_Expression_Plots/Trend Within Ast (DE Only).pdf"), height = 6, width = 8)  
  for (j in use.genes) {
    print(j)
    print(plot.gene.devTrend(j, check.trend = TRUE))
  }
  dev.off()
  
  
  use.genes <- intersect(res.final$Gene[which(res.final$HitPermissive)],
                       rownames(astro_dev_pb))

  pdf(file = paste0("Gene_Expression_Plots/Trend Within Ast.pdf"), height = 6, width = 8)  
  for (j in use.genes) {
    print(j)
    print(plot.gene.devTrend(j, check.trend = FALSE))
  }
  dev.off()
  

  
## For genes that mark astrocytes
  pb$log2CPM <- pb$Exp
  pb$log2CPM <- apply(pb$log2CPM, 2, function(x) x / (sum(x) / 10^6))
  pb$log2CPM <- log2(pb$log2CPM + 1)
  
    # function
  plot.gene.marker <- function(gene, check.mark = TRUE) {
    if (check.mark) {
      m <- match(gene, pb_lm_allSamps$Gene)
      mark.stats <- paste0("log2fc=", round(pb_lm_allSamps$log2fc[m], 2), ", FDR=", signif(pb_lm_allSamps$FDR[m], 1))
      ylab <- paste0(gene, " (", mark.stats, ")")
    } else {
      ylab <- gene
    }  
    
    
    b <- data.frame(Exp = pb$log2CPM[gene,],
                      Stage = pb$Meta$Stage,
                      Ct = pb$Meta$Celltype,
                      Individual = pb$Meta$Individual)
    b <- b[which(pb$Meta$nCells >= 10),]
    b$Stage <- factor(b$Stage,levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
      
    ggplot(b, aes(x = Ct, y = Exp, colour = Stage)) +
      # geom_violin(scale = "width", alpha = 0.5) +
      geom_quasirandom(size = 2, alpha = 1, dodge.width = 0.75, width = 0.1) +
      stat_summary(fun = "median", shape = "-", colour = "black", size = 2, show.legend = FALSE) +
      theme_bw() +
      facet_wrap(~Ct, scales = "free_x", nrow = 1, strip.position = "bottom") +
      scale_color_carto_d(palette = "Temps") +
      theme(panel.border = invis, axis.line.y = element_line(),
            panel.grid.major.x = invis, axis.title.x = invis, axis.text.x = invis, axis.ticks.x = invis) +
      scale_y_continuous(limits = c(0,NA)) +
      labs(y = paste0(ylab, "\nPseudobulk Expression (log2(CPM+1))"))
  }
  
  
  
  use.genes <- intersect(res.final$Gene[which(res.final$HitPermissive)],
                         pb_lm_allSamps$Gene[which(pb_lm_allSamps$AstMarker)])
  
  pdf(file = paste0("Gene_Expression_Plots/Trend Across Cell-types (Ast Markers Only).pdf"), height = 3, width = 8)  
  for (j in use.genes) {
    print(j)
    print(plot.gene.marker(j, check.mark = TRUE))
  }
  dev.off()
  
  # and for genes not marking ast
  use.genes <- intersect(res.final$Gene[which(res.final$HitPermissive)],
                           pb_lm_allSamps$Gene)
    
    pdf(file = paste0("Gene_Expression_Plots/Trend Across Cell-types.pdf"), height = 3, width = 8)  
    for (j in use.genes) {
      print(j)
      print(plot.gene.marker(j, check.mark = TRUE))
    }
    dev.off()
  
    