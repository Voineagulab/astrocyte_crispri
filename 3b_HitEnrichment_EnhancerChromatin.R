## This script annotates hit enhancers by overlap to "interesting" coordinates from various resources, including:

  # enhancer annotations across cell-types

  # regions of different accessibility in various brain disorders

  # ctcf

  # TADs

  # transcription factors

  # ttseq

  # etc

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/")
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
  library(gplots)
  library(ggbeeswarm)
  library(ggsci)

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
  loj <- function(a, b, out) {
    call <- paste("intersectBed",
                "-a", a,
                "-b", b,
                "-loj", 
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
  
  
## Commonly used paths
  # candidate coordinates
  nha_dir_38 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed" 
  nha_dir_19 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19.bed" 
  
  # and for speed, read in existing annotations
  # load("Final.rda")

################################################################################################################################ #
## Cell-type specificity: coverage at our peaks ----
  
    
## For these analyses, I will be quantifying the activity of each of our peaks using publically-available (hg19) BigWig files
## This negates the issue of having to overlap binary peak calls
    
    
## Load required libraries
    library(Rsamtools)
    library(rtracklayer)
    
    
## First, create a windowed version of the NHA peaks
  ## Read in
    x <- read.delim(nha_dir_19, header = FALSE)
    
  ## Extend by 1000bp 
    window <- 1000
    x$V2 <- x$V2 - (window/2)
    x$V3 <- x$V3 + (window/2)
    
  ## Save
    nha_dir_19_w1000 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19_window1000.bed"
    write.bed(x, nha_dir_19_w1000)
    
    
## Process dataset 1: Herring 2022
  # this is snATAC-seq from the human brain across maturation (foetal to early adulthood)
    
  ## Directories
    bed_nha_hg19_file <- import.bed(nha_dir_19)
    bed_nha_hg19_w1000_file <- import.bed(nha_dir_19_w1000)
    bw_herring <- list.files("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig/") # this is hg19
    bw_herring <- bw_herring[grep("bigwig", bw_herring)]
    bw_herring <- paste0("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig/", bw_herring)
    

  ## Calculate coverage in each of the snATAC-seq bigwigs
    covHerring <- list()
    
    covHerring$window0 <- lapply(bw_herring, function(x) { # requires ~20s
      print(x)  
      
      # get coverage
      y <- import(x, selection = bed_nha_hg19_file, as = "NumericList")
      z <- sapply(y, sum) # this is the sum of coverages across the bases, and is not normalised to peak width. this makes comparisons between peaks difficult
      
      # collect peak ids
      names(z) <- bed_nha_hg19_file$name[-grep("chrX", bed_nha_hg19_file$name)] # $name is the GRCh38 coordinate! this notes that anything on ChrX dropped from the analyses as the bigwig lacks them
      
      return(z)
    })
    
    covHerring$window1000 <- lapply(bw_herring, function(x) { # requires ~20s
      print(x)  
      
      # get coverage
      y <- import(x, selection = bed_nha_hg19_w1000_file, as = "NumericList")
      z <- sapply(y, sum) # this is the sum of coverages across the bases, and is not normalised to peak width. this makes comparisons between peaks difficult
      
      # collect peak ids
      names(z) <- bed_nha_hg19_file$name[-grep("chrX", bed_nha_hg19_file$name)] # $name is the GRCh38 coordinate! this notes that anything on ChrX dropped from the analyses as the bigwig lacks them
      
      return(z)
    })
    
    
  ## Process the coverage file
    covHerring <- lapply(covHerring, function(x) {
      names(x) <- splitter(bw_herring, "BigWig\\/", 2) %>% sub("_atac_insertions.bigwig", "", .)
      x <- do.call("cbind", x)
      
      # get the enhancer id
      m <- match(rownames(x), guides$TargetCoord)
      rownames(x) <- guides$TargetID[m]
      
      # rename columns
      colnames(x) <- gsub("_der", "", colnames(x)) %>%
        gsub("L5_6_", "L56_", .) %>% 
        gsub("L2_3_", "L23_", .)
      
      neuro_cols <- grep("^L|^MGE|^CGE", colnames(x))
      colnames(x)[neuro_cols] <- paste0("Neuro_", colnames(x)[neuro_cols])
      
      # reorder by stage
      colnames(x) <- gsub("Fetal", "1Fetal", colnames(x)) %>%
        gsub("Neonatal", "2Neonatal", .) %>% 
        gsub("Infancy", "3Infancy", .) %>%
        gsub("Childhood", "4Childhood", .) %>% 
        gsub("Adolescence", "5Adolescence", .) %>% 
        gsub("Adult", "6Adult", .) 
      
      x <- x[,order(colnames(x))]
      
      # return
      return(x)
      
    })
    
    
  ## Save
    save(covHerring, file = "Coverage/Herring - Coverage Matrices.rda") 
    write.csv(covHerring$window0, "Coverage/Herring - Coverage Matrix (Window 0).csv")
    write.csv(covHerring$window1000, "Coverage/Herring - Coverage Matrix (Window 1000).csv")
  
    
## Process dataset 2: Nott 2019
  # these are antibody-sorted nuclei from the human brain, derived from samples from childhood and adolescence (4-18) 
  
  ## Directories
    bw_nott <- list.files("../../../../PublicData/Nott2019/BigWigs/") # this is hg19
    bw_nott <- bw_nott[-grep("README", bw_nott)]
    bw_nott <- paste0("../../../../PublicData/Nott2019/BigWigs/", bw_nott)
    

  ## Calculate coverage in each of the snATAC-seq bigwigs
    covNott <- list()
    
    covNott$window0 <- lapply(bw_nott, function(x) { # requires ~20s
      print(x)  
      
      # get coverage
      y <- import(x, selection = bed_nha_hg19_file, as = "NumericList")
      z <- sapply(y, sum) # this is the sum of coverages across the bases, and is not normalised to peak width. this makes comparisons between peaks difficult
      
      # collect peak ids
      names(z) <- bed_nha_hg19_file$name # $name is the GRCh38 coordinate! 
      
      return(z)
    })
    
    covNott$window1000 <- lapply(bw_nott, function(x) { # requires ~20s
      print(x)  
      
      # get coverage
      y <- import(x, selection = bed_nha_hg19_w1000_file, as = "NumericList")
      z <- sapply(y, sum) # this is the sum of coverages across the bases, and is not normalised to peak width. this makes comparisons between peaks difficult
      
      # collect peak ids
      names(z) <- bed_nha_hg19_file$name # $name is the GRCh38 coordinate! 
      
      return(z)
    })
    
  ## Process
     covNott <- lapply(covNott, function(x) {
      names(x) <- splitter(bw_nott, "human_", 2) %>% 
        splitter(., "_epilepsy", 1) %>% 
        gsub("nuclei", "", .)
      x <- do.call("cbind", x)
      
      # get the enhancer id
      m <- match(rownames(x), guides$TargetCoord)
      rownames(x) <- guides$TargetID[m]
      
      # rename markers to cell-types
      colnames(x) <- gsub("LHX2", "Astro", colnames(x)) %>%
        gsub("NEUN", "Neuro", .) %>% 
        gsub("OLIG2", "Oligo", .) %>%
        gsub("PU1", "Micro", .) 
      
      # return
      return(x)
      
    })
     
 
  ## Save
    save(covNott, file = "Coverage/Nott - Coverage Matrices.rda") 
    write.csv(covNott$window0, "Coverage/Nott - Coverage Matrix (Window 0).csv")
    write.csv(covNott$window1000, "Coverage/Nott - Coverage Matrix (Window 1000).csv")
    
    
## Compare coverages!
  
  ## First, window versus non-window
    ## Herring
      x <- melt(covHerring$window0)
      y <- melt(covHerring$window1000)
      z <- data.frame(Enh = x$Var1,
                      Sample = x$Var2,
                      w0 = log2(x$value + 1),
                      w1000 = log2(y$value + 1))
      z$Celltype <- splitter(z$Sample, "_", 1)
      z$Stage <- gsub("Neuro_", "", z$Sample) %>% splitter(., "_", 2) %>% substr(., 2, 1000) 
      z$Stage <- factor(z$Stage, levels = c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))
      
      lim <- max(c(z$w0, z$w1000))
      
      pdf(file = "Coverage/Herring - Window Comparison.pdf", height = 10, width = 10)
      ggplot(z, aes(x = w0, y = w1000)) +
        geom_point() +
        facet_grid(Stage~Celltype) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
        scale_x_continuous(limits = c(0, lim), breaks = c(0, 5, 10)) +
        scale_y_continuous(limits = c(0, lim), breaks = c(0, 5, 10)) +
        theme_bw() +
        theme() +
        labs(x = "log2 Sum Coverage at Peak", y = "log2 Sum Coverage at Peak ±500bp")
      dev.off()
      
    ## Nott
      x <- melt(covNott$window0)
      y <- melt(covNott$window1000)
      z <- data.frame(Enh = x$Var1,
                      Sample = x$Var2,
                      w0 = log2(x$value + 1),
                      w1000 = log2(y$value + 1))
      z$Celltype <- splitter(z$Sample, "_", 1)
      z$Protocol <- splitter(z$Sample, "_", 2) 
    
      lim <- max(c(z$w0, z$w1000))
      
      pdf(file = "Coverage/Nott - Window Comparison.pdf", height = 5, width = 10)
      ggplot(z, aes(x = w0, y = w1000)) +
        geom_point() +
        facet_grid(Protocol~Celltype) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
        scale_x_continuous(limits = c(0, lim), breaks = c(0, 5, 10)) +
        scale_y_continuous(limits = c(0, lim), breaks = c(0, 5, 10)) +
        theme_bw() +
        theme() +
        labs(x = "log2 Sum Coverage at Peak", y = "log2 Sum Coverage at Peak ±500bp")
      dev.off()
    
  ## And ATAC versus H3K27ac
      # plot
      y <- melt(covNott$window1000)
      y$Celltype <- splitter(y$Var2, "_", 1)
      y$Protocol <- splitter(y$Var2, "_", 2) 
      # y <- y[,3:5]
      z <- dcast(y[,-2], Celltype+Var1~Protocol)
      z$atac <- log2(z$atac + 1)
      z$H3K27ac <- log2(z$H3K27ac + 1)
      
      lim <- summary(c(z$atac, z$H3K27ac))[c(1,6)]
      
      pdf(file = "Coverage/Nott - Protocol Comparison.pdf", height = 5, width = 5)
      ggplot(z, aes(x = atac, y = H3K27ac)) +
        geom_point() +
        facet_wrap(~Celltype) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
        scale_x_continuous(limits = lim, breaks = c(0, 5, 10)) +
        scale_y_continuous(limits = lim, breaks = c(0, 5, 10)) +
        theme_bw() +
        theme() +
        labs(x = "ATAC\nlog2 Sum Coverage at Peak±500bp", y = "H3K27ac\nlog2 Sum Coverage at Peak ±500bp")
      dev.off()
    
    
  ## And Herring versus Nott
    # bind
    x <- covHerring$window1000
    colnames(x) <- paste0("H_", colnames(x))
    y <- covNott$window1000[rownames(x),]
    colnames(y) <- paste0("N_", colnames(y))
    x <- cbind(x, y)
    
    # correlate
    cor <- cor(x, method = "s")
    
    # heatmap
    csc <- splitter(colnames(cor), "_", 1) %>% 
      gsub("H", "firebrick1", .) %>% 
      gsub("N", "dodgerblue4", .)
  
    pdf(file = "Coverage/Publication Comparison (Spearman).pdf", height = 10, width = 10)
    heatmap.2(cor, trace = "none", col = carto_pal(n = 7, name = "Tropic"), margins = c(10, 10), symbreaks = FALSE, ColSideColors = csc)
    dev.off()
    
    # compare the depth
    # z <- colnames(x)
    # z <- data.frame(Sample = z,
    #                 Study = splitter(z, "_", 1),
    #                 Ct = splitter(z, "_", 2),
    #                 MedianCoverage = apply(x, 2, median))
    # ggplot(z, aes(x = Ct, y = MedianCoverage, colour = Study)) +
    #   geom_quasirandom()
    
    z <- melt(x)
    z$Study <- splitter(z$Var2, "_", 1)
    z$Ct <- splitter(z$Var2, "_", 2)
    z$Protocol <- "ATAC"; z$Protocol[grep("H3K27ac", z$Var2)] <- "H3K27ac"
    z$value <- (z$value + 1)
    
    pdf(file = "Coverage/Publication Comparison (Distribution).pdf", height = 3.5, width = 5)
    ggplot(z, aes(x = Ct, y = value, linetype = Protocol, fill = Study)) +
      # geom_quasirandom(dodge.width = 0.5) +
      geom_violin(position = position_dodge(width = 0.75), width = 0.7, scale = "width", draw_quantiles = 0.5, colour = "black") +
      scale_y_continuous(trans = "log2", breaks = c(0, 8, 32, 128, 512, 2048, 8192, 32768)) +
      theme_bw() +
      scale_fill_manual(values = c("firebrick1", "dodgerblue4")) +
      labs(x = "Celltype", y = "Coverage at Autosomal Enhancers (log2 axis)")
    dev.off()
    
    
    
## Visualise
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
                        RowSideColors = rsc, ColSideColors = csc, col = carto_pal(7, "Geyser"), main = ttl)
      
    }
    
    # run
    pdf(file = "Coverage/Herring - Heatmaps.pdf", height = 10, width = 10)
    heatmap.scaled(covHerring$window0, do.scale = TRUE, ttl = "Window 0, Row-wise Scaling")
    heatmap.scaled(covHerring$window0, do.scale = FALSE, ttl = "Window 0, No Scaling")
    heatmap.scaled(covHerring$window1000, do.scale = TRUE, ttl = "Window 1000, Row-wise Scaling")
    heatmap.scaled(covHerring$window1000, do.scale = FALSE, ttl = "Window 1000, No Scaling")
    dev.off()
     
    atac <- grep("atac", colnames(covNott$window0))
    
    pdf(file = "Coverage/Nott - Heatmap ATAC.pdf", height = 10, width = 6)
    heatmap.scaled(covNott$window0[,atac], do.scale = TRUE, ttl = "Window 0, Row-wise Scaling")
    heatmap.scaled(covNott$window0[,atac], do.scale = FALSE, ttl = "Window 0, No Scaling")
    heatmap.scaled(covNott$window1000[,atac], do.scale = TRUE, ttl = "Window 1000, Row-wise Scaling")
    heatmap.scaled(covNott$window1000[,atac], do.scale = FALSE, ttl = "Window 1000, No Scaling")
    dev.off()
    
    pdf(file = "Coverage/Nott - Heatmap H3K27ac.pdf", height = 10, width = 6)
    heatmap.scaled(covNott$window0[,-atac], do.scale = TRUE, ttl = "Window 0, Row-wise Scaling")
    heatmap.scaled(covNott$window0[,-atac], do.scale = FALSE, ttl = "Window 0, No Scaling")
    heatmap.scaled(covNott$window1000[,-atac], do.scale = TRUE, ttl = "Window 1000, Row-wise Scaling")
    heatmap.scaled(covNott$window1000[,-atac], do.scale = FALSE, ttl = "Window 1000, No Scaling")
    dev.off()
    
    
## Temporal trends within astrocytes
  ## Get data
    x <- covHerring$window1000
    x <- x[,grep("Ast", colnames(x))]
    colnames(x) <- splitter(colnames(x), "_", 2) %>% substr(2, 100)
    
  ## A visualisation
    # scale
    y <- apply(x, 1, function(y) y / max(y) )
    y <- t(y)
    y <- y[-which(is.na(y[,1])),]
    
    # plot
    pdf(file = "Coverage/Temporal - Heatmap Scaled.pdf", height = 10, width = 6)
    heatmap.2(y, 
              trace = "none", col = viridis_pal()(10),
              Colv = FALSE, dendrogram = "row", margins = c(10,5))
    dev.off()
    
  ## Statistics
    # as there is but one replicate per stage, calculate specificity using fold-change
    calc.fc <- function(covMx) {
     
       out <- apply(covMx, 1, function(a) {
        
        # highest coverage stage
        w <- which.max(a)
        
        # output
        out <- data.frame(MaxID = NaN,
                          Coverage = a[w],
                          RatioVsMedian = a[w] / median(a[-w]),
                          RatioVsRank2 = a[w] / max(a[-w]),
                          RatioVsLast = a[w] / min(a[-w]))
        
        out <- round(out, 2)
        out$MaxID <- colnames(covMx)[w]
        return(out)
        
      })
      
      out <- do.call("rbind", out)
      
    }
    
    temporal.fc <- calc.fc(x)
    
    # call as stage specific if fc to second highest stage is > 2
    temporal.fc$Specific <- temporal.fc$RatioVsRank2 >= 2
    
  ## Save
    write.csv(temporal.fc, file = "Coverage/Temporal - Foldchanges.csv")
    
    
## Cell-type trends in Nott
  ## Calculate in Nott
    x <- covNott$window1000

    # as there is but one replicate per stage, calculate specificity using fold-change
    
    # calculate for atac and h3k27ac separately
    ct.fc.atac <- calc.fc(x[,grep("atac", colnames(x))])
    ct.fc.atac$Specific <- ct.fc.atac$RatioVsRank2 > 2
    colnames(ct.fc.atac) <- paste0("ATAC_", colnames(ct.fc.atac))
    
    ct.fc.H3K27ac <- calc.fc(x[,grep("H3K27ac", colnames(x))])
    ct.fc.H3K27ac$Specific <- ct.fc.H3K27ac$RatioVsRank2 > 2
    colnames(ct.fc.H3K27ac) <- paste0("H3K27ac_", colnames(ct.fc.H3K27ac))
    
  ## Cell-type trends in Herring
    x <- covHerring$window1000
    y <- list()
    
    for (j in c("Astro", "Neuro", "Micro", "Oligo")) {
      y[[j]] <- rowMeans(x[,grep(j, colnames(x))])
    }
    
    y <- do.call("cbind", y)
    
    ct.fc.herring <- calc.fc(y)
    ct.fc.herring$Specific <- ct.fc.herring$RatioVsRank2 > 2
    colnames(ct.fc.herring) <- paste0("Herring_", colnames(ct.fc.herring))
    
  ## Combine and assess concordance
    ct.fc <- cbind(ct.fc.atac, ct.fc.H3K27ac)
    m <- match(rownames(ct.fc), rownames(ct.fc.herring))
    ct.fc <- cbind(ct.fc, ct.fc.herring[m,])
    
    ct.fc$ATAC_MaxID <- splitter(ct.fc$ATAC_MaxID, "_", 1)
    ct.fc$H3K27ac_MaxID <- splitter(ct.fc$H3K27ac_MaxID, "_", 1)
    
    ct.fc <- ct.fc[,grep("MaxID|VsRank2|Specific", colnames(ct.fc))]
    
    ct.fc$Enh <- rownames(ct.fc)
    ct.fc <- relocate(ct.fc, "Enh") 
    

  ## Save
    write.csv(ct.fc, file = "Coverage/Cell-type-specificity - Foldchanges.csv")
    
    
# ## Compare hit enhancers to non-hit peaks
#   ## Nott
#     x <- covNott$window1000 %>% as.data.frame() 
#     x$Hit <- rownames(x) %in% hit.enh
#     x <- melt(x)
#     x$Celltype <- splitter(x$variable, "_", 1)
#     x$Protocol <- splitter(x$variable, "_", 2)
# 
#     pdf(file = "Coverage/Nott - Hit Coverage.pdf", height = 3, width = 8)
#     ggplot(x, aes(x = Celltype, y = value, fill = Hit)) +
#       # geom_quasirandom(dodge.width = 0.5) +
#       geom_violin(position = position_dodge(width = 0.75), width = 0.7, scale = "width", draw_quantiles = 0.5, colour = "black", alpha = 0.75) +
#       scale_y_continuous(trans = "log2", breaks = c(0, 8, 32, 128, 512, 2048, 8192, 32768)) +
#       facet_wrap(~Protocol) +
#       theme_bw() +
#       scale_fill_manual(values = c("firebrick1", "dodgerblue4")) +
#       labs(x = "Celltype", y = "Coverage at Peak (log2 axis)")
#     dev.off()
#     
#   ## Herring
#     x <- covHerring$window1000 %>% as.data.frame() 
#     x$Hit <- rownames(x) %in% hit.enh
#     x <- melt(x)
#     x$Celltype <- splitter(x$variable, "_", 1)
#     x$Stage <- gsub("Neuro_", "", x$variable) %>%  splitter("_", 2) %>% substr(2, 100)
# 
#     pdf(file = "Coverage/Herring - Hit Coverage.pdf", height = 6, width = 10)
#     ggplot(x, aes(x = Celltype, y = value, fill = Hit)) +
#       # geom_quasirandom(dodge.width = 0.5) +
#       geom_violin(position = position_dodge(width = 0.75), width = 0.7, scale = "width", draw_quantiles = 0.5, colour = "black", alpha = 0.75) +
#       scale_y_continuous(trans = "log2", breaks = c(0, 8, 32, 128, 512, 2048, 8192)) +
#       facet_wrap(~Stage) +
#       theme_bw() +
#       scale_fill_manual(values = c("firebrick1", "dodgerblue4")) +
#       labs(x = "Celltype", y = "Coverage at Peak (log2 axis)")
#     dev.off()
#     
#     
#   ## Stats
#     covStats <- list()
#     
#     # herring
#     covStats$Herring <- apply(covHerring$window1000, 2, function(x) {
#       h <- rownames(covHerring$window1000) %in% hit.enh
# 
#       t <- t.test(log2(x + 1) ~  h)
#       w <- wilcox.test(log2(x + 1) ~  h)
#       out <- data.frame(MeanCandidate = t$estimate[1], MeanHitEnh = t$estimate[2], Ttest_P = t$p.value, Wtest_P = w$p.value)
#     })
#     covStats$Herring <- do.call("rbind", covStats$Herring)
#     
#     # nott
#     covStats$Nott <- apply(covNott$window1000, 2, function(x) {
#       h <- rownames(covNott$window1000) %in% hit.enh
# 
#       t <- t.test(log2(x + 1) ~  h)
#       w <- wilcox.test(log2(x + 1) ~  h)
#       out <- data.frame(MeanCandidate = t$estimate[1], MeanHitEnh = t$estimate[2], Ttest_P = t$p.value, Wtest_P = w$p.value)
#     })
#     covStats$Nott <- do.call("rbind", covStats$Nott)
#     
#     # combine
#     covStats <- do.call("rbind", covStats)
#     covStats$Study <- splitter(rownames(covStats), "\\.", 1)
#     
#     # save
#     write.csv(covStats, file = "Coverage/Hit Coverage Stats.csv")
    
    
## Output
  # ct.fc <- read.csv("Coverage/Cell-type-specificity - Foldchanges.csv", row.names = 1)
    
  # annotation table 
  candidate.annot$AstSpecific_NottAtac <- candidate.annot$Enh %in% ct.fc$Enh[which(ct.fc$ATAC_Specific & ct.fc$ATAC_MaxID == "Astro")]  
  candidate.annot$AstSpecific_NottH3K27ac <- candidate.annot$Enh %in% ct.fc$Enh[which(ct.fc$H3K27ac_Specific & ct.fc$H3K27ac_MaxID == "Astro")]  
  candidate.annot$AstSpecific_Herring <- candidate.annot$Enh %in% ct.fc$Enh[which(ct.fc$Herring_Specific & ct.fc$Herring_MaxID == "Astro")]  

  # enrichment table
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_NottAtac", rbind = FALSE)
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_NottH3K27ac", rbind = TRUE)
  candidate.enrich <- run.EnhancerFisher(data.column = "AstSpecific_Herring", rbind = TRUE)

  
  
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
  
  ## Regions of human-accelerated evolution rate from Doan 2016 (pooled from 5 past studies)
    doan_data_hg19 <- read_xlsx("../../../../PublicData/Doan2016/1-s2.0-S0092867416311692-mmc2.xlsx") %>% as.data.frame()
    colnames(doan_data_hg19) <- c("Chr", "Start", "End", "Class", "ReferenceGene", "ID_1", "ID_2", "ChromHMM", "TargetGene")
    doan_data_hg19 <- doan_data_hg19[,c("Chr", "Start", "End", "Class", "TargetGene")]    
    doan_dir_in <- "../../../../PublicData/Doan2016/Doan_hg19.bed"
    doan_dir_out <- "Doan.bed"
    write.bed(doan_data_hg19, doan_dir_in)

  ## Regions of human-accelerated evolution from Keough 2023 in the Zoonomia consortium
    # keough_data <- read_xlsx("../../../../PublicData/Zoonomia/Keough_HARs/science.abm1696_table_s1.xlsx", sheet = "zooHARs") %>% as.data.frame()
    # keough_data <- keough_data[,c("chrom", "start", "end", "simple_name")]    
    # keough_dir_in <- "../../../../PublicData/Zoonomia/Keough_HARs/zooHARs.bed"
    # keough_dir_out <- "Keough.bed"
    # write.bed(keough_data, keough_dir_in)
    
    # removed as no intersect
    
  ## Regions of human-specific insertions from Kronenberg 2018 (and further processed in Keough 2023)
    # kron_data <- read_xlsx("../../../../PublicData/Zoonomia/Keough_HARs/science.abm1696_table_s3.xlsx", sheet = "disruption_scores_hsSVs") %>% as.data.frame()
    # kron_data <- kron_data[which(kron_data$type == "hsIns"),] # 12k insertions of 17k total
    # kron_data <- kron_data[,c("chrom", "start", "end")]
    # kron_dir_in <- "../../../../PublicData/Zoonomia/Keough_HARs/hsSVs.bed"
    # kron_dir_out <- "hsSVs.bed"
    # write.bed(kron_data, kron_dir_in)
    
    # removed as no intersect
    

## Intersect
  # doan suffices with a binary left outer join
  loj(nha_dir_19, doan_dir_in, doan_dir_out)
  
  # full intersect information for Vermunt
  wawb(nha_dir_38, vermunt_dir_in, vermunt_dir_out)
  
## Read in and wrangle
  ## Doan
    # read in
    doan_intersect <- read.bed(doan_dir_out)
    doan_intersect <- doan_intersect[which(doan_intersect$V8 != -1),] # just the one...
    
    # add to annotation
    candidate.annot$HumanAccelerated_Doan2016 <- candidate.annot$Coord %in% doan_intersect$V4
    candidate.enrich <- run.EnhancerFisher(data.column = "HumanAccelerated_Doan2016", rbind = TRUE)
    
  ## Vermunt
    vermunt_intersect <- read.bed(vermunt_dir_out)
    
    # add to annotation
    for (j in unique(vermunt_intersect$V9)) {
      x <- vermunt_intersect[which(vermunt_intersect$V9 == j),]
      y <- paste0("HumanH3K27ac_", j, "_Vermunt2016")
      
      candidate.annot[,y] <- candidate.annot$Coord %in% x$V4
      candidate.enrich <- run.EnhancerFisher(data.column = y, rbind = TRUE) 
    }
    

################################################################################################################################ #
## Ageing ----
    
  
## Nativio 2018
  # h4k16ac in the temporal cortex
  # statistically controlled for neuronal proportion (e.g. 10% of neuronal-proportion associated peaks removed)
  # compares ages ~50 to ~70
  # hg 19
    
## Already processed when downloaded
    
## Overlap
  loj(nha_dir_19, 
      "../../../../PublicData/Nativio2018_H4K16ac_Ageing_AD/GSE84618_Age-Regulated_H4K16ac_Changes.hg19.bed",
      "Nativio_AgeRegulated.bed")
    
  loj(nha_dir_19, 
    "../../../../PublicData/Nativio2018_H4K16ac_Ageing_AD/GSE84618_Age-Disregulated_H4K16ac_Changes.hg19.bed",
    "Nativio_AgeDysregulated.bed")
  

## Read in
  x <- list.files()
  x <- x[grep("Nativio", x)]
  nativio <- lapply(x, read.bed)
  names(nativio) <- gsub(".bed", "", x) %>% splitter("_", 2)
  nativio <- do.call("rbind", nativio)
  nativio <- nativio[which(nativio$V8 != -1),]
  nativio$Condition <- splitter(rownames(nativio), "\\.", 1)
  
## Annotate
  candidate.annot$Ageing_Nativio2018 <- candidate.annot$Coord %in% nativio$V4[grep("Age", nativio$Condition)]
  candidate.enrich <- run.EnhancerFisher(data.column = "Ageing_Nativio2018", rbind = TRUE) 
  
    
        
################################################################################################################################ #
## Brain disorder-associated peaks ----

  
## Alzheimer's (Nativio 2018)
  # simply use the above
  # bulk temporal cortex H4K16ac
  candidate.annot$AD_Nativio2018 <- candidate.annot$Coord %in% nativio$V4[grep("AD", nativio$Condition)]
  candidate.enrich <- run.EnhancerFisher(data.column = "AD_Nativio2018", rbind = TRUE) 
  
  
## Alzheimer's (Morabito 2021)
  # snATAC 
  
  ## Read in 
    morabito_data <- read_xlsx("../../../../PublicData/Morabito2021_AD_snATAC_snRNA/41588_2021_894_MOESM9_ESM.xlsx", skip = 2)
    
    # filter
    morabito_data <- morabito_data[which(morabito_data$cell_type == "ASC"),] # astrocytes
    morabito_data <- morabito_data[which(morabito_data$p_val_adj < 0.05),] # and significant
    
    # convert to bed
    x <- morabito_data$Peak
    morabito_data <- data.frame(Chr = splitter(x, ":", 1),
                                Start = splitter(splitter(x, ":", 2), "-", 1),
                                End = splitter(splitter(x, ":", 2), "-", 2),
                                PAdj = morabito_data$p_val_adj,
                                logFC = morabito_data$avg_logFC)
    
    morabito_dir_in <- "../../../../PublicData/Morabito2021_AD_snATAC_snRNA/AD_Ast_DiffAtac.bed"
    morabito_dir_out <- "Morabito_AD_Ast_DiffAtac.bed"
    
    write.bed(morabito_data, morabito_dir_in)
    
  ## Intersect
    wawb(nha_dir_38, morabito_dir_in, morabito_dir_out)
    
  ## Analyse
    morabito_overlap <- read.bed(morabito_dir_out)
    candidate.annot$AD_Morabito2021 <- candidate.annot$Coord %in% morabito_overlap$V4
    candidate.enrich <- run.EnhancerFisher(data.column = "AD_Morabito2021", rbind = TRUE) 
    
    
## Alzheimer's (Bendl 2022)
  # ATAC from NeuN- sorted temporal and entorhinal cortex
  
  ## Intersect
    wawb(nha_dir_38, "../../../../PublicData/Bendl2022_AD_ATAC/Data_peaks_glia.bed", "Bendl.bed")
    
  ## Process
    bendl <- read.bed("Bendl.bed")
    
    # get annotation of differential atac regions
    bendl_hits <- read.delim("../../../../PublicData/Bendl2022_AD_ATAC/Data_peaks_diff_glia.tsv")
    bendl_hits <- bendl_hits[,grep("_all_", colnames(bendl_hits))] # results from pooling both brain regions
    colnames(bendl_hits) <- c("Braak", "AD", "ClinDemRating", "Plaque")
    
    x <- lapply(colnames(bendl_hits), function (x) {
      rownames(bendl_hits)[which(bendl_hits[,x] == 1)]
    })
    
    names(x) <- colnames(bendl_hits)
    bendl_hits <- x
    
  ## Annotate candidate peaks
    for (j in names(bendl_hits)) {
      
      label <- paste0("AD_Bendl2022_", j)
      x <- bendl[which(bendl$V8 %in% bendl_hits[[j]]),]
      
      if (nrow(x) == 0) next
      
      candidate.annot[,label] <- candidate.annot$Coord %in% x$V4
      candidate.enrich <- run.EnhancerFisher(data.column = label, rbind = TRUE) 
    }
    
  
## Glioblastoma superenhancers
  # from Xu et al, 2021, Science
  # h3k27ac from a 50-100 GBM samples
  # this table focus on reported superenhancers
  # hg38
  
  ## Read in
    gbm_dir_base <- "../../../../PublicData/Xu2021_GBM_H3K27ac/abd4676_table_s1.xlsx"
    gbm_dir_in <- "../../../../PublicData/Xu2021_GBM_H3K27ac/GBM_Tissue_Superenhancers.bed"
    gbm_dir_out <- "GBM_Superenhancers.bed"
  
    # superenhancer coordinates
    gbm_se <- read_xlsx(gbm_dir_base, sheet = "Table S1h", skip = 1)
    
    # metadata annotation
    gbm_meta_tissue <- read_xlsx(gbm_dir_base, sheet = "Table S1a", skip = 1)
    gbm_meta_tissue <- gbm_meta_tissue[which(gbm_meta_tissue$`Pathology/Grade` == "GBM"),] # gbm samples only
    gbm_meta_tissue <- gbm_meta_tissue[which(gbm_meta_tissue$`GIE subtypes` != "NA"),] # samples with a GIE subtype
    gbm_meta_tissue <- gbm_meta_tissue[,c("Sample ID", "GIE subtypes")]
    colnames(gbm_meta_tissue) <- c("Sample", "Subtype")
    
    gbm_meta_cells <- read_xlsx(gbm_dir_base, sheet = "Table S1b", skip = 1)
    gbm_meta_cells <- gbm_meta_cells[grep("GBM", gbm_meta_cells$`Sample type`),]
    gbm_meta_cells <- data.frame(Sample = gbm_meta_cells$`Sample ID`,
                                 Subtype = "Cellline")
    
    gbm_meta <- rbind(gbm_meta_tissue, gbm_meta_cells) # combined
    
  ## Annotate SEs by their type
    gbm_se <- gbm_se[which(gbm_se$Sample %in% gbm_meta$Sample),]
    gbm_se <- gbm_se[,c("Chr", "Start", "Stop", "Sample")]
    gbm_se$Type <- gbm_meta$Subtype[match(gbm_se$Sample, gbm_meta$Sample)] # no cells, just tissue?
    
    gbm_meta <- gbm_meta[which(gbm_meta$Sample %in% gbm_se$Sample),]
    
    write.bed(gbm_se, gbm_dir_in)
    
  ## Intersect
    # i note that this intersection contains SEs called in each individual GBM, rather than a merged list
    # however, each candidate peak is assessed for having at least one 
    wawb(nha_dir_38, gbm_dir_in, gbm_dir_out)
    
  ## Analyse
    gbm_overlap <- read.bed(gbm_dir_out)
    
    for (j in unique(gbm_overlap$V9)) { # for each subtype
      # label
      label <- paste0("GBM_", j, "_Xu2021")
      
      # count overlaps
      x <- gbm_overlap$V4[which(gbm_overlap$V9 == j)] # candidate peaks which overlap a SE of the given subtype j
      x <- table(x)
      
      # filter to peaks overlapping SEs found in >= 50 of samples
      min_n <- length(which(gbm_meta$Subtype == j)) * 0.5
      x <- names(x)[which(x >= min_n)]
      
      # annotate 
      candidate.annot[,label] <- candidate.annot$Coord %in% x
      candidate.enrich <- run.EnhancerFisher(data.column = label, rbind = TRUE) 
    }
    
## Glioblastoma coverage calls (Xu 2021)
    
  # drop because it is only GBM, with no control...
    
  # # per the above resource, but using bigwigs to call coverage directly in our preaks
  # 
  # ## Create windowed peaks for hg38  
  #   x <- read.delim(nha_dir_38, header = FALSE)
  #   window <- 1000
  #   x$V2 <- x$V2 - (window/2)
  #   x$V3 <- x$V3 + (window/2)
  #   
  #   # save
  #   nha_dir_38_w1000 <- "../../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg38_window1000.bed"
  #   write.bed(x, nha_dir_38_w1000)
  #   
  # ## Directories
  #   bed_nha_hg38_w1000_file <- import.bed(nha_dir_38_w1000)
  #   xu_dir <- "../../../../PublicData/Xu2021_GBM_H3K27ac/Bigwigs/"
  #   bw_xu <- list.files(xu_dir) 
  #   bw_xu <- bw_xu[grep("bw", bw_xu)]
  #   bw_xu <- paste0(xu_dir, bw_xu)
  #   
  # 
  # ## Get metadata
  #   xu_meta_tissue <- read_xlsx(gbm_dir_base, sheet = "Table S1a", skip = 1)
  #   xu_meta_tissue <- xu_meta_tissue[which(xu_meta_tissue$`Pathology/Grade` == "GBM"),] # gbm samples only
  #   xu_meta_tissue <- xu_meta_tissue[which(xu_meta_tissue$`GIE subtypes` != "NA"),] # samples with a GIE subtype
  #   xu_meta_tissue <- xu_meta_tissue[,c("Sample ID", "GIE subtypes")]
  #   colnames(xu_meta_tissue) <- c("Sample", "Subtype")
  #   
  #   xu_meta_cells <- read_xlsx(gbm_dir_base, sheet = "Table S1b", skip = 1)
  #   xu_meta_cells <- xu_meta_cells[grep("GBM", xu_meta_cells$`Sample type`),]
  #   xu_meta_cells <- data.frame(Sample = xu_meta_cells$`Sample ID`,
  #                                Subtype = "Cellline")
  #   
  #   xu_meta <- rbind(xu_meta_tissue, xu_meta_cells); rm(xu_meta_tissue, xu_meta_cells) # combined
  #   
  # ## Calculate coverage in each of the H3K27ac bigwigs
  #   # get
  #   covXu_w1000 <-  lapply(bw_xu, function(x) { # requires ~20s
  #     print(x)  
  #     
  #     # get coverage
  #     y <- import(x, selection = bed_nha_hg38_w1000_file, as = "NumericList")
  #     z <- sapply(y, sum) # this is the sum of coverages across the bases, and is not normalised to peak width. this makes comparisons between peaks difficult
  #     
  #     # collect peak ids
  #     names(z) <- bed_nha_hg38_w1000_file$name # $name is the GRCh38 coordinate of the peak
  #     
  #     return(z)
  #   })
  #   
  #   # combine
  #   names(covXu_w1000) <- splitter(bw_xu, "_", 4)
  #   covXu_w1000 <- do.call("cbind", covXu_w1000)
  #     
  #   # ids
  #   m <- match(rownames(covXu_w1000), guides$TargetCoord)
  #   rownames(covXu_w1000) <- guides$TargetID[m]
  #     

## Glioblastoma stem cells vs non-stem cells (Tome-Garcia 2018)
  # by sorting GBM cells for ability to bind EGF
  # atac seq performed
  # only upregulated peaks reported
    
  ## Directories
    tome_dir <- "../../../../PublicData/TomeGarcia2018_GBM_ATAC/41467_2018_6258_MOESM4_ESM.xlsx"
    tome_dir_in <- "../../../../PublicData/TomeGarcia2018_GBM_ATAC/GBM_StemCells_Up.bed"
    tome_dir_out <- "GBM_StemCells.bed"
    
  ## Read in and process raw data
    tome_data <- read_xlsx(tome_dir, sheet = "ALL GSC-Specific UP peaks")
    tome_data <- tome_data[,c("Chr", "Start", "End")]
    write.bed(tome_data, tome_dir_in)
    
  ## Intersect
    wawb(nha_dir_19, tome_dir_in, tome_dir_out)
    
  ## Analyse
    tome_overlap <- read.bed(tome_dir_out)
    
    candidate.annot$GBM_StemCellUp_TomeGarcia2018 <- candidate.annot$Coord %in% tome_overlap$V4
    candidate.enrich <- run.EnhancerFisher(data.column = "GBM_StemCellUp_TomeGarcia2018", rbind = TRUE)
    
  
    
    
## ASD differential acetylation (Sun 2016)
  # h3k27ac profiling of the asd and ctl pfc, tc, and cb
  # hg19
    
  ## Directories
    sun_dir <- "../../../../PublicData/Sun2016_ASD_H3K27ac/Sun2016_ST4.xlsx"
    sun_dir_in <- "../../../../PublicData/Sun2016_ASD_H3K27ac/Bed.bed"
    sun_dir_out <- "ASD_H3K27ac.bed"
    
  ## Process
    sun_names <- excel_sheets(sun_dir)
    sun_names <- sun_names[grep("PFC|TC", sun_names)]
    sun_data <- lapply(sun_names, function(x) { read_xlsx(sun_dir, sheet = x) })
    names(sun_data) <- gsub(" ", "_", sun_names)
    sun_data <- do.call("rbind", sun_data)
    sun_data$Condition <- splitter(rownames(sun_data), "\\.", 1)
    sun_data <- sun_data[,c("Chrom", "Start", "End", "Condition")]
    write.bed(sun_data, sun_dir_in)
    
  ## Intersect
    wawb(nha_dir_19, sun_dir_in, sun_dir_out)
    
  ## Analyse
    sun_overlap <- read.bed(sun_dir_out)
    candidate.annot$ASD_Sun2016 <- candidate.annot$Coord %in% sun_overlap$V4
    candidate.enrich <- run.EnhancerFisher(data.column = "ASD_Sun2016", rbind = TRUE)   
    
    
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
  candidate.annot$KnownEnh_Yao2022 <- candidate.annot$Coord %in% known$V4
  candidate.enrich <- run.EnhancerFisher(data.column = "KnownEnh_Yao2022", rbind = TRUE)
    
   
################################################################################################################################ #
## Transcription factor footprints ----

    
## These results were generated using the TOBIAS algorithm, which looks for motifs within footprints in the ATAC data
  # or, rather, it calculates the inverse: which motifs have footprints?
    
    
## To begin, filter data to only motifs with all components expressed
  ## First: what TFs are used as input from JASPAR?
    # read in JASPAR
    tob_TF <- read.table("../../../../PublicData/TF_BindingSites/JASPAR2022_human_TF.txt")
    tob_TF$V1 <- sub(".*\\.","",tob_TF$V1)
    tob_TF$V2 <- sub("(.*)::(.*)","\\2",tob_TF$V1)
    tob_TF$V1 <- sub("(.*)::(.*)","\\1",tob_TF$V1)
    tob_TF <- unique(tob_TF)
    
    # get EnsID  
    library(biomaRt)
    getENSID <- function(symbols, attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", aggregate = TRUE) {
      mart <- useMart("ensembl", "hsapiens_gene_ensembl")
      mart <- useDataset("hsapiens_gene_ensembl",mart)
      ENSB <- getBM(attributes = attributes, filters = filters ,values = symbols,
                    mart = mart, uniqueRows = TRUE)
      
      if (aggregate == TRUE) {
        ENSB <- aggregate(ENSB$ensembl_gene_id, by = list(ENSB[,filters]), FUN = function (x) { paste0(x, collapse = ",") } )
      }
      
      return(ENSB)
    }
    
    ENSB <- getENSID(union(tob_TF$V1, tob_TF$V2 ))

    tob_TF <- merge(tob_TF, ENSB, all.x = TRUE, by.x = "V1", by.y = "Group.1")
    tob_TF <- merge(tob_TF, ENSB, all.x = TRUE, by.x = "V2", by.y = "Group.1")
    colnames(tob_TF) <- c("TF.1", "TF.2", "ENS.1", "ENS.2")

    tob_TF$concat_names <- tob_TF$TF.2
    tob_TF[tob_TF$TF.2 != tob_TF$TF.1, ]$concat_names <- paste0(tob_TF[tob_TF$TF.2 != tob_TF$TF.1, ]$TF.2, tob_TF[tob_TF$TF.2 != tob_TF$TF.1, ]$TF.1)

  ## Second: which genes are expressed?
    # rnaseq
    load("../../../../PreprocessedData/Bulk_ATAC_RNAseq/RESULTS/RNAseq/geneData.rda")
    genes.exp.rnaseq <- geneData$geneRPKM[which(geneData$geneRPKM$NHA > 1), c("Symbol", "EnsID")]
    genes.exp.rnaseq <- c(genes.exp.rnaseq$Symbol, genes.exp.rnaseq$EnsID) %>% as.character() %>% unique()
    
    # scrnaseq
    genes.exp.10x <- rowMeans(nha@assays$RNA@data) 
    genes.exp.10x <- genes.exp.10x[which(genes.exp.10x > 2^-6)] %>% names()
    
    # combine
    genes.exp <- union(genes.exp.10x, genes.exp.rnaseq)
    
  ## Combine 1 and 2
    tob_TF$TF.1_Exp <- (tob_TF$TF.1 %in% genes.exp) | (tob_TF$ENS.1 %in% genes.exp)
    tob_TF$TF.2_Exp <- (tob_TF$TF.2 %in% genes.exp) | (tob_TF$ENS.2 %in% genes.exp)
    
    tob_TF_exp <- tob_TF$concat_names[which(tob_TF$TF.1_Exp & tob_TF$TF.2_Exp)] # 369. It is 342 when using only symbol
    
    
## Now, read in TOBIAs output 
    ## Setup
    tob <- list()
    
    tob_baseMx <- candidate.annot[,1:4] # z is a dataframe that will be filled with output from Tobias
    tob_baseMx$nMotif <- NA # number of motifs for this peak
    tob_baseMx[,tob_TF_exp] <- NA # count columns for each expressed TF
    
  ## Read in
    # bound TFs
    tob$bound <- read.bed("../../../../EnhancerPredictionModels/Results/Tobias/Footprint/BINDetect/bound_overlaps.bed")
    
    # unbound TFs
    tob$unbound <- read.bed("../../../../EnhancerPredictionModels/Results/Tobias/Footprint/BINDetect/unbound_overlaps.bed")
    
    # process
    tob <- lapply(tob, function(x) {
      # filter columns
      m <- match(x$V4, candidate.annot$Coord)
      x$Enh <- candidate.annot$Enh[m]
      x$TF <- splitter(x$V8, "_", 1) %>% splitter("\\.", 3)
      x <- x[,c("Enh", "TF")]
      
      # add annotation
      y <- tob_baseMx
      for (j in tob_TF_exp) {
        
        
        if (any(x$TF == j)) { # if the tf has any peaks
          
          z <- x[which(x$TF == j),] # filter to peaks-motif pairs
          # y[,j] <-y$Enh %in% z$Enh # old code: annotate the peaks that have it as TRUE, regardless of how many
          z$Enh <- factor(z$Enh, levels = y$Enh)
          y[,j] <- table(z$Enh) %>% as.numeric() # this counts 
          
        
          } else {
            
          # y[,j] <- NA # rather than FALSE
          y[,j] <- 0 
          
        }
        
        
      }
      
      y$nMotif <- rowSums(y[,tob_TF_exp], na.rm = TRUE)
      
      # finish
      return(y)
      
    })
  
  write.csv(tob$bound, file = "TF/Bound Counts.csv")
  write.csv(tob$unbound, file = "TF/Unbound Counts.csv")  

  
## Explore the bound hits
  ## Are any TFs enriched within hits?
    # setup  
    candidate.tf.enrich <- candidate.enrich[0,] # retains only structure and column names
    
    # run
    for (j in tob_TF_exp) {
      candidate.tf.enrich <- run.EnhancerFisher(data = tob$bound, data.column = j, rbind = TRUE, rbind.to = candidate.tf.enrich)
    }
  
    # fdr correct
    candidate.tf.enrich$FDR <- p.adjust(candidate.tf.enrich$p, method = "fdr")
    candidate.tf.enrich <- candidate.tf.enrich[order(candidate.tf.enrich$p),]
    
    # save
    write.csv(candidate.tf.enrich, file = "TF/Bound Enrichment.csv")
  
  ## Is the total number of bound TFs different within hits?
    p <- tob$bound
    p$Call <- factor(p$Hit)
    levels(p$Call) <- c("ns Peak", "Hit Enhancer")
    
    pdf(file = "TF/Bound Sum Across Peak.pdf", height = 3, width = 3)
    ggplot(p, aes(y = nMotif, x = Call)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.2) +
      theme_bw() +
      theme(axis.title.x = invis, panel.grid.major.x = invis, panel.border = invis,
            axis.line.y = element_line()) +
      labs(y = "Number of Bound TF Motifs")
    dev.off()
  
## Explore motif binding ratios
  # ultimately, TOBIAS calculates TFBS based on motifs
  # then categorises as bound if there is an overlapping footprint, and unbound if not
   
  ## Are a greater percentage of motifs bound than unbound in hits?
    p <- data.frame(Hit = tob$bound$Hit,
                    Tested = tob$bound$Tested,
                    nMotif_Bound = tob$bound$nMotif,
                    nMotif_Unbound = tob$unbound$nMotif)
    p <- p[which(p$Tested),]
    p$BoundVsUnbound <- p$nMotif_Bound / (p$nMotif_Bound + p$nMotif_Unbound)
    p$Call <- factor(p$Hit)
    levels(p$Call) <- c("ns Peak", "Hit Enhancer")
    
    pdf(file = "TF/Ratio BoundvsUnbound - Across Peak Violins.pdf", height = 3, width = 3)
    ggplot(p, aes(x = Call, y = BoundVsUnbound)) +
      geom_violin(scale = "width") +
      geom_quasirandom(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = "+", colour = "red", size = 10) +
      theme_bw() +
      theme(axis.title.x = invis, panel.grid.major.x = invis, panel.border = invis,
            axis.line.y = element_line()) +
      labs(y = "Fraction of TF Motifs in Footprints")
    dev.off()
    
  ## For each motif, is is more frequently bound in hit than non-hit enhancers?
    motif_BoundVsUnbound <- list()
    
    for (j in tob_TF_exp) {
      # get motif bound and unbound counts
      unbnd_ns <- sum(tob$unbound[,j] > 0 & !tob$unbound$Hit)
      unbnd_hit <- sum(tob$unbound[,j] > 0 & tob$unbound$Hit)
      
      bnd_ns <- sum(tob$bound[,j] > 0 & !tob$bound$Hit)
      bnd_hit <- sum(tob$bound[,j] > 0 & tob$bound$Hit)
      
      # run fisher test
      tab <- data.frame(Unbound = c(unbnd_ns, unbnd_hit),
                        Bound = c(bnd_ns, bnd_hit),
                        row.names = c("ns Peak", "Hit Enhancer"))
      
      f <- fisher.test(tab)
      
      # output
      motif_BoundVsUnbound[[j]] <- data.frame(Unbound_ns = unbnd_ns,
                                              Unbound_hit = unbnd_hit,
                                              Bound_ns = bnd_ns,
                                              Bound_hit = bnd_hit,
                                              p = f$p.value,
                                              OR = f$estimate,
                                              Lower = f$conf.int[1],
                                              Upper = f$conf.int[2])
      
    }
    
    motif_BoundVsUnbound <- do.call("rbind", motif_BoundVsUnbound)
    motif_BoundVsUnbound$FDR <- p.adjust(motif_BoundVsUnbound$p, method = "fdr")
    motif_BoundVsUnbound <- motif_BoundVsUnbound[order(motif_BoundVsUnbound$p),]
    
    write.csv(motif_BoundVsUnbound, file = "TF/Ratio BoundvsUnbound - Enrichment.csv")
    
      
## Motif footprint plots
  ## Setup
    footprints <- list(bound = read.bed("../../../../EnhancerPredictionModels/Results/Tobias/Footprint/BINDetect/bound_overlaps.bed"),
               unbound = read.bed("../../../../EnhancerPredictionModels/Results/Tobias/Footprint/BINDetect/unbound_overlaps.bed"))
    footprints <- do.call("rbind", footprints)
    footprints$Footprint <- splitter(rownames(footprints), "\\.", 1)
    m <- match(footprints$V4, candidate.annot$Coord)
    footprints$Enh <- candidate.annot$Enh[m]
    footprints$Hit <- footprints$Enh %in% hit.enh
    footprints <- footprints[,-c(1:4)] # remove coordinates for the peak, not the motif
    colnames(footprints)[1:4] <- c("Motif_Chr", "Motif_Start", "Motif_End", "Motif")
    footprints$Motif <- splitter(footprints$Motif, "_", 1) %>% splitter("\\.", 3)

  ## Function to plot
    plot.footprint <- function(set) {
      ## Write out beds to scratch directory
        dir.hit <- "../../Scratchspace/Footprint_Set_Hit.bed"
        dir.ns <- "../../Scratchspace/Footprint_Set_ns.bed"
        dir.out <- "../../Scratchspace/Footprint_Aggregation.txt"
        
        set_hit <- which(set & footprints$Hit)
        set_ns <- which(set & !(footprints$Hit))
        
        write.bed(footprints[set_hit,], dir.hit)
        write.bed(footprints[set_ns,], dir.ns)
      
      ## Get aggregated atac signal using TOBIAS
        system("bash TF/Tobias_PlotAggregate.sh", intern = FALSE, wait = TRUE, ignore.stdout = TRUE) 
      
        # the shell code looks as follows:
        # cd /mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/
        # source /home/rna2/PROGRAMS/miniconda3/etc/profile.d/conda.sh
        # conda activate tobias
        # 
        # TOBIAS PlotAggregate \
        # --TFBS ../../Scratchspace/Footprint_Set_Hit.bed ../../Scratchspace/Footprint_Set_ns.bed \
        # --signals ../../../../EnhancerPredictionModels/Results/Tobias/Footprint/NHA_ATAC_S3.filtered_corrected.bw \
        # --output-txt ../../Scratchspace/Footprint_Aggregation.txt 
      
      ## Create plot in R
        # read in and process
        x <- read.table(dir.out, skip = 2, sep = "\t")
        x <- data.frame(Hit = do.call("c", strsplit(x$V3[1], ",")),
                        Ns = do.call("c", strsplit(x$V3[2], ",")))
        x$Position <- (1:nrow(x)) - (round(nrow(x) / 2))
        x <- melt(x, id.vars = "Position")
        x$value <- as.numeric(x$value)
        
        levels(x$variable) <- c(paste0("Hit Enhancer (", length(set_hit), ")"), 
                                paste0("ns Peak (", length(set_ns), ")"))
        
        x <- x[-which(abs(x$Position)  > 45),]
        
        
        # plot
        ggplot(x, aes(x = Position, y = value, colour = variable, group = variable)) +
          geom_line(alpha = 0.8) +
          theme_bw() +
          theme(panel.grid = invis, panel.border = invis, axis.line = element_line(), legend.title = invis,
                legend.position = "bottom") +
          scale_colour_lancet() +
          scale_x_continuous(expand = c(0,0)) +
          labs(y = "ATAC Signal", x = "Position relative to motif centre (bp)")
    }
    
  ## Make versions
    # KLF5
    pdf(file = "TF/Aggregate Plot - KLF5.pdf", height = 3, width = 3.5)
    plot.footprint(footprints$Motif == "KLF5") 
    plot.footprint(footprints$Motif == "KLF5" & footprints$Footprint == "bound") + labs(x = "Position relative to bound motif centre (bp)")
    dev.off()
    
    # JUN
    pdf(file = "TF/Aggregate Plot - JUN.pdf", height = 3, width = 3.5)
    plot.footprint(footprints$Motif == "JUN")
    plot.footprint(footprints$Motif == "JUN" & footprints$Footprint == "bound") + labs(x = "Position relative to bound motif centre (bp)")
    dev.off()
    
    # all bound motifs
    pdf(file = "TF/Aggregate Plot - All Motifs.pdf", height = 3, width = 3.5)
    plot.footprint(footprints$Footprint == "bound" | footprints$Footprint == "unbound") 
    plot.footprint(footprints$Footprint == "bound") + labs(x = "Position relative to bound motif centre (bp)")
    dev.off()
        
     
 
    


  
################################################################################################################################ #
## Transcription factor ChIP ----
    
## Compare to an alternative TF resource: Remap
  # uses ChIP data from the ENCODE consortium, across a range of cell-types
  
## Directories
  remap.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/TF_BindingSites/remap2022_nr_macs2_hg38_v1_0.bed"
  remap.out <- "TF/Remap_Overlaps.bed"


## Intersect
  call <- paste("intersectBed",
                "-a", nha_dir_38,
                "-b", remap.in,
                "-wb", 
                ">", remap.out)

  system(call, intern = FALSE, wait = TRUE) 
  
  
## Read in and wrangle
  remap <- read.bed(remap.out)
  remap$V4 <- sub("_", ":", remap$V4) %>% sub("_", "-", .)
  m <- match(remap$V4, guides$TargetCoord)
  
  remap <- data.frame(Enh = guides$TargetID[m],
                      Enh.Coord = guides$TargetCoord[m],
                      TF = splitter(remap$V8, ":", 1),
                      TF.Expressed = ".",
                      TF.Coord = paste0(remap$V5, ":", remap$V6, "-", remap$V7),
                      Celltypes = splitter(remap$V8, ":", 2))
  
## Filter TFs
  # is the TF expressed?
  remap$TF.Expressed <- remap$TF %in% genes.exp
  remap <- remap[which(remap$TF.Expressed),]
  
  # to those in the TOBIAS analysis
  remap$Tobias_TFtested <- remap$TF %in% tob_TF_exp
  
  
## Categorise peaks as: LevelA (in Ast); LevelB (2+ brain); LevelC (1 brain and 1+ other); LevelD (1 brain); LevelE (2+ other); and LevelF (1 non-brain)
  brain.ct <- c("astrocyte", "neural-progenitor", "glioma", "tibial-nerve", "primary-glioblastoma", "glioblastoma", "NB-1643", "hippocampus",
                "cortical-interneuron", "neuron-progenitor", "neuroepithelilal-cells", "neuron", "nerve", "NPC", "dopaminergic-neuron",
                "neural", "brain-prefrontal-cortex", "SH-SY5Y", "SK-N-BE2-C", "HNPC", "SK-N-MC")
  
  split <- strsplit(remap$Celltypes, ",")
  nBrain <- sapply(split, function(x) length(which(x %in% brain.ct)))
  nNonBrain <- sapply(split, function(x) length(which(!(x %in% brain.ct))))
  
  remap$LvlA_Ast <- grepl("astrocyte", remap$Celltypes)
  remap$LvlB_BrainRobust <- (nBrain >= 2)
  remap$LvlB_BrainSingle <- (nBrain == 1)
  remap$LvlO_OtherRobust <- (nBrain == 0 & nNonBrain >= 2)
  
## Save
  write.csv(remap, file = "TF/Remap Processed.csv", row.names = FALSE)
  
  
## Convert to a wide format
  # focus on: expressed in one brain cell-type, or 2+ non-brain cell-types
  x <- remap[which(remap$LvlO_OtherRobust | remap$LvlB_BrainSingle | remap$LvlB_BrainRobust | remap$LvlA_Ast),]
  # x <- remap[which(remap$LvlO_OtherRobust | remap$LvlB_BrainSingle),] 
  
  remap_wide <- candidate.annot[,1:4] # z is a dataframe that will be filled with output from Tobias
  remap_wide[,unique(x$TF)] <- NA # count columns for each expressed TF
  
  
  # fill the matrix
  for (j in unique(x$TF)) {
    
    z <- x[which(x$TF == j),] # filter to peaks-motif pairs
    z$Enh <- factor(z$Enh, levels = remap_wide$Enh)
    remap_wide[,j] <- table(z$Enh) %>% as.numeric() # this counts 
    
  }
  
## Stats
  remap.tf.enrich <- candidate.enrich[0,] # retains only structure and column names
  
  # run
  for (j in unique(x$TF)) {
    remap.tf.enrich <- run.EnhancerFisher(data = remap_wide, data.column = j, rbind = TRUE, rbind.to = remap.tf.enrich)
  }
  
  # fdr correct
  remap.tf.enrich$FDR <- p.adjust(remap.tf.enrich$p, method = "fdr")
  remap.tf.enrich <- remap.tf.enrich[order(remap.tf.enrich$p),]
  
  # save
  write.csv(remap.tf.enrich, file = "TF/Remap - Enrichment.csv")
    
    
## Compare to TOBIAS
  # which tfs are found in both datasets?
  common_TFs <- intersect(colnames(remap_wide), tob_TF_exp)
  
  # fisher test of overlap
  tf_consistency <- data.frame(Tobias_Bound_TRUE = NA,
                      Tobias_Unbound_TRUE = NA,
                      Remap_TRUE = NA,
                      Jaccard_Bound = NA,
                      Jaccard_Unbound = NA,
                      OR_Bound = NA,
                      OR_Unbound = NA,
                      p_Bound = NA,
                      p_Unbound = NA)[0,]
  
  # run fisher test
  subset <- which(remap_wide$Tested)
  
  for (j in common_TFs) {
    
    if (j == "ZBED1") next
    
    total <- length(subset)
    x <- tob$bound[,j] > 0
    y <- tob$unbound[,j] > 0
    z <- remap_wide[,j] > 0

    # run stats
    fx <- table(x, z) %>% fisher.test()
    fy <- table(y, z) %>% fisher.test()

    out <- data.frame(Tobias_Bound_TRUE = sum(x),
                      Tobias_Unbound_TRUE = sum(y),
                      Remap_TRUE = sum(z),
                      Jaccard_Bound = sum(x & z) / sum(x | z), 
                      Jaccard_Unbound = sum(y & z) / sum(y | z), 
                      OR_Bound = fx$estimate,
                      OR_Unbound = fy$estimate,
                      p_Bound = fx$p.value,
                      p_Unbound = fy$p.value,
                      row.names = j)
    
    tf_consistency[j,] <- out
    
  }
  
  tf_consistency$Sig_Bound <- p.adjust(tf_consistency$p_Bound, "fdr") < 0.05
  tf_consistency$Sig_Unbound <- p.adjust(tf_consistency$p_Unbound, "fdr") < 0.05
  
  write.csv(tf_consistency, file = "TF/Annotation Comparison Across Methods.csv")
  
  
  ## Plot
    # what is the jaccard
    pdf(file = "TF/Annotation Consistency - Jaccard.pdf", height = 2, width = 2)
    boxplot(tf_consistency$Jaccard_Bound)
    dev.off()
    
    # odds ratio
    tf_consistency$Cut <- cut(tf_consistency$OR_Bound, c(-1, 1, 2, 5, 10, Inf), include.lowest = FALSE)
    levels(tf_consistency$Cut) <- c("<1", "1-2", "2-5", "5-10", "10+")
    
    pdf(file = "TF/Annotation Consistency - Odds Ratio.pdf", height = 3, width = 4)
    ggplot(tf_consistency, aes(x = Cut)) +
      geom_bar() +
      theme_bw() +
      labs (x = "Consistency in ReMap and TOBIAS TF Annotation\nOdds Ratio Bin", y = "Number of TFs") +
      scale_y_continuous(expand = c(0,0))
    dev.off()
      
  # ## And compare calls
  #   common_TFs <- intersect(colnames(remap_wide), tob_TF_exp)
  #   
  #   x <- remap.tf.enrich[common_TFs,]
  #   y <- candidate.tf.enrich[common_TFs,]
  #   z <- data.frame(TF = common_TFs, 
  #                   OR_T = y$OR, 
  #                   OR_R = x$OR,
  #                   P_T = -log10(y$p),
  #                   P_R = -log10(x$p),
  #                   Sig_T = y$FDR < 0.001,
  #                   Sig_R = x$FDR < 0.001)
  #   
  #   table(z$Sig_T, z$Sig_R) %>% fisher.test()
  # 
  #   cor(z[,-1])
  
    
  
################################################################################################################################ #
## Save ----
  
write.csv(candidate.enrich, file = "Final - Enrichments.csv")
write.csv(candidate.annot, file = "Final - Annotation Logical.csv", row.names = FALSE)
save(candidate.enrich, candidate.annot, file = "Final.rda")
# write.csv(annotation.detailed, file = "Final - Annotation Detailed.csv", row.names = FALSE)
  
