## This script explores our TTseq data

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

  
################################################################################################################################ #
## FANTOM5 ----
  
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
    
  ## Enhancer expression across celltypes
    # fantom.exp <- read.delim("../../../PublicData/FANTOM5/F5.hg38.enhancers.expression.tpm.matrix") # one row per enhancer id, one column per sample 
    # fantom.exp <- fantom.exp[,fantom.sampsAst]
    # fantom.exp <- rowMeans(fantom.exp) # mean expression
    
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
   
  ## There are 20 screen peaks that overlap >1 (i.e. 2) F5 enhancers
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
    
  ## Plot
    p <- apply(fantom.overlap[,-c(1:2)], 2, function(x) {
      tab <- table(fantom.overlap$Hit, x > 0) # note this coerces the logical, thus is valid
      tab <- tab / rowSums(tab)
      
      tab <- as.data.frame(tab)
      colnames(tab) <- c("Hit", "FANTOM5", "Freq")
      
      return(tab)
      
    })
    
    p <- do.call("rbind", p)
    p$Sample <- splitter(rownames(p), "\\.", 1)
    p$Sample <- gsub("Used", "", p$Sample) %>% 
      gsub("FANTOM5", "Any Sample", .) %>%
      gsub("Ast", "Astrocytes", .)
    p <- p[-which(p$FANTOM5 == "FALSE"),]
    levels(p$Hit) <- c("ns Peak", "Hit Enhancer")
    
    pdf(file = "FANTOM5/Overlap Frequency.pdf", height = 3.5, width = 4)
    ggplot(p, aes(x = Sample, y = Freq, fill = Hit)) +
      geom_col(position = "dodge", colour = "black", width = 0.75) +
      theme_bw() +
      scale_fill_lancet() +
      scale_y_continuous(limits = c(0,0.8), expand = c(0,0)) +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis,
            legend.position = c(0.8, 0.7)) +
      guides(fill = guide_legend(title = "Screen Peak")) +
      labs(y = "Fraction of Screen Peaks Overlapping\nFANTOM5 Transcribed Enhancers", x = "FANTOM5 Sample Set")
    dev.off()
  
################################################################################################################################ #
## TTseq Feature Count at Enhancers ----
  
## Prepare input files
  # convert peak bed into saf
  bed <- read.bed(nha_dir_38)
  colnames(bed) <- c("chr", "start", "end", "id")
  bed$id <- convert.id(bed$id)
  
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
  
  
## Run feature counts
  ## Setup functions and directories
    # bam_tt <- "/mnt/Data0/PROJECTS/GWP/rnaseq_GOK10844/results-circ/star_salmon/GOK10844A1.markdup.sorted.bam"
    bam_tt <- "/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/TTseq/NFcore_RNAseq_Basic/star_salmon/GOK10844A1.markdup.sorted.bam"
    # bam_rnaseq <- "/mnt/Data0/PROJECTS/GWP/rnaseq_GOK10844/results-circ/star_salmon/GOK10856A1.markdup.sorted.bam"
    bam_rnaseq <- "/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/TTseq/NFcore_RNAseq_Basic/star_salmon/GOK10856A1.markdup.sorted.bam"
    
    run_fc <- function(saf, saf_dir, stranded = 1) {
      
      write.table(saf, saf_dir, sep = "\t", row.names = FALSE, quote = FALSE)
      
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
    counts$Unstranded <- run_fc(saf$Unstranded, "TTseq/SAF_Unstranded_w1000.saf", stranded = 0)
    counts$Pos <- run_fc(saf$Pos, "TTseq/SAF_Pos_w1000.saf", stranded = 1)
    counts$Neg <- run_fc(saf$Neg, "TTseq/SAF_Neg_w1000.saf", stranded = 1)
    
  ## Save
    save(counts, file = "TTseq/FeatureCounts.rda")
    # load("TTseq/FeatureCounts.rda")
    
  ## Explore stats
    seqStats <- counts$Unstranded$stat
    colnames(seqStats) <- c("Status", "TTseq", "RNAseq")
    libSize_M <- colSums(seqStats[,-1]) / 10 ^ 6
    
  
## Process
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
    x$TPM_TTseq_Total <- x$TTseq_Total / libSize_M["TTseq"] # note that no need to correct for length, as the size is 1kb
    x$TPM_TTseq_Pos <- x$TTseq_Pos / libSize_M["TTseq"]
    x$TPM_TTseq_Neg <- x$TTseq_Neg / libSize_M["TTseq"]
    x$TPM_RNAseq_Total <- x$RNAseq_Total / libSize_M["RNAseq"] 
    x$TPM_RNAseq_Pos <- x$RNAseq_Pos / libSize_M["RNAseq"] 
    x$TPM_RNAseq_Neg <- x$RNAseq_Neg / libSize_M["RNAseq"] 
 
  ## Comparisons across libraries
    # strand bias
    x$FractionPos_TTseq <- x$TTseq_Pos / (x$TTseq_Total)
    x$FractionPos_RNAseq <- x$RNAseq_Pos / (x$RNAseq_Total)
    x$FractionPos_TTseq[which(is.na(x$FractionPos_TTseq))] <- 0 # this and line below: as NA was returned when both strands 0, set to 0
    x$FractionPos_RNAseq[which(is.na(x$FractionPos_RNAseq))] <- 0
    
    # TTseq to RNAseq ratio
    x$Ratio_TTversusRNA <- (x$TTseq_Total + 0.01) / (x$RNAseq_Total + 0.01)
    x$PassThresh_TT <- x$TTseq_Total >= 3
    x$PassThresh_RNA <- x$RNAseq_Total >= 3
    x$TT_Enrich <- (x$Ratio_TTversusRNA > 1) & (x$PassThresh_TT)

  ## Save
    featCounts <- x
    write.csv(featCounts, file = "TTseq/Results Table.csv")
    
    
## Plots
  ## Distribution of expression
    # CPM
    # x <- exon_counts_TPM
    # x$Hit <- rownames(x) %in% res.final$Gene[which(res.final$HitPermissive)]
    # x <- x[,c("Hit", "TTseq", "RNAseq")]
    # x$Type <- "Gene"
    
    p <- data.frame(Hit = featCounts$Hit,
                    TTseq = (featCounts$TPM_TTseq_Total),
                    RNAseq = (featCounts$TPM_RNAseq_Tota))
    lim <- max(c(p$TTseq, p$RNAseq))
    
    pdf(file = "TTseq/TPM Scatterplot.pdf", height = 3, width = 3.5)
    ggplot(p, aes(x = log2(RNAseq+(2^-5)), y = log2(TTseq+(2^-5)), colour = Hit)) +
      geom_point() +
      scale_colour_lancet() +
      theme_bw() +
      # facet_wrap(~Hit) +
      # geom_smooth(method = "lm") +
      theme(panel.border = invis, axis.line = element_line()) +
      labs() +
      geom_abline(intercept = 0, slope = 1, linetype = 2) #+
      # scale_y_continuous(expand = c(0,0), limits = c(0, lim)) +
      # scale_x_continuous(expand = c(0,0), limits = c(0, lim)) 
    # , trans = pseudo_log_trans(sigma = 0.1)
    dev.off()
      
    
    p <- melt(p)
    pdf(file = "TTseq/TPM.pdf", height = 3, width = 3)
    ggplot(p, aes(x = variable, colour = Hit, fill = Hit, y = value)) +
      geom_violin(scale = "width", fill = "white", position = position_dodge(width = 0.7), width = 0.7, alpha = 0.2) +
      geom_quasirandom(dodge.width = 0.7, alpha = 0.2, size = 1) +
      stat_summary(fun = mean, geom = "point", shape = "-", colour = "black", size = 8, position = position_dodge(width = 0.7)) +
      # stat_summary(fun = mean, geom = "point", shape = "+", colour = "white", size = 4, position = position_dodge(width = 0.7)) +
      theme_bw() +
      scale_fill_lancet() +
      scale_colour_lancet() +
      scale_y_continuous(expand = c(0,0)) +
      theme(axis.title.x = invis, panel.grid.major.x = invis, panel.border = invis,
            axis.line.y = element_line()) +
      labs(y = "TPM")
    dev.off()
      
  # ## Expression thresholds
  #   thresh <- c(0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 1)
  #   p <- list()
  #   for (j in thresh) {
  #     tab <- table(featCounts$Hit, featCounts$TPM_TTseq_Total >= j)
  #     frac <- tab / rowSums(tab)
  #     frac <- frac[,2]
  #     
  #     fish <- fisher.test(tab)
  #     or <- signif(fish$estimate, 2)
  #     pval <- fish$p.value %>% formatC(format = "e", digits = 1)
  #     
  #     
  #     
  #     i <- paste0("Thresh_", j)
  #     p[[i]] <- data.frame(TPM = j,
  #                          Frac_ns = frac["FALSE"],
  #                          Frac_Hit = frac["TRUE"],
  #                          OddsRatio = or,
  #                          Lower = fish$conf.int[1],
  #                          Upper = fish$conf.int[2],
  #                          p = pval)
  #   }
  #   
  #   p <- do.call("rbind", p)
  #   p$TPM <- factor(p$TPM)
  #   
  #   q <- melt(p[,1:3], id.vars = "TPM")
  #   pdf(file = "TTseq/Threshold Comparison.pdf", height = 4, width = 6)
  #   ggplot(q, aes(x = TPM, fill = variable, y = value)) +
  #     geom_col(alpha = 0.8, colour = "black", width = 0.7, position = "dodge") +
  #     theme_bw() +
  #     labs(y = "Fraction of Peaks Above TPM Threshold",
  #          x = "TTseq TPM Threshold") +
  #     scale_fill_lancet() +
  #     scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  #     theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.8, 0.7))
  #   
  #   
  #   p$p_bin <- cut(as.numeric(p$p), c(0,1e-10, 1e-5, 1e-2, 1))
  #   levels(p$p_bin) <- c("< 1e-10", "< 1e-5", "< 1e-2", "ns")
  #   p_colours <- carto_pal(7, "Geyser")[(c(7,6,4,1))]
  #   names(p_colours) <- levels(p$p_bin)
  # 
  #   
  #    ggplot(p, aes(x = TPM, y = OddsRatio, size = Frac_Hit*100, ymin = Lower, ymax = Upper, fill = p_bin)) +
  #       # geom_segment(mapping = aes(xend = Study, y = 0, yend = OR), size = 0.5) +
  #       geom_errorbar(width = 0.2, size = 0.5, colour = "black", alpha = 0.5) +
  #       geom_point(shape = 21, colour = "black") +
  #       scale_size_continuous(limits = c(0,90), breaks = c(0, 10, 25, 40, 60, 90)) +
  #       scale_fill_manual(values = p_colours, limits = names(p_colours)) +
  #       theme_bw() +
  #       scale_y_continuous(limits = c(0, 10), expand = c(0,0)) +
  #       geom_hline(yintercept = 1, linetype = 2) +
  #       theme(legend.position = "right", legend.box = "horizontal", panel.border = invis, 
  #             axis.line.y = element_line(), axis.ticks.x = invis, panel.grid.major.x = invis) +
  #       labs(y = "Odds Ratio", x = "TTseq TPM Threshold")
  #  dev.off()
   
  ## Expression thresholds, but taking strand into account
    # thresh <- c(0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 1)
    thresh <- c(2, 3, 5, 10, 25)
    p <- list()
    for (j in thresh) {
      
      # classify reads passing threshold on each strand
      strands <- apply(featCounts[,c("TTseq_Pos", "TTseq_Neg")], 1, function(x) {
        length(which(x >= j))
      })
      
      strands <- factor(strands)
      
      # tabulate versus hits
      tab <- table(featCounts$Hit, strands)
      frac <- tab / rowSums(tab)
      # frac <- frac[,2]
      
      ## Statistics
        # on the full matrix
        fish_all <- fisher.test(tab)
        pval_all <- fish_all$p.value %>% formatC(format = "e", digits = 1)
      
        # rate of bidirectionality 
        fish_bi <- table(featCounts$Hit, strands == 2) %>% fisher.test()
        pval_bi <- fish_bi$p.value %>% formatC(format = "e", digits = 1)
        or_bi <- signif(fish_bi$estimate, 2)
        or_bi_lower <- fish_bi$conf.int[1]
        or_bi_upper <- fish_bi$conf.int[2]
        
        # rate of bidirectional | unidirectional
        fish_uni <- table(featCounts$Hit, strands != 0) %>% fisher.test()
        pval_uni <- fish_uni$p.value %>% formatC(format = "e", digits = 1)
        or_uni <- signif(fish_uni$estimate, 2)
        or_uni_lower <- fish_uni$conf.int[1]
        or_uni_upper <- fish_uni$conf.int[2]
      
      # output
      i <- paste0("Thresh_", j)
      p[[i]] <- data.frame(Reads = j,
                           Uni_ns = frac["FALSE", "1"],
                           Uni_Hit = frac["TRUE", "1"],
                           Bi_ns = frac["FALSE", "2"],
                           Bi_Hit = frac["TRUE", "2"],
                           p_all = pval_all,
                           p_uni = pval_uni,
                           or_uni = or_uni,
                           or_uni_lower = or_uni_lower,
                           or_uni_upper = or_uni_upper,
                           p_bi = pval_bi,
                           or_bi = or_bi,
                           or_bi_lower = or_bi_lower,
                           or_bi_upper = or_bi_upper)
    }
    
    p <- do.call("rbind", p)
    p$Reads <- factor(p$Reads)
    
    q <- melt(p[,1:5], id.vars = "Reads")
    q$Hit <- splitter(q$variable, "_", 2)
    q$Hit <- factor(q$Hit, levels = c("Hit", "ns"))
    levels(q$Hit) <- c("Hit Enhancer", "ns Peak")
    q$Strand <- splitter(q$variable, "_", 1)
    q$Strand <- factor(q$Strand, levels = c("Uni", "Bi"))
    levels(q$Strand) <- c("One", "Two")
   
    # fills <- carto_pal(7, "Tropic")[c(3,1,5,7)]
    # fills <- carto_pal(7, "Geyser")[c(2,1,6,7)]
    # fills <- carto_pal(7, "Earth")[c(2,1,6,7)]
    fills <- c(carto_pal(7, "Geyser")[c(6,7)], carto_pal(7, "Earth")[c(6,7)])
    
    pdf(file = "TTseq/Threshold Comparison (Stranded).pdf", height = 4, width = 6)
    ggplot(q, aes(x = Hit, fill = interaction(Strand, Hit), y = value)) +
      geom_col(alpha = 0.8, width = 0.9, position = "stack") +
      theme_bw() +
      facet_wrap(~Reads, scales = "free_x", nrow = 1, strip.position = "bottom") +
      labs(y = "Fraction of peaks above threshold", x = "TTseq read threshold") +
      scale_fill_manual(values = fills, labels = c("Hit Enhancer / Unidirectional", 
                                                   "Hit Enhancer / Bidirectional", 
                                                   "ns Peak / Unidirectional", 
                                                   "ns Peak / Bidirectional")) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      guides(fill = guide_legend(ncol = 2)) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = c(0.65, 0.95),
            panel.grid = invis, axis.text.x = invis, axis.ticks.x = invis, 
            strip.background = element_rect(fill = "white", linetype = 0),
            legend.title = invis)
   dev.off()
    

  p_bins <- c(0, 1e-15, 1e-10, 1e-5, 1)
  q <- p
  q$p_uni <- cut(as.numeric(q$p_uni), p_bins)
  q$p_bi <- cut(as.numeric(q$p_bi), p_bins)
  levels(q$p_uni) <- levels(q$p_bi) <- c("< 1e-15", "< 1e-10", "< 1e-5", "ns")
  # p_colours <- carto_pal(7, "Temps")[(c(7,6,4,1))]
  # p_colours <- pal_lancet()(9)[c(1,2,3,9)]
  p_colours <- c(carto_pal(2, "TealRose"), carto_pal(2, "Earth")[1], "black")
  names(p_colours) <- levels(q$p_bi)
  q <- q[,-c(2,4,6)]
  
  g <- grep("bi", colnames(q), ignore.case = TRUE)
  q <- list(Uni = q[,c(-g)],
            Bi = q[,c(1, g)])
  q <- lapply(q, function(x) {
    colnames(x) <- c("Reads", "Fraction_Hit", "P", "OR", "Lower", "Upper")
    return(x)
  })
  q <- do.call("rbind", q)
  q$Strand <- splitter(rownames(q), "\\.", 1) %>% paste0(., "directional")
  
  # q <- melt(q, id.vars = "Reads")
  
  pdf(file = "TTseq/Threshold Comparison (Stranded) Odds Ratio.pdf", height = 4, width = 6)
  ggplot(q, aes(x = Reads, y = OR, size = Fraction_Hit*100, colour = P, ymin = Lower, ymax = Upper, fill = P, shape = Strand, linetype = Strand)) +
    # geom_segment(mapping = aes(xend = Study, y = 0, yend = OR), size = 0.5) +
    geom_errorbar(width = 0.2, size = 0.5, colour = "black", alpha = 0.5, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_size_continuous(limits = c(1,50), breaks = c(1, 15, 30, 50)) +
    scale_shape_manual(values = c(23, 21)) +
    scale_linetype_manual(values = c("solid", "solid")) +
    scale_fill_manual(values = p_colours, limits = names(p_colours)) +
    scale_colour_manual(values = p_colours, limits = names(p_colours)) +
    theme_bw() +
    # scale_y_continuous(limits = c(NA, NA), expand = c(0,0), trans = pseudo_log_trans(), breaks = c(0.1, 0.5, 1, 2, 4, 8, 16, 64, 256)) +
    scale_y_continuous(limits = c(0, 11.5), expand = c(0,0)) +
    geom_hline(yintercept = 1, linetype = 2) +
    theme(legend.position = "right", panel.border = invis, 
          axis.line.y = element_line(), axis.ticks.x = invis, panel.grid.major.x = invis) +
    labs(y = "Enrichment of Hits (Odds Ratio)", x = "TTseq read threshold")
  dev.off()
    
  ## TTseq enrichment of enhancers
    p <- table(featCounts$Hit, featCounts$TT_Enrich)
    fish <- fisher.test(p)
    or <- signif(fish$estimate, 2)
    pval <- fish$p.value %>% formatC(format = "e", digits = 1)
    p <- p / rowSums(p)
    p <- as.data.frame(p)
    colnames(p) <- c("Hit", "TTseq", "Freq")
    p <- p[which(p$TTseq == "TRUE"),]
    
    pdf(file = "TTseq/TTseq Enrichment.pdf", height = 3, width = 2.5)
    ggplot(p, aes(x = Hit, y = Freq, fill = Hit)) +
      geom_col(alpha = 0.8, colour = "black", width = 0.7) +
      theme_bw() +
      labs(y = "Fraction of Peaks with TTseq\nEnrichment (Counts >= 3)", 
           x = paste0("Peak is Active Enhancer\nOR=", or, ", p=", pval)) +
      scale_fill_lancet() +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      theme(panel.border = invis, axis.line.y = element_line(), legend.position = "none")
    dev.off()
    
      
    
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

  
  
################################################################################################################################ #
## TTseq QC ----
  
    
## For this, you shall compare exonic versus intronic counts


## Exonic counts
  gtf <- "/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/cellranger_refdata-gex-GRCh38-2020-A/genes/genes.gtf"
  exon_counts <- featureCounts(c(bam_tt, bam_rnaseq),
                                 annot.ext = gtf, isGTFAnnotationFile = TRUE,
                                 useMetaFeatures = TRUE,
                                 GTF.attrType = "gene_name",
                                 allowMultiOverlap = TRUE,
                                 isPairedEnd = TRUE, # self-explanatory
                                 nthreads = 8, # self-explanatory
                                 strandSpecific = 2, # reverse-stranded, suitable for Tru-seq-based libraries
                                 minMQS = 10, # minimum quality score at either end of the pair
                                 checkFragLength = FALSE,
                                 countMultiMappingReads = FALSE,
                                 requireBothEndsMapped = FALSE,
                                 countChimericFragments = FALSE)

    save(exon_counts, file = "TTseq/FeatureCounts_Exons.rda")

    exon_counts$counts <- as.data.frame(exon_counts$counts)
    colnames(exon_counts$counts) <- c("TTseq", "RNAseq")
    exon_counts$TPM <- exon_counts$counts
    exonicLengths <- exon_counts$annotation$Length / 1000

    exon_counts$TPM$TTseq <- exon_counts$TPM$TTseq / exonicLengths / libSize_M["TTseq"]
    exon_counts$TPM$RNAseq <- exon_counts$TPM$RNAseq / exonicLengths / libSize_M["RNAseq"]

    exon_counts_TPM <- exon_counts$TPM
    exon_counts_TPM <- exon_counts_TPM[unique(res.final$Gene),] # only genes tested in the screen


## Intron counts
  ## Make bed of intron coordinates
    w <- exon_counts$annotation
    # w <- w[which(w$GeneID %in% res.final$Gene),]

    intron_bed <- apply(w, 1, function(x) {

      ## Wrangle into one row per exon
        names(x) <- colnames(exon_counts$annotation)
        print(x["GeneID"])
        y <- data.frame(GeneID = x["GeneID"],
                        Chr = splitter(x["Chr"], ";", 1),
                        Start = do.call("c", strsplit(x["Start"], ";")),
                        End = do.call("c", strsplit(x["End"], ";")),
                        Strand = splitter(x["Strand"], ";", 1))

      ## Make bed files
        # the first bed: the entire gene body
        body.start <- min(as.numeric(y$Start))
        body.end <- max(as.numeric(y$End))
        bed.body <- data.frame(Chr = unique(y$Chr),
                               Start = body.start,
                               End = body.end)
        write.bed(bed.body, "../Scratchspace/TTseq_geneBody.bed")

        # second bed: each exon
        bed.exon <- y[,c(2,3,4,1,5)]
        write.bed(bed.exon, "../Scratchspace/TTseq_geneExon.bed")


      ## Intersect using bedtools, to subtract the exons from the entire gene body
        call <- paste("subtractBed",
                  "-a", "../Scratchspace/TTseq_geneBody.bed",
                  "-b", "../Scratchspace/TTseq_geneExon.bed",
                  ">", "../Scratchspace/TTseq_geneIntron.bed")

        system(call, intern = FALSE, wait = TRUE)

      ## Read in
        z <- tryCatch(read.bed("../Scratchspace/TTseq_geneIntron.bed"),
                      error = function(x) NA )

        # return NA is the file is empty (when exons tile the entire gene body)
        if (is.na(z)) return(z)

        # otherwise, process and return
        colnames(z) <- c("Chr", "Start", "End")
        z$ID <- unique(y$GeneID)
        z$Strand <- unique(y$Strand)
        return(z)

    })

    intron_bed <- do.call("rbind", intron_bed)
    intron_bed <- intron_bed[-which(is.na(intron_bed$Chr)),]

  ## Make saf from bed
    intron_saf <- data.frame(GeneID = intron_bed$ID,
                      Chr = intron_bed$Chr,
                      Start = intron_bed$Start,
                      End = intron_bed$End,
                      Strand = intron_bed$Strand)
    write.table(intron_saf, "../Scratchspace/TTseq_geneIntron.saf", sep = "\t", row.names = FALSE, quote = FALSE)

  ## Feature counts
    intron_counts <- featureCounts(c(bam_tt, bam_rnaseq),
                                 annot.ext = "../Scratchspace/TTseq_geneIntron.saf", isGTFAnnotationFile = FALSE,
                                 useMetaFeatures = TRUE,
                                 GTF.attrType = "gene_name",
                                 allowMultiOverlap = TRUE,
                                 isPairedEnd = TRUE, # self-explanatory
                                 nthreads = 8, # self-explanatory
                                 strandSpecific = 2,
                                 minMQS = 10, # minimum quality score at either end of the pair
                                 checkFragLength = FALSE,
                                 countMultiMappingReads = FALSE,
                                 requireBothEndsMapped = FALSE,
                                 countChimericFragments = FALSE)

    save(intron_counts, file = "TTseq/FeatureCounts_Introns.rda")

    intron_counts$counts <- as.data.frame(intron_counts$counts)
    colnames(intron_counts$counts) <- c("TTseq", "RNAseq")
    intron_counts$TPM <- intron_counts$counts
    intronicLengths <- intron_counts$annotation$Length / 1000

    intron_counts$TPM$TTseq <- intron_counts$TPM$TTseq / intronicLengths / libSize_M["TTseq"]
    intron_counts$TPM$RNAseq <- intron_counts$TPM$RNAseq / intronicLengths / libSize_M["RNAseq"]

    intron_counts_TPM <- intron_counts$TPM


## Plot
   
  # scatterplot ttseq and rnaseq expression for each gene for exons and introns
  p <- intron_counts$counts
  # p <- intron_counts_TPM
  p <- p[-which(apply(p, 1, max) == 0),] # remove zeros in both
  p <- p[which(rownames(p) %in% (res.final$Gene)),]
  # p <- log2(p + (2^-10))
  p <- log2(p+1)
  #
  p$Density <- get_density(p$TTseq, p$RNAseq, n = 100)

  pdf(file = "TTseq/QC - Intronic Counts.pdf", height = 3.5, width = 4.5)
  ggplot(p, aes(x = RNAseq, y = TTseq, colour = Density)) +
    geom_point() +
    theme_bw() +
    scale_colour_viridis_c() +
    geom_abline(linetype = 2, colour = "red", size = 1) +
    labs(x = "RNASeq Intronic Counts\n(log2, offset = 1)",
         y = "TTseq Intronic Counts\n(log2, offset = 1)") +
    theme(panel.border = invis)
  dev.off()

  # q <- melt(p[,-3])
  # ggplot(q, aes(x = variable, y = value)) +
  #   geom_violin(scale = "width") +
  #   geom_boxplot(outlier.shape = NA, width = 0.3)

# ## Compare exon counts to intron counts ratio
#   # get counts
#   common <- intersect(exon_counts$annotation$GeneID, intron_counts$annotation$GeneID)
#   p <- data.frame(Exon_TTseq = exon_counts$TPM[common, "TTseq"],
#                   Exon_RNAseq = exon_counts$TPM[common, "RNAseq"],
#                   Intron_TTseq = intron_counts$TPM[common, "TTseq"],
#                   Intron_RNAseq = intron_counts$TPM[common, "RNAseq"],
#                   row.names = common)
# 
#   # threshold
#   p <- p[which(rownames(p) %in% res.final$Gene),] # to tested genes
#   # p <- p[which(apply(p, 1, min) >= 5),] # at least 5 counts
#   # p <- p[which(apply(p, 1, min) >= 0.1),] # at least 0.5 TPM
#   p <- log2(p + 2^-4)
# 
#   
#   # plot
#   ggplot(p, aes(x = Intron_RNAseq, y = Intron_TTseq)) +
#     geom_point()
#   
#   # wrangle
#   p$Ratio_TTseq <- p$Intron_TTseq / p$Exon_TTseq
#   p$Ratio_RNAseq <- p$Intron_RNAseq / p$Exon_RNAseq
#   
#   

  
  