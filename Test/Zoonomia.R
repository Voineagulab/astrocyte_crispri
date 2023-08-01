## Please load the start of script 3b (HitEnrichment Chromatin) first




## Conservation in the Zoonomia consortium
  # PhyloP scores measure evolutionary conservation at individual alignment sites. 
  # Interpretations of the scores are compared to the evolution that is expected under neutral drift.
  # Positive scores — Measure conservation, which is slower evolution than expected, at sites that are predicted to be conserved.
  # Negative scores — Measure acceleration, which is faster evolution than expected, at sites that are predicted to be fast-evolving.
  # The absolute values of phyloP scores represent –log p-values under a null hypothesis of neutral evolution
    
  ## Read in scores for each enhancer
    bed_nha_hg38_file <- import.bed(nha_dir_38)
    phylop_dir <- "../../../../PublicData/Zoonomia/241-mammalian-2020v2.bigWig" 
    phylop <- import(phylop_dir, selection = bed_nha_hg38_file, as = "NumericList")
    names(phylop) <- bed_nha_hg38_file$name
    
  ## Average at each peak
    phylop_mean <- sapply(phylop, mean)
    p <- data.frame(candidate.annot[,1:4], Mean = phylop_mean)
    wilcox.test(p$Mean ~ p$Hit) # p = 0.93
      

  ## Number of ultraconserved stretches
    # an ultraconserved stretch is on which is significantly conserved (Phylop > 2, e.g.), for 20+ bp
    

  ## Read in scores for KLF5
    

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
    klf5_motifs <- footprints[which(footprints$Motif == "KLF5"),]
    klf5_motifs$ID <- paste0(klf5_motifs$Enh, "_", klf5_motifs$Footprint, "_", 1:nrow(klf5_motifs))
    klf5_motifs <- klf5_motifs[,c(1:3, 8)]
    write.bed(klf5_motifs, "Zoonomia_KLF5_in.bed")
    bed_klf5_file <- import.bed("Zoonomia_KLF5_in.bed")
    
    klf5_phylop <- import(phylop_dir, selection = bed_klf5_file, as = "NumericList")
    names(klf5_phylop) <- bed_klf5_file$name
    
    # plot
    x <- do.call("rbind",klf5_phylop) %>% as.data.frame()
    colnames(x) <- paste0("Pos_", 1:ncol(x))
    x$Enh <- splitter(rownames(x), "_", 1)
    x$Bound <- splitter(rownames(x), "_", 2)
    x$Hit <- x$Enh %in% hit.enh
    
    y <- melt(x)
    # y$variable <- sub("Pos_", "", y$variable) %>% as.numeric()
    # y <- y[which(y$Bound == "bound"),]
    
    ggplot(y, aes(x = variable, y = value, fill = Bound)) +
       # geom_line(alpha = 0.8) +
      geom_violin(draw_quantiles = 0.5, width = 0.7, scale = "width") +
      # stat_summary(fun = mean, geom = "point", shape = 21, colour = "black", size = 10) +
      theme_bw() +
      facet_wrap(~Hit) +
      theme(panel.grid = invis, panel.border = invis, axis.line = element_line(), legend.title = invis,
            legend.position = "bottom") +
      scale_colour_lancet() +
      geom_hline(yintercept = 0) +
      scale_y_continuous(expand = c(0,0), breaks = -5:5) +
      labs(y = "ATAC Signal", x = "Position relative to motif centre (bp)")