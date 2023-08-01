## Recalculating fold-change!

## Loop
  pb <- list()
  g <- unique(res.final$Gene)

  # for each enhancer
  a <- Sys.time()
  for (j in unique(res.final$Enh)) {
    
    print(j)
    
    ## Get pseudobulk of targeted and non-targeted cells
      targeted <- which(nha@meta.data[,j])
      
      # raw counts
      e <- list()
      e$targeted <- rowSums(nha@assays$RNA@counts[,targeted])
      e$nonTargeted <- rowSums(nha@assays$RNA@counts[,-targeted])
      
      # vst-stabilised counts
      e$targeted_VST <- rowSums(nha@assays$VST@counts[,targeted])
      e$nonTargeted_VST <- rowSums(nha@assays$VST@counts[,-targeted])
      
    # normalise expression to cpm
    e <- sapply(e, function(x) {
      x <- x / (sum(x) / (10^6))
      x <- x[g]
    })
    
    # filter to tested genes (after CPM)
    e <- as.data.frame(e[g,])
    colnames(e) <- paste0("CPM_", colnames(e))
    
    
    # fold-change 
    offset <- 1
    e$log2fc <- log2(e$CPM_targeted + offset) - log2(e$CPM_nonTargeted + offset)
    e$log2fc_vst <- log2(e$CPM_targeted_VST + offset) - log2(e$CPM_nonTargeted_VST + offset)
    
    # add pair test information
    e$Pair <- paste0(j, "_", rownames(e))
    e <- relocate(e, Pair)
    e$Tested <- e$Pair %in% res.final$Pair
    
    
    # output
    pb[[j]] <- e
  }
  
  b <- Sys.time()
  b-a

## Save
  save(pb, file = "Fold-changes Pseudobulk.rda")
    

## Get
  x <- do.call("rbind", pb)
  x <- x[which(x$Tested),]
  m <- match(x$Pair, res.final$Pair)
  x$logfc.old <- res.final$logfc[m]
  x$logfc.oldVST <- res.final$logfc.vst[m]
  x$Expression.Seurat <- res.final$Gene.Exp[m]
  x$Hit <- res.final$HitPermissive[m]
  write.csv(x, file = "Fold-changes Pseudobulk Filtered.csv")

## Compare  
  p <- x
  p$Hit <- factor(p$Hit, levels = c("FALSE", "TRUE"))
  levels(p$Hit) <- c("Ns", "Hit")
  
  pdf(file = "Fold-changes Pseudobulk Comparison.pdf", height = 4, width = 8)
  
## Compare expression values
  
  ggplot(p, aes(x = Expression.Seurat, y = log(CPM_nonTargeted_VST+1))) +
    geom_point() +
    facet_wrap(~Hit) +
    theme_bw() +
    labs(x = "Mean Expression (Seurat)", y = "Mean Expression (Pseudobulk log1p(CPM))")
  
  cor(p$CPM_nonTargeted_VST, p$Expression.Seurat, method = "s")
  
## Compare fold-changes
  
  ggplot(p, aes(x = logfc.oldVST, y = log2fc_vst)) +
    geom_point() +
    facet_wrap(~Hit) +
    theme_bw() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "firebrick1") +
    labs(x = "EGP Fold-change (Seurat)", y = "EGP Fold-change (Pseudobulk log2(CPM))")
  
  cor(p$logfc.oldVST, p$log2fc_vst)
  cor(p$logfc.oldVST, p$log2fc_vst, method = "s")
  
  dev.off()
  