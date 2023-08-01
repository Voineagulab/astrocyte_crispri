################################################################################################################################ #
## Setup ----


## Load
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/1_Processing/")
  load("../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv")
  load("Temp Metadata.rda")
  meta <- cbind(meta, nha@reductions$umap@cell.embeddings)
  
## The list
  filt.cells <- filt.guides <- list()


################################################################################################################################ #
## Functions ----
  
  
## Evaluation function: apply to metadata to determine which cells pass/fail and their respective properties
  # for each proposed filter, you need: 
    # number failing
    # number failing per library
    # retained guide cell assignments
    # failures vs. kept on ribo, malat1, mito, nUMI, moi
    # which barcodes fail

  thresh.cells <- function(w, label) { # where w is a vector of cells to remove
    l <- list()
    
    l$FiltBarc <- rownames(meta)[w]
    
    
    
    
    l$Stats <- data.frame(RemovedN = length(w),
                          # per.lib,
                          MeanUMI = mean(meta$UMI_Total[w]),
                          MeanMOI = mean(meta$MOI[w]),
                          MeanRibo = mean(meta$Ribo_Pct[w]),
                          MeanMALAT1 = mean(meta$MALAT1_Pct[w]),
                          MeanMito = mean(meta$Mito_Pct[w]),
                          AssignmentsLost = sum(meta$MOI[w]),
                          row.names = label)
    
    per.lib <- table(substr(l$FiltBarc, 1, 5))
    lib.total <- table(meta$Library)        
    m <- match(names(per.lib), names(lib.total))
    
                
                 
    l$FractionRemovedPerLib <- (per.lib / lib.total[m]) %>% round(3)
    

    return(l)
  }
  
  
## Plotting function: promoter suppression. to be applied to each cell.
  pos.guides <- guides$GuideID[which(guides$TargetCat == "Promoter")]
  thresh <- 3
  keep.cells <- colnames(guide.exp)[(colSums(guide.exp >= thresh) > 0)]
  check.fraction <- apply(guide.exp, 2, function(x) {
    x[which(x < thresh)] <- 0
    x <- x / sum(x)
  })
  check.fraction <- as.data.frame(check.fraction)
  
  plot.pos <- function(label) {
    for (j in pos.guides) {
      print(j)
      
      gene <- guides$TargetID[which(guides$GuideID == j)]
      
      if (!(gene %in% rownames(nha))) next
      if (gene == "EIF3A") next
      
      # basic dataframe with binary guide assignment
      p <- data.frame(Exp = nha@assays$VST@data[gene, keep.cells],
                      GuideUMI = t(guide.exp[j,keep.cells] >= thresh),
                      GuideFraction = t(check.fraction[j,keep.cells]),
                      GuideDominant = t(check.fraction[j,keep.cells] >= 0.02),
                      row.names = keep.cells)
      p[,label] <- factor(rownames(p) %in% filt.cells[[label]]$FiltBarc)
      colnames(p) <- c("Exp", "GuideUMI", "GuideFraction", "GuideDominant", label)
      levels(p[,label]) <- paste0(levels(p[,label]), " (n=", table(p[,label]), ")")
      p$GuideUMI <- factor(p$GuideUMI)
      levels(p$GuideUMI) <- paste0(levels(p$GuideUMI), " (n=", table(p$GuideUMI), ")")
      
      # plot
      print(ggplot(p, aes_string(y = "Exp", x = "GuideUMI", shape = "GuideUMI", fill = label, colour = label, colour = label)) +
              geom_violin(draw_quantiles = 0.5, width = 0.7, scale = "width", colour = "black") +
              geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.4, fill = "black") +
              scale_shape_manual(values = c(NA, 21)) +
              theme_bw() +
              scale_fill_carto_d(palette = "Earth") +
              scale_colour_carto_d(palette = "Earth") +
              scale_y_continuous() +
              labs(y = paste0("Normalised ", gene, " Expression"), x = paste0("Has ", j)) +
              theme(panel.border = invis, axis.line.y = element_line(),
                    axis.ticks.x = invis, panel.grid.major.x = invis, legend.position = "none"))
    }
  }
  
  
## Plotting function: UMAP
  umap.plot <- function(label) {
    x <- meta
    x[,label] <- rownames(x) %in% filt.cells[[label]]$FiltBarc
    
    x <- x[order(x[,label], decreasing = FALSE),]
    title <- paste0(label, " (n=", length(which(x[,label])), ")")
    ggplot(x, aes_string(x = "UMAP_1", y = "UMAP_2", colour = label, alpha = label)) +
      geom_point(size = 1) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) +
      NoLegend() +
      scale_colour_manual(values = c("grey80", "darkorange1")) +
      scale_alpha_manual(values = c(0.25, 1)) +
      labs(title = title)
  }
  
## Statistics function: positive control suppression
  calc.suppression <- function(label) {
    s <- list()
    for (j in pos.guides) {
      print(j)
      
      gene <- guides$TargetID[which(guides$GuideID == j)]
      
      if (!(gene %in% rownames(nha))) next
      if (gene == "EIF3A") next
      
      # basic dataframe with binary guide assignment
      p <- data.frame(Exp = nha@assays$VST@data[gene, keep.cells],
                      GuideUMI = t(guide.exp[j,keep.cells] >= thresh),
                      GuideFraction = t(check.fraction[j,keep.cells]),
                      GuideDominant = t(check.fraction[j,keep.cells] >= 0.02),
                      row.names = keep.cells)
      p[,label] <- (rownames(p) %in% filt.cells[[label]]$FiltBarc)
      colnames(p) <- c("Exp", "GuideUMI", "GuideFraction", "GuideDominant", label)
      
      # stats
      s[[j]] <- data.frame(Bg = median(p$Exp[which(!(p$GuideUMI))]),
                           PosNotLabel = median(p$Exp[which((p$GuideUMI & !(p[,label])))]),
                           DomNotLabel = median(p$Exp[which((p$GuideDominant & !(p[,label])))]),
                           PosAndLabel = median(p$Exp[which((p$GuideUMI & p[,label]))]),
                           DomAndLabel = median(p$Exp[which((p$GuideDominant & p[,label]))]),
                           n = paste(length(which(p$GuideUMI)),
                                     length(which(p$GuideDominant)),
                                     length(which(p$GuideUMI & p[,label])),
                                     length(which(p$GuideDominant & p[,label])), sep = "/"))
      
      colnames(s[[j]]) <- gsub("Label", label, colnames(s[[j]]))
    }
    ss <- do.call("rbind", s)
    ss[,1:5] <- round(ss[,1:5], 2)
    return(ss)
  }
  

  
################################################################################################################################ #
## Get lists of filtered cells ----
  
filt.cells <- list()
  
## With no filtering
  filt.cells$NoFilt <- thresh.cells(1:nrow(meta), "NoFilt")

## MOI > 20  
  filt.cells$HighMOI <- thresh.cells(which(meta$MOI > 20), "HighMOI")
  
## NHA_8
  filt.cells$NHA_8 <- thresh.cells(which(meta$Library == "NHA_8"), "NHA_8")
  
## DropletQC
  filt.cells$Empty <- thresh.cells(which(meta$DropletQC == "empty_droplet"), "Empty")
  filt.cells$Damaged <- thresh.cells(which(meta$DropletQC == "damaged_cell"), "Damaged")
  # filt.cells$EmptyOrDamaged <- thresh.cells(which(meta$DropletQC != "cell"), "EmptyOrDamaged")
  
## Relative to the median UMI of MOI1 cells
  # per library, > 99th percentile of MOI1 cells
  w <- list()
  for (j in unique(meta$Library)) {
    k <- meta[which(meta$Library == j),]
    m <- quantile(k$UMI_Total[k$MOI == 1], 0.99)
    w[[j]] <- rownames(k)[which(k$UMI_Total > (m))]
  }
  w <- do.call("c", w)
  w <- which(rownames(meta) %in% w)
  filt.cells$Quantile99 <- thresh.cells(w, "Quantile99")
  
################################################################################################################################ #
## Combine stats ----
  
## Combine stats
  filt.stats <- sapply(filt.cells, function(x) x$Stats) %>% t()
  
  
  
################################################################################################################################ #
## UMAP Plots ----
  
filts <- names(filt.cells)[-1]
  
plots <- lapply(filts, umap.plot)
pdf(file = "Final/Filtering - UMAP.pdf", height = 5, width = 8)
plot_grid(plotlist = plots, ncol = 3)  
dev.off()

################################################################################################################################ #
## Positive control suppression violins ----

  for (j in filts) {
    pdf(file = paste0("Final/Filtering Violins (", j, ").pdf"), height = 3, width = 4)
    print(plot.pos(j))
    dev.off()
  }


################################################################################################################################ #
## Positive control suppression median ----

filts <- names(filt.cells)[-1]

## Run
  filt.suppression <- lapply(filts, calc.suppression)
  names(filt.suppression) <- filts
  
## Plot  
  plots <- lapply(filts, plot.suppression)
  pdf(file = "Final/Filtering - Suppression.pdf", height = 5, width = 8)
  plot_grid(plotlist = plots, nrow = 2)  
  dev.off()
  
################################################################################################################################ #
## Save ----  

  save(filt.cells, filt.stats, filt.suppression, file = "Final/Filtering.rda")
  
  
  

  
  