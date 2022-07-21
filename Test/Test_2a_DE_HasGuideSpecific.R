## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Seurat)
  library(sceptre)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(rcartocolor)

## Load
  source("../../Scripts/Functions.R")
  load("../../Data/Preprocessed/NHA Pooled.rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist.csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]

## Data information
  samples <- c(paste0("NHA_", c(1:5, 7:8)))
  sample.colours <- carto_pal(7, "Bold")
  names(sample.colours) <- samples
  
  pos <- guides$TargetID[which(guides$TargetCat == "Promoter")]
  pos <- unique(pos)
  enh <- guides$TargetID[which(guides$TargetCat == "Enh")]
  enh <- unique(enh)
  
  wilcox <- sceptre <- lp <- list()

## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()


################################################################################################################################ #
## Run Sceptre ----

## Load
  load("../2_DE/Sceptre NHA Input.rda")
  # load("../2_DE/Sceptre.rda")
  
## Make guide-gene test matrix
    sceptre.pairs.enh.guide <- list()
    for (j in guides$GuideID[which(guides$TargetCat == "Enh")]) {
      print(j)
      
      # get enh coordinates
      m <- match(j, guides$GuideID)
      w <- guides$TargetCoord[m]
      
      
      # find nearby tss to the enh
      nearby <- find.nearby.tss(query = w, expand.by = 10^6, reduced.output = TRUE)
      nearby <- intersect(nearby, rownames(nha))
      sceptre.pairs.enh.guide[[j]] <- data.frame(gene_id = nearby, gRNA_id = j, pair_type = "candidate") # recommended name for pair-type
    }
    
    sceptre.pairs.enh.guide <- do.call("rbind", sceptre.pairs.enh.guide)
  
 
## Run
    sceptre.direction <- "both"
    
    # sceptre$enh.guide <- list()
    
    use <- guides$GuideID[which(guides$TargetCat == "Enh")]
    use <- use[-which(use %in% names(sceptre$enh.guide))]
    
    
    for (j in use) {
      k <- which(sceptre.pairs.enh.guide$gRNA_id == j)
      
      if (substr(j, 1, 6) == "Enh809") {
          sceptre$enh.guide[[j]] <- paste0("Skipped ", j, " manually")
          next
      }
      
      
      if (length(k) == 1) { # throws an error in this case due to it not being a matrix, so bind a random row and remove post-analysis 
        too.few <- TRUE
        k <- c(k, sample((1:nrow(sceptre.pairs.enh.guide)[-k]), 9, replace = FALSE )) # so pick 9 rows that are not k
        k <- sceptre.pairs.enh.guide[k,]  
        k$gRNA_id <- j
        k <- unique(k) # in case, by chance, the gene is duplicated 
      } else {
        too.few <- FALSE
        k <- sceptre.pairs.enh.guide[k,]  
      }
      
      # skip if too few guides, as sceptre auto-removes
      if (sum(sceptre.guide[j,]) < 30) {
        sceptre$enh.guide[[j]] <- "<30 guides, auto-removed"
        next
      }
      
      ## A vector to filter out cells with guides targeting the same enhancer
        t <- splitter(j, "_", 1)
        replicate.guides <- guides$GuideID[which(guides$TargetID == t & guides$GuideID != j)]
        
        if (length(replicate.guides) == 1) { 
          remove <- which(sceptre.guide[replicate.guides,] == 1 & sceptre.guide[j,] != 1) # has a different guide targeting the same enhancer, but not in combination with the given guide
        } else {
          remove <- apply(sceptre.guide[replicate.guides,], 2, function(x) any(x==1)) # where a cell has a different guide targeting the same enhancer...
          remove <- which(remove & sceptre.guide[j,] != 1) # ... but not the same enhancer
          
        }
        
      
  
      # run sceptre
      a <- Sys.time()
      sceptre$enh.guide[[j]] <- run_sceptre_high_moi(gene_matrix = sceptre.gene[k$gene_id,-remove], 
                                                     combined_perturbation_matrix = sceptre.guide[,-remove], # rather than sceptre.perturb
                                                     covariate_matrix = sceptre.covar[-remove,],
                                                     gene_gRNA_group_pairs = k,
                                                     side = sceptre.direction,
                                                     B = 500)
      b <- Sys.time()
      
      # when there are too few genes, remove the extraneous gene runs
      if (too.few) {
        m <- sceptre.pairs.enh.guide$gene_id[which(sceptre.pairs.enh.guide$gRNA_id == j)]
        w <- which(sceptre$enh.guide[[j]]$gene_id == m)
        if (length(w) == 0) {
          sceptre$enh.guide[[j]] <- "No Output"  
        } else {
          sceptre$enh.guide[[j]] <- sceptre$enh.guide[[j]][w,]  
        }
        
      }
      
      print(paste0(j, ": ", b-a))
      gc()
    }
  
    
## Save results
  save(sceptre, file = "Sceptre (with missing 2 single guides, redone with removing same target, May 7th).rda")
    

################################################################################################################################ #
## Analyse SCEPTRE ----
  
  
## Create dataframe, based on SCEPTRE calls
  gs <- do.call("rbind", sceptre$enh.guide) # for "guide-specific
  gs$Pair <- paste0(gs$gRNA_id, "_", gs$gene_id) 
  gs$Pair2 <- paste0(splitter(gs$gRNA_id, "_", 1), "_", gs$gene_id) # annotates with the enhancer-gene pair, rather than just guide-gene pair
  
## Rearrange
  gs <- gs[,c(2,1,6,7,5,4)]
  
  colnames(gs) <- c("Enh", "Gene", "Pair", "Pair2", "Z", "P")  
  
## Add n
  n <- rowSums(sceptre.guide)
  m <- match(gs$Enh, names(n))
  gs$N <- n[m]
  
  
## Filter
  gs <- gs[-which(gs$Enh == "<30 guides, auto-removed"), ]
  gs <- gs[-grep("Enh809", gs$Enh), ]
  

  
## Add distance and mean expression
  ## Distance
    # guide coord
    m <- match(gs$Enh, guides$GuideID)
    gs$Enh.Pos <- guides$TargetCoord[m]
    
    # gene coord
    m <- match(gs$Gene, geneInfo$Symbol)
    gs$Gene.Pos <- geneInfo$ID[m]
   
    # distance from enh midpoint to tss
    gs$Distance <- apply(gs, 1, function(y) {
      e.start <- as.numeric(splitter(splitter(y[8], "-", 1), ":", 2))
      e.end <- as.numeric(splitter(y[8], "-", 2))
      e.mid <- round((e.start + e.end) / 2)
      
      g.start <- as.numeric(splitter(splitter(y[9], "-", 1), ":", 2))
      
      dist <- abs(g.start - e.mid)
      return(dist)
    })
    
    # nearest gene
    min.dist <- aggregate(Distance~Enh, gs, min, na.rm = TRUE)
    m <- match(gs$Enh, min.dist$Enh)
    gs$Nearest <- min.dist$Distance[m]
    gs$Nearest <- gs$Distance == gs$Nearest
    
    # bin distance
    gs$Distance.Bin <- cut(gs$Distance, c(2000,10000,50000,100000,500000, max(gs$Distance) + 1)) 
    levels(gs$Distance.Bin) <- c("2-10", "10-50", "50-100", "100-500", "500-1000")
    
   

## Annotate with mean expression
  # expression normalised to library size via Seurat
  mean <- rowMeans(nha@assays$RNA@data[rownames(nha@assays$RNA@data) %in% gs$Gene,]) # rna assay, not sct
  m <- match (gs$Gene, names(mean))
  gs$Gene.Exp <- mean[m]
  
  # expgssion corrected for covariates via VST, and normalised to library size
  mean <- rowMeans(nha@assays$SCT@data[rownames(nha@assays$SCT@data) %in% gs$Gene,]) # rna assay, not sct
  m <- match (gs$Gene, names(mean))
  gs$Gene.ExpVST <- mean[m]
  
  # expression below threshold
  exp.thresh <- 2^-6
  gs$High.Expression <- gs$Gene.Exp > exp.thresh
  high <- which(gs$High.Expression)
  
## Numericise
  gs$Z <- as.numeric(gs$Z)
  gs$P <- as.numeric(gs$P)

## Hits
  gs$FDR <- p.adjust(gs$P, method = "fdr")
  gs$Hit <- gs$FDR < 0.1

## Add information on pooled enhancer analyses
  load("../2_DE/Enhancers - All Results Summary.rda")
  
  # match
  
  
  # add Z, hit, and FDR
  
  
  m <- match(gs$Pair2, res$Pair)
  gs$Pair2.Z <- res$Sceptre.Z[m]
  gs$Pair2.P <- res$Sceptre.P[m]
  gs$Pair2.FC <- res$logfc.vst[m]
  
  gs$Pair2.Hit <- gs$Pair2 %in% res$Pair[which(res$Hit)]
  
  
## Threshold
  gs.filt <- gs[which(gs$High.Expression & gs$Distance < 500000),]
  gs.filt$FDR <- p.adjust(gs.filt$P, method = "fdr")
  gs.filt$Hit <- gs.filt$FDR < 0.1
  
## Save
  save(gs, gs.filt, file = "Guide-specific - Dataframes.rda")
  write.csv(gs, file = "Guide-specific - Dataframe Unthresholded.csv")
  
################################################################################################################################ #
## Compare guide-level stats across guides, and in relation to the significance of the pooled analyses ----
  
## Setup dataframe
  any <- unique(gs.filt$Pair2[which(gs.filt$Hit | gs.filt$Pair2.Hit)]) # this collects all enh-gene pairs where at least one guide or the pool is significant. it differs from the previous "x" in that it also has some non-significant guide-gene based on another single guide reaching significance without the pool doing so
  x <- gs.filt[which(gs.filt$Pair2 %in% any),]
  
  x$Significance <- "."
  x$Significance[which(x$Pair2.Hit & !(x$Hit))] <- "Pooled Only"
  x$Significance[which(x$Hit & !(x$Pair2.Hit))] <- "Guide Only"
  x$Significance[which(x$Pair2.Hit & x$Hit)] <- "Guide & Pooled"
  x$Significance[which(!(x$Pair2.Hit | x$Hit))] <- "Another Guide Only"
      

## Plot -log10(FC)
  pdf(file = "Guide-specific - Scatterplots of Pooled Versus Single.pdf", height = 4, width = 5)
  ggplot(x, aes(x = -log10(Pair2.P), y = -log10(P), colour = Significance)) +
    geom_point() +
    theme_bw() +
    theme(panel.border = invis, axis.line = element_line()) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 17)) + 
    scale_x_continuous(expand = c(0,0), limits = c(0, 17)) +
    labs(x = "Pooled DE", y = "Guide DE", title = "-log10(P)") +
    scale_color_lancet()
  
  
## Plot Z
  ggplot(x, aes(x = Pair2.Z, y = Z, colour = Significance)) +
    geom_point() +
    theme_bw() +
    theme(panel.border = invis) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
  # scale_y_continuous(expand = c(0,0), limits = c(0, 17)) +
  # scale_x_continuous(expand = c(0,0), limits = c(0, 17)) +
    labs(x = "Pooled DE", y = "Guide DE", title = "Z") +
    scale_color_lancet()
  dev.off()  

## Replicablility across guides: the fraction of guides reaching significance
  z <- table(x$Pair2, x$Hit) %>% as.data.frame() %>% dcast(Var1 ~ Var2)
  colnames(z) <- c("Pair", "NS", "Sig")
  z$nGuide <- rowSums(z[,-1]) %>% paste0(., " Guide(s) Total")
  z$Sig <- factor(z$Sig)
  levels(z$Sig) <- c("Pool\nOnly", 1, 2, 3, 4, 5)
  z$PooledHit <- z$Pair %in% res$Pair[which(res$Hit)]

  pdf(file = "Guide-specific - Replicability Across Guides.pdf", height = 3, width = 8)
  ggplot(z, aes(x = Sig, fill = PooledHit)) +
    geom_bar(colour = "black") +
    facet_wrap(~nGuide, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme() + 
    labs(y = "Count of Enhancers", x = "Number of Significant Guides") +
    scale_y_continuous(expand = c(0,0))
    dev.off()

## Hit as a function of n
    pdf(file = "Guide-specific - Effect of nCells.pdf", height = 4, width = maxw)
    ggplot(x, aes(x = Significance, colour = Significance, y = N)) +
      geom_jitter() +
      stat_summary(fun.min = function(z) { quantile(z,0.25) }, 
                   fun.max = function(z) { quantile(z,0.75) }, 
                   fun = median, 
                   colour = "black") +
      theme_bw() +
      scale_color_lancet() +
      theme(panel.border = invis, axis.line.y = element_line()) +
      labs(y = "Number of Cells Per Guide")
    dev.off()

  
    z <- table(x$Pair2, x$Hit) %>% as.data.frame() %>% dcast(Var1 ~ Var2)
    z <- z[which(z$Var1 %in% x$Pair2[grep("Pooled", x$Significance)]),]
    z$Fraction <- z[,3] / (z[,3] + z[,2])
    table(z$Fraction >= 0.5)    
    
################################################################################################################################ #
## For enhancers that are significant with but a single guide, visualise ----
  

## Identities
  use <- x$Pair2[which(x$Significance == "Guide Only")] # no need to call unique

## Function
  violin.any.enh <- function(gene, enhs, neg = FALSE) {
    # get guides for the enhs
    pattern <- paste0(enhs, "_") %>% paste(collapse = "|")
    use <- guides[grep(pattern, guides$GuideID),]
    n <- nrow(use)

    # name guides
    names <- splitter(use$GuideID, "_chr", 1)


    p <- data.frame(VST = nha@assays$SCT@data[gene,],
                    Neg = nha@meta.data$TransductionPool == "Neg",
                    Category = "NTC",
                    nha@meta.data[,use$GuideID])

    guide.cols <- 4:(n+3)
    colnames(p)[guide.cols] <- names


    # categorise
    for (j in names) p$Category[which(p[,j])] <- j
    p$Category[which(rowSums(p[,guide.cols]) > 1)] <- "Multiple"

    if (neg) {
      p$Category[which(p$Neg)] <- "Pure Negative"
      p$Category <- factor(p$Category, levels = c("NTC", "Pure Negative", names, "Multiple"))

    } else {
      p$Category <- factor(p$Category, levels = c("NTC", names, "Multiple"))
    }

    # finally, add n
    levels(p$Category) <- paste0(levels(p$Category), "\nn=", table(p$Category))

    # colour by category
    p$Colour <- "_NTC"
    for (j in unique(use$TargetID)) {
      p$Colour[grep(paste0(j, "_"), p$Category)] <- j
    }

    enh.colours <- ggsci::pal_locuszoom()(length(unique(p$Colour)))
    names(enh.colours) <- levels(p$Colour)


    # drop cells categorised with "multiple"
    p <- p[-which(p$Category == levels(p$Category)[length(levels(p$Category))]),]

    # plot
    ggplot(p, aes(x = Category, y = log2(exp(VST)), fill = Colour, colour = Colour)) +
      geom_violin(colour = "black", scale = "width", width = 0.8, draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
      # geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2, position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = enh.colours) +
      scale_colour_manual(values = enh.colours) +
      # scale_colour_manual(values = enh.colours) +
      theme_bw() +
      labs(y = paste0(gene, " log2-Normalised Expression")) +
      scale_y_continuous(expand = c(0,0), limits = ylim) +
      # annotate("text", x = 1.5, y = max(p$Exp * 0.85), label = label, size = 2) +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid.major.x = invis,
            axis.ticks.x = invis, legend.position = "none",
            axis.title.x = invis, axis.text.x = element_text(hjust = 1, vjust = 0.3, angle = 90))
  }

## Apply
  pdf(file = "Guide-specific - Violins for Guide-Specific Pairs.pdf", height = 3, width = 4)
  for (j in use) {
    print(j)
    k <- grep(j, x$Pair2) %>% intersect(., which(x$Hit))
    label <- paste0(splitter(x$Enh[k], "_", 2), ": Z = ", round(x$Z[k], 2), ", P = ", signif(x$P[k], 2))
    
    print(violin.any.enh(splitter(j, "_", 2), splitter(j, "_", 1)) + labs(title = label))
  }
  dev.off()

    