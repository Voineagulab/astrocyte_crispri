## Load
  rm(list = ls())
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
  load("../../Data/Preprocessed/NHA Pooled (Final).rda")
  exp <- nha@assays$RNA@counts
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv")
  load("../../Data/Preprocessed/Guide Expression Matrices.rda")
  # load("Temp Metadata.rda")
  meta <- nha@meta.data
  # meta <- cbind(meta, nha@reductions$umap@cell.embeddings)
  library(tidyverse)
  library(sceptre)
  
  source("../../Scripts/Functions.R")

  
## Setup
  # information
  samples <- c(paste0("NHA_", c(1:5, 7:8)))
  targ.pos <- guides$TargetID[which(guides$TargetCat == "Promoter")] %>% unique()
  targ.enh <- guides$TargetID[which(guides$TargetCat == "Enh")] %>% unique()
  targ.gde <- guides$GuideID[which(guides$TargetCat == "Enh")] %>% unique()
  targ.neg <- guides$GuideID[which(guides$TargetCat == "Negative")] %>% unique()
#   
#   # group and categorise guides
#     groups <- data.frame(gRNA_id = guides$GuideID, 
#                          gRNA_group = guides$TargetID, 
#                          gRNA_type = factor(guides$TargetCat))
# 
#     levels(groups$gRNA_type) <- c("enh_target", "non_target", "tss_target")
# 
#   # enhancer pairs
#     pairs.enh <- list()
#     for (j in targ.enh) {
#       print(j)
#       
#       # get enh coordinates
#       m <- match(j, guides$TargetID)
#       w <- guides$TargetCoord[m]
#       
#       
#       # find nearby tss to the enh
#       nearby <- find.nearby.tss(query = w, expand.by = 10^6, reduced.output = TRUE)
#       nearby <- intersect(nearby, rownames(nha))
#       pairs.enh[[j]] <- data.frame(gene_id = nearby, gRNA_id = j, pair_type = "candidate") # recommended name for pair-type
#     }
#     
#     pairs.enh <- do.call("rbind", pairs.enh)
#   
#     # quick filter to highly expressed genes
#     mean <- rowMeans(nha@assays$RNA@data) # rna assay, not sct
#     mean <- names(mean)[which(mean > 2^-6)]
#     pairs.enh.filt <- pairs.enh[which(pairs.enh$gene_id %in% mean),]
#     
#   # guide assignments
#     frac02 <- guide.frac > 0.02
#     frac02 <- as.data.frame(frac02)
#     frac02 <- frac02[,which(colSums(frac02) > 0)]
#     use02 <- which(colnames(nha) %in% colnames(frac02))
#     exp02 <- exp[,use02] # keep it as a sparse matrix
#     guide02 <- apply(frac02, 2, as.numeric)  # I note that this also transposes teh data, which is necessary
#     rownames(guide02) <- rownames(frac02)
#     guide02 <- as(guide02, "sparseMatrix")
#     guide02 <- guide02[,colnames(exp02)]
#     pool02 <- combine_perturbations(perturbation_matrix = guide02,
#                                                 gRNA_groups_table = groups)
#   
# ## Step 0: Define cells
#   # use cells with any guide
#   
#     
#     
#     binary.guide <- meta[,guides$GuideID[which(guides$TargetCat == "Enh")]] 
#     binary.guide <- apply(binary.guide, 2, as.numeric) %>% t() %>% as("sparseMatrix")
#     colnames(binary.guide) <- rownames(meta)
#     binary.guide.pool <- combine_perturbations(perturbation_matrix = binary.guide,
#                                                gRNA_groups_table = groups)
#     
#     
# 
#     
# 
#  
  
################################################################################################################################ #
## Version 1: The negative control is based on matching guides ----  

  
  ## Use the code in test sceptre 
  
    meta02 <- meta[use02,]
  
## Function
  # g <- rownames(pool02)[1]
  
  libs <- list()
  # for (j in samples) { libs[[j]] <- (meta$Library == j)}
  for (j in samples) { libs[[j]] <- (meta02$Library == j)}
  
  run.cooper <- function(g) {
    time.start <- Sys.time()
    
    # index targeting cells
    # i <- (meta[,g])
    i <- pool02[g,] == 1
    n <- list()
    
    # pseudobulk targeting cells
    pb.target <- list()
    for (j in samples) {
      use <- which(libs[[j]] & i)
      n[[paste0(j, "_Target")]] <- length(use)
      pb.target[[j]] <- rowSums(exp02[,use])
    }
    
    pb.target <- do.call("cbind", pb.target) %>% as.data.frame()
    
    # pseudobulk non-targeting cells
    pb.nontarget <- list()
    for (j in samples) {
      # print(j)
      # targeting cells
      avoid <- which(libs[[j]] & i)
      
      # other guides in targeting cells
      other.guides <- list()
      for (k in avoid) {
        other.guides[[paste0("Cell_", k)]] <- targ.enh[which(as.logical(pool02[targ.enh,k]))]
      }
      
      other.guides <- do.call("c", other.guides) %>% unique()
      
      # cells expressing other guides
      a <- pool02[other.guides, which(libs[[j]])]
      b <- colSums(a) > 0
      b <- names(b)[which(b)]
      
      # cells with the above guides
      use <- b[-which(b %in% rownames(meta)[avoid])] 
      pb.nontarget[[j]] <- rowSums(exp[,use])
      n[[paste0(j, "_NonTarget")]] <- length(use)
    }
    
    pb.nontarget <- do.call("cbind", pb.nontarget) %>% as.data.frame()
    
    # combined pseudobulk matrix
    pb <- cbind(pb.target, pb.nontarget)
    colnames(pb) <- paste0(colnames(pb), c(rep("_Target", 7), rep("_NonTarget", 7)))
      
    
    ## DE
      # setup
      p <- 0.05
      m <- data.frame(Library = c(samples, samples),
                      Group = c(rep("Targeting", 7), rep("Nontargeting", 7)),
                      row.names = colnames(pb))
      formula <- "~ Library + Group"
      
      # run DESeq2
      dds <- DESeqDataSetFromMatrix(countData = pb, colData = m, design = as.formula(formula)) 
      
      dds <- tryCatch(DESeq(dds), error = function(e) e)
      
      if (inherits(dds, "error")) {
        output <- list(Exp = pb, n = do.call("c", n), DE = dds$message)
        time.stop <- Sys.time()
        print(paste0("Error output after ", time.stop - time.start))
        return(output)
      }
      
      dds <- results(dds, contrast = c("Group", "Targeting", "Nontargeting"), alpha = 0.05)
      res <- as.data.frame(dds@listData); rownames(res) <- rownames(dds)

    # return
    output <- list(Exp = pb, n = do.call("c", n), DE = res)
    time.stop <- Sys.time()
    print(paste0("Completed in ", time.stop - time.start))
    return(output)
  }
  
  
## Run
  
  cooper <- list()
  for (g in targ.enh) {
    print(g)
    cooper[[g]] <- run.cooper(g)
  }
  
  save(cooper, file = "Cooper DE Method.rda")
  
## Analyse
  ## How many failed
    in.error <- sapply(cooper, function(x) class(x$DE))
    table(in.error)    
    e <- which(in.error == "character")
    errors <- sapply(cooper[e], function(x) x$DE)
    table(errors) 
    
      # in sum: 48 enhancers have the error "every gene contains at least one zero, cannot compute log geometric means"
  
  ## Bind data and analyse
    # bind
    de <- lapply(cooper, function(x) {
      if (class(x$DE) == "character") {
        x <- data.frame(Mean = "Error",
                        log2fc = "Error",
                        p = "Error")
      } else {
        x <- x$DE
        x <- data.frame(Mean = x$baseMean,
                        log2fc = x$log2FoldChange,
                        p = x$pvalue,
                        Gene = rownames(x))
      }
      
      x <- x[which(x$Gene %in% mean),] # mean was made in the temp sceptre script
      
      return(x)
    })
    
  de <- do.call("rbind", de)
  de$Enh <- splitter(rownames(de), "\\.", 1)
  de$Pair <- paste0(de$Enh, "_", de$Gene)

  old <- read.csv("../2_DE/Archive_V1/Enhancers - All Results Summary Filtered.csv")
  
  d <- de[which(de$Pair %in% old$Pair),]
  d$FDR <- p.adjust(d$p, method = "fdr")
  d$Hit <- d$FDR < 0.1
  m <- match(d$Pair, old$Pair)
  d$OldP <- old$Sceptre.P[m]
  d$OldFDR <- old$Sceptre.FDR[m]
  d$OldHit <- old$Hit[m]
  d$OldFC <- old$logfc.vst[m]
  table(d$Hit, d$OldHit)
  d$Cat <- "Neither"
  d$Cat[which(d$Hit & !(d$OldHit))] <- "Cooper Only"
  d$Cat[which(d$OldHit & !(d$Hit))] <- "SCEPTRE Only"
  d$Cat[which(d$Hit & d$OldHit)] <- "Both"
  
  # Visualisation
    ggplot(d, aes(x = -log10(OldP), y = -log10(p), colour = Cat)) +
      geom_point() +
      theme_bw() +
      labs(x = "SCEPTRE -log10(P)", y = "Cooper -log10(P)") +
      scale_x_continuous(limits = c(0, max(-log10(d$p))))
    
      scale_colour_manual(values = c("black", "red")) 
        
      scale_alpha_manual(values = c(0.2, 1))
    
    pdf(file = "Cooper FC.pdf", height = 3.5, width = 5)
    ggplot(d, aes(x = Cat, y = log2fc)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0)
    dev.off()
    
################################################################################################################################ #
## Version 2: The negative control is derived from PC space ----  

## Determine KNN for each cell    
  # pc02 <- nha@reductions$pca@cell.embeddings[use02,]
  # 
  knn <- list()
  for (j in samples) {
    w <- nha@assays$RNA@counts[,intersect(use02, which(nha$Library == j))] %>% as.matrix()
    # pc <- princomp(w[nha@assays$RNA@var.features,], cor = FALSE)
    pc <- RunPCA(w[nha@assays$RNA@var.features,], verbose = FALSE)

    x <- FindNeighbors(object = pc@cell.embeddings[,1:25], # same as clustering
                       k.param = 50, # calculate many, filtering can occur later
                       return.neighbor = TRUE,
                       verbose = FALSE)
    x <- data.frame(x@nn.idx[,-1],
                    row.names = x@cell.names)

    knn[[j]] <- x
  }
 
    
## Function
  run.cooper.local <- function(g, neighbours = 10) {
    time.start <- Sys.time()
    
    # index targeting cells
    # i <- (meta[,g])
    i <- pool02[g,] == 1
    n <- list()
    
    # pseudobulk targeting cells
    pb.target <- list()
    for (j in samples) {
      use <- which(libs[[j]] & i)
      n[[paste0(j, "_Target")]] <- length(use)
      pb.target[[j]] <- rowSums(exp02[,use])
    }
    
    pb.target <- do.call("cbind", pb.target) %>% as.data.frame()
    
    # pseudobulk non-targeting cells
    pb.nontarget <- list()
    for (j in samples) {
      # print(j)
      # targeting cells
      avoid <- which(libs[[j]] & i)
      
      # neighbours to targeting cells
      use <- knn[[j]][names(avoid),1:neighbours] #### THIS IS EMPTY?
      use <- c(use) %>% do.call("c", .) %>% unique()
      if (any(use %in% avoid)) use <- use[-which(use %in% avoid)]
      
      # pseudobulk the above cells
      pb.nontarget[[j]] <- rowSums(exp[,use])
      n[[paste0(j, "_NonTarget")]] <- length(use)
    }
    
    pb.nontarget <- do.call("cbind", pb.nontarget) %>% as.data.frame()
    
    # combined pseudobulk matrix
    pb <- cbind(pb.target, pb.nontarget)
    colnames(pb) <- paste0(colnames(pb), c(rep("_Target", 7), rep("_NonTarget", 7)))
      
    
    ## DE
      # setup
      p <- 0.05
      m <- data.frame(Library = c(samples, samples),
                      Group = c(rep("Targeting", 7), rep("Nontargeting", 7)),
                      row.names = colnames(pb))
      formula <- "~ Library + Group"
      
      # run DESeq2
      dds <- DESeqDataSetFromMatrix(countData = pb, colData = m, design = as.formula(formula)) 
      
      dds <- tryCatch(DESeq(dds), error = function(e) e)
      
      if (inherits(dds, "error")) {
        output <- list(Exp = pb, n = do.call("c", n), DE = dds$message)
        time.stop <- Sys.time()
        print(paste0("Error output after ", time.stop - time.start))
        return(output)
      }
      
      dds <- results(dds, contrast = c("Group", "Targeting", "Nontargeting"), alpha = 0.05)
      res <- as.data.frame(dds@listData); rownames(res) <- rownames(dds)

    # return
    output <- list(Exp = pb, n = do.call("c", n), DE = res)
    time.stop <- Sys.time()
    print(paste0("Completed in ", time.stop - time.start))
    return(output)
  }
  
## Run
  local.cooper.v2 <- list()
  for (g in targ.enh) {
    print(g)
    local.cooper.v2[[g]] <- run.cooper.local(g)
  }
  
  save(local.cooper.v2, file = "Cooper (Local V2, New PCA) DE Method.rda")
  
  
# ## Process
#   in.error <- sapply(local.cooper.v2, function(x) class(x$DE))
#   table(in.error)    
#   e <- which(in.error == "character")
#   errors <- sapply(cooper[e], function(x) x$DE)
#   table(errors) 
#     
#       # in sum: 48 enhancers have the error "every gene contains at least one zero, cannot compute log geometric means"
  
  ## Bind data and analyse
    # bind
    localpb <- lapply(local.cooper.v2, function(x) {
      if (class(x$DE) == "character") {
        x <- data.frame(Mean = "Error",
                        log2fc = "Error",
                        p = "Error")
      } else {
        x <- x$DE
        x <- data.frame(Mean = x$baseMean,
                        log2fc = x$log2FoldChange,
                        p = x$pvalue,
                        Gene = rownames(x))
      }

      x <- x[which(x$Gene %in% mean),] # mean was made in the temp sceptre script

      return(x)
    })

  localpb <- do.call("rbind", localpb)
  localpb$Enh <- splitter(rownames(localpb), "\\.", 1)
  localpb$Pair <- paste0(localpb$Enh, "_", localpb$Gene)

  # old <- read.csv("../2_DE/Archive_V1/Enhancers - All Results Summary Filtered.csv")
  res <- read.csv("../2_DE/Enh/Results Summary.csv")

  d <- localpb[which(localpb$Pair %in% res$Pair),]
  d$FDR <- p.adjust(d$p, method = "fdr")
  d$Hit <- d$FDR < 0.1
  m <- match(d$Pair, old$Pair)
  d$OldP <- old$Sceptre.P[m]
  d$OldFDR <- old$Sceptre.FDR[m]
  d$OldHit <- old$Hit[m]
  d$OldFC <- old$logfc.vst[m]
  table(d$Hit, d$OldHit)
  d$Cat <- "Neither"
  d$Cat[which(d$Hit & !(d$OldHit))] <- "Cooper Only"
  d$Cat[which(d$OldHit & !(d$Hit))] <- "SCEPTRE Only"
  d$Cat[which(d$Hit & d$OldHit)] <- "Both"
  # 
  # # pdf(file = "Cooper FC.pdf", height = 3.5, width = 5)
  #   ggplot(d, aes(x = Cat, y = log2fc)) +
  #     geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
  #     geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
  #     theme_bw() +
  #     geom_hline(yintercept = 0)
  #   dev.off()
    
    
################################################################################################################################ #
## Version 3: ----  
    
## PC performed separately per libray
## PC does no control for cell size
## DE: for each target cell, compare sign versus neighbours %>% each library is a fraction of hit cells %>% t.test of libraries versus zero
    
# pc02 <- nha@reductions$pca@cell.embeddings[use02,]

knn <- list()
  for (j in samples) {
    print(j)
    w <- nha@assays$RNA@counts[,intersect(use02, which(nha$Library == j))] %>% as.matrix()
    # pc <- princomp(w[nha@assays$RNA@var.features,], cor = FALSE)
    pc <- RunPCA(w[nha@assays$RNA@var.features,], verbose = FALSE)

    x <- FindNeighbors(object = pc@cell.embeddings[,1:25], # same as clustering
                       k.param = 50, # calculate many, filtering can occur later
                       return.neighbor = TRUE,
                       verbose = FALSE)
    x <- data.frame(x@nn.idx[,-1],
                    row.names = x@cell.names)

    knn[[j]] <- x
  }
    
    
## For speed optimisation, convert knn indices to cellnames
  # knn <- lapply(knn, function(z) {
  #   ids <- rownames(z)
  # })
    
    library(matrixStats)
    
    
  # run.cooper.signed.v1 <- function(g, neighbours = 20, percentile = 0.75) {
  #   time.start <- Sys.time()
  #   
  #   # thresholds
  #   # down.thresh <- neighbours * ((100-percentile)/100) # the number of cells it must be less than
  #   # up.thresh <- neighbours * ((percentile)/100) # the number of cells it must be greater than
  #   genes <- sceptre.pairs.enh$gene_id[which(sceptre.pairs.enh$gRNA_id == g)]
  #   e <- as.matrix(exp[genes,use02])
  #   
  #   
  #   # index targeting cells
  #   # i <- (meta[,g])
  #   i <- pool02[g,] == 1
  #   n <- list()
  #   
  #   # for each target cell
  #   pb <- list()
  #   
  #   for (j in samples) {
  #     # print(j)
  #     use <- which(libs[[j]] & i)
  #     pb[[j]] <- list()
  #     n[[j]] <- length(use)
  #     
  #     
  #     baseMx <- matrix(nrow = length(genes), ncol = length(use), dimnames = list(genes, names(use)))
  #     # baseMx[] <- 0
  #     
  #     for (k in names(use)) {
  #       # print(k)
  #       
  #       # expression in hit cell
  #       x.hit <- e[,k]
  #       
  #       # expression in neighbours
  #       other.cells <- knn[[j]][k,1:neighbours] %>% as.numeric()
  #       other.cells <- rownames(knn[[j]])[other.cells]
  #       x.other <- e[,other.cells] 
  #       
  #       qnt <- rowQuantiles(x.other, probs = c(1-percentile, percentile))
  #       
  #       baseMx[x.hit < qnt[,1],k] <- -1
  #       baseMx[x.hit > qnt[,2],k] <- 1
  #       
  #       
  #       
  #     }
  #     
  #     
  #     pb[[j]] <- rowSums(baseMx, na.rm = TRUE) / ncol(baseMx)
  #       
  #   }
  #   
  #   pb <- do.call("cbind", pb)
  #   
  #   res <- apply(pb, 1, function(z) t.test(z, alternative = "two.sided"))
  #   res <- lapply(res, function(x) as.numeric(x[c("p.value", "estimate")] )) 
  #   res <- do.call("rbind", res)
  #   colnames(res) <- c("P", "Mean")
  #   
  #   output <- list(pb = pb,
  #                  n = do.call("c", n),
  #                  res = res)
  # }
  # 
  # cooper.signed <- list()
  # a <- Sys.time()
  # for(j in targ.enh) {
  #   print(j)
  #   cooper.signed[[j]] <- run.cooper.signed(g = j, percentile = 0.75, neighbours = 20)
  # }
  # b <- Sys.time()
  # b - a  
  # 
  # save(cooper.signed, file = "Cooper (Signed) (Incomplete).rda")
  # 
  # s <- lapply(cooper.signed, function(x) as.data.frame(x$res))
  # s <- do.call("rbind", s)
  # s$Gene <- splitter(rownames(s), "\\.", 2)
  # s$Enh <- splitter(rownames(s), "\\.", 1)
  # s$Pair <- paste0(s$Enh, "_", s$Gene)
  # 
  
## Redo, but make every cell a comparison vs all its neig
  run.cooper.signed.v2 <- function(g, neighbours = 10, genes = "auto", use.pool = TRUE) {
    time.start <- Sys.time()
    
    # thresholds
    if (genes == "auto") genes <- c(sceptre.pairs.enh$gene_id[which(sceptre.pairs.enh$gRNA_id == g)], "dCas9-KRAB-T2A-BLAST-WPRE")
    if (genes == "all") genes <- c(res$Gene, sceptre.pairs.pos$gene_id) %>% unique()
    e <- as.matrix(exp[genes,use02])
    
    # index targeting cells
    # i <- (meta[,g])
    if (use.pool) {
      i <- pool02[g,] == 1
    } else {
      i <- guide02[g,] == 1
    }
    # n <- list()
    
    # for each target cell
    pb <- list()
    
    
    for (j in samples) {
      # print(j)
      use <- which(libs[[j]] & i)
      pb[[j]] <- list()

      baseMx <- matrix(nrow = length(genes), ncol = length(use), dimnames = list(genes, names(use)))
      # baseMx[] <- 0
      
      for (k in names(use)) {
        # print(k)
        
        # expression in hit cell
        x.hit <- e[,k]
        
        # expression in neighbours
        other.cells <- knn[[j]][k,1:neighbours] %>% as.numeric()
        other.cells <- rownames(knn[[j]])[other.cells]
        x.other <- e[,other.cells] 

        baseMx[,k] <- apply(x.other, 2, function(x.neighbour) sign(x.hit - x.neighbour)) %>% rowSums()
                
      }
      pb[[j]] <- baseMx / neighbours
    }

    # combine data    
    cellLevel <- do.call("cbind", pb) # retain all cells as separate
    libLevel <- sapply(pb, rowMeans) # average cells within each library
    
    ## DE
      # cell-level (note: ignore it)
      # cov <- sapply(strsplit(colnames(cellLevel), "_"), "[", 2) %>% factor()
      
    
      # library-level
      res.lib <- apply(libLevel, 1, function(z) {
        z <- t.test(z, alternative = "two.sided")
        z <- as.numeric(z[c("p.value", "estimate")] )
      })
      
      res.lib <- t(res.lib) %>% as.data.frame()
    colnames(res.lib) <- c("Lib.P", "Lib.Mean")
    

      # pb <- do.call("cbind", pb)
    
    
    
    output <- list(pb = libLevel,
                   # n = do.call("c", n),
                   res = res.lib)
  }
  
  cooper.signed.v2 <- list()
  
  for(j in targ.enh) { # ~20min
    print(j)
    cooper.signed.v2[[j]] <- run.cooper.signed.v2(g = j, neighbours = 10)
  }

  save(cooper.signed.v2, file = "Cooper (Signed V2).rda")
  
  res <- read.csv("../2_DE/Enh/Results Summary.csv")
  
  cooper.signed.v2 <- lapply(cooper.signed.v2, function(x) cbind(x$res, x$pb))
  cooper.signed.v2 <- do.call("rbind", cooper.signed.v2)
  colnames(cooper.signed.v2)[1:2] <- c("P", "Mean") 
  cooper.signed.v2 <- cooper.signed.v2[-grep("dCas9", rownames(cooper.signed.v2)),]  # this gene was added as a dummy to prevent matrices being coerced to vectors
  cooper.signed.v2$Gene <- splitter(rownames(cooper.signed.v2), "\\.", 2)
  cooper.signed.v2$Enh <- splitter(rownames(cooper.signed.v2), "\\.", 1)
  cooper.signed.v2$Pair <- paste0(cooper.signed.v2$Enh, "_", cooper.signed.v2$Gene)
  cooper.signed.v2 <- cooper.signed.v2[,c(10:12, 1:9)]  
  cooper.signed.v2 <- cooper.signed.v2[which(cooper.signed.v2$Pair %in% res$Pair),]  
  save(cooper.signed.v2, file = "Cooper (Signed V2, Processed).rda")

  # cooper.signed.v2$FDR <- p.adjust(cooper.signed.v2$P, method = "fdr")  
  # cooper.signed.v2$Hit <- cooper.signed.v2$FDR < 0.1  
  
  ## Quick run on neg controls
    cooper.signed.pos <- list()
    for(j in targ.pos) { # ~20min
      print(j)
      a <- Sys.time()
      if (j == "EIF3A") next
      cooper.signed.pos[[j]] <- run.cooper.signed.v2(g = j, neighbours = 10, genes = "all")
      b <- Sys.time()
      b - a
    }
  
  # look at top
    cooper.signed.pos.processed <- lapply(cooper.signed.pos, function(x) {
      x$res$FDR <- p.adjust(x$res$Lib.P, method = "fdr")
      x$res$Hit <- x$res$FDR < 0.1
      return(x$res)
    })
    
    cooper.signed.pos.processed <- do.call("rbind", cooper.signed.pos.processed)
    cooper.signed.pos.processed <- cooper.signed.pos.processed[order(cooper.signed.pos.processed$Lib.P),]
    cooper.signed.pos.processed$Gene <- splitter(rownames(cooper.signed.pos.processed), "\\.", 2)
    cooper.signed.pos.processed$Guide <- splitter(rownames(cooper.signed.pos.processed), "\\.", 1)
    cooper.signed.pos.processed$Paired <- cooper.signed.pos.processed$Gene == cooper.signed.pos.processed$Guide
  
  ## Neg guides
    cooper.signed.neg <- list()
    for(j in targ.neg) { # ~20min
      print(j)
      a <- Sys.time()
      if (j == "Neg_136") next
      cooper.signed.neg[[j]] <- run.cooper.signed.v2(g = j, neighbours = 10, genes = "all", use.pool = FALSE)
      b <- Sys.time()
      b - a
    }
    
    cooper.signed.neg.processed <- lapply(cooper.signed.neg, function(x) {
      x$res$FDR <- p.adjust(x$res$Lib.P, method = "fdr")
      x$res$Hit <- x$res$FDR < 0.1
      return(x$res)
    })
    
    cooper.signed.neg.processed <- do.call("rbind", cooper.signed.neg.processed)
    cooper.signed.neg.processed <- cooper.signed.neg.processed[order(cooper.signed.neg.processed$Lib.P),]
    cooper.signed.neg.processed$Gene <- splitter(rownames(cooper.signed.neg.processed), "\\.", 2)
    cooper.signed.neg.processed$Guide <- splitter(rownames(cooper.signed.neg.processed), "\\.", 1)

    save(cooper.signed.neg, cooper.signed.pos, file = "Cooper (Signed V2, Controls).rda")  
    
  ## Let's look at thresholding on effect size?
    load("Cooper (Signed V2, Processed).rda")
    load("Cooper (Signed V2, Controls).rda")
    
    cse <- cooper.signed.v2
    csn <- cooper.signed.neg.processed
    csp <- cooper.signed.pos.processed    

    # pos
    x <- csp[csp$Paired,]
    table((x$Lib.P < (0.05/125))) # 26 are false!
    
    load("../2_DE/Pos/SCEPTRE Output.rda")
    m <- match(x$Gene, de.pos$gene_id)
    x$Sceptre <- de.pos$p_value[m]
    thresh <- -log10(0.05 / 125)
    pdf(file = "Cooper (Signed) Pos.pdf", height = 4, width = 4)
    qplot(-log10(x$Lib.P), -log10(x$Sceptre)) +
      scale_y_continuous(limits = c(0, 17)) +
      scale_x_continuous(limits = c(0, 17)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_hline(yintercept = thresh, linetype = 2, colour = "red", alpha = 0.5) +
      geom_vline(xintercept = thresh, linetype = 2, colour = "red", alpha = 0.5) +
      labs(x = "Neighbour Signed Perturbation -log10(P)", y = "SCEPTRE -log10(P)")
    dev.off()
    
    # neg
    pdf(file = "Cooper (Signed) Neg.pdf", height = 4, width = 4)
    make_qq_plot(csn$Lib.P[sample(1:nrow(csn), 100000)])
    dev.off()
    
    # enh
    res <- read.csv("../2_DE/Enh/Results Summary.csv")
    cse <-  cse[which(cse$Pair %in% res$Pair),]   
    cse$FDR <- p.adjust(cse$P, method = "fdr")
    cse$Hit <- cse$FDR < 0.1
    m <- match(cse$Pair, res$Pair)
    cse$SceptreHit <- res$Hit[m]
    cse <- cse[order(cse$Sceptre),]

    m <- (cse$P[which(cse$Hit)]) %>% max() %>% -log10(.)
    
    pdf(file = "Cooper (Signed) Volcano.pdf", height = 3, width = 4)
    ggplot(cse, aes(x = Mean, y = -log10(P), colour = SceptreHit)) +
      geom_point() +
      labs() +
      scale_x_continuous(limits = c(-1,1)) +
      theme_bw() +
      geom_hline(yintercept = m, linetype = 2, colour = "red") +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") 
  dev.off()    
  
################################################################################################################################ #
## Version 4: Cooper + DESeq2, neighbours are those with the same MOI ----  
  
run.cooper.v4 <- function(g) {
    time.start <- Sys.time()
    
    # index targeting cells
    # i <- (meta[,g])
    i <- pool02[g,] == 1
    n <- list()
    
    # pseudobulk targeting cells
    pb.target <- list()
    for (j in samples) {
      use <- which(libs[[j]] & i)
      n[[paste0(j, "_Target")]] <- length(use)
      pb.target[[j]] <- rowSums(exp02[,use])
    }
    
    pb.target <- do.call("cbind", pb.target) %>% as.data.frame()
    
    # pseudobulk non-targeting cells
    pb.nontarget <- list()
    for (j in samples) {
      # print(j)
      # targeting cells' moi
      get.mois <- meta02$MOI[which(libs[[j]] & i)]
      
      
      # other guides in targeting cells
      con.cells <- list()
      for (k in get.mois) {
        check <- which(meta02$MOI == k & meta02$Library == j & !(i))
        nSamp <- min(10, length(check))
        con.cells[[paste0("Cell_", k)]] <- sample(rownames(meta02)[check], nSamp, FALSE) 
      }
      
      con.cells <- do.call("c", con.cells) %>% unique()
      
      # cells expressing other guides
      # a <- pool02[other.guides, which(libs[[j]])]
      # b <- colSums(a) > 0
      # b <- names(b)[which(b)]
      
      # cells with the above guides
      # use <- b[-which(b %in% rownames(meta)[avoid])] 
      pb.nontarget[[j]] <- rowSums(exp[,con.cells])
      n[[paste0(j, "_NonTarget")]] <- length(con.cells)
    }
    
    pb.nontarget <- do.call("cbind", pb.nontarget) %>% as.data.frame()
    
    # combined pseudobulk matrix
    pb <- cbind(pb.target, pb.nontarget)
    colnames(pb) <- paste0(colnames(pb), c(rep("_Target", 7), rep("_NonTarget", 7)))
      
    
    ## DE
      # setup
      p <- 0.05
      m <- data.frame(Library = c(samples, samples),
                      Group = c(rep("Targeting", 7), rep("Nontargeting", 7)),
                      row.names = colnames(pb))
      formula <- "~ Library + Group"
      
      # run DESeq2
      dds <- DESeqDataSetFromMatrix(countData = pb, colData = m, design = as.formula(formula)) 
      
      dds <- tryCatch(DESeq(dds), error = function(e) e)
      
      if (inherits(dds, "error")) {
        output <- list(Exp = pb, n = do.call("c", n), DE = dds$message)
        time.stop <- Sys.time()
        print(paste0("Error output after ", time.stop - time.start))
        return(output)
      }
      
      dds <- results(dds, contrast = c("Group", "Targeting", "Nontargeting"), alpha = 0.05)
      res <- as.data.frame(dds@listData); rownames(res) <- rownames(dds)

    # return
    output <- list(Exp = pb, n = do.call("c", n), DE = res)
    time.stop <- Sys.time()
    print(paste0("Completed in ", time.stop - time.start))
    return(output)
  }
  
  
  cooper.v4 <- list()
  for (g in targ.enh) {
    print(g)
    cooper.v4[[g]] <- run.cooper.v4(g)
  }
  
  save(cooper.v4, file = "Cooper (V4).rda")
 
## QC
 in.error <- sapply(cooper.v4, function(x) class(x$DE))
  table(in.error)    

  ## Bind data and analyse
    # bind
    cooper.v4.de <- lapply(cooper.v4, function(x) {
      if (class(x$DE) == "character") {
        x <- data.frame(Mean = "Error",
                        log2fc = "Error",
                        p = "Error")
      } else {
        x <- x$DE
        x <- data.frame(Mean = x$baseMean,
                        log2fc = x$log2FoldChange,
                        p = x$pvalue,
                        Gene = rownames(x))
      }
      
      x <- x[which(x$Gene %in% mean),] # mean was made in the temp sceptre script
      
      return(x)
    })
    
  cooper.v4.de <- do.call("rbind", cooper.v4.de)
  cooper.v4.de$Enh <- splitter(rownames(cooper.v4.de), "\\.", 1)
  cooper.v4.de$Pair <- paste0(cooper.v4.de$Enh, "_", cooper.v4.de$Gene)
  
  # read in old data to filter it
  old <- read.csv("../2_DE/Enh/Results Summary.csv")
  
  d <- cooper.v4.de[which(cooper.v4.de$Pair %in% old$Pair),]
  d$FDR <- p.adjust(d$p, method = "fdr")
  d$Hit <- d$FDR < 0.1
  
  m <- match(d$Pair, old$Pair)
  d$OldP <- old$P[m]
  d$OldFDR <- old$FDR[m]
  d$OldHit <- old$Hit[m]
  d$OldFC <- old$logfc.vst[m]
  table(d$Hit, d$OldHit)
  d$Cat <- "Neither"
  d$Cat[which(d$Hit & !(d$OldHit))] <- "Cooper V4 Only"
  d$Cat[which(d$OldHit & !(d$Hit))] <- "SCEPTRE Only"
  d$Cat[which(d$Hit & d$OldHit)] <- "Both"
  
  # Visualisation
  pdf(file = "Cooper (V4).pdf", height = 5, width = 5)
    ggplot(d, aes(x = -log10(OldP), y = -log10(p), colour = Cat)) +
      geom_point() +
      theme_bw() +
      labs(x = "SCEPTRE -log10(P)", y = "Cooper -log10(P)") +
      scale_x_continuous(limits = c(0, max(-log10(d$p))))
      scale_colour_manual(values = c("black", "red")) 
  
    
    ggplot(d, aes(x = Cat, y = log2fc)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0)
    dev.off()
  
  
################################################################################################################################ #
## Version 5: Signed Cooper, neighbours are those with the same MOI ----  
  
  
cpm02 <- apply(exp[unique(c(res$Gene, sceptre.pairs.pos$gene_id, "dCas9-KRAB-T2A-BLAST-WPRE")),], 2, function(x) {
  libSize <- sum(x) 
  x <- x / (libSize / 10^6)
  return(x)
})
    
run.cooper.signed.v5 <- function(g, neighbours = 10, use.pool = TRUE, moi.range = c(0,1000), genes = "auto", use.cpm = FALSE) {
    time.start <- Sys.time()
    
    # thresholds
    if (genes == "auto") genes <- c(sceptre.pairs.enh$gene_id[which(sceptre.pairs.enh$gRNA_id == g)], "dCas9-KRAB-T2A-BLAST-WPRE")
    if (genes == "all") genes <- c(res$Gene, sceptre.pairs.pos$gene_id) %>% unique()
    # genes <- c(sceptre.pairs.enh$gene_id[which(sceptre.pairs.enh$gRNA_id == g)], "dCas9-KRAB-T2A-BLAST-WPRE")
    
    if (use.cpm) {
      genes <- genes[which(genes %in% rownames(cpm02))]
      e <- as.matrix(cpm02[genes,use02])
    } else {
      e <- as.matrix(exp[genes,use02])
    }
      
    # index targeting cells
    # i <- (meta[,g])
    if (use.pool) {
      i <- pool02[g,] == 1
    } else {
      i <- guide02[g,] == 1
    }
    # n <- list()
    
    # for each target cell
    pb <- list()
    
    
    for (j in samples) {
      # print(j)
      use <- which(libs[[j]] & i)
      pb[[j]] <- list()

      baseMx <- matrix(nrow = length(genes), ncol = length(use), dimnames = list(genes, names(use)))
      # baseMx[] <- 0
      
      for (k in names(use)) {
        # print(k)
        
        # expression in hit cell
        x.hit <- e[,k]
        
        # expression in neighbours
        # other.cells <- knn[[j]][k,1:neighbours] %>% as.numeric()
        # other.cells <- rownames(knn[[j]])[other.cells]
        # x.other <- e[,other.cells] 
        m <- meta02[k,"MOI"] # MOI of target cells
        other.cells <- which(meta02$MOI == m & meta02$Library == j & !(rownames(meta02) %in% (use))) # cells from same lib, same MOI, without guide.
        other.cells <- rownames(meta02)[other.cells]
        
        # check if cell within MOI range
        if (m < moi.range[1] | m > moi.range[2]) {
          baseMx[,k] <- NA # filtered on MOI!
        }
        
        if (length(other.cells) >= neighbours) {
          other.cells <- sample(other.cells, neighbours, FALSE) 
          x.other <- e[,other.cells]
          baseMx[,k] <- apply(x.other, 2, function(x.neighbour) sign(x.hit - x.neighbour)) %>% rowSums()
        } else {
          baseMx[,k] <- NA # if too few cells, ignore!
        }
        

        
                
      }
      pb[[j]] <- baseMx / neighbours
    }

    # combine data    
    cellLevel <- do.call("cbind", pb) # retain all cells as separate
    libLevel <- sapply(pb, rowMeans, na.rm = TRUE) # average cells within each library
    
    ## DE
      # cell-level (note: ignore it)
      # cov <- sapply(strsplit(colnames(cellLevel), "_"), "[", 2) %>% factor()
      
    
      # library-level
      res.lib <- apply(libLevel, 1, function(z) {
        z <- t.test(z, alternative = "two.sided")
        z <- as.numeric(z[c("p.value", "estimate")] )
      })
      
      res.lib <- t(res.lib) %>% as.data.frame()
      colnames(res.lib) <- c("Lib.P", "Lib.Mean")
    

      # pb <- do.call("cbind", pb)
    
    
    
    output <- list(pb = libLevel,
                   # n = do.call("c", n),
                   res = res.lib)
}
  
  cooper.signed.v5 <- list()
  
  for(j in targ.enh) { # ~20min
    print(j)
    cooper.signed.v5[[j]] <- run.cooper.signed.v5(g = j, neighbours = 10)
  }
  
  save(cooper.signed.v5, file = "Cooper (Signed V5.rda)")
  
  y <- lapply(cooper.signed.v5, function(x) cbind(x$res, x$pb))
  y <- do.call("rbind", y)
  
  
## And with MOI < 5
  cooper.signed.v5.lowmoi <- list()
  
  for(j in targ.enh) { # ~20min
    print(j)
    cooper.signed.v5.lowmoi[[j]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(0,5))
  }
  save(cooper.signed.v5.lowmoi, file = "Cooper (Signed V5 LowMOI.rda)")
  
  z <- lapply(cooper.signed.v5.lowmoi, function(x) cbind(x$res, x$pb))
  z <- do.call("rbind", z)
  
  
## And with MOI 10-20 
  cooper.signed.v5.highmoi <- list()
  
  for(j in targ.enh) { # ~20min
    print(j)
    cooper.signed.v5.highmoi[[j]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(10,20))
  }
  save(cooper.signed.v5.highmoi, file = "Cooper (Signed V5 HighMOI.rda)")
  
  w <- lapply(cooper.signed.v5.highmoi, function(x) cbind(x$res, x$pb))
  w <- do.call("rbind", z)
  
  
## Compare three methods (plus SCEPTRE)
  v5.all <- lapply(cooper.signed.v5, function(x) cbind(x$res, x$pb)) %>% do.call("rbind", .)
  v5.low <- lapply(cooper.signed.v5.lowmoi, function(x) cbind(x$res, x$pb)) %>% do.call("rbind", .)
  v5.high <- lapply(cooper.signed.v5.highmoi, function(x) cbind(x$res, x$pb)) %>% do.call("rbind", .)
  
  v5 <- data.frame(v5.all[,1:2], v5.low[,1:2], v5.high[,1:2])
  colnames(v5) <- c("P.All", "Mean.All", "P.Low", "Mean.Low", "P.High", "Mean.High")
  
  rownames(v5) <- sub("\\.", "_", rownames(v5))
  v5$Gene <- splitter(rownames(v5), "_", 2)
  v5$Enh <- splitter(rownames(v5), "_", 1)
  v5$Pair <- paste0(v5$Enh, "_", v5$Gene)
  
  v5 <- v5[which(v5$Pair %in% old$Pair),]
  m <- match(v5$Pair, old$Pair)
  v5$SCEPTRE <- old$P[m]  

  v5$FDR.all <- p.adjust(v5$P.All, method = "fdr")  
  v5$FDR.low <- p.adjust(v5$P.Low, method = "fdr")  
  v5$FDR.high <- p.adjust(v5$P.High, method = "fdr")  
  v5$FDR.SCEPTRE <- p.adjust(v5$SCEPTRE, method = "fdr")  
  
  # Categorise
  p <- melt(v5[,c(1,3,5,10),], id.vars = "SCEPTRE")
  p$SCEPTRE <- -log10(p$SCEPTRE)
  p$value <- -log10(p$value)
  levels(p$variable) <- c("All Cells", "MOI 1-5 Only", "MOI 10-20 Only")
  
  pdf(file = "Cooper (Signed V5) vs SCEPTRE.pdf", height = 3, width = 7)
  ggplot(p, aes(x = SCEPTRE, y = value)) +
    geom_point() +
    facet_wrap(~variable) +
    theme_bw() +
    scale_y_continuous(limits = c(0,16)) +
    labs(x = "SCEPTRE -log10P", y = "Signed MOI Control Test -log10P")
dev.off()  


## Quick check of negative controls...
  cooper.signed.v5.neg <- list()
  for(j in targ.neg) { # ~20min
    print(j)
    a <- Sys.time()
    if (j == "Neg_136") next
    cooper.signed.v5.neg[[j]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(0,1000), genes = "all", use.pool = FALSE)
    # cooper.signed.v5.neg[[paste0("Low_", j)]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(0,5), genes = "all", use.pool = FALSE)
    # cooper.signed.v5.neg[[paste0("High_", j)]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(10,20), genes = "all", use.pool = FALSE)
    b <- Sys.time()
    b - a
  }
  
  cooper.signed.v5.neg.cpm <- list()
  for(j in targ.neg) { # ~20min
    print(j)
    a <- Sys.time()
    if (j == "Neg_136") next
    cooper.signed.v5.neg.cpm[[j]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(0,1000), genes = "all", use.pool = FALSE, use.cpm = TRUE)
    # cooper.signed.v5.neg[[paste0("Low_", j)]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(0,5), genes = "all", use.pool = FALSE)
    # cooper.signed.v5.neg[[paste0("High_", j)]] <- run.cooper.signed.v5(g = j, neighbours = 10, moi.range = c(10,20), genes = "all", use.pool = FALSE)
    b <- Sys.time()
    print(b - a)
  
  
  
  save(cooper.signed.v5.neg, cooper.signed.v5.neg.cpm, file = "Cooper (Signed V5 Neg).rda")
  
  
  x <- lapply(cooper.signed.v5.neg.cpm, function(x) cbind(x$res, x$pb))
  x <- do.call("rbind", x)
  y <- lapply(cooper.signed.v5.neg, function(x) cbind(x$res, x$pb))
  y <- do.call("rbind", y)
  
  p <- data.frame(Counts = y$Lib.Mean, CPM = x$Lib.Mean, Run = rownames(x))
  p$Gene <- splitter(p$Run, "\\.", 2)
  
  # mean.exp <- rowMeans(cpm02)
  m <- match(p$Gene, names(mean.exp))
  p$Exp <- mean.exp[m]
  p$ExpBin <- cut(p$Exp, c(0,10, 25, 50, 100, 500, 1000000))
  levels(p$ExpBin) <- c("<10cpm", "10-25cpm", "25-50cpm", "50-100cpm", "100-500cpm", ">500cpm", "NA")
  
  # y$MOI.filter <- splitter(rownames(y), "_", 1)
  
  # pdf(file = "Cooper (Signed V5) Negative Control Guides.pdf", height = 3, width = 5)
  # ggplot(y, aes(x = ".", y = Lib.Mean)) +
  #   geom_boxplot() +
  #   geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  #   theme_bw() +
  #   theme(panel.border = element_blank(), axis.line.y = element_line()) +
  #   scale_y_continuous(limits = c(-1, 1))
  # dev.off()  
  # 

  

## Scatterplot cpm vs count on p and libmean
  # q <- melt(p[,c(1,2,4)])
  p$Density <- get_density(p$Counts, p$CPM)
  pdf(file = "Cooper (Signed V5) (Negative) (CPM vs Count.pdf", height = 4, width = 5)
  ggplot(p[sample(1:nrow(p), 100000),], aes(x = Counts, y = CPM, colour = Density)) +
    geom_point() +
    theme_bw() +
    theme() +
    scale_colour_viridis_c() +
    labs(x = "Lib Mean Using Counts", y = "Lib Mean Using CPM")
  dev.off()
    
  
## Stratify lib mean by expression
  
  pdf(file = "Cooper (Signed V5) (Negative) (Lib Mean Across Expression Bins).pdf", height = 4, width = 5)
  q <- melt(p[,c(1,2,6)])
  ggplot(q, aes(x = variable, fill = ExpBin, y = value)) +
    # geom_boxplot() +
    geom_violin() +
    theme_bw() +
    labs(x = "Expression Normalisation As Input", y = "Library Mean")
  
  q <- melt(p[,c(1,2,5)], id.vars = "Exp")
  ggplot(q[sample(1:nrow(q), 100000),], aes(x = Exp, y = value)) +
    geom_point() +
    facet_wrap(~variable) +
    # scale_colour_manual(values = c("grey1", "grey2")) +
    geom_smooth(colour = "red", linetype = 2) +
    theme_bw() +
    geom_hline(yintercept = 0, colour = "grey50") +
    scale_x_continuous(trans = "log2") +
    labs(x = "Log2 CPM", y = "Library Mean")
  dev.off()
## Scatterplot mean exp vs lib mean

  
## Rerun enh on cpm
  cooper.signed.v5.cpm <- list()
  sceptre.pairs.enh <- pairs.enh.filt
  
  for(j in targ.enh[-which(targ.enh %in% names(cooper.signed.v5.cpm))]) { # ~20min
    print(j)
    if (j %in% c("Enh326", "Enh527", "Enh528", "Enh689", "Enh690", "Enh839", "Enh891", "Enh921", "Enh953", "Enh679", "Enh713", "Enh784", "Enh840")) next
    cooper.signed.v5.cpm[[j]] <- run.cooper.signed.v5(g = j, neighbours = 10, use.cpm = TRUE, genes = "auto")
  }
  
  save(cooper.signed.v5.cpm, file = "Cooper (Signed V5 CPM.rda)")
  
  y <- lapply(cooper.signed.v5.cpm, function(x) cbind(x$res, x$pb))
  y <- do.call("rbind", y)
    
  y$Pair <- sub("\\.", "_", rownames(y))
  y <- y[which(y$Pair %in% res$Pair),]
  
  y$Lib.PAdj <- p.adjust(y$Lib.P, method = "fdr")
  table(y$Lib.PAdj < 0.1)
  z <- y$Pair[which(y$Lib.PAdj < 0.1)]
  
  table(res.final$HitCategory[which(res$Pair %in% z)])
  
  m <- match(y$Pair, res.final$Pair)
  z <- data.frame(Cooper = -log10(y$Lib.P), SCEPTRE = -log10(res.final$P.SCEPTRE[m]), N50 = -log10(res.final$P.N50[m]), Hit = res.final$HitCategory[m])
  z <- melt(z, id.vars = c("Cooper", "Hit"))
  
  pdf(file = "Cooper (Signed V5 CPM) - Versus Other Algs.pdf", height = 3.5, width = 6)
  ggplot(z, aes(x = Cooper, y = value, colour = Hit)) +
    geom_point() +
    facet_wrap(~variable) +
    theme_bw() +
    scale_colour_manual(values = sig.colours2) +
    geom_vline(xintercept = -log10(3.915283e-04), linetype = 2) +
    labs(x = "Cooper -log10(P)", y = "Other Alg -log10(P)")
  dev.off()
  
  y$Hit <- y$Lib.PAdj < 0.1
  pdf(file = "Cooper (Signed V5 CPM) - Enh Effect Size Boxplots.pdf", height = 4, width = 3)
  ggplot(y, aes(x = Hit, y = Lib.Mean)) +
    geom_jitter(width = 0.1) +
    geom_violin(scale = "width", alpha = 0.2, draw_quantiles = c(0.5)) +
    theme_bw() +
    geom_hline(yintercept = 0, colour = "red") +
    scale_y_continuous(limits = c(-1, 1))
  dev.off()
    