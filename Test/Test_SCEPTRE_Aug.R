################################################################################################################################ #
## Setup ----

## Load
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/1_Processing/")
  load("../../Data/Preprocessed/NHA Pooled (Final).rda")
  exp <- nha@assays$RNA@counts
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv")
  load("../../Data/Preprocessed/Guide Expression Matrices.rda")
  # load("Temp Metadata.rda")
  meta <- nha@meta.data
  # meta <- cbind(meta, nha@reductions$umap@cell.embeddings)
  library(tidyverse)
  library(sceptre)

# ## Fraction matrix
#   umi.thresh <- 3
#   guide.frac <- apply(guide.exp.raw, 2, function(x) {
#       x[which(x < umi.thresh)] <- 0
#       x <- x / sum(x)
#     })
#   guide.frac <- as.data.frame(guide.frac)

## Guide groupings
  targ.pos <- guides$TargetID[which(guides$TargetCat == "Promoter")] %>% unique()
  targ.enh <- guides$TargetID[which(guides$TargetCat == "Enh")] %>% unique()
  targ.neg <- guides$GuideID[which(guides$TargetCat == "Negative")] %>% unique()
  
  ## Step 4: Group and categorise guides
    groups <- data.frame(gRNA_id = guides$GuideID, 
                         gRNA_group = guides$TargetID, 
                         gRNA_type = factor(guides$TargetCat))
    
    # for negative control guides, create "groups" for a non-targeting and AAVS1-targeting guide
    for (j in 1:125) {
      m <- match(paste0("Neg_", j), groups$gRNA_id)
      groups$gRNA_group[c(m, m+125)] <- paste0("NegPair_", j)
      
    }
    
    for (j in 1:25) {
      m <- match(paste0("Neg_E_", j), groups$gRNA_id)
      groups$gRNA_group[c(m, m+25)] <- paste0("NegPair_E_", j)
      
    }
    
    levels(groups$gRNA_type) <- c("enh_target", "non_target", "tss_target")
  
  
    
  ## Enhancers
    pairs.enh <- list()
    for (j in targ.enh) {
      print(j)
      
      # get enh coordinates
      m <- match(j, guides$TargetID)
      w <- guides$TargetCoord[m]
      
      
      # find nearby tss to the enh
      nearby <- find.nearby.tss(query = w, expand.by = 10^6, reduced.output = TRUE)
      nearby <- intersect(nearby, rownames(nha))
      pairs.enh[[j]] <- data.frame(gene_id = nearby, gRNA_id = j, pair_type = "candidate") # recommended name for pair-type
    }
    
    pairs.enh <- do.call("rbind", pairs.enh)
  
    # quick filter to highly expressed genes
    mean <- rowMeans(nha@assays$RNA@data) # rna assay, not sct
    mean <- names(mean)[which(mean > 2^-6)]
    pairs.enh.filt <- pairs.enh[which(pairs.enh$gene_id %in% mean),]
    
## Pos
sceptre.pairs.pos <- data.frame(gene_id = targ.pos,
                                    gRNA_id = targ.pos,
                                    pair_type = "positive_control")
  
    # two genes are not present in the expression data; this is because a synonym is used...  
    sceptre.pairs.pos$gene_id[sceptre.pairs.pos$gRNA_id == "FAM96B"] <- "CIAO2B" # synonym
    sceptre.pairs.pos$gene_id[sceptre.pairs.pos$gRNA_id == "SEPT11"] <- "SEPTIN11" # synonym
    
## Direction to test
    direction <- "both"
  
## Clean
    rm(nha)
    
################################################################################################################################ #
## Test with a guide fraction of 0.02 ----


frac02 <- guide.frac > 0.02
frac02 <- as.data.frame(frac02)
frac02 <- frac02[,which(colSums(frac02) > 0)]


## Step 0: Define cells
  # use cells with any guide
  use02 <- which(colnames(nha) %in% colnames(frac02))

## Step 1: Get expression matrix
  exp02 <- exp[,use02] # keep it as a sparse matrix
    
## Step 2: gRNA binary expression matrix
  # guide <- as.data.frame(nha@meta.data[use02,guide.cols])
  guide02 <- apply(frac02, 2, as.numeric)  # I note that this also transposes teh data, which is necessary
  rownames(guide02) <- rownames(frac02)
  guide02 <- as(guide02, "sparseMatrix")
  guide02 <- guide02[,colnames(exp02)]
    

      
    
## Step 3: now also pool perturbations!
  pool02 <- combine_perturbations(perturbation_matrix = guide02,
                                                gRNA_groups_table = groups)
  
  
## Step 3: Covariate matrix
  covar02 <- data.frame(UMI_count = log(colSums(exp02)),
                            gRNA_count = log(nha$MOI)[use02],
                            # gRNA_count = log(colSums(pool02)),
                            Mito_P = meta$Mito_Pct[use02],
                            Batch = meta$Library[use02],
                            C_Cycle = meta$Cycle_Seurat_S[use02])


    
## Run loop (requires 1-2min per chunk on our server)
  chunkSize <- 100
  nChunks <- ceiling(nrow(pairs.enh.filt) / chunkSize)
  
  res02 <- list()
  
  for (j in 1:nChunks) {
    gc()
    
    # setup
    a <- Sys.time()
    run.name <- paste0("Chunk_", j)
    
    # get the current chunk of a given chunkSize
    chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
    chunk <- chunk[which(chunk < nrow(pairs.enh.filt))] # truncate if the range is outside the full dataframe
    chunk <- pairs.enh.filt[chunk,]
    
    # run sceptre
    
    res02[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp02[chunk$gene_id,],
                                                         combined_perturbation_matrix = pool02,
                                                         covariate_matrix = covar02,
                                                         gene_gRNA_group_pairs = chunk,
                                                         side = direction,
                                                         B = 1000)
    
    # sign off
    b <- Sys.time()
    print(paste0(run.name, ": ", b-a))
    save(res02, file = "Temp SCEPTRE.rda")
    
  }

## Save  
  res02 <- do.call("rbind", res02)
  save(res02, file = "../2_DE/Temp3 (Dom2 SCEPTRE).rda")
  
## Positive
  res02.pos <- run_sceptre_high_moi(gene_matrix = exp02[sceptre.pairs.pos$gene_id,],
                                                         combined_perturbation_matrix = pool02,
                                                         covariate_matrix = covar02,
                                                         gene_gRNA_group_pairs = sceptre.pairs.pos,
                                                         side = direction,
                                                         B = 1000)
  
  
## Run with alternative covariates!
  covar02.alt <- data.frame(UMI_count = log(colSums(exp02)),
                            # gRNA_count = log(nha$MOI)[use02],
                            gRNA_count = log(colSums(pool02)),
                            Mito_P = meta$Mito_Pct[use02],
                            Batch = meta$Library[use02],
                            Tricycle_1 = meta$Cycle_Tricycle_PC1[use02],
                            Tricycle_2 = meta$Cycle_Tricycle_PC2[use02],
                            SPhase = meta$Cycle_Seurat_S[use02],
                            NuclearFraction = meta$NuclearFraction[use02],
                            Ribo_P = meta$Ribo_Pct[use02])
  
  chunkSize <- 100
  nChunks <- ceiling(nrow(pairs.enh.filt) / chunkSize)
  
  covar.test <- list()
  covar.test$PlusRiboTriNF <- covar.test$PlusRiboTri <- covar.test$PlusNF <- covar.test$PlusRibo <- covar.test$PlusTri <- covar.test$Base <- list()
  covar.use <- covar.test
  covar.use$Base <- c("UMI_count", "gRNA_count", "Mito_P", "Batch", "SPhase")
  covar.use$PlusTri <- c("UMI_count", "gRNA_count", "Mito_P", "Batch", "Tricycle_1")
  covar.use$PlusRibo <- c("UMI_count", "gRNA_count", "Mito_P", "Batch", "SPhase", "Ribo_P")
  covar.use$PlusNF <- c("UMI_count", "gRNA_count", "Mito_P", "Batch", "SPhase", "NuclearFraction")
  covar.use$PlusRiboTri <- c("UMI_count", "gRNA_count", "Mito_P", "Batch", "Tricycle_1",  "Ribo_P")
  covar.use$PlusRiboTriNF <- c("UMI_count", "gRNA_count", "Mito_P", "Batch", "Tricycle_1",  "Ribo_P", "NuclearFraction")
  
  for (j in 1:nChunks) {
    gc()
    
    # setup
    a <- Sys.time()
    run.name <- paste0("Chunk_", j)
    
    # get the current chunk of a given chunkSize
    chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
    chunk <- chunk[which(chunk < nrow(pairs.enh.filt))] # truncate if the range is outside the full dataframe
    chunk <- pairs.enh.filt[chunk,]
    
    # run sceptre
    # covar.test$Base[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp02[chunk$gene_id,],
    #                                                      combined_perturbation_matrix = pool02,
    #                                                      covariate_matrix = covar02.alt[,covar.use$Base],
    #                                                      gene_gRNA_group_pairs = chunk,
    #                                                      side = direction,
    #                                                      B = 1000)
    # gc()
    # covar.test$PlusTri[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp02[chunk$gene_id,],
    #                                                      combined_perturbation_matrix = pool02,
    #                                                      covariate_matrix = covar02.alt[,covar.use$PlusTri],
    #                                                      gene_gRNA_group_pairs = chunk,
    #                                                      side = direction,
    #                                                      B = 1000)
    # gc()
    # covar.test$PlusRibo[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp02[chunk$gene_id,],
    #                                                      combined_perturbation_matrix = pool02,
    #                                                      covariate_matrix = covar02.alt[,covar.use$PlusRibo],
    #                                                      gene_gRNA_group_pairs = chunk,
    #                                                      side = direction,
    #                                                      B = 1000)
    # gc()
    covar.test$PlusNF[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp02[chunk$gene_id,],
                                                         combined_perturbation_matrix = pool02,
                                                         covariate_matrix = covar02.alt[,covar.use$PlusNF],
                                                         gene_gRNA_group_pairs = chunk,
                                                         side = direction,
                                                         B = 1000)
    gc()
    covar.test$PlusRiboTri[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp02[chunk$gene_id,],
                                                         combined_perturbation_matrix = pool02,
                                                         covariate_matrix = covar02.alt[,covar.use$PlusRiboTri],
                                                         gene_gRNA_group_pairs = chunk,
                                                         side = direction,
                                                         B = 1000)
    gc()
    covar.test$PlusRiboTriNF[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp02[chunk$gene_id,],
                                                         combined_perturbation_matrix = pool02,
                                                         covariate_matrix = covar02.alt[,covar.use$PlusRiboTriNF],
                                                         gene_gRNA_group_pairs = chunk,
                                                         side = direction,
                                                         B = 1000)
    gc()

    # sign off
    b <- Sys.time()
    print(paste0(run.name, ": ", b-a))
    save(covar.test, file = "Temp SCEPTRE CovarTest.rda")
    
  }

  # x <- res02.covar
  # res02.covar <- do.call("rbind", res02.covar)
  
  covar.test <- lapply(covar.test, function(x) do.call("rbind", x))
  save(covar.test, file = "../2_DE/Temp6 (covariate eval).rda")
  
  y <- do.call("cbind", covar.test)
  
  y <- y[,c(1:3, seq(4, 30, 5), seq(5, 30, 5))]
  y$Base.pair_type <- paste0(y$Base.gRNA_id, "_", y$Base.gene_id)
  colnames(y)[1:3] <- c("Gene", "Enh", "Pair") 
  colnames(y) <- gsub("\\.p_value", "_P", colnames(y)) %>% gsub("\\.z_value", "_z", .) %>% gsub("Plus", "", .)
  
  # filter
  y <- y[which(y$Pair %in% old$Pair),] # where "old" is loaded some way below, as a version filtered for distance
  m <- match(y$Pair, old$Pair)
  y$Distance <- old$Distance[m]
  y$OldHit <- old$Hit[m]
   
  # fdr
  runs <- gsub("Plus", "", names(covar.use))
  
  for (j in runs) {
    y[,paste0(j, "_FDR")] <- p.adjust(y[,paste0(j, "_P")], method = "fdr")
    y[,paste0(j, "_Hit")] <- y[,paste0(j, "_FDR")] < 0.1
  }
  
  # overlap runs
  overlaps <- apply(y[,paste0(runs, "_Hit")], 2, function(x) {
    apply(y[,paste0(runs, "_Hit")], 2, function(z) { length(which(x & z)) })
  })
  
  overlaps <- as.data.frame(overlaps)
  overlaps$Model <-rownames(overlaps) <- colnames(overlaps) <- gsub("_Hit", "", rownames(overlaps)) 
  
  overlaps <- melt(overlaps)
  
  levels(overlaps$variable)
  overlaps$Model <- factor(overlaps$Model, levels = rev(levels(overlaps$variable)))
  
  pdf(file = "../2_DE/Temp6 - Heatmap.pdf", height = 5, width = 5)
  ggplot(overlaps, aes(x = variable, y = Model, label = value)) +
    geom_tile(fill = "white", colour = "grey50") +
    geom_text() +
    theme_bw() +
    theme(panel.border = invis, panel.grid = invis, axis.title = invis)
  dev.off()
  
  # compare the gain and loss for the base to the most complex model
  gain <- y[which(y$RiboTriNF_Hit & !(y$Base_Hit)),] # 5, 4 down, distance 110
  common <- y[which(y$RiboTriNF_Hit & (y$Base_Hit)),] 
  lost <- y[which(!(y$RiboTriNF_Hit) & (y$Base_Hit)),] # 4, 3 down, distance 303
  
  # compare gain and loss from the original analyses (with base model and no dom thresh and no cell filtering)
  gain <- y[which(y$RiboTriNF_Hit & !(y$OldHit)),] # 14, 8 down, distance 150
  common <- y[which(y$RiboTriNF_Hit & (y$OldHit)),] # 81, 75 down, distance 27
  lost <- y[which(!(y$RiboTriNF_Hit) & (y$OldHit)),] # 23, 16 down, distance 120
    # interesting losses: atxn10; syt14
  
## Plots for pondering... Are these sensible to regress?
  ## Versus PCA
    # collect pca informatin
    pca <- nha@reductions$pca@cell.embeddings[,1:20] # top 20 principal components
    pctvar <- (nha@reductions$pca@stdev)^2 
    pctvar <- pctvar / sum(pctvar)
    pctvar <- pctvar[1:20]
    pctvar <- round(pctvar * 100, 1)
    names(pctvar) <- colnames(pca)
    
    # correlate
    cor <- cor(covar02.alt[,-4], pca[use02,]) %>% as.data.frame()
    cor$Covariate <- rownames(cor)
    cor <- melt(cor)
    levels(cor$variable) <- paste0(levels(cor$variable), " (", pctvar,"%)") %>% gsub("PC_", "", .)
    cor$Bin <- cut(abs(cor$value), c(1, 0.5, 0.2, 0.1, 0.05, 0))
  
    # plot
    pdf(file = "Final/Covariates - Versus PCA.pdf", height = 5.5, width = 9.5)
    ggplot(cor, aes(x = variable, y = Covariate, fill = Bin, label = round(value*100,1))) +
      geom_tile() +
      geom_text(size = 3) +
      scale_fill_carto_d(palette = "Geyser") +
      theme_bw() +
      labs(x = "Principal Component") +
      guides(fill = guide_legend(title = "Abs. Pearson")) +
      theme(panel.border = invis, panel.grid = invis, axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    dev.off()
      
 
  
################################################################################################################################ #
## Test with a guide fraction of 0.05 ----


frac05 <- guide.frac > 0.05
frac05 <- as.data.frame(frac05)
frac05 <- frac05[,which(colSums(frac05) > 0)]


## Step 0: Define cells
  # use cells with any guide
  use05 <- which(colnames(exp) %in% colnames(frac05))

## Step 1: Get expression matrix
  exp05 <- exp[,use05] # keep it as a sparse matrix
    
## Step 2: gRNA binary expression matrix
  # guide <- as.data.frame(nha@meta.data[use05,guide.cols])
  guide05 <- apply(frac05, 2, as.numeric)  # I note that this also transposes teh data, which is necessary
  rownames(guide05) <- rownames(frac05)
  guide05 <- as(guide05, "sparseMatrix")
  guide05 <- guide05[,colnames(exp05)]
    
    
## Step 3: now also pool perturbations!
  pool05 <- combine_perturbations(perturbation_matrix = guide05,
                                                gRNA_groups_table = groups)
## Step 4: Covariate matrix
  covar05 <- data.frame(UMI_count = log(colSums(exp05)),
                              # gRNA_count = log(meta$MOI)[use05],
                              gRNA_count = log(colSums(pool05)),
                              Mito_P = meta$Mito_Pct[use05],
                              Batch = meta$Library[use05],
                              C_Cycle = meta$Cycle_Seurat_S[use05])
    


    
## Run loop (requires 1-2min per chunk on our server)
  chunkSize <- 100
  nChunks <- ceiling(nrow(pairs.enh.filt) / chunkSize)
  
  res05 <- list()
  
  for (j in 1:nChunks) {
    gc()
    
    # setup
    a <- Sys.time()
    run.name <- paste0("Chunk_", j)
    
    # get the current chunk of a given chunkSize
    chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
    chunk <- chunk[which(chunk < nrow(pairs.enh.filt))] # truncate if the range is outside the full dataframe
    chunk <- pairs.enh.filt[chunk,]
    
    # run sceptre
    
    res05[[run.name]] <- run_sceptre_high_moi(gene_matrix = exp05[chunk$gene_id,],
                                                         combined_perturbation_matrix = pool05,
                                                         covariate_matrix = covar05,
                                                         gene_gRNA_group_pairs = chunk,
                                                         side = direction,
                                                         B = 1000)
    
    # sign off
    b <- Sys.time()
    print(paste0(run.name, ": ", b-a))
    save(res05, file = "Temp SCEPTRE.rda")
    
  }

## Save  
  res05 <- do.call("rbind", res05)
  save(res05, file = "../2_DE/Temp4 (Dom5 SCEPTRE).rda")
  
 res05.pos <- run_sceptre_high_moi(gene_matrix = exp05[sceptre.pairs.pos$gene_id,],
                                                         combined_perturbation_matrix = pool05,
                                                         covariate_matrix = covar05,
                                                         gene_gRNA_group_pairs = sceptre.pairs.pos,
                                                         side = direction,
                                                         B = 1000)
 

################################################################################################################################ #
## Test with a guide fraction > 0 ----
 
 frac00 <- guide.frac > 0.00
frac00 <- as.data.frame(frac00)
frac00 <- frac00[,which(colSums(frac00) > 0)]


## Step 0: Define cells
  # use cells with any guide
  use00 <- which(colnames(exp) %in% colnames(frac00))

## Step 1: Get expression matrix
  exp00 <- exp[,use00] # keep it as a sparse matrix
    
## Step 2: gRNA binary expression matrix
  # guide <- as.data.frame(nha@meta.data[use00,guide.cols])
  guide00 <- apply(frac00, 2, as.numeric)  # I note that this also transposes teh data, which is necessary
  rownames(guide00) <- rownames(frac00)
  guide00 <- as(guide00, "sparseMatrix")
  guide00 <- guide00[,colnames(exp00)]
    
    
## Step 3: now also pool perturbations!
  pool00 <- combine_perturbations(perturbation_matrix = guide00,
                                                gRNA_groups_table = groups)
## Step 4: Covariate matrix
  covar00 <- data.frame(UMI_count = log(colSums(exp00)),
                              # gRNA_count = log(meta$MOI)[use00],
                              gRNA_count = log(colSums(pool00)),
                              Mito_P = meta$Mito_Pct[use00],
                              Batch = meta$Library[use00],
                              C_Cycle = meta$Cycle_Seurat_S[use00])
    


  
 res00.pos <- run_sceptre_high_moi(gene_matrix = exp00[sceptre.pairs.pos$gene_id,],
                                                         combined_perturbation_matrix = pool00,
                                                         covariate_matrix = covar00,
                                                         gene_gRNA_group_pairs = sceptre.pairs.pos,
                                                         side = direction,
                                                         B = 1000)
  
 
################################################################################################################################ #
## Test with cell filtering for high UMI and damaged cells ----


## Step 0: Define cells
  # cells with high UMI
  w <- list()
  for (j in unique(meta$Library)) {
    k <- meta[which(meta$Library == j),]
    m <- quantile(k$UMI_Total[k$MOI == 1], 0.99)
    w[[j]] <- rownames(k)[-which(k$UMI_Total > (m))]
  }
  w <- do.call("c", w)
  
  # damaged cells
  w <- intersect(w, rownames(meta)[which(meta$DropletQC != "damaged_cell")])
  
  ## Define guides
    # use cells with any guide
    fracCellRm <- guide.exp >= 3
    fracCellRm <- fracCellRm[,w]
    fracCellRm <- as.data.frame(fracCellRm)
    fracCellRm <- fracCellRm[,which(colSums(fracCellRm) > 0)]

    useCellRm <- which(colnames(exp) %in% colnames(fracCellRm))
    
## Step 1: Get expression matrix
  expCellRm <- exp[,useCellRm] # keep it as a sparse matrix
    
## Step 2: gRNA binary expression matrix
  # guide <- as.data.frame(nha@meta.data[useCellRm,guide.cols])
  guideCellRm <- apply(fracCellRm, 2, as.numeric)  # I note that this also transposes teh data, which is necessary
  rownames(guideCellRm) <- rownames(fracCellRm)
  guideCellRm <- as(guideCellRm, "sparseMatrix")
  guideCellRm <- guideCellRm[,colnames(expCellRm)]
    
    
## Step 3: now also pool perturbations!
  poolCellRm <- combine_perturbations(perturbation_matrix = guideCellRm,
                                                gRNA_groups_table = groups)
## Step 4: Covariate matrix
  covarCellRm <- data.frame(UMI_count = log(colSums(expCellRm)),
                              # gRNA_count = log(meta$MOI)[useCellRm],
                              gRNA_count = log(colSums(poolCellRm)),
                              Mito_P = meta$Mito_Pct[useCellRm],
                              Batch = meta$Library[useCellRm],
                              C_Cycle = meta$Cycle_Seurat_S[useCellRm])
    


    
## Run loop (requires 1-2min per chunk on our server)
  chunkSize <- 100
  nChunks <- ceiling(nrow(pairs.enh.filt) / chunkSize)
  
  resCellRm <- list()
  
  for (j in 1:nChunks) {
    gc()
    
    # setup
    a <- Sys.time()
    run.name <- paste0("Chunk_", j)
    
    # get the current chunk of a given chunkSize
    chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
    chunk <- chunk[which(chunk < nrow(pairs.enh.filt))] # truncate if the range is outside the full dataframe
    chunk <- pairs.enh.filt[chunk,]
    
    # run sceptre
    
    resCellRm[[run.name]] <- run_sceptre_high_moi(gene_matrix = expCellRm[chunk$gene_id,],
                                                         combined_perturbation_matrix = poolCellRm,
                                                         covariate_matrix = covarCellRm,
                                                         gene_gRNA_group_pairs = chunk,
                                                         side = direction,
                                                         B = 1000)
    
    # sign off
    b <- Sys.time()
    print(paste0(run.name, ": ", b-a))
    save(resCellRm, file = "Temp SCEPTRE.rda")
    
  }

## Save  
  resCellRm <- do.call("rbind", resCellRm)
  save(resCellRm, file = "../2_DE/Temp5 CellRM SCEPTRE).rda")
  
  
################################################################################################################################ #
## Evaluate ----
  
  
## Thresholds
  exp.thresh <- 2^-6
  p.thresh <- 0.1

  

  
## Load
  old <- read.csv("../2_DE/Archive_V1/Enhancers - All Results Summary Filtered.csv")
  
# s <- data.frame(Gene = res02$gene_id,
#                 Enh = res02$gRNA_id,
#                 Pair = paste0(res02$gRNA_id, "_", res02$gene_id),
#                 Dom2_P = ".",
#                 Dom2_Z = ".",
#                 Dom5_P = res05$p_value,
#                 Dom5_Z = res05$z_value,
#                 CellRm_P = resCellRm$p_value,
#                 CellRm_Z = resCellRm$z_value)
  
  s <- data.frame(Gene = res02$gene_id,
                Enh = res02$gRNA_id,
                Pair = paste0(res02$gRNA_id, "_", res02$gene_id))
  
  m <- match(s$Pair, paste0(res05$gRNA_id, "_", res05$gene_id))
  s <- cbind(s,  
             Dom2_P = res02$p_value,
                Dom2_Z = res02$z_value,
                Dom5_P = res05$p_value[m],
                Dom5_Z = res05$z_value[m],
                CellRm_P = resCellRm$p_value[m],
                CellRm_Z = resCellRm$z_value[m])
  

s <- s[which(s$Pair %in% old$Pair),]
m <- match(s$Pair, old$Pair)
s$Base_P <- old$Sceptre.P[m]
s$Base_Z <- old$Sceptre.Z[m]
s$Distance <- old$Distance[m]

s$Dom2_FDR <- p.adjust(s$Dom2_P, method = "fdr")
s$Dom5_FDR <- p.adjust(s$Dom5_P, method = "fdr")
s$CellRm_FDR <- p.adjust(s$CellRm_P, method = "fdr")
s$Base_FDR <- p.adjust(s$Base_P, method = "fdr")

table(s$Dom2_FDR < 0.1) # 104
table(s$Dom5_FDR < 0.1) # 64
table(s$CellRm_FDR < 0.1) # 99
table(s$Base_FDR < 0.1) # 104

table(s$Dom2_FDR < 0.1, s$Base_FDR < 0.1) # 104
table(s$Dom5_FDR < 0.1, s$Base_FDR < 0.1)
table(s$CellRm_FDR < 0.1, s$Base_FDR < 0.1)
table(s$Dom2_FDR < 0.1, s$Dom5_FDR < 0.1)


gain <- s[which(s$Dom2_FDR < 0.1 & s$Base_FDR > 0.1),] # 15, 9 down, distance 235
common <- s[which(s$Dom2_FDR < 0.1 & s$Base_FDR < 0.1),] # 79
lost <- s[which(s$Dom2_FDR > 0.1 & s$Base_FDR < 0.1),] # 25, 18 down, distance 114


gain <- s[which(s$Dom2_FDR < 0.1 & s$CellRm_FDR > 0.1),] # 15, 9 down, distance 235
common <- s[which(s$Dom2_FDR < 0.1 & s$CellRm_FDR < 0.1),] # 79
lost <- s[which(s$Dom2_FDR > 0.1 & s$CellRm_FDR < 0.1),] # 20, 18 down, distance 80
# 
# dom2_hits <- s[which(s$Dom2_FDR < 0.1),]
# base_hits <- s[which(s$Base_FDR < 0.1),]
# 
# table(sign(abs(common$Dom2_Z) - abs(common$Base_Z)))
# boxplot(abs(common$Dom2_Z) - abs(common$Base_Z))

# gain <- s[which(s$CellRm_FDR < 0.1 & s$Base_FDR > 0.1),]
# common <- s[which(s$CellRm_FDR < 0.1 & s$Base_FDR < 0.1),]
# lost <- s[which(s$CellRm_FDR > 0.1 & s$Base_FDR < 0.1),]




## For positive controls
pos.z <- data.frame(
  Target = res00.pos$gene_id,
  Z_Frac00 = res00.pos$z_value,
  Z_Frac02 = res02.pos$z_value,
  Z_Frac05 = res05.pos$z_value,
  P_Frac00 = res00.pos$p_value,
  P_Frac02 = res02.pos$p_value,
  P_Frac05 = res05.pos$p_value,
  Adj_Frac00 = p.adjust(res00.pos$p_value, method = "bonferroni"),
  Adj_Frac02 = p.adjust(res02.pos$p_value, method = "bonferroni"),
  Adj_Frac05 = p.adjust(res05.pos$p_value, method = "bonferroni"))

  

pA <- ggplot(pos.z, aes(x = Z_Frac00, y = Z_Frac02)) +
  geom_point() +
  geom_abline(linetype = 2, slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(y = "Negative Binomial Z Statistic\nDom02", x = "Negative Binomial Z Statistic\nAll assignments") +
  theme(panel.border = invis, axis.line = element_line())

pB <- ggplot(pos.z, aes(x = Z_Frac00, y = Z_Frac05)) +
  geom_point() +
  geom_abline(linetype = 2, slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(y = "Negative Binomial Z Statistic\nDom05", x = "Negative Binomial Z Statistic\nAll assignments") +
  theme(panel.border = invis, axis.line = element_line())

pC <- ggplot(pos.z, aes(x = Z_Frac02, y = Z_Frac05)) +
  geom_point() +
  geom_abline(linetype = 2, slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(x = "Negative Binomial Z Statistic\nDom02", y = "Negative Binomial Z Statistic\nDom05") +
  theme(panel.border = invis, axis.line = element_line())  

pdf(file = "../2_DE/Pos/Comparison of Dom on Z.pdf", height = 3, width = 8)
plot_grid(pA, pB, pC, nrow = 1)
dev.off()

p <- data.frame(Z_02_00 = pos.z$Z_Frac02 - pos.z$Z_Frac00,
                Z_05_00 = pos.z$Z_Frac05 - pos.z$Z_Frac00,
                Z_02_05 = pos.z$Z_Frac02 - pos.z$Z_Frac05,
                Target = pos.z$Target)

p <- melt(p)
levels(p$variable) <- c("02 versus No Filt", "05 versus No Filt", "02 versus 05")
pdf(file = "../2_DE/Pos/Comparison of Dom on Z, Boxplot.pdf", height = 4, width = 4)
ggplot(p, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  labs(y = "Difference in NB Z-statistic\n(Negative Means More Significant)") +
    theme(axis.title.x = invis)


ggplot(p, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(y = "Difference in NB Z-statistic\n(Negative Means More Significant)") +
    theme(axis.title.x = invis)
dev.off()

## Compare p-value and significance
pA <- ggplot(pos.z, aes(x = -log10(P_Frac00), y = -log10(P_Frac02))) +
  geom_point() +
  geom_abline(linetype = 2, slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(y = "Sceptre P\nDom02", x = "Sceptre P\nAll assignments") +
  theme(panel.border = invis, axis.line = element_line())

pB <- ggplot(pos.z, aes(x = -log10(P_Frac00), y = -log10(P_Frac05))) +
  geom_point() +
  geom_abline(linetype = 2, slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(y = "Sceptre P\nDom05", x = "Sceptre P\nAll assignments") +
  theme(panel.border = invis, axis.line = element_line())

pC <- ggplot(pos.z, aes(x = -log10(P_Frac02), y = -log10(P_Frac05))) +
  geom_point() +
  geom_abline(linetype = 2, slope = 1, intercept = 0, colour = "red") +
  theme_bw() +
  labs(x = "Sceptre P\nDom02", y = "Sceptre P\nDom05") +
  theme(panel.border = invis, axis.line = element_line())  

pdf(file = "../2_DE/Pos/Comparison of Dom on P.pdf", height = 3, width = 8)
plot_grid(pA, pB, pC, nrow = 1)
dev.off()  

p <- data.frame(Z_02_00 = -log10(pos.z$P_Frac02) - -log10(pos.z$P_Frac00),
                Z_05_00 = -log10(pos.z$P_Frac05) - -log10(pos.z$P_Frac00),
                Z_02_05 = -log10(pos.z$P_Frac02) - -log10(pos.z$P_Frac05),
                Target = pos.z$Target)

p <- melt(p)
levels(p$variable) <- c("02 versus No Filt", "05 versus No Filt", "02 versus 05")
pdf(file = "../2_DE/Pos/Comparison of Dom on Z, Boxplot.pdf", height = 4, width = 4)
ggplot(p, aes(x = variable, y = value)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  labs(y = "Difference in NB Z-statistic\n(Negative Means More Significant)") +
    theme(axis.title.x = invis)


ggplot(p, aes(x = variable, y = value)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(y = "Difference in NB Z-statistic\n(Negative Means More Significant)") +
    theme(axis.title.x = invis)
dev.off()




  