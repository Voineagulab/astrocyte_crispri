## This script runs differential expression

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(Seurat)
  library(sceptre)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(rcartocolor)
  library(tidyverse)

## Load
  source("../../Scripts/Functions.R")
  load("../../Data/Preprocessed/NHA Pooled (Final).rda")
  load("Sceptre Input Files.rda")
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)


## Data information
  targ.pos <- guides$TargetID[which(guides$TargetCat == "Promoter")] %>% unique()
  targ.enh <- guides$TargetID[which(guides$TargetCat == "Enh")] %>% unique()
  targ.neg <- guides$GuideID[which(guides$TargetCat == "Negative")] %>% unique()
  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()


################################################################################################################################ #
## Prepare objects for sceptre ----
  
## Per the tutorial at https://katsevich-lab.github.io/sceptre/articles/using_sceptre_v2.html, this is a 7 step process
  
## Step 0: Define cells
  # use cells with any guide
  sceptre.use <- which(nha$AnyGuide)
  # sceptre.use <- which(nha$AnyGuide & nha$MOI <= 20) # change on 1/8/2022

## Step 1: Get expression matrix
  sceptre.exp <- nha@assays$RNA@counts[,sceptre.use] # keep it as a sparse matrix

## Step 2: gRNA binary expression matrix
  guide.cols <- grep("Enh|Pos|Neg", colnames(nha@meta.data))
  sceptre.guide <- as.data.frame(nha@meta.data[sceptre.use,guide.cols])
  x <- apply(sceptre.guide, 1, as.numeric)  # I note that this also transposes teh data, which is necessary
  rownames(x) <- colnames(sceptre.guide)
  sceptre.guide <- as(x, "sparseMatrix")

## Step 3: Covariate matrix
    sceptre.covar <- data.frame(UMI_count = log(colSums(sceptre.exp)),
                                gRNA_count = log(nha$MOI[sceptre.use]), # note: cannot use colSums(sceptre.guide), as this double counts the supermajority of guides but not all...
                                Mito_P = nha$Mito_Pct[sceptre.use],
                                Batch = nha$Library[sceptre.use],
                                C_Cycle = nha$Cycle_Tricycle_PC1[sceptre.use],
                                NuclearFraction = nha$NuclearFraction[sceptre.use],
                                Ribo_P = nha$Ribo_Pct[sceptre.use])

## Step 4: Group and categorise guides
  sceptre.groups <- data.frame(gRNA_id = guides$GuideID,
                               gRNA_group = guides$TargetID,
                               gRNA_type = factor(guides$TargetCat))

  # for negative control guides, create "groups" for a non-targeting and AAVS1-targeting guide
  for (j in 1:125) {
    m <- match(paste0("Neg_", j), sceptre.groups$gRNA_id)
    sceptre.groups$gRNA_group[c(m, m+125)] <- paste0("NegPair_", j)

  }

  for (j in 1:25) {
    m <- match(paste0("Neg_E_", j), sceptre.groups$gRNA_id)
    sceptre.groups$gRNA_group[c(m, m+25)] <- paste0("NegPair_E_", j)

  }

  levels(sceptre.groups$gRNA_type) <- c("enh_target", "non_target", "tss_target")

  # now also pool perturbations!
  sceptre.guide.pooled <- combine_perturbations(perturbation_matrix = sceptre.guide,
                                                gRNA_groups_table = sceptre.groups)


## Step 5: define pairs to test
  ## Positive controls
    # p <- targ.pos[which(targ.pos %in% rownames(sceptre.exp))]
    sceptre.pairs.pos <- data.frame(gene_id = targ.pos,
                                    gRNA_id = targ.pos,
                                    pair_type = "positive_control")

    # two genes are not present in the expression data; this is because a synonym is used...
    sceptre.pairs.pos$gene_id[sceptre.pairs.pos$gRNA_id == "FAM96B"] <- "CIAO2B" # synonym
    sceptre.pairs.pos$gene_id[sceptre.pairs.pos$gRNA_id == "SEPT11"] <- "SEPTIN11" # synonym


  ## Enhancers
    sceptre.pairs.enh <- list()
    for (j in targ.enh) {
      print(j)

      # get enh coordinates
      m <- match(j, guides$TargetID)
      w <- guides$TargetCoord[m]


      # find nearby tss to the enh
      nearby <- find.nearby.tss(query = w, expand.by = 10^6, reduced.output = TRUE)
      nearby <- intersect(nearby, rownames(nha))
      sceptre.pairs.enh[[j]] <- data.frame(gene_id = nearby, gRNA_id = j, pair_type = "candidate") # recommended name for pair-type
    }

    sceptre.pairs.enh <- do.call("rbind", sceptre.pairs.enh)


## Step 6: direction to test
    sceptre.direction <- "both"

## Step X: save!
    save(sceptre.covar,
         sceptre.direction,
         sceptre.exp,
         sceptre.groups,
         sceptre.guide,
         sceptre.guide.pooled,
         sceptre.pairs.enh,
         sceptre.pairs.pos,
         file = "Sceptre Input Files.rda")
    

    
################################################################################################################################ #
## Run sceptre: Positive controls ----
    
## Run
  de.pos <- run_sceptre_high_moi(gene_matrix = sceptre.exp,
                                 combined_perturbation_matrix = sceptre.guide.pooled,
                                 covariate_matrix = sceptre.covar,
                                 gene_gRNA_group_pairs = sceptre.pairs.pos,
                                 side = sceptre.direction,
                                 B = 1000)
  

  save(de.pos, file  = "Pos/SCEPTRE Output.rda")
  
  
    
################################################################################################################################ #
## Run sceptre: Enhancers ----
    

## Define list
  de.enh <- list()

## Chunk 
  chunkSize <- 200
  nChunks <- ceiling(nrow(sceptre.pairs.enh) / chunkSize)
    
## Run loop (requires 1-2min per chunk on our server)
  for (j in 1:nChunks) {
    
    # setup
    a <- Sys.time()
    run.name <- paste0("Chunk_", j)
    
    # get the current chunk of a given chunkSize
    chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
    chunk <- chunk[which(chunk < nrow(sceptre.pairs.enh))] # truncate if the range is outside the full dataframe
    chunk <- sceptre.pairs.enh[chunk,]
    
    # run sceptre
    de.enh[[run.name]] <- run_sceptre_high_moi(gene_matrix = sceptre.exp[chunk$gene_id,],
                                                         combined_perturbation_matrix = sceptre.guide.pooled,
                                                         covariate_matrix = sceptre.covar,
                                                         gene_gRNA_group_pairs = chunk,
                                                         side = sceptre.direction,
                                                         B = 1000)
    
    # sign off
    b <- Sys.time()
    gc()
    print(paste0(run.name, ": ", b-a))
    
  }

## Save  
  de.enh <- do.call("rbind", de.enh)
  save(de.enh, file = "Enh/SCEPTRE Output.rda")


################################################################################################################################ #
## Run sceptre: Enhancers, but guide-specific rather than pooled at each target enhancer ----
  
## Generate a guide-gene correspondence for use with sceptre
  sceptre.pairs.enh.guidelvl <- list()
  for (j in targ.enh) {
    w <- sceptre.pairs.enh[which(sceptre.pairs.enh$gRNA_id == j),]
    x <- guides$GuideID[which(guides$TargetID == j)]
    y <- lapply(x, function(k) {
      z <- w
      z$gRNA_id <- k
      return(z)
    })
    
    z <- do.call("rbind", y)
    sceptre.pairs.enh.guidelvl[[j]] <- z
  }
  
  sceptre.pairs.enh.guidelvl <- do.call("rbind", sceptre.pairs.enh.guidelvl)
  
## Define list
  de.enh.guidelvl <- list()

## Chunk into runs of 100 tests
  chunkSize <- 200
  nChunks <- ceiling(nrow(sceptre.pairs.enh.guidelvl) / chunkSize)
    
## Run loop (requires 1-2min per chunk on our server)
  for (j in (length(de.enh.guidelvl) + 1):nChunks) {
    
    # setup
    a <- Sys.time()
    run.name <- paste0("Chunk_", j)
    # if (j %in% c(231, 274)) de.enh.guidelvl[[run.name]] <- NA
    
    # get the current chunk of a given chunkSize
    chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
    chunk <- chunk[which(chunk < nrow(sceptre.pairs.enh.guidelvl))] # truncate if the range is outside the full dataframe
    # k <- which(sceptre.pairs.enh$gRNA_id == j) # in this loop, only test genes within 1mb of the enhancer
    chunk <- sceptre.pairs.enh.guidelvl[chunk,]
    
    # run sceptre
    de.enh.guidelvl[[run.name]] <- run_sceptre_high_moi(gene_matrix = sceptre.exp[chunk$gene_id,],
                                                        combined_perturbation_matrix = sceptre.guide,
                                                        covariate_matrix = sceptre.covar,
                                                        gene_gRNA_group_pairs = chunk,
                                                        side = sceptre.direction,
                                                        B = 1000)
    
    # sign off
    b <- Sys.time()
    gc()
    print(paste0(run.name, ": ", b-a))
    
  }
  
  de.enh.guidelvl <- do.call("rbind", de.enh.guidelvl)
  save(de.enh.guidelvl, file = "Enh/Guide-level/Guide-level SCEPTRE.rda")
     
################################################################################################################################ #
## Run sceptre: Negative controls ----
     
## 
  
## Using the 50 negative controls in the Enh transduction pool
  de.negE <- list()

  ## Pair each guide with all genes tested in the enhancer pool
    x <- read.csv("Enh/Results Summary.csv")
    used.genes <- unique(c(x$Gene))

    sceptre.pairs.negE <- data.frame(gene_id = rep(used.genes, each = 50),
                                     gRNA_id = guides$GuideID[grep("Neg_E_", guides$GuideID)],
                                     pair_type = "non_target")


  ## Chunk into runs of 200 tests
    chunkSize <- 200
    nChunks <- ceiling(nrow(sceptre.pairs.negE) / chunkSize)
    save.points <- seq(0, nChunks, 100)
    
  ## Run loop (requires 1-2min per chunk on our server)
    for (j in 1:nChunks) {

      # setup
      a <- Sys.time()
      run.name <- paste0("Chunk_", j)

      # get the current chunk of a given chunkSize
      chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
      chunk <- chunk[which(chunk < nrow(sceptre.pairs.negE))] # truncate if the range is outside the full dataframe
      # k <- which(sceptre.pairs.enh$gRNA_id == j) # in this loop, only test genes within 1mb of the enhancer
      chunk <- sceptre.pairs.negE[chunk,]

      # run sceptre

      de.negE[[run.name]] <- run_sceptre_high_moi(gene_matrix = sceptre.exp[chunk$gene_id,],
                                                  combined_perturbation_matrix = sceptre.guide, # note: not sceptre.guide.pooled, as we don't want to pool the negE guides!
                                                  covariate_matrix = sceptre.covar,
                                                  gene_gRNA_group_pairs = chunk,
                                                  side = sceptre.direction,
                                                  B = 1000)

      # sign off
      b <- Sys.time()
      gc()
      print(paste0(run.name, ": ", b-a))
      if (j %in% save.points) save(de.negE, file = paste0("Temp NegE Guide-level (", j, ").rda"))

    }

  ## Save  
    de.negE <- do.call("rbind", de.negE)
    save(de.negE, file = "Neg/SCEPTRE Output (NegE, Guide-level).rda")


## Using the 250 negative controls in the Neg transduction pool, at the guide-level
    sceptre.pairs.neg.guidelevel <- data.frame(gene_id = rep(used.genes, each = 250),
                                     gRNA_id = guides$GuideID[grep("AAVS1|NoTarget", guides$TargetID)],
                                     pair_type = "non_target")

    # b <- do.call("rbind", sceptre.pairs.neg)


    # run
    
    de.neg.guidelevel <- list()
    chunkSize <- 300 # just to ensure there's > 1 guide in the matrix...
    nChunks <- ceiling(nrow(sceptre.pairs.neg.guidelevel) / chunkSize)
    save.points <- seq(0, nChunks, 100)

  ## Run loop (requires 1-2min per chunk on our server)
    # for (j in 1:nChunks) {
    for (j in (length(de.neg.guidelevel)+1):nChunks) {

      # setup
      a <- Sys.time()
      run.name <- paste0("Chunk_", j)

      # get the current chunk of a given chunkSize
      chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
      chunk <- chunk[which(chunk < nrow(sceptre.pairs.neg.guidelevel))] # truncate if the range is outside the full dataframe
      chunk <- sceptre.pairs.neg.guidelevel[chunk,]

      # run sceptre

      de.neg.guidelevel[[run.name]] <- run_sceptre_high_moi(gene_matrix = sceptre.exp[chunk$gene_id,],
                                                  combined_perturbation_matrix = sceptre.guide, 
                                                  covariate_matrix = sceptre.covar,
                                                  gene_gRNA_group_pairs = chunk,
                                                  side = sceptre.direction,
                                                  B = 1000)
      b <- Sys.time()

      print(paste0(j, ": ", b-a))
      if (j %in% save.points) save(de.neg.guidelevel, file = paste0("Temp Neg Guidelevel (", j, ").rda"))
      gc()
    }


  ## Save results
    de.neg.guidelevel <- do.call("rbind", de.neg.guidelevel)
    save(de.neg.guidelevel, file =  "Neg/SCEPTRE Output (Neg, Guide-level).rda")

  
  