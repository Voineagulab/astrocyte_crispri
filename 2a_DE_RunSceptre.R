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
  
## Step 1: Get expression matrix
  sceptre.exp <- nha@assays$RNA@counts[,sceptre.use] # keep it as a sparse matrix
    
## Step 2: gRNA binary expression matrix
  guide.cols <- grep("Enh|_A$|_B$|Neg", colnames(nha@meta.data))
  sceptre.guide <- as.data.frame(nha@meta.data[sceptre.use,guide.cols])
  x <- apply(sceptre.guide, 1, as.numeric)  # I note that this also transposes teh data, which is necessary
  rownames(x) <- colnames(sceptre.guide)
  sceptre.guide <- as(x, "sparseMatrix")
    
## Step 3: Covariate matrix
  sceptre.covar <- data.frame(UMI_count = log(colSums(sceptre.exp)),
                              gRNA_count = log(colSums(sceptre.guide)),
                              Mito_P = nha$Mito_Pct[sceptre.use],
                              Batch = nha$Library[sceptre.use],
                              C_Cycle = nha$Cycle_Seurat_S[sceptre.use])
      
    
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
## Load from here ----
    
  rm(nha)
  load("Sceptre Input Files.rda")
    
################################################################################################################################ #
## Run sceptre: Positive controls ----
    
a <- Sys.time()

de.pos <- run_sceptre_high_moi(gene_matrix = sceptre.exp,
                               combined_perturbation_matrix = sceptre.guide.pooled,
                               covariate_matrix = sceptre.covar,
                               gene_gRNA_group_pairs = sceptre.pairs.pos,
                               side = sceptre.direction,
                               B = 1000)

b <- Sys.time() # 2.53 minutes!!

save(de.pos, file  = "Pos/SCEPTRE Output.rda")

    
################################################################################################################################ #
## Run sceptre: Enhancers ----
    

## Define list
  de.enh <- list()

## Chunk into runs of 100 tests
  chunkSize <- 100
  nChunks <- ceiling(nrow(sceptre.pairs.enh) / chunkSize)
    
## Run loop (requires 1-2min per chunk on our server)
  for (j in 1:nChunks) {
    
    # setup
    a <- Sys.time()
    run.name <- paste0("Chunk_", j)
    
    # get the current chunk of a given chunkSize
    chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
    chunk <- chunk[which(chunk < nrow(sceptre.pairs.enh))] # truncate if the range is outside the full dataframe
    # k <- which(sceptre.pairs.enh$gRNA_id == j) # in this loop, only test genes within 1mb of the enhancer
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
    print(paste0(run.name, ": ", b-a))
    
  }

## Save  
  de.enh <- do.call("rbind", de.enh)
  save(de.enh, file = "Enh/SCEPTRE Output.rda")


   
################################################################################################################################ #
## Run sceptre: Negative controls ----
     
  
## Using the 50 negative controls in the Enh transduction pool
  de.negE <- list()
  
  ## Pair each guide with all genes tested in the enhancer pool
    used.genes <- unique(c(de.enh$gene_id, de.pos$Gene))

    sceptre.pairs.negE <- data.frame(gene_id = rep(used.genes, each = 50),
                                     gRNA_id = guides$GuideID[grep("Neg_E_", guides$GuideID)],
                                     pair_type = "non_target")    
    
    
  ## Chunk into runs of 200 tests
    chunkSize <- 200
    nChunks <- ceiling(nrow(sceptre.pairs.negE) / chunkSize)
    
  ## Run loop (requires 1-2min per chunk on our server)
    for (j in 1:nChunks) {
      
      # setup
      a <- Sys.time()
      run.name <- paste0("Chunk_", j)
      
      # get the current chunk of a given chunkSize
      chunk <- ((chunkSize*(j-1)) + 1) : (chunkSize*j) # get a range of rows totaling chunkSize
      chunk <- chunk[which(chunk < nrow(sceptre.pairs.negE))] # truncate if the range is outside the full dataframe
      # k <- which(sceptre.pairs.enh$gRNA_id == j) # in this loop, only test genes within 1mb of the enhancer
      chunk <- sceptre.pairs.negE[chunk,uni]
      
      # run sceptre
      
      de.negE[[run.name]] <- run_sceptre_high_moi(gene_matrix = sceptre.exp[chunk$gene_id,],
                                                  combined_perturbation_matrix = sceptre.guide, # note: not sceptre.guide.pooled, as we don't want to pool the negE guides!
                                                  covariate_matrix = sceptre.covar,
                                                  gene_gRNA_group_pairs = chunk,
                                                  side = sceptre.direction,
                                                  B = 1000)
      
      # sign off
      gc()
      b <- Sys.time()
      print(paste0(run.name, ": ", b-a))
      
    }

## Save  
    x <- de.negE
  de.negE <- do.call("rbind", de.negE)
  save(de.negE, file = "Neg/SCEPTRE Output (NegE).rda")

    
#     n <- ceiling(length(used.genes) / 100)
#     for (j in 1:n) { # this chunks sceptre in runs with 125 guide pairs for 100 genes
#       range <- ((100*(j-1)) + 1) : (100*j) # get 100 genes
#       range <- range[which(range < length(used.genes))]
#       
#       sceptre.pairs.neg[[paste0("Genes", min(range), "_", max(range))]] <- data.frame(gene_id = rep(used.genes[range], each = 125),
#                                                                                       gRNA_id = paste0("NegPair_", 1:125),
#                                                                                       pair_type = "non_target")
#     }
#     
#     b <- do.call("rbind", sceptre.pairs.neg)
#     
#     
#     # run
#     sceptre$neg <- list()
#      for (j in names(sceptre.pairs.neg)) {
#       k <- sceptre.pairs.neg[[j]]
#       
#       # run sceptre
#       a <- Sys.time()
#       sceptre$neg[[j]] <- run_sceptre_high_moi(gene_matrix = sceptre.exp[unique(k$gene_id),],
#                                                combined_perturbation_matrix = sceptre.guide.pooled,
#                                                covariate_matrix = sceptre.covar,
#                                                gene_gRNA_group_pairs = k,
#                                                side = sceptre.direction,
#                                                B = 1000)
#       b <- Sys.time()
#     
#       print(paste0(j, ": ", b-a))
#       save(sceptre, file = "Temp.rda")s
#       gc()
#     }
#     
#     
# ## Save results
#   save(de, file = "Sceptre.rda")
#     

  

  