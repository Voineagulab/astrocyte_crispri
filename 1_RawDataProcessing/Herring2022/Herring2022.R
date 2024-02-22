## In this script, we convert the scanpy processed data produced by Herring et al 2022 
## to Seurat. Then, these data are pseudobulked by stage, individual, and cell-type.

################################################################################################################################ #
## Setup ----

setwd("/mnt/Data0/PROJECTS/CROPSeq/PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/")
library(anndata)
library(SeuratDisk)
library(Seurat)
source("../../../../FullScale/Scripts/Functions.R")

## The direct reading in of the h5ad via SeuratDisk returns the following error:
  # Error in dfgroup[[tname]][] : 
  # object of type 'environment' is not subsettable

## Thus I will first use anndata to load it in
  scanpy <- read_h5ad("Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad")
  
  
################################################################################################################################ #
## Get counts ----

  
## Collect the counts as a sparse matrix
  a <- (scanpy$X)
  
  
  # this is a dgr, whilst Seurat would like a dgc
  # use code from: https://rdrr.io/github/barkasn/nbHelpers/
  transpose_dgRMatrix <- function(inmat) {
    if(class(inmat) != 'dgRMatrix')
        stop('inmat is not of class dgRMatrix')
    out <- new('dgCMatrix',
               i=inmat@j,
               p=inmat@p,
               x=inmat@x,
               Dim=rev(inmat@Dim),
               Dimnames=rev(inmat@Dimnames)
               )
    out
}      
  
  a <- transpose_dgRMatrix(a)
  

## Read in metadata
  # meta <- read.csv("Processed_data_RNA-all_BCs-meta-data.csv", row.names = 1)
  meta <- read.csv("Processed_data_RNA-all_BCs-meta-data_updated.csv", row.names = 1)

## Create seurat object
  brainMat <- CreateSeuratObject(counts = a,
                                 project = "Herring2022",
                                 meta.data = meta,
                                 min.cells = 0,
                                 min.features = 0)
  
## Save
  save(brainMat, file = "Processed/scRNAseq_UMIs_Seurat.rda")
  
## (Oh, and save the UMAP coordinates)
  write.csv(scanpy$obsm$X_umap, file = "Processed/UMAP_coordinates.csv")
  # write.csv(scanpy$obsm$X_pca, file = "Processed/PCA_coordinates.csv")
  

  
################################################################################################################################ #
## Pseudobulk ----
  
  
## Load, if necessary
  # load("Processed/scRNAseq_UMIs_Seurat.rda")
  
  
## Setup
  # setup
  pb <- list()
  inds <- unique(brainMat$batch)
  # ct.type <- c("PN", "IN", "Astro", "Micro", "Oligo", "OPC", "Vas")
  ct.type <- c("PN", "IN", "Astro", "Micro", "Oligo", "OPC")
  
  ct.annot <- brainMat$cell_type
  w <- which(ct.annot == "Non-Neu")
  ct.annot[w] <- brainMat$major_clust[w]
  
 
  ## Collect expression profiles
    pb$Exp <- list()
    pb$Meta <- list()
    
  
    for (j in ct.type) {
      
      print(j)
      
      for (k in inds) {
        # id of loop
        id <- paste0(k, "-", j)
        
        # index of cells
        w <- which( (brainMat$batch == k) & (ct.annot == j) )
        nonzero <- length(w) != 0
        
        
        # get expression data
        if (nonzero) {
          pb$Exp[[id]] <- rowSums(brainMat@assays$RNA@counts[,w])  
        } else {
          pb$Exp[[id]] <- NA
        }
        
        # note the metadata 
        pb$Meta[[id]] <- data.frame(FinalID = NA,
                                    OriginalID = id,
                                    Celltype = j,
                                    Individual = k,
                                    Stage = NA,
                                    nCells = length(w),
                                    nUMI = sum(pb$Exp[[id]]),
                                    Chemistry = NA)
        
      }
    }
  
## Bind
  pb$Exp <- do.call("cbind", pb$Exp)
  pb$Meta <- do.call("rbind", pb$Meta)
    
  
## Fill out the metadata
  # stage
  m <- match(pb$Meta$Individual, brainMat$batch) # not an issue that there are many matches, as it takes the first
  pb$Meta$Stage <- brainMat$stage_id[m]
  
  # chemistry
  pb$Meta$Chemistry <- splitter(pb$Meta$OriginalID, "-", 1) %>% splitter(., "_", 3) %>% toupper()
  
  # celltype
  pb$Meta$Celltype <- gsub("PN", "Exc", pb$Meta$Celltype)
  pb$Meta$Celltype <- gsub("IN", "Inh", pb$Meta$Celltype)
  
  # individual
  pb$Meta$Individual <- splitter(pb$Meta$Individual, "_", 1)
  
  # final id
  pb$Meta$FinalID <- rownames(pb$Meta) <- paste(pb$Meta$Celltype, pb$Meta$Stage, pb$Meta$Individual, sep = "_")
  
  
## Rename sample ids in expression, as easier to parse
  colnames(pb$Exp) <- pb$Meta$FinalID
  
## Add gene metadata
  genes <- read.csv("Processed_data_RNA-all_genes-meta-data.csv")
  pb$GeneInfo <- genes
  
## Save
  save(pb, file = "Processed/Pseudobulk_byGJS.rda")
  
  
################################################################################################################################ #
## Pseudobulk, 2023-16-10 ----

## This pseudobulk aims to match the ATACseq pooling process
## No individuals
## Two inh subtypes
## Three exc subtypes
## Once subtypes are pseudobulk, normalise to depth, then average the exc and inh subtypes
  
  
  
## Setup
  # setup
  pb2 <- list()

  ct.type <- list(Astro = "Astro",
                  Oligo = c("Oligo", "OPC"),
                  Micro = "Micro",
                  Inh_CGE = c("ID2", "VIP", "LAMP5_NOS1", "CGE_dev"),
                  Inh_MGE = c("SST", "PV", "PV_SCUBE3", "MGE_dev"),
                  Exc_L23 = c("L2-3_CUX2", "PN_dev"),
                  Exc_L4 = c("L4_RORB", "PN_dev"),
                  Exc_L56 = c("L5-6_THEMIS", "L5-6_TLE4", "PN_dev"))
  
  # note that PN_dev appears in three different exc subtypes. here, these share the same progenitor pool (i.e., foetal and neonatal stages), as was used in the original paper for trajectory analyses
  # however, given that thesesubtypes will be averaged, it makes little difference to the final expression value
  
  stages <- unique(brainMat$stage_id)[c(5,1,3,2,4,6)]
 
  ## Collect expression profiles
    pb2$Exp <- list()
    pb2$Meta <- list()
    
  
    for (j in names(ct.type)) {
      
      for (k in stages) {
        
        # id of loop
        id <- paste0(j, "_", k)
        print(id)
        
        # index of cells
        w <- which( (brainMat$stage_id == k) & (brainMat$major_clust %in% ct.type[[j]]) )
        nonzero <- length(w) != 0

        # get expression data
        if (nonzero) {
          pb2$Exp[[id]] <- rowSums(brainMat@assays$RNA@counts[,w])  
        } else {
          pb2$Exp[[id]] <- NA
        }
        
        # note the metadata 
        pb2$Meta[[id]] <- data.frame(OriginalID = id,
                                     Celltype = j,
                                     Stage = k,
                                     nCells = length(w),
                                     nUMI = sum(pb2$Exp[[id]]))
        
      }
    }
  
## Bind
  pb2$Exp <- do.call("cbind", pb2$Exp)
  pb2$Meta <- do.call("rbind", pb2$Meta)
  
## Add gene metadata
  genes <- read.csv("Processed_data_RNA-all_genes-meta-data.csv")
  pb2$GeneInfo <- genes
  
## Normalise to cpm
  pb2$Cpm <- apply(pb2$Exp, 2, function(x) {
    x / (sum(x) / 10^6)
  })
  
## Combine subtypes
  pb2$Final <- pb2$Cpm[,grepl("^Astro|^Oligo|^Micro", colnames(pb2$Cpm))] %>% as.data.frame()
  
  for (j in c("Exc", "Inh")) {
    
    for (k in stages) {
      
      id <- paste0(j, "_", k)
      
      which_j <- grep(paste0("^", j), colnames(pb2$Cpm))
      which_k <- grep(paste0("_", k), colnames(pb2$Cpm))
      w <- intersect(which_j, which_k)
      
      pb2$Final[,id] <- rowMeans(pb2$Cpm[,w]) # the mean of depth-normalised xpression across the subtypes
      
    }
  }
  
  # in effect, this is the non-weighted average expression. it does not account for differences in subtype abundance.
  # this mimics the extent of possible normalisation in the snATAC
  
## Save
  # save(pb, file = "Processed/Pseudobulk_byGJS.rda")
  save(pb2, file = "Processed/Pseudobulk2_ATACesquePooling.rda")
  