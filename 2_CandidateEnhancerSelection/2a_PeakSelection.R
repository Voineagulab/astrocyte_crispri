## Herein I explore and annotate open chromatin regions (OCRs) across numerous brain datasets

################################################################################################################################ #
## Setup ----

## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/")
options(stringsAsFactors = FALSE)


## Libraries
library(Seurat)


## Functions
  # convert counts to cpm
  make.cpm <- function(counts) {
    lib.size <- colSums(counts) # total reads per sample
    lib.size <- 10^6 / lib.size # normalise to per million
    cpm <- counts; cpm[,] <- 0
  
    for(j in colnames(counts)) {
      cpm[,j] <- counts[,j] * lib.size[j]
    }
    return(cpm)
  }
  
  # what percentage of peaks are above a given expression threshold?
  check.thresh <- function(cpm, thresh) {
    apply(cpm, 2, function(x) length(which(x > thresh)) / length(x))
  }
  
  # get the indices of peaks crossing a given threshold
  threshold <- function(cpm, n = round(ncol(cpm) / 2), thresh) {
    keep <- apply(cpm, 1, function(x) length(which(x > thresh))) # at least half of all samples are above the given threshold
    names(keep) <- rownames(cpm)
    keep <- keep[which(keep >= n)]
    return(keep)
  }
  
  # write a bed file
  write.bed <- function(x, dir) {
    write.table(x, dir, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
  
  
## Key parameters!
  cpm.threshold <- 2 # this is the cpm threshold used to filter peaks
  extend.genes.by <- 2000 # genes in gencode V32 will have their coordinates expanded by this amount to either side
  extended.gencode.dir <- paste0("/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/PublicData_forLibrarySelectionOnly/Other/gencodeV32_extend", extend.genes.by, ".bed")
  
  bed.intersect.dir <- "/Applications/bedtools2/bin/intersectBed" # location of intersect bed
  bed.annotate.dir <- "/Applications/bedtools2/bin/annotateBed" # location of annotate bed
  bed.getfasta.dir <- "/Applications/bedtools2/bin/fastaFromBed" # location of getFasta (equivalent names)
  
  inhouse.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/Bulk_ATAC_RNAseq/GOK8505-GOK8837/GOK8837A1/Mapping/MACS2/NHA_ATAC_S3.filtered.BAM_peaks.bed"
  
################################################################################################################################ #
## Load in data, and ready for intersections! ----

  
## Whilst I am not process this in, I will note some statistics of it: in-house atac-seq from cultured astrocytes
  GOK8873A1 <- read.table(inhouse.dir, sep = "\t")
  colnames(GOK8873A1) <- c("chr",	"start", "end",	"length",	"abs_summit",	"pileup",	"-LOG10(pvalue)",	"fold_enrichment",	"-LOG10(qvalue)",	"name")

  # get annotation for highly-expressed peaks
  high.expression <- GOK8873A1
  high.expression$cpm <- (high.expression$pileup * 10^6) / sum(high.expression$pileup)
  high.expression <- high.expression[which(high.expression$cpm > 2),] # from 260,193 to 140,813
  high.expression <- paste0(high.expression$chr, high.expression$start, high.expression$end)  
    
## Data from Song et al. 2019: atac-seq from cultured astrocytes
  # load
  song <- read.table("PublicData_forLibrarySelectionOnly/CulturedCells/Song2019/astrocyte.atac.seq.peaks (1).narrowPeak", sep = "\t", header = TRUE)

  # the above file contains a bed-like format. 
  # prep for conversion to hg38
  song[,4] <- paste0("hg19_", song$chrom, "-", song$chromStart, "-", song$chromEnd) # note that the LiftOver tool despises ":"
  write.table(song[,1:4], file = "PublicData_forLibrarySelectionOnly/CulturedCells/Song2019/astrocyte_hg19.bed", 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE) 
  
  # use liftover online tool to create hg38

  
## Data from Trevino et al 2020: atac-seq from organoids
  # load  
  trev.counts <- read.table("PublicData_forLibrarySelectionOnly/CulturedCells/Trevino2020/GSE132403_Data_04_merged_peak_counts_updated2.tsv", sep = "\t", header = TRUE, row.names = 1) # in this are peak labels and their expression
  trev.peaks <- read.table("PublicData_forLibrarySelectionOnly/CulturedCells/Trevino2020/GSE132403_Data_02_merged_peaks_updated.tsv", sep = "\t", header = TRUE) # in this are the peak labels and their corresponding coordinates
  
  # understand the meta data
  trev.meta <- strsplit(colnames(trev.counts), "_")
  trev.meta <- as.data.frame(do.call("rbind", trev.meta)) # i note that this fails to cleanly bind the data, as some rows have four entries instead of five. however, these samples are not ones in which we're interested
  rownames(trev.meta) <- colnames(trev.counts)
  
  # filter data to glial samples, from human corticospheroids
  keep <- which(trev.meta$V2 == "glial" & trev.meta$V5 == "hCS")
  trev.meta <- trev.meta[keep,]
  trev.counts <- trev.counts[,keep]
  
  # determine a threshold for "highly-expressed"
  trev.cpm <- make.cpm(trev.counts)
  
  check.thresh(trev.cpm, 1)
  check.thresh(trev.cpm, 2) # I'll use this...
  check.thresh(trev.cpm, 5)
  
  # apply threshold
  keep <- threshold(trev.cpm, thresh = cpm.threshold)
  
  # create a bed file of coordinates for these peaks
  trev.bed <- trev.peaks[which(trev.peaks$peak %in% names(keep)),]
  write.table(trev.bed, file = "PublicData_forLibrarySelectionOnly/CulturedCells/Trevino2020/Glia_hCS_cpm2_hg38.bed", 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE) # I gleefully note that it's already hg38
  
  
## From Rizzardi 2019, PFC/NAc and ATAC Seq, split into NeuN+/-, n=22 adults
  # load
  rizz <- read.table("PublicData_forLibrarySelectionOnly/AdultBrain/Rizzardi2019/GSE96614_flow_sorted_brain.ATAC-seq_counts.txt", sep = "\t", header = TRUE)
  
  # split into two matrices: bed peaks and expression
  rizz.peaks <- rizz[,1:3]
  rizz.counts <- rizz[,-c(1:3)]
  
  # collect metadata
  rizz.meta <- strsplit(colnames(rizz.counts), split = "\\.")
  rizz.meta <- as.data.frame(do.call("rbind", rizz.meta))
  
  # filter data to glial (NeuN- samples) from the PFC
  keep <- which(rizz.meta$V2 == "BA9" & rizz.meta$V3 == "neg")
  rizz.counts <- rizz.counts[,keep]
  
  # determine a threshold for peaks
  rizz.cpm <- make.cpm(rizz.counts)
  
  check.thresh(rizz.cpm, 1)
  check.thresh(rizz.cpm, 2) # I'll use this... mostly for consistency!
  check.thresh(rizz.cpm, 5)
  
  # apply threshold!
  keep <- threshold(rizz.cpm, thresh = cpm.threshold)
  
  # create a bed file of coordinates 
  rizz.bed <- rizz.peaks[names(keep),]
  rizz.bed$id <- paste0("hg19_", rizz.bed$seqnames, "-", rizz.bed$start, "-", rizz.bed$end) # such that, when liftover occurs, a record of the original remains
  write.table(rizz.bed, file = "PublicData_forLibrarySelectionOnly/AdultBrain/Rizzardi2019/NeuN-_cpm2_BA_hg19.bed", 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # liftover using online tool
  
  
## From Fullard 2018, NeuN+/- ATAC seq from 14 regions in 5 adult humans
  # load
  full.counts <- read.table("PublicData_forLibrarySelectionOnly/AdultBrain/Fullard2018/counts.tsv", sep = "\t")
  full.peaks <- read.table("PublicData_forLibrarySelectionOnly/AdultBrain/Fullard2018/peak_hg38.bed", sep = "\t")
  
  # collect metadata
  full.meta <- strsplit(colnames(full.counts), "_")
  full.meta <- as.data.frame(do.call("rbind", full.meta))
  
  # filter to glial samples from the frontal cortex: DLPFC, OFC, and VLPFC 
  keep <- which(full.meta$V2 == "G" & full.meta$V3 %in% c("DLPFC", "VLPFC", "OFC"))
  full.counts <- full.counts[,keep]
  full.meta <- full.meta[keep,]
  
  # threshold
  full.cpm <- make.cpm(full.counts)
  
  check.thresh(full.cpm, 1)
  check.thresh(full.cpm, 2) # stick with this, as it returns ~100K peaks, on a similar scale to the other datasets. 
  check.thresh(full.cpm, 5)
  
  # apply threshold
  keep <- threshold(full.cpm, thresh = cpm.threshold)
  full.bed <- full.peaks[which(full.peaks$V4 %in% names(keep)),]
  write.table(full.bed, file = "PublicData_forLibrarySelectionOnly/AdultBrain/Fullard2018/NeuN-_cpm2_FC_hg38.bed", 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  
## From Rajaran 2018, of HiC from cultured glia (differentiated from iPSC)
  # load in
  raja <- read.table("PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_blocks.bed", sep = "\t")
  
  # modify and save
  raja$V1 <- paste0("chr", raja$V1)
  write.table(raja, file = "PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_blocks_modified.bed", 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # the above file was lifted over to hg38 using UCSC
  
  # intersect with gencodeV32 genes. this is a shell call, whose output must later be loaded in
  system("bedtools intersect -a PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38.bed -b PublicData_forLibrarySelectionOnly/Other/gencode32_bed.txt -bed -wa -wb > PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38_geneOverlap.bed",
         intern = FALSE, wait = TRUE)
  
  # load
  raja.genes <- read.table("PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38_geneOverlap.bed", sep = "\t")
  raja.genes <- raja.genes[,c(1,2,3,7)]
  raja.genes$V7 <- sapply(strsplit(raja.genes$V7, "\\."), "[", 1)

  ## Next: load in SC data, aiming to collect genes expressed at > 0.3...
  
    # this section is commented out, due to computational run time. the relevant output is loaded as part of the script, however
  
    # load("../Pilot/Data/seuratObjects/allBarcodes_Annotated_iSeq.rda")
    # load(file = "../Pilot/Results/GJS getProtospacers.rda")
    # obj <- obj[1:2]
    # 
    # obj$C1$has.guide <- (colnames(obj$C1) %in% hits.count$C1$cell_barcode)
    # obj$C2$has.guide <- (colnames(obj$C2) %in% hits.count$C2$cell_barcode)
    # 
    # obj <- lapply(obj, function(x) {
    #   x$Assignment <- "None"
    #   x$Assignment[which(x$cellranger.call & !(x$has.guide))] <- "CR Only"
    #   x$Assignment[which(!(x$cellranger.call) & x$has.guide)] <- "PS Only"
    #   x$Assignment[which(x$cellranger.call & x$has.guide)] <- "CR+PS"
    #   return(x)
    # })  
    # 
    # for (j in names(obj)) {
    #   for (k in levels(as.factor(hits.count[[j]]$protospacer))) {
    #     use <- hits.count[[j]][which(hits.count[[j]]$protospacer == k),]
    #     use$cell_barcode <- as.character(use$cell_barcode)
    #     pure <- use$cell_barcode[which(use$pure)]
    #     dom <- use$cell_barcode[which(use$dom)]
    #     minor <- use$cell_barcode[which(use$minor)]
    #     
    #     # the default is absent
    #     obj[[j]]@meta.data[,k] <- "Absent"
    #     
    #     # annotate other situations
    #     obj[[j]]@meta.data[pure,k] <- "Pure"
    #     obj[[j]]@meta.data[dom,k] <- "Dominant"
    #     obj[[j]]@meta.data[minor,k] <- "Minor"
    #   }
    # }
    # 
    # for (j in names(obj)) {
    #   neg <- obj[[j]]@meta.data
    #   neg <- neg[,grep("neg", colnames(neg))]
    #   
    #   any.pure <- apply(neg, 1, function(x) length(grep("Pure", x)) == 1)
    #   any.dom <- apply(neg, 1, function(x) length(grep("Dominant", x)) == 1)
    #   obj[[j]]$neg.pure <- any.pure
    #   obj[[j]]$neg.dom <- any.dom
    #   obj[[j]]$neg.pure.dom <- any.pure | any.dom
    # }
    # 
    # obj.mg <- lapply(obj, function(x) subset(x, subset = `has.guide` & `cellranger.call`))
    # obj.mg <- merge(x = obj.mg$C1, y = obj.mg$C2, add.cell.ids = c("C1", "C2"), project = "SC")
    # obj.mg <- subset(obj.mg, subset = nCount_RNA > 10000)
    # obj.mg <- NormalizeData(object = obj.mg, normalization.method = "LogNormalize", scale.factor = 10000) # standard parameters for Seurat
    # rm(obj); gc()
    # 
    # mean.exp <- rowMeans(obj.mg@assays$RNA@data[,which(obj.mg$neg.pure.dom)])
    # mean.exp <- mean.exp[order(mean.exp, decreasing = TRUE)]
    # exp.genes <- names(mean.exp)[which(mean.exp > 0.3)]
    # write.csv(mean.exp, file = "PublicData_forLibrarySelectionOnly/Other/SC_Ast_meanExp.csv")
  
  mean.exp <- read.csv("PublicData_forLibrarySelectionOnly/Other/SC_Ast_meanExp.csv")
  exp.genes <- mean.exp$X[which(mean.exp$x > 0.3)]
  
  # ...filter Rajaran's TADs to those with expressed genes!
  library(EnsDb.Hsapiens.v79)
  library(ensembldb)
  exp.txid <- ensembldb::select(EnsDb.Hsapiens.v79, keys = exp.genes, keytype = "SYMBOL", columns = c("SYMBOL", "TXID"))
  raja.genes <- raja.genes[which(raja.genes$V7 %in% exp.txid$TXID),]
  dup <- duplicated(raja.genes[,1:3])
  raja.genes <- raja.genes[-which(dup),]
  write.table(raja.genes, file = "PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38_geneOverlap_exp0.3.bed",
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  
## Processing Gencode V32. I'm simply extending the coordinates of each gene by (extend.genes.by) to either side
  gencode <- read.table("PublicData_forLibrarySelectionOnly/Other/gencode32_bed.txt")
  gencode <- gencode[,1:6]
    
  gencode$V2 <- gencode$V2 - extend.genes.by
  gencode$V3 <- gencode$V3 + extend.genes.by
  
  negative <- which(gencode$V2 < 0)
  gencode$V2[negative] <- 0

  
  write.table(gencode, file = extended.gencode.dir, 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # note that these extended gene bodies were not used for the TAD analyses

  

################################################################################################################################ #
## Intersect! Compare your NHA ATAC-seq to the above annotations ----
  
  
## Annotate in-house data with overlaps to other datasets!
  grand.annotation.dir <- "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/GOK8873A1_annotated.bed"
  
  call.grand.annotation <- paste("bedtools annotate",
                                  "-i", inhouse.dir,
                                  "-files", # after this comes a list of beds with which you will intersect
                                  extended.gencode.dir,
                                  "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/PublicData_forLibrarySelectionOnly/AdultBrain/Fullard2018/NeuN-_hg38.bed", # oversight: not actually using the 2cpm 
                                  "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/PublicData_forLibrarySelectionOnly/AdultBrain/Rizzardi2019/NeuN-_cpm2_BA_hg38.bed",
                                  "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Song2019/astrocyte_hg38.bed", 
                                  "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Trevino2020/Glia_hCS_cpm2_hg38.bed",  
                                  "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38_geneOverlap_exp0.3.bed",
                                  "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/PublicData_forLibrarySelectionOnly/AdultBrain/PsychENCODE_PEC_enhancers_hg38.bed",
                                  ">", grand.annotation.dir,
                                  sep = " ")
  
  system(call.grand.annotation, intern = FALSE, wait = TRUE)
  
## Read in the above annotation
  # read in
  annotated <- read.table(grand.annotation.dir, sep = "\t")
  colnames(annotated) <- c(colnames(GOK8873A1), "GENCODE32", "Full.Brain", "Rizz.Brain", "Song.Cult", "Trev.HCS", "TAD", "PE")
  
  # convert to a binary overlap, i.e. if overlap > 0.00 is TRUE
  annotated[,11:17] <- annotated[,11:17] > 0
  
## Annotate with the previously calculated cpm
  annotated$cpm <- (annotated$pileup * 10^6) / sum(annotated$pileup)
  annotated$high.cpm <- paste0(annotated$chr, annotated$start, annotated$end) %in% high.expression
  
## Annotate with the number of replications in four datasets
  annotated$nRep <- apply(annotated[,12:15], 1, sum)
  
## Add an informative name to column 4
  colnames(annotated)[4] <- "id"
  annotated$id <- paste0(annotated$chr, ":", annotated$start, "-", annotated$end)
    

## Save!
  # final list
  final.peaks <- annotated[which(annotated$GENCODE32 == FALSE & 
                                   annotated$TAD &
                                   annotated$nRep >= 2 & 
                                   annotated$PE & 
                                   annotated$high.cpm),]
  
  ## All annotations...
    write.table(annotated, file = "Results/Peaks_Annotated.bed",
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    
    
    write.csv(annotated, file = "Results/Peaks_Annotated.csv", quote = FALSE, row.names = FALSE)
    
  ## Filtered peaks...
    write.table(final.peaks[,1:4], file = "Results/Final_List/NHA_Peaks.bed",
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    write.csv(final.peaks, file = "Results/Final_List/NHA_Peaks_Annotated.csv", quote = FALSE, row.names = FALSE)
    
    
  ## Liftover
    ## Liftover to hg19
      lift.bed(chain = chainFile_38to19, bed.in = "Results/Final_List/NHA_Peaks.bed", bed.out = "Results/Final_List/NHA_Peaks_hg19.bed",
               isList = FALSE)
      
      
      # fix the liftover, in which one enhancer (chr1:17531728-17532439) is split into two rows, running from X-Y and Y-Z, rather than X-Z in one row
      x <- read.delim("Results/Final_List/NHA_Peaks_hg19.bed", header =  FALSE)
      colnames(x) <- c("chr", "start", "end", "id", "score", "strand") # just dummies
      split.peak <- x$id[duplicated(x$id)]
      split.peak <- grep(split.peak, x$id)
      x$start[split.peak] <- min(x$start[split.peak])
      x$end[split.peak] <- max(x$end[split.peak])
      x <- unique(x) # in effect, this code sets each side of the peak to be its full width, then removes the duplicated row this creates
      write.bed(x, "Results/Final_List/NHA_Peaks_hg19.bed")
      
    ## And to mouse
      lift.bed(chain = chainFile_38toMm39, bed.in = "Results/Final_List/NHA_Peaks.bed", bed.out = "Results/Final_List/NHA_Peaks_mm39.bed",
               isList = FALSE)
      
      lift.bed(chain = chainFile_38toMm10, bed.in = "Results/Final_List/NHA_Peaks.bed", bed.out = "Results/Final_List/NHA_Peaks_mm10.bed",
               isList = FALSE)
      
      # and concatenate the split enhancers
      combine.lifted <- function(x) {
        
        y <- list()
        
        for (j in x$V4) {
          z <- x[which(x$V4 == j),]
          
          if (nrow(z) == 1) {
            y[[j]] <- z
            next
          } 
          
          if (length(unique(z$V1)) == 1) { 
            z$V2 <- min(z$V2)
            z$V3 <- max(z$V3)
            z <- unique(z)
            
            y[[j]] <- z
          } else { # if all peaks not on same chromosome, skip it
            next 
          }
          
          
        }
        
        y <- do.call("rbind", y)
        
        return(y)
        
      }
      
      x <- read.bed("Results/Final_List/NHA_Peaks_mm39.bed") %>% combine.lifted()
      write.bed(x, "Results/Final_List/NHA_Peaks_mm39.bed") 
      
      y <- read.bed("Results/Final_List/NHA_Peaks_mm10.bed") %>% combine.lifted()
      write.bed(y, "Results/Final_List/NHA_Peaks_mm10.bed") 
      
  
    
################################################################################################################################ #
## Getting DNA sequences ----
    
## This step is necessary for sgRNA design, being used as input for the webtool
      
## Split peaks into segments!
  ## Get short peaks (< 200 bp)
    short <- 200
    # final.peaks.short <- final.peaks[which((final.peaks$end - final.peaks$start) <= short),]
    
  ## Expand these to 200bp
    peak.expander <- function(x, short = short) {
      # restrict to peaks < short
      x <- x[which((x$end - x$start) <= short),]
      
      # get the midpoint of each peak
      widths <- x$end - x$start
      midpoints <- round(x$start + (widths/2))
      
      # expand the peak to be of width short centred around the midpoint
      x$start <- midpoints - (short / 2)
      x$end <- midpoints + (short / 2)
      
      # rename the peak
      x$id <- paste0(x$id, "-exp-", x$chr, ":", x$start, "-", x$end)
      
      return(x)
    }
    
    final.peaks.short <- peak.expander(final.peaks, 200)
    
  
  ## Get long peaks, and segment them!
    # function
    peak.segmenter <- function(x = final.peaks, short = short) {
      x <- x[which((x$end - x$start) > short),] # get just the long peaks (201 or greater)
      
      output <- list()
      for (j in 1:nrow(x)) {
        y <- x[j,] # storing the row of interest in a new variable
        
        # determine the number of divisions
        l <- y$end - y$start 
        n.segs <- ceiling(l / short) + 1 
        seg.length <- l / n.segs
        
        # split the existing peak into n.segs peaks or seg.length
        z <- vector("list", length = n.segs)
        
        for(k in 1:n.segs) {
          z[[k]] <- y
          z[[k]]$start <- y$start + ((k-1)* seg.length)
          z[[k]]$end <- z[[k]]$start + seg.length
          
          z[[k]]$start <- round(z[[k]]$start)
          z[[k]]$end <- round(z[[k]]$end)
          
          z[[k]]$id <- paste0(z[[k]]$id, "-seg", k, "-", z[[k]]$chr, ":", z[[k]]$start, "-", z[[k]]$end)
        }
        
        z <- do.call("rbind", z)
        # z$id <- paste0(z$chr, ":", z$start, "-", z$end)
        
        output[[j]] <- z
      }
      output <- do.call("rbind", output)
      return(output)
    }
    
    # run function
    # final.peaks.long <- peak.segmenter(x = final.peaks, short = short)
    final.peaks.long <- final.peaks[which((final.peaks$end - final.peaks$start) > short),]
    
    
    
  ## Save
  write.table(final.peaks.short, file = "Results/FinalPeaksShort.bed",
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  write.table(final.peaks.long, file = "Results/FinalPeaksLong.bed",
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  
    
    
    
    
## Get fasta sequences
  hg38 <- "/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/Genome/UCSC/hg38.fa"
  candidate.dna.dir <- "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/CandidateDNA_"
  
  call.getFasta.short <- paste0("bedtools getfasta", 
                          " -name", # use the name field and coordinates for the FASTA header
                          " -fi ", hg38, # fasta of the genome
                          " -bed /mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/FinalPeaksShort.bed", # the bed file of the final peaks
                          " -fo ", candidate.dna.dir, "short.fa")
  
  
  system(call.getFasta.short, intern = FALSE, wait = TRUE)
  
  call.getFasta.long <- paste0("bedtools getfasta", # bed annotate
                          " -name", # use the name field and coordinates for the FASTA header
                          " -fi ", hg38, # fasta of the genome
                          " -bed /mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/FinalPeaksLong.bed", # the bed file of the final peaks
                          " -fo ", candidate.dna.dir, "long.fa")
  
  
  system(call.getFasta.long, intern = FALSE, wait = TRUE)
  
      
  
################################################################################################################################ #
## Validation of peaks using external datasets ----
  
  
## Here, we're looking at two external datasets (Nott et al 2019 and ENCODE) to validate our candidate selection as being enhancers / in open chromatin
  
## Process input files
  
  ## ENCODE
    encode.file <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/ENCODE/ENCODE_CREs_050422_ucscDownload"
  
    # read in 
    encode <- read.delim(encode.file)
  
    # partition by ucscLabel
    table(encode$ucscLabel)
    for(j in unique(encode$ucscLabel)) {
      x <- encode[which(encode$ucscLabel == j),]
      write.bed(x[,1:4], dir = paste0(encode.file, ".", j, ".bed"))
    }
    
  ## Nott
    nott.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Nott2019/"
    nott.files <- list.files(nott.dir)[-c(1:2)]
    nott.files <- nott.files[-length(nott.files)]
    nott.files <- nott.files[-grep("sorted", nott.files)]
    # nott.files <- paste0(nott.dir, nott.files)
    
    # convert to hg38
    for (j in nott.files) {
      
      lift.bed(chainFile = chainFile_19to38,
               bed.in = paste0(nott.dir, j), 
               bed.out = paste0(nott.dir, "hg38/", j),
              isList = FALSE)
    }

## Run
  # directories
  nott.in <- paste0(nott.dir, "hg38/", nott.files)  
  encode.in <- paste0(encode.file, ".", unique(encode$ucscLabel), ".bed")
  validation.out <- "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Candidate_Validation/BED_Annotations_ValidationSets.bed"
  validation.in <- paste0(c(nott.in, encode.in), collapse = " ")
  
  # run  
  call <- paste("annotateBed",
                "-i", inhouse.dir,
                "-files", validation.in,
                ">", validation.out)

  system(call, intern = FALSE, wait = TRUE) 
  
## Read in and wrangle
  # read
  vld <- read.table(validation.out, sep = "\t")
  vld <- vld[,-c(4:9)]
  
  # rename columns
  annot.cols <- 5:21
  colnames(vld)[-annot.cols] <- c("chr", "start", "end", "peak")
  
  x <- gsub("/mnt/Data0/PROJECTS/CROPSeq/PublicData/", "", validation.in) %>%
    gsub(".bed", "", .) %>%
    gsub("Nott2019/hg38/", "Nott_", .) %>%
    gsub("/PROCESSED/liftOver_hg38/Astrocytes.hg38.super", "", .) %>%
    gsub("/ENCODE_CREs_050422_ucscDownload.", "_", .) %>%
    gsub("optimal_peak.", "", .) %>%
    gsub("IDR_ENCODE.", "", .) %>%
    gsub("Cell2013_", "", .)
    
  x <- strsplit(x, " ")[[1]]
  
  colnames(vld)[annot.cols] <- x
  
  # binarise
  vld[annot.cols] <- vld[annot.cols] > 0
  
  # add annotation for candidate and hit peaks
  vld$peak <- paste0(vld$chr, ":", vld$start, "-", vld$end)
  vld$Candidate <- vld$peak %in% guides$TargetCoord[which(guides$TargetCat == "Enh")]
  
  # save
  write.csv(vld, file = "Results/Candidate_Validation/Annotations.csv")
  
## Analyse
  x <- apply(vld[,annot.cols], 2, table) %>% t() %>% as.data.frame
  colnames(x) <- c("NoOverlap", "Overlap")

  x$OverlapInCandidates <- apply(vld[which(vld$Candidate),annot.cols], 2, function(y) length(which(y)))

  x$FisherCanditates <- apply(vld[,annot.cols], 2, function(y) fisher.test(table(y, vld$Candidate))$estimate)
  x$FisherCanditatesP <- apply(vld[,annot.cols], 2, function(y) fisher.test(table(y, vld$Candidate))$p.value)
  
  write.csv(x, file = "Results/Candidate_Validation/Statistics.csv")
  
  