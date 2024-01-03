############
##
##
##Annotation Functions
##
##
###########

#########

#####
#Add Tobias Bound/Unbound enhancer level 
library(biomaRt)

addnumTSSbetween <- function (all_rf_factors, geneinfo) {
  geneinfo$TSS <- sub(".*:","",geneinfo$TSS)
  numTSSEnhGene <- apply(all_rf_factors,1, function(pair) {
    chr_genes <- geneinfo[geneinfo$Chr == pair["Gene.chr"],]
    midpoint <- as.numeric(pair["Enh.Midpoint"])
    TSS <- as.numeric(pair["TSS"])
    if (TSS > midpoint) {
      out <- sum(midpoint <  chr_genes[,"TSS"] & chr_genes[,"TSS"] < TSS )
    } else {
      out <- sum(midpoint >  chr_genes[,"TSS"] & chr_genes[,"TSS"] > TSS)
    }
    return(out)
  })
  return(numTSSEnhGene)
}



addExpTobiasBound<- function(df, exp_bound_overlaps, exp_unbound_overlaps) {
  bound <- read.table(exp_bound_overlaps)
  unbound <- read.table(exp_unbound_overlaps)
  bound_counts <- TFmatrix(bound,outnumeric = TRUE)
  unbound_counts <- TFmatrix(unbound,outnumeric = TRUE)
  #Add total bound/unbound
  ubdf <- data.frame(Enh = unbound_counts$Group.1, Tobias.Exp_Unbound_TF_counts = rowSums(unbound_counts[,2:ncol(unbound_counts)]))
  bdf <- data.frame(Enh = bound_counts$Group.1, Tobias.Exp_Bound_TF_counts =  rowSums(bound_counts[,2:ncol(bound_counts)]))
  colnames(unbound_counts) <- paste0("Tobias.Unbound_", colnames(unbound_counts))
  colnames(bound_counts) <- paste0("Tobias.Bound_", colnames(bound_counts))
  df <- merge(df, bdf, all.x = T)
  df <- merge(df, ubdf, all.x = T)
  df <- merge(df, bound_counts,by.x = "Enh", by.y = "Tobias.Bound_Group.1", all.x = T) #bound
  df <- merge(df, unbound_counts, by.x = "Enh", by.y = "Tobias.Unbound_Group.1", all.x = T) #Gasp_unbound
  df[is.na(df)] <- 0
  return(df)
}




##
addGeneStability <- function(EGPs, geneinfo, ensembl) {
  housekeeping_genes <- read.table("../PublicData/HousekeepingGenes/Eisenberg2013_Genes.txt")
  housekeeping_genes <- data.frame(Gene = housekeeping_genes$V1)
  syns <- getBM(attributes = c("external_synonym", "hgnc_symbol") , filters = c("external_synonym"),
                values = housekeeping_genes$Gene, mart = ensembl) 
  colnames(syns) <- c("Gene", "hgnc_symbol") #check if synoymn is there
  housekeeping_genes <- merge(housekeeping_genes, syns, all.x = T)
  housekeeping_genes <- merge(housekeeping_genes, geneinfo[,c("Gene", "EnsID")], all.x = T)
  EGPs$Gene.Housekeeping <- EGPs$Gene %in% housekeeping_genes$Gene | toupper(EGPs$Gene) %in% housekeeping_genes$hgnc_symbol | EGPs$EnsID %in% housekeeping_genes$EnsID
  
  #This section is quite convoluted 
  #because there are so many different gtfs being used by the ENCODE paper
  lin_hk <- read_excel("../PublicData/HousekeepingGenes/Lin2019_SingleCell.xlsx")
  lin_hk <- as.data.frame(lin_hk[c("GeneSymbol","Stability index")])
  syns <- getBM(attributes = c("external_synonym", "hgnc_symbol") , filters = c("external_synonym"),
                values = lin_hk$GeneSymbol, mart = ensembl) 
  colnames(syns) <- c("Gene", "hgnc_symbol")
  colnames(lin_hk) <- c("Gene", "Gene.StabilityIndex")
  lin_hk <- merge(lin_hk, syns, all.x = T)
  lin_hk <- merge(lin_hk, geneinfo[,c("Gene", "EnsID")], all.x = T)
  #merge on Gene or hgnc if gene stability not there
  
  EGPs <- merge(EGPs, unique(lin_hk[,c("Gene", "Gene.StabilityIndex")]), all.x = T)
  
  lin_hk_hgnc <- lin_hk[match(unique(lin_hk$hgnc_symbol),lin_hk$hgnc_symbol),c("hgnc_symbol", "Gene.StabilityIndex")]
  NAcols<- merge(EGPs[is.na(EGPs$Gene.StabilityIndex),colnames(EGPs) != "Gene.StabilityIndex"], 
           unique(lin_hk_hgnc[,c("hgnc_symbol", "Gene.StabilityIndex")]),
           by.x = "Gene", by.y = "hgnc_symbol", all.x = T)
  EGPs <- rbind(EGPs[!is.na(EGPs$Gene.StabilityIndex),], NAcols)
  
  #None match on EnsID but will leave in just in case
  lin_hk <- lin_hk[match(unique(lin_hk$EnsID),lin_hk$EnsID),c("EnsID", "Gene.StabilityIndex")]
  NAcols<- merge(EGPs[is.na(EGPs$Gene.StabilityIndex),colnames(EGPs) != "Gene.StabilityIndex"], 
                 unique(lin_hk[,c("EnsID", "Gene.StabilityIndex")]),
                 all.x = T)
  EGPs <- rbind(EGPs[!is.na(EGPs$Gene.StabilityIndex),], NAcols)
  #EGPs[is.na(EGPs$Gene.StabilityIndex),"Gene.StabilityIndex"] <- 0
  return(EGPs)
}

addFeatureCounts <- function (EGPs, feature.file, try.genecol = F) {
  Exp <- read.csv(feature.file)
  if (ncol(Exp) != 5){
    warning("Not in expected format feature not added")
    return(EGPs)
  }
  colnames(Exp) <- c("Gene", "EnsID",paste0("Gene.",colnames(Exp)[3:5]))
  unique(EGPs[! EGPs$EnsID %in%  Exp$EnsID,c("EnsID")])
  Exp_Gene <- unique(Exp[,c("Gene",colnames(Exp)[3:5])]) #replace with match?
  
  Exp_EnsID <- unique(Exp[,c("EnsID",colnames(Exp)[3:5])])
  #Try merge on ENSID first
  EGPs <- merge(EGPs, Exp_EnsID, all.x = T)
  if (try.genecol){
    missing_exp <- is.na(EGPs[,colnames(EGPs) == colnames(Exp_EnsID)[2]])
    EGPs_missing <- EGPs[missing_exp,! colnames(EGPs) %in% colnames(Exp_EnsID)[2:4]]
    EGPs <- EGPs[!missing_exp,]
    EGPs_missing <- merge(EGPs_missing, Exp_Gene, all.x = T)
    EGPs <- rbind(EGPs, EGPs_missing)
  }
  return(EGPs)
}

getGeneNearest <- function (EGPs, geneInfo) { #
  nearest_Gene <- list()
  nearest_EnsID <- list()
  for (enh  in unique(EGPs$Enh)) {
    # enhancer coordinates
    enh.mid <- unique(EGPs[EGPs$Enh == enh,"Enh.midpoint"])
    enh.chr <- paste0("chr",unique(EGPs[EGPs$Enh == enh,"Enh.chr"]))
    chr_genes <- geneInfo[geneInfo$Chr == enh.chr,]
    nearest_Gene[[enh]] <- chr_genes[which.min(abs(enh.mid - chr_genes$TSS)),"Gene"] #gene.col
    nearest_EnsID[[enh]] <- chr_genes[which.min(abs(enh.mid - chr_genes$TSS)), "EnsID"]
  }
  nearest_EnsID <- do.call("c", nearest_EnsID)
  nearest_Gene <- do.call("c", nearest_Gene)
  EGPs$Gene.Nearest <- FALSE
  for (i in 1:nrow(EGPs)) {
    EGPs$Gene.Nearest[i] <- EGPs[,"EnsID"][i] == nearest_EnsID[names(nearest_EnsID) == EGPs$Enh[i]] | EGPs[,"Gene"][i] %in% nearest_Gene[names(nearest_Gene) == EGPs$Enh[i]]
  }
  return(EGPs)
}
#################
## Setup
# bed file with regions to quantify coverage in
getTTfromBigWig <- function (bed , bwloc ) {
  imported_bed <- import.bed(bed)
  #TODO Some regions are split up should handle this
  chain <- import.chain("../PublicData/hg19ToHg38.over.chain")
  imported_bed <- liftOver(imported_bed, chain = chain)
  #duplicates <- unlist(lapply(imported_bed, function(x) length(x) > 1))
  #imported_bed[duplicates,] <- lapply(imported_bed[duplicates,],)
  imported_bed <- unlist(imported_bed)
  for (name in unique(imported_bed$name)) { #chang
    start <- start(imported_bed[imported_bed$name == name])
    end <- end(imported_bed[imported_bed$name == name])
    start(imported_bed[imported_bed$name == name]) <- min(start)
    end(imported_bed[imported_bed$name == name]) <- max(end)
  }
  imported_bed <- imported_bed[match(unique(imported_bed$name), imported_bed$name)]
  # bigwig fileliest
  bw_filelist <- list.files(bwloc, pattern = ".bw", full.names = T) # this is hg19
  ## Calculate coverage in each of the bigwigs
  coverage <- lapply(bw_filelist, function(bigWigString) { # requires ~20s
    print(bigWigString)  
    # get coverage
    y <- import(bigWigString, selection = imported_bed, as = "NumericList") 
    # get a summary of the coverage for each region in the bed
    z <- sapply(y, sum) # this is the sum of coverages across the bases, and is not normalised to peak width. this makes comparisons between peaks difficult
    # collect peak ids
    names(z) <- imported_bed$name
    return(z)
  })
  TT_minus <- data.frame(Enh = names(coverage[[1]]), TTseq_minus = coverage[[1]])
  TT_plus <- data.frame(Enh = names(coverage[[2]]), TTseq_plus = coverage[[2]])
  TT_all <- merge(TT_minus, TT_plus)
  TT_all$TTseq.TTseq_Total <- TT_all$TTseq_minus + TT_all$TTseq_plus
  return(TT_all)
}

#
makeWindowBed <- function (EGPs, midpoint.col = "Enh.midpoint",window = 500) {
  bed<-data.frame(Enh.chr = EGPs$Enh.chr,
                  Enh.start = EGPs[,midpoint.col] - window,
                  Enh.end = EGPs[,midpoint.col] + window,
                  Enh = EGPs$Enh )
  return(bed)
}

#.tab file format
#name - name field from bed, which should be unique
#size - size of bed (sum of exon sizes 
#covered - # bases within exons covered by bigWig
#                      sum - sum of values over all bases covered
#                    mean0 - average over bases with non-covered bases counting as zeroes
#                    mean
addChipBigWigs <- function (df, filenames, folder) {
 for (file in filenames) {
  averages <- read.table(paste0(folder, file))
  name=sub(".*(ENCFF.*).tab", "\\1", file)
  name <- accession2Name(name)
  averages <- averages[,c(1,6)]
  colnames(averages) <- c("Enh",name)
  df <- merge(df, averages, by = "Enh") #only getting the mean atm
 }
 return(df)
}

addTobiasBound <- function (df, bound_file, unbound_file, merge_col = "Enh.Pos", expressed = T, exp.file = "Results/Tobias/ExpressedTFs/Expressed_TFs.csv") {
  bound_overlaps <- read.table(bound_file)
  unbound_overlaps <- read.table(unbound_file)
  unbound_matrix <- TFmatrix(unbound_overlaps, outnumeric = T)
  bound_matrix <- TFmatrix(bound_overlaps, outnumeric = T)
  
  if (expressed == T) {
    homo_TFs <- read.csv(exp.file)
    bound_matrix <- bound_matrix[,colnames(bound_matrix) %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names | colnames(bound_matrix) == "Group.1"]
    unbound_matrix <- unbound_matrix[,colnames(unbound_matrix) %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names | colnames(unbound_matrix) == "Group.1"]
  }
  colnames(unbound_matrix) <- paste0("Tobias.Unbound_", colnames(unbound_matrix))
  colnames(bound_matrix) <- paste0("Tobias.Bound_", colnames(bound_matrix))
  #Remove bound matrix
  #all_rf_factors <- all_rf_factors[, ! colnames(all_rf_factors) %in% colnames(all_rf_factors)[which(colnames(all_rf_factors) == "Tobias.Bound_ALX3"):which(colnames(all_rf_factors) =="Tobias.Unbound_ZSCAN4")]]
  df <- merge(df, bound_matrix,by.x = merge_col, by.y = "Tobias.Bound_Group.1", all.x = T)
  df <- merge(df, unbound_matrix, by.x = merge_col, by.y = "Tobias.Unbound_Group.1", all.x = T)
  #The only columns with NA value are Tobias values which should all be 0
  colnames(df)[which(apply(df, MARGIN = 2, function (x) {sum(is.na(x))}) > 0)]
  df[is.na(df)] <- 0 
  if (expressed == T) {
    df$Tobias.Exp_Bound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Bound")])
    df$Tobias.Exp_Unbound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Unbound")])
  } else {
    df$Tobias.Bound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Bound")])
    df$Tobias.Unbound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Unbound")])
  }
  return(df)
}

sameTADs <- function (all_rf_factors ,tads) {
  same_tad <- apply(tads[tads$chr %in% unique(all_rf_factors$Enh.chr),], 1, function (tad) {
    chr_EGPs <- all_rf_factors[all_rf_factors$Enh.chr == tad["chr"],]
    tad_res <- apply(chr_EGPs, 1, function(EGP) {
      if (as.numeric(EGP["Enh.start"]) <= as.numeric(tad["end"]) & as.numeric(EGP["Enh.end"]) >= as.numeric(tad["start"])) {
        if (as.numeric(EGP["TSS"]) > as.numeric(tad["start"]) & as.numeric(EGP["TSS"]) < as.numeric(tad["end"])) {
          return(EGP['Pair'])
        } 
      }
      return("")
    })
    return(tad_res)
  }) 
  same_tad <- as.data.frame(table(unname(unlist(same_tad))))
  colnames(same_tad) <- c("Pair", "EGP.Shared_TADs")
  return(same_tad)
}

inTAD <- function  (id, chr,start, end, tads) {
  in_tad <- apply(tads[tads$chr %in% unique(chr),], 1, function (tad) {
    chr_EGPs <- which(chr == tad["chr"])
    subset_ids <- id[chr_EGPs]
    return(subset_ids[start[chr_EGPs] <= as.numeric(tad["end"]) & end[chr_EGPs] >= as.numeric(tad["start"])])
  })
  return(unique(unlist(in_tad)))
}
NarrowPeaksHeader <- c("chr", "start", "end", "name", "Peakscore", "strand", "signalValue", "pValue", "qValue", "peak")


