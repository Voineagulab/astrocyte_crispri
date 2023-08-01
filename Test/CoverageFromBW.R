## Setup
  # bed file with regions to quantify coverage in
  imported_bed <- import.bed(nha_dir_19)
  
  # bigwig fileliest
  bw_filelist <- list.files("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig/") # this is hg19
  bw_filelist <- bw_filelist[grep("bigwig", bw_filelist)]
  bw_filelist <- paste0("../../../../PublicData/snRNAseq/Herring2022_MaturationBrain/ATAC/BigWig/", bw_herring)
    

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