setwd("/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/TTseq/NFcore_RNAseq_Basic/")
samps <- read.csv("../TT-seq_2022_NHA_SH-SY5Y.csv")
keep.samps <- samps$LIMSID[grep("NHA_TT", samps$Sample)]

## Read in stringtie transcripts
  # id
  j <- keep.samps[1]
  
  # read in
  x <- read.table(paste0("star_salmon/stringtie/", j, ".transcripts.gtf"), sep = "\t", header = FALSE)
  colnames(x) <- c("Chr", "Database", "Type", "Start", "End", "V6", "Strand", "V8", "AllInfo")
  
  y <- x[which(x$Database == "StringTie"), c("Chr", "Start", "End", "Database", "Type", "AllInfo")]
  y$Chr <- paste0("chr", y$Chr)
  
  # convert to bed
  write.bed(y, "../../../FullScale/Results/Scratchspace/Stringtie.bed")

  # intersect
  wawb("../../../FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed", "../../../FullScale/Results/Scratchspace/Stringtie.bed", "../../../FullScale/Results/Scratchspace/Stringtie_Intersect.bed")
    
  # read in
  z <- read.bed("../../../FullScale/Results/Scratchspace/Stringtie_Intersect.bed")
  