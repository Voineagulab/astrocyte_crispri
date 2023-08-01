## Directories
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
  screened <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks_hg19lift.bed"
  screened <- read.delim(screened, header = FALSE)
  screened$ID <- paste(screened$V1, screened$V2, screened$V3, sep = "_")
  nha.all <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/Bulk_ATAC_RNAseq/RESULTS/ATACseq/NHA.allPeaks.hg19.bed"
  allPeaks <- read.delim(nha.all, header = FALSE) 
  allPeaks$ID <- paste(allPeaks$V1, allPeaks$V2, allPeaks$V3, sep = "_")
  

## Process Cooper data
  cooper <- read.csv("../../../PublicData/Cooper2022_science.abi8654_data_s1.csv", skip = 4, header = TRUE)
  cooper.bed <- data.frame(chr = cooper$chr, start = cooper$pos, end = cooper$pos, id = cooper$MPRA_ID)
  write.bed(cooper.bed, dir = "../../../PublicData/Cooper2022.bed")

## Intersect at 0kb window
  call <- paste("windowBed",
                 "-w", "0",
                "-a", nha.all,
                "-b", "../../../PublicData/Cooper2022.bed",
                # "-loj", # Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B
                ">", "Cooper_0kb.bed")
  system(call, intern = FALSE, wait = TRUE) 
  
  window0 <- read.delim("Cooper_0kb.bed", header = FALSE)  
  
## Intersect at 1kb window
  call <- paste("windowBed",
                 "-w", "1000",
                "-a", nha.all,
                "-b", "../../../PublicData/Cooper2022.bed",
                # "-loj", # Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B
                ">", "Cooper_1kb.bed")
  system(call, intern = FALSE, wait = TRUE) 
  
  window1 <- read.delim("Cooper_1kb.bed", header = FALSE)  
  window1$ID <- paste(window1$V1, window1$V2, window1$V3, sep = "_")
  window1$Screen <- window1$ID %in% screened$V4
  