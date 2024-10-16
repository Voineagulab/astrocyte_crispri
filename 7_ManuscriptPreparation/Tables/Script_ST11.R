setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable1_ATAC/")
options(stringsAsFactors = FALSE)
source("../../../FullScale/Scripts/Functions.R")


## 11C: scRNAseq statistics
  ## Read in
    samps <- paste0("NHA_", c(1:5,7:8)) 
    qc <- lapply(samps, function(x) {
      x <- read.csv(paste0("../../../PreprocessedData/FullScale/FinalCount/", x, "/outs/metrics_summary.csv"))
    })
  
    qc <- do.call("rbind", qc)
    
    rownames(qc) <- samps
    qc <- t(qc)
    qc <- as.data.frame(qc)
    
  ## Remove stats on cells, as these are derived from the default CellRanger caller rather than EmptyDrops
    remRows <- c(1,2,3,17,19)
    qc <- qc[-remRows,]
  
  ## Add information on the heminested PCR
    # shell code to collect sequencing depth
    # necessary as these arenot processed using normal pipelines
    libSize_hnPCR <- list()

    for (j in samps) {
      print(j)
      
      # get paths for FASTQs
      wd <- paste0("/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/FullScale/GOK10211_hnPCR_2/FASTQ/", j)
      f <- list.files(wd, full.names = TRUE)
      # f <- f[grep("_R1_", f)] # this gets R1 only
      
      # count number of lines
      wc <- sapply(f, function(x) {
        call <- paste("wc -l", x)
        system(call, intern = TRUE, wait = TRUE)
      })
      
      # sum across lanes
      libSize_hnPCR[[j]] <- splitter(wc, " ", 1) %>% as.numeric() %>% sum()
      
    }

      libSize_hnPCR <- do.call("c", libSize_hnPCR)
      
      qc["Heminested.PCR.Reads",] <- libSize_hnPCR
      
  ## Save
    write.csv(qc, file = "11C_scRNAseq.csv", row.names = TRUE)
    
    
 
    