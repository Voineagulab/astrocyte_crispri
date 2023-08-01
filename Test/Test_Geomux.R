
## IN

  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
  load(file = "../1_Processing/GuideAssignment/GuideAssignments.rda")
  h <- lapply(g, function(x) {
    # filter rows
    x <- x[which(x$umi_count >= 3 & x$cr_call & !(is.na(x$GuideID))),]
    
    # filter columns
    x <- x[,c("cell", "GuideID", "umi_count")]
    
    # return
    return(x)
  })
  
  for (j in names(h)) write.table(h[[j]], file = paste0("Geomux_Input_", j, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  
  
  
  ## In shell
  # cd /mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/
  # conda activate geomux
  

  geomux -i Geomux_Input_NHA_1.txt -o Geomux_Output_NHA_1.txt -u 3
  geomux -i Geomux_Input_NHA_2.txt -o Geomux_Output_NHA_2.txt -u 3
  geomux -i Geomux_Input_NHA_3.txt -o Geomux_Output_NHA_3.txt -u 3
  geomux -i Geomux_Input_NHA_4.txt -o Geomux_Output_NHA_4.txt -u 3
  geomux -i Geomux_Input_NHA_5.txt -o Geomux_Output_NHA_5.txt -u 3
  geomux -i Geomux_Input_NHA_7.txt -o Geomux_Output_NHA_7.txt -u 3
  geomux -i Geomux_Input_NHA_8.txt -o Geomux_Output_NHA_8.txt -u 3
  
  
  a <- read.table("Geomux_Output_NHA_1.txt", sep = "\t", header = TRUE)
  