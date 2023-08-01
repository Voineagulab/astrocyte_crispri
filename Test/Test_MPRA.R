################################################################################################################################ #
## Hansen ----

## Directories
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
  nha.all <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/Bulk_ATAC_RNAseq/RESULTS/ATACseq/NHA.allPeaks.bed"
  hansen.activator <-"/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hansen2022_ATACSTARR/GSE181317_GM12878_active_regions.bed"
  hansen.silencer <-"/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hansen2022_ATACSTARR/GSE181317_GM12878_silent_regions.bed"
  hansen.neither <-"/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hansen2022_ATACSTARR/GSE181317_GM12878_accessible-peaks_genrich.narrowPeak"
  hansen.out <- "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Hansen2022_NHA.bed"
  
## Categories of screened hits
  enh.activator <- res$Enh.Pos[which(res$Hit & (res$Z < 1))]
  enh.silencer <- res$Enh.Pos[which(res$Hit & (res$Z > 1))]
  enh.neutral <- res$Enh.Pos[-which(res$Enh.Pos %in% enh.activator | res$Enh.Pos %in% enh.silencer)] %>% unique()

## Overlap
  call <- paste("intersectBed",
                  "-a", nha.all,
                  "-b", hansen.activator, hansen.silencer,
                  "-loj", # Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B
                  ">", hansen.out)
  
  system(call, intern = FALSE, wait = TRUE)
    
## Read in for screen
  han <- read.delim(hansen.out, header = FALSE)
  han$ID <- paste0(han$V1, ":", han$V2, "-", han$V3)
  han$Hit <- "No"
  han$Hit[which(han$ID %in% enh.activator)] <- "Screen Enhancer"
  han$Hit[which(han$ID %in% enh.silencer)] <- "Screen Silencer"
  han$Hit[which(han$ID %in% enh.neutral)] <- "Screen Neutral"
  
  a <- unique(han$ID[which(han$V7 == 1)])
  s <- unique(han$ID[which(han$V7 == 2)])
  n <- unique(han$ID[which(han$V7 == 0)])
  
  table(enh.activator %in% a)
  table(enh.activator %in% s)
  
  table(enh.silencer %in% a)
  table(enh.silencer %in% s)

## Now look across all NHA peaks
  hansen.all  <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Hansen2022_ATACSTARR/GSE181317_GM12878_accessible-peaks_genrich.narrowPeak"
  hansen.all.out <- "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Hansen2022_NHA_all.bed"
  
    call <- paste("intersectBed",
                  "-a", nha.all,
                  "-b", hansen.all,
                  "-loj", # Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B
                  ">", hansen.all.out)
  
  system(call, intern = FALSE, wait = TRUE)
  
  tested <- read.delim(hansen.all.out, header = FALSE)
  tested <- tested[which(tested$V8 != "-1"),]
  tested <- paste0(tested$V1, ":", tested$V2, "-", tested$V3)
  
  x <- data.frame(Peak = unique(han$ID))
  x$MPRA <- x$Peak %in% tested
  x$MPRA.Effect <- "MPRA Absent"
  x$MPRA.Effect[which(x$MPRA)] <- "MPRA ns"
  x$MPRA.Effect[which(x$Peak %in% a)] <- "MPRA Activator"
  x$MPRA.Effect[which(x$Peak %in% s)] <- "MPRA Silencer"
  x$MPRA.Effect[which(x$Peak %in% s & x$Peak %in% a)] <- "MPRA A+S"
  
  x$Screen <- "Screen Absent"
  x$Screen[which(x$Peak %in% enh.activator)] <- "Screen Activator"
  x$Screen[which(x$Peak %in% enh.silencer)] <- "Screen Silencer"
  x$Screen[which(x$Peak %in% enh.neutral)] <- "Screen ns"
  
  x <- x[which(x$MPRA),]
  table(x$MPRA.Effect, x$Screen)
  
  # x$
  # x <- unique(han$ID)
  # names(x) <- x
  # x[which((x %in% a) & !(x %in% s))] <- "A"
  # x[which(!(x %in% a) & (x %in% s))] <- "S"
  # x[which((x %in% a) & (x %in% s))] <- "AS"
  # 
  # x <- data.frame(Peak = names(x),
  #                 Cat = x)
  # x$Screen <- "Not Screened"
  # x$Screen[which(x$Peak %in% enh.activator)] <- "Screened Activator"
  # x$Screen[which(x$Peak %in% enh.silencer)] <- "Screened Silencer"
  # x$Screen[which(x$Peak %in% enh.neutral)] <- "Screened Neutral"
  # 
  # 

## Check Manolis Kellis
hidra.in <- "Hidra.txt"
nha.all.19 <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/Bulk_ATAC_RNAseq/RESULTS/ATACseq/NHA.allPeaks.hg19.bed"
hidra.out <- "Hidra_out.txt"

call <- paste("intersectBed",
                "-a", nha.all.19,
                "-b", hidra.in,
                "-c", 
                ">", hidra.out)

system(call, intern = FALSE, wait = TRUE)
hidra <- read.delim(hidra.out, header = FALSE)
m <- match(hidra$V4, han$V4)
hidra$Hit <- han$Hit[m]

bg <- read.delim("GSE104001_HiDRA_counts_per_fragmentgroup.txt")
bg2 <- data.frame(chr = splitter(bg$FragmentGroupPosition_NumUniqFragments, ":", 1),
                  start = splitter(splitter(bg$FragmentGroupPosition_NumUniqFragments, "-", 1), ":", 2),
                  end = splitter(bg$FragmentGroupPosition_NumUniqFragments, "-", 2))
bg2$end <- splitter(bg2$end, "_", 1)
write.bed(bg2, "GSE104001_HiDRA_counts_per_fragmentgroup.bed")

call <- paste("intersectBed",
                "-a", nha.all.19,
                "-b", "GSE104001_HiDRA_counts_per_fragmentgroup.bed",
                "-c", 
                ">", "GSE104001_HiDRA_counts_per_fragmentgroup_intersected.bed")

system(call, intern = FALSE, wait = TRUE)

bg3 <- read.delim("GSE104001_HiDRA_counts_per_fragmentgroup_intersected.bed", header = FALSE)
