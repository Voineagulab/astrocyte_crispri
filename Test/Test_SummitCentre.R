## Read in MAC2 data
macs <- read.delim("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/Bulk_ATAC_RNAseq/GOK8505-GOK8837/GOK8837A1/Mapping/MACS2/NHA_ATAC_S3.filtered.BAM_peaks_reformatedGJS.bed.txt", header = FALSE)
macs$id <- paste(macs$V1, macs$V2, macs$V3, sep = "_")

## Restrict to our peaks
screened <- read.delim("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed", header = FALSE)
macs <- macs[which(macs$id %in% screened$V4),c(1,2,3,5,11)]
colnames(macs) <- c("chr", "start", "end", "summit", "id")

## Function to create window of size x around peak peak (i.e. absolute summit of peak)
find.peak.centre.window <- function(data = macs, window = 150) {
  # get summits
  summits <- data$summit
  
  # draw window
  dist <- round(window / 2)
  centred.window <- data.frame(chr = data$chr,
                               start = data$summit - dist,
                               end = data$summit + dist,
                               id = data$id)
  
  # output
  return(centred.window)
}

w <- 128
x <- find.peak.centre.window(window = w)
centres.dir <- "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Peak_Centres_259bp.bed"
write.bed(x, dir = centres.dir)

## Map to SNPs
  db.dir <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/dbSNP153_280422.bed"
  centre.out <- paste0("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Peak_Centres_", w, "bp_SNPs.bed")
  
## Run
  call <- paste("windowBed",
                "-a", centres.dir,
                "-b", db.dir,
                "-w", "0", # window of 1000bp
                ">", centre.out)
  
  system(call, intern = FALSE, wait = TRUE) 

## Read
  centre.snps <- read.delim(centre.out, header = FALSE)
  colnames(centre.snps) <- c("SummitChr", "SummitStart", "SummitEnd", "EnhCoord", "snpChr", "snpStart", "snpEnd", "snpID")
  
## Match to Enh name
  g <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  g$TargetCoord <- gsub(":", "_", g$TargetCoord)
  g$TargetCoord <- gsub("-", "_", g$TargetCoord)
  m <- match(centre.snps$V4, g$TargetCoord)  
  centre.snps$EnhId <- g$TargetID[m]

## Save
  write.csv(centre.snps, file = paste0("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Peak_Centres_", w, "bp_SNPs.csv"))
  