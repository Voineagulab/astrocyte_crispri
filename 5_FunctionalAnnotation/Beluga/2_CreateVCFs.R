################
##
##
## Create VCFs for all possible sequence variants at each nucleotide position for the 145 functional enhancers
## @author: Sam Bagot
## @date: 24-10-23
################

#create .vcf for beluga
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/")
source("Scripts/Header_functions.R")
library("gtools")

enhancers <- read.csv("Data/Enhancer_factors.csv")
fasta <- read.table("Data/PeakData/PeaksFasta.tab")
colnames(fasta) <- c("Enh", "fasta")
fasta$fasta <- toupper(fasta$fasta)
#Process fasta &#Create VCF
#results <- data.frame(Enh = sub(">","",fasta$V1[1:nrow(fasta) %% 2 == 1 ]), fasta = fasta$V1[1:nrow(fasta) %% 2 == 0 ])
results <- merge(results, enhancers[,c("Enh","HitPermissive")])
results <- results[results$HitPermissive > 0,]
split <- str_split(results$fasta, "")
names(split) <- paste0(results$Enh, "_")
split<- unlist(split)

vcf <- data.frame(Enh = sub("_.*", "", names(split)), Enh.Pos = as.numeric(sub(".*_", "", names(split))) - 1, ref = split)
vcf <- merge(vcf, enhancers[,c("Enh", "Enh.chr", "Enh.start")])
vcf$Pos <- vcf$Enh.Pos + vcf$Enh.start
vcfs <- list(vcf, vcf, vcf, vcf)
names(vcfs) <- c("G", "C", "A" , "T")
vcfs[["G"]]$alt <- "G"
vcfs[["C"]]$alt <- "C"
vcfs[["A"]]$alt <- "A"
vcfs[["T"]]$alt <- "T"
vcfs<- do.call("rbind",vcfs)
vcfs <- vcfs[vcfs$ref != vcfs$alt,]
vcfs$name <- paste0(vcfs$Enh,"_", vcfs$Enh.Pos)
vcfs <- vcfs[,c("Enh.chr", "Pos", "name", "ref", "alt")]
vcfs <- vcfs[mixedorder(vcfs$name),]


#Save vcf into 10000 SNP segments
savevcfs <- function(vcfs, name = "") {
  dir.create(paste0("Results/Beluga/HitEnhVCFs/", name))
  if (name != "") {
    name <- paste0("_", name)
  }
  write.table(vcfs, file = paste0("Results/Beluga/HitEnhVCFs/AllAlleles_Hitenhancers",name,".vcf"), 
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  split_vcfs <- split(vcfs, (seq(nrow(vcfs))-1) %/% 10000) 
  for (split in names(split_vcfs) ) {
    file.name = paste0("Results/Beluga/HitEnhVCFs/", sub("_(.*)","\\1",name),"/",split, name,".vcf")
    write.table(split_vcfs[[split]] ,file = file.name, sep = "\t", row.names = F, col.names = F, quote = F)
  }
}

#
head(vcfs)
dir.create(paste0("Results/Beluga/HitEnhVCFs/", name))

#Don't use the identity
savevcfs(vcfs, "identity")
message("Warning fixing positional error")
#This is the correct one
vcfs$Pos <- vcfs$Pos + 1
savevcfs(vcfs, "posplus1")

#vcfs$Pos <- vcfs$Pos - 2
#savevcfs(vcfs, "posmin1")

