####################
#
# BelugaVariants
# @author: Sam Bagot
# @date: 25-10-23
#Process the Beluga outputs for the 145 functional enhancers. 
####################

setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/")
source("Scripts/Header_functions.R")
library("gtools")
enhancers <- read.csv("Data/Enhancer_factors.csv")

#Ran VCFs on Beluga webtool at this address https://hb.flatironinstitute.org/sei/ 
#using our vcf but adding 1 to the position so it matches their reference
#unzip all files from beluga variants
unzipfiles <- function(unzip = F) {
  if (unzip == F) return("Skipping Unzip")
  setwd("Results/Beluga/BelugaVariants/Webtool_Variants/")
  files <- list.files(pattern = "tar.gz")
  for (file in files) {
    system(paste("tar -xf", file))
  }
}
unzipfiles(unzip =F)


#Get all disease variants
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/")
files <- list.files("Results/Beluga/BelugaVariants/Webtool_Variants",pattern = "VARIANT_dis.tsv", full.names =T)
all_res <- lapply(files, function(x){read.table(x, header = T)})
all_res <- do.call("rbind",all_res)
colnames(all_res)[colnames(all_res) == "max_score"] <- "DIS_max_score"

files <- list.files("Results/Beluga/BelugaVariants/Webtool_Variants",pattern = "vcf", full.names =T)
all_vcf <- lapply(files, function(x){read.table(x)})
all_vcf <- do.call("rbind",all_vcf)
colnames(all_vcf) <- c("chrom", "end", "enh_pos", "ref", "alt")


all_res$chrom <- sub("chr","",all_res$chrom)
all_res <- merge(unique(all_vcf[,c("chrom", "end", "enh_pos", "ref")]), all_res)
all_res$Enh <- sub("_.*","",all_res$enh_pos)
all_res$Enh <- factor(all_res$Enh, levels = mixedsort(unique(all_res$Enh )))
all_res$Enh.pos <- as.numeric(sub(".*_","",all_res$enh_pos))
all_res <- merge(all_res, enhancers[,c("Enh", "Enh.size")])
all_res <- all_res[order(all_res$Enh),]
#Add Escore
all_res <- all_res[order(all_res$DIS_max_score, decreasing = T),]
all_res$AllVariants_Escore <-1 - rank(all_res$DIS_max_score, ties.method = "average") / nrow(all_res)
all_res$AllVariants_Log10Escore <- -log10(all_res$AllVariants_Escore)
all_res$DIS_max_score_p <- 10^-all_res$DIS_max_score


sum(all_res$DIS_max_score > -log10(0.05))
sum(p.adjust(all_res$DIS_max_score_p, "fdr") < 0.05)


#show that there appears to be a hard limit on the maximum disease score
which(all_res$DIS_max_score == max(all_res$DIS_max_score))
#check if having a max score increase the maxZ (no difference)
max_score_enhs <- as.character(unique(all_res[all_res$DIS_max_score == max(all_res$DIS_max_score),"Enh"]))
t.test(enhancers[enhancers$Enh %in% max_score_enhs,"Z"], enhancers[! enhancers$Enh %in% max_score_enhs & enhancers$HitPermissive > 0 ,"Z"])
#Apparently 89% of our variants are significant e-scores???
sum(all_res$DIS_max_score > -log10(0.05)) / nrow(all_res)
table(all_res$AllVariants_Log10Escore > -log10(0.05), all_res$DIS_max_score / 3 > -log10(0.05))
fisher.test(all_res$AllVariants_Log10Escore > -log10(0.05), all_res$DIS_max_score > -log10(0.05))

#Max score is biased (since 3 tests per postition)
table(all_res$AllVariants_Log10Escore > -log10(0.05), all_res$DIS_max_score / 3 > -log10(0.05), dnn = c("my E sig", "E / 3 sig"))
fisher.test(all_res$AllVariants_Log10Escore > -log10(0.05), all_res$DIS_max_score / 3 > -log10(0.05))

cor(all_res$AllVariants_Log10Escore, all_res$DIS_max_score)


beluga_knownvariants <-  read.csv("Results/Beluga/KnownVariants/BelugaDiseaseScores.csv")
colnames(beluga_knownvariants)[colnames(beluga_knownvariants) == "CHROM"] <- "chrom"
beluga_knownvariants <- beluga_knownvariants[,c("end","chrom","Enh","KnownVariants_DIS_max_score", "SNP","KnownVariants_Escore")]
all_res <- merge(all_res, beluga_knownvariants, all.x = T)

sum(all_res$max_score == all_res$max_score_known_variants)

#Cutoff for max_score based on our e_values 
knownvariants_E_5pt_cutoff <- min(beluga_knownvariants[beluga_knownvariants$KnownVariants_Escore <= 0.05,"KnownVariants_DIS_max_score"])
all_res$max_score_above_cutoff <- all_res$DIS_max_score > knownvariants_E_5pt_cutoff

write.csv(all_res, file = "Results/Beluga/BelugaVariants/BelugaVariants.csv", row.names = F)
