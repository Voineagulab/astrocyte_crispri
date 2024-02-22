################
##
##
## Beluga Known variants VCF analysis
## @author: Sam Bagot
## @date: 20-10-23
################

#Beluga predicts the effect of individual variants on Disease
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/")
source("Scripts/Header_functions.R")
enhancers <- read.csv("Data/Enhancer_factors.csv")
library("ggsignif")
library("data.table")

########
##Get known variants and write to vcf for upload to Beluga
########
#Read snp2peak subset to SNPS within 1kb window
snp2peak <- read.csv("../FullScale/Results/3_HitEnrichment/Variants/Final - SNP-to-Peak List.csv")
snp2peak <- snp2peak[snp2peak$DistanceFromCentre <= 1000 & snp2peak$DistanceFromCentre >=- 1000,]

#This changes to long (based on synonyms) so we can catch all SNPs possible
snp2peak$Synonym <- paste0(snp2peak$SNP, "|", snp2peak$Synonym )
snp2peak <- data.table(snp2peak)
expand_syns <- as.data.frame(unique(snp2peak[ , list( Synonym = unlist( strsplit( Synonym , "\\|" ) ) ) , by = SNP ]))
snp2peak <- as.data.frame(snp2peak)
snp2peak <- merge(snp2peak[,colnames(snp2peak) != "Synonym"], expand_syns)
snp2peak$CHROM <- sub(":.*","",snp2peak$Peak)


#Gets the Alleles from the vcf
call <- paste("vcftools --gzvcf /home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/dbSNPs/GRCh38_latest_dbSNP_all.vcf.gz", #OverlappingSNPs.txt was replaced by Gavins version
              "--recode --snps Results/Beluga/KnownVariants/OverlappingSNPs_GJS_20231108.txt --out Results/Beluga/KnownVariants/OverlappingSNPs" )
#"conda activate mosdepth" (or any environment with vcftools installed)
system(call) #this doesn't work I just ran it in terminal

overlap_vcf <- read.table("Results/Beluga/KnownVariants/OverlappingSNPs.recode.vcf")[3:5]
colnames(overlap_vcf) <- c("Synonym" ,"REF", "ALT") #Changed SNP to Synonym
out_vcf <- merge(unique(snp2peak[,c("Synonym", "SNP.pos", "CHROM")]), overlap_vcf)
out_vcf <- out_vcf[match(unique(out_vcf$SNP), out_vcf$SNP),c("CHROM", "SNP.pos","Synonym", "REF", "ALT")]
out_vcf$ALT <- sub(",.*","",out_vcf$ALT) #Remove all except first reference
write.table(out_vcf, file = "Results/Beluga/KnownVariants/OverlappingSNPs.recode.formatted.vcf", quote = F, row.names = F, col.names = F, sep = "\t")
#It appears that there is a one BP out of position
#Upload THIS FILE TO BELUGA
out_vcf$SNP.pos <- out_vcf$SNP.pos + 1
write.table(out_vcf, file = "Results/Beluga/KnownVariants/OverlappingSNPs.recode.formatted_posplus1.vcf", quote = F, row.names = F, col.names = F, sep = "\t")

#Links to download the Sei/Beluga Results
#New Sei with updated SNPs and pos +1 
#https://hb.flatironinstitute.org/deepsea/jobs/3de40fcf-ddb9-4ca9-8699-37ac1c69b521/
#New Beluga with updated SNPs and pos +1 
#https://hb.flatironinstitute.org/deepsea/jobs/746aa96f-3894-4344-b085-b8b1b81b08c2/

beluga_results<- read.table("Results/Beluga/KnownVariants/BelugaResults/746aa96f-3894-4344-b085-b8b1b81b08c2_OverlappingSNPs.recode.formatted_posplus1_VARIANT_dis.tsv", sep = "\t", header = T)
colnames(beluga_results) <- c("CHROM", "SNP.pos", "end", "KnownVariants_DIS_max_score")
beluga_results <- merge(beluga_results, unique(snp2peak[,c("SNP","CHROM", "SNP.pos","Enh", "Hit")]))
nrow(beluga_results)
beluga_results <- beluga_results[order(beluga_results$KnownVariants_DIS_max_score, decreasing = T),]
beluga_results$KnownVariants_Escore <-  1- rank(beluga_results$KnownVariants_DIS_max_score, ties.method = "average") / nrow(beluga_results)
beluga_results$KnownVariants_Log10Escore <- -log10(beluga_results$KnownVariants_Escore)
beluga_results$CHROM <- sub("chr","",beluga_results$CHROM)
write.csv(beluga_results, "Results/Beluga/KnownVariants/BelugaDiseaseScores.csv", row.names = F)


max_Dis <- aggregate(beluga_results$KnownVariants_DIS_max_score, by = list(beluga_results$Enh), max)
colnames(max_Dis) <- c("Enh", "Beluga.MaxDisScore")
write.csv(max_Dis,"Results/Beluga/KnownVariants/MaxBelugaDiseaseScores.csv", row.names = F)

#Statss on beluga max score
ft.dis <-fisher.test(beluga_results$Hit,beluga_results$KnownVariants_DIS_max_score > - log10(0.05)) #fisher.test(beluga_results$KnownVariants_DIS_max_score > 2, beluga_results$Hit)
t.dis <- t.test(beluga_results[beluga_results$Hit,]$KnownVariants_DIS_max_score,beluga_results[! beluga_results$Hit,]$KnownVariants_DIS_max_score)
ft.dis
t.dis
wt.e_score <- wilcox.test(beluga_results[beluga_results$Hit,]$KnownVariants_Escore,beluga_results[! beluga_results$Hit,]$KnownVariants_Escore)
wt.DIS_score <- wilcox.test(beluga_results[beluga_results$Hit,]$KnownVariants_DIS_max_score,beluga_results[! beluga_results$Hit,]$KnownVariants_DIS_max_score)


beluga_results <- beluga_results[order(beluga_results$Hit, decreasing = T),]
beluga_results_unique <- beluga_results[match(unique(beluga_results$SNP), beluga_results$SNP),]

tests.list <- list()
tests.list[["ft.dis.uniq"]] <-fisher.test(beluga_results_unique$Hit,beluga_results_unique$KnownVariants_DIS_max_score > - log10(0.05)) #fisher.test(beluga_results$KnownVariants_DIS_max_score > 2, beluga_results$Hit)
tests.list[["ft.e.uniq"]] <-fisher.test(beluga_results_unique$Hit,beluga_results_unique$KnownVariants_Escore < 0.05) #fisher.test(beluga_results$KnownVariants_DIS_max_score > 2, beluga_results$Hit)
tests.list[["t.dis.uniq"]] <- t.test(beluga_results_unique$KnownVariants_DIS_max_score ~ beluga_results_unique$Hit)
#Fisher TEst  score < 0.05 &  hits
saveRDS(tests.list, "Results/Beluga/KnownVariants/Unique_SNPs_TestResults.RDS")

