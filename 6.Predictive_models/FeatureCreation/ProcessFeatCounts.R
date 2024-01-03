##############
#
# Process Feature counts 
# @author: Gavin Sutton & Sam Bagot
# @date: Unknown
#
###############

## Load data produced in FullScale/Scripts/4_EnhancerTranscription.R
# see the section on TTseq_QC (line ~600 for how this was generated)
#run feature counts
setwd("/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/")
library("Rsubread")  
bam_astro_rnaseq <- "Data/GeneData/Ref_AstroRNAseq/ENCFF205MGD.bam"
gtf <- "/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/cellranger_refdata-gex-GRCh38-2020-A/genes/genes.gtf"
exon_counts <- featureCounts(bam_astro_rnaseq,
                             annot.ext = gtf, isGTFAnnotationFile = TRUE,
                             useMetaFeatures = TRUE,
                             GTF.attrType = "gene_name",
                             allowMultiOverlap = TRUE,
                             isPairedEnd = TRUE, # self-explanatory
                             nthreads = 8, # self-explanatory
                             strandSpecific = 2, # stranded reverse
                             minMQS = 10, # minimum quality score at either end of the pair
                             checkFragLength = FALSE,
                             countMultiMappingReads = FALSE,
                             requireBothEndsMapped = FALSE,
                             countChimericFragments = FALSE)
saveRDS(exon_counts, "Data/GeneData/Ref_AstroRNAseq/FeatureCounts_Exons.rds")
#######
#Repeat process for K562 (not making a function as only running feature counts)
bam_k562_rnaseq <- "Data/K562Data/RNAseq/ENCFF950KXS.bam"
exon_counts <- featureCounts(bam_k562_rnaseq,
                             annot.ext = gtf, isGTFAnnotationFile = TRUE,
                             useMetaFeatures = TRUE,
                             GTF.attrType = "gene_name",
                             allowMultiOverlap = TRUE,
                             isPairedEnd = TRUE, # self-explanatory
                             nthreads = 8, # self-explanatory
                             strandSpecific = 2, # stranded reverse
                             minMQS = 10, # minimum quality score at either end of the pair
                             checkFragLength = FALSE,
                             countMultiMappingReads = FALSE,
                             requireBothEndsMapped = FALSE,
                             countChimericFragments = FALSE)
saveRDS(exon_counts, "Data/K562Data/RNAseq/FeatureCounts_Exons.rds")


#Gavin's workflow
load("../FullScale/Results/4_EnhancerTranscription/TTseq/FeatureCounts_Exons.rda")
x <- exon_counts
## Wrangle
feature_counts <- as.data.frame(exon_counts$counts)
colnames(feature_counts) <- c("TTseq_Counts", "RNAseq_Counts")
# make CPM
libSize <- sum(exon_counts$stat$GOK10856A1.markdup.sorted.bam)
feature_counts$RNAseq_CPM <- feature_counts$RNAseq_Counts / (libSize / 10^6)
feature_counts$TTseq_CPM <- feature_counts$TTseq_Counts / (libSize / 10^6)
# make RPKM
feature_counts$RNAseq_RPKM <- feature_counts$RNAseq_CPM / (exon_counts$annotation$Length / 1000)
feature_counts$TTseq_RPKM <- feature_counts$TTseq_CPM / (exon_counts$annotation$Length / 1000)
## Save
write.csv(feature_counts, file = "Data/GeneData/RNAseq_FeatCounts.csv", row.names = F)

######## (Sam Bagot) Get feature counts
#Feature counts for new  data
process_featureCounts<- function (rdsfile, colname = "Ref_RNAseq_Counts") {
  featureCounts_refernce <- readRDS(rdsfile)
  fc <- as.data.frame(featureCounts_refernce$counts)
  colnames(fc) <- c(colname)
  libSize <- sum(featureCounts_refernce$stat[,2]) #this is only correct if only mapping one Bam
  fc$Ref_RNAseq_CPM <- fc$Ref_RNAseq_Counts / (libSize / 10^6)
  fc$Ref_RNAseq_RPKM <- fc$Ref_RNAseq_CPM / (featureCounts_refernce$annotation$Length / 1000)
  return(fc)
}
fc <- process_featureCounts("Data/GeneData/Ref_AstroRNAseq/FeatureCounts_Exons.rds")
write.csv(fc, file = "Data/GeneData/Ref_RNAseq_FeatCounts.csv", row.names = F)

#K562 gene expression
fc <- process_featureCounts("Data/K562Data/RNAseq/FeatureCounts_Exons.rds")
geneinfo <- read.table("../FullScale/Data/Whitelists/GeneInfo.txt")
fc$Gene <- rownames(fc)
fc <- merge(geneinfo[,c( "EnsID","Symbol")], fc,  by.y = "Gene", by.x = "Symbol", all.y = T)
colnames(fc)[colnames(fc) == "Symbol"] <- "TargetGene"
write.csv(fc, file = "Data/K562Data/RNAseq/Ref_RNAseq_FeatCounts.csv", row.names = F)


  