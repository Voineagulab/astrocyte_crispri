# @author: Irina Voineagu + Sam Bagot
# @date: 18-11-23

beluga=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Beluga/BelugaVariants/BelugaVariants.csv")
makebedGraph <- function(bed, filename) {
  bed$chrom=paste0("chr", bed$chrom)
  bed$Id=paste(bed$chrom, bed$pos, bed$end, sep="_")
  t=as.data.frame(table(bed$Id))
  bed=bed[which(bed$Id%in%t$Var1[t$Freq==1]), c(1:4)]
  write.table(bed, filename,
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


makebedGraph(beluga[, c("chrom", "pos", "end", "DIS_max_score")], 
             "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Beluga/BelugaVariants/BelugaVariants.bedGraph")
makebedGraph(beluga[, c("chrom", "pos", "end", "AllVariants_Log10Escore")], 
             "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Beluga/BelugaVariants/BelugaVariants_Escore.bedGraph")
#> length(unique(bed$Id))
#[1] 51782
#> dim(bed)
#[1] 51792 
#There are 10 positions that are duplicated. Is it because of overlapping enhancers?
#Also, strangely the disease scores are different for identical positions:must check with Sam
# chrom     pos     end DIS_max_score                   Id
# 6668  chr1 8890764 8890765     0.9625695 chr1_8890764_8890765
# 6669  chr1 8890764 8890765     0.6123095 chr1_8890764_8890765
# Will remove those as they are ambiguous