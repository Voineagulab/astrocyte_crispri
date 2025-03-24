library(data.table)

rm(list=ls())
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.Results.rda")
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.DEP.rda")
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
# subset for powered EG pairs
power=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/PowerCalculationSuppTable.csv")
power=power[power$WellPowered015, ]
keep=union(which(res$HitPermissive==TRUE), which(res$Pair%in%power$Pair))
res=res[keep,]

genes=data.frame(unique(res$Gene), FALSE)
colnames(genes)=c("Symbol", "Hit")
genes$Hit[which(genes$Symbol%in%res$Gene[res$HitPermissive])]=TRUE

res$DE_AD=res$Gene%in%c(rosmap.results$FourGroup$sig.genes$Symbol, rosmap.results$DE$sig.genes$Symbol)

table(res$HitPermissive, res$DE_AD)
# FALSE TRUE
# FALSE  1053 1115
# TRUE     58  100
fisher.test(res$HitPermissive, res$DE_AD)
# Fisher's Exact Test for Count Data
# 
# data:  res$HitPermissive and res$DE_AD
# p-value = 0.003888
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.152770 2.316272
# sample estimates:
# odds ratio 
#   1.627922 

###### Core regulatory circuitries https://genome.cshlp.org/content/26/3/385.long 
crc=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/SaintAndre_CRC/Extended_RC_Astro.csv")
genes$CRC=genes$Symbol%in%crc$Astrocytes
fisher.test(genes$Hit, genes$CRC)

#Fisher's Exact Test for Count Data

# data:  genes$Hit and genes$CRC
# p-value = 2.719e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.350832 8.009876
# sample estimates:
# odds ratio 
#   4.399298 

