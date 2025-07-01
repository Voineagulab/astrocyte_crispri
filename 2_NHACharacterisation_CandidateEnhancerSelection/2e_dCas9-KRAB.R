rm(list=ls())
library(pheatmap)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/NHAcharacterisation_resubmission/RNA_seq/star_salmon")
## read in RNA-seq data
data=read.delim("salmon.merged.gene_tpm.tsv", sep="\t", header=TRUE)
## read in the results of the CRISPRi screen
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")

#Subset the data for NHA samples
data=data[, c(1,2,grep("NHA_", colnames(data)))]

# Genes tested in the CRISPRi screen
t=which(data$gene_name%in%res$Gene)
# Genes with significant effects in the CRISPRi screen
s=which(data$gene_name%in%res$Gene[res$HitPermissive])

# RNA-seq from NHAs
x=data$NHA_2019_GRE13789A15_22KTVVLT3
# RNA-seq from NHAs transduced with dCas9-KRAB
y=data$NHA_2019_dCas9_KRAB_GRE13789A13_22KTVVLT3

#Scatterplot shown in Extended Data Figure 1
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/dCas9-KRAB.pdf")
plot(log2(y) ~ log2(x),xlab="NHA:log2(TPM)", ylab="NHA-dCas9-KRAB:log2(TPM)",
     pch=20, col="grey", 
     sub=paste0("r=",round(cor(x,y, method="p"),2), " , rho=", round(cor(x,y, method="s"),2) ," ; blue:tested , red: significant") 
     )
points(log2(y[t]) ~ log2(x[t]),
     pch=20, col="midnightblue")

points(log2(y[s]) ~ log2(x[s]),
       pch=20, col="indianred")
dev.off()
