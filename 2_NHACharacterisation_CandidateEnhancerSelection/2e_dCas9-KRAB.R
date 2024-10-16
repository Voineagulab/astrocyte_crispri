rm(list=ls())
library(pheatmap)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/NHAcharacterisation_resubmission/RNA_seq/star_salmon")
data=read.delim("salmon.merged.gene_tpm.tsv", sep="\t", header=TRUE)
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
data=data[, c(1,2,grep("NHA_", colnames(data)))]
pheatmap(cor(data[, -c(1,2)], method="s"))

t=which(data$gene_name%in%res$Gene)
s=which(data$gene_name%in%res$Gene[res$HitPermissive])

x=data$NHA_2019_GRE13789A15_22KTVVLT3
y=data$NHA_2019_dCas9_KRAB_GRE13789A13_22KTVVLT3
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


plot(log2(y/x) ~ log2(x),xlab="NHA:log2(TPM)", ylab="log2FC NHA-dCas9-KRAB/NHA",
     pch=20, col="grey", sub=paste0("rho=",round(cor(x,y, method="s"),2)))

points(log2(y[t]/x[t]) ~ log2(x[t]),
       pch=20, col="midnightblue", sub=round(cor(x,y, method="s"),2))

points(log2(y[s]/x[s]) ~ log2(x[s]),
       pch=20, col="indianred", sub=round(cor(x,y, method="s"),2))
