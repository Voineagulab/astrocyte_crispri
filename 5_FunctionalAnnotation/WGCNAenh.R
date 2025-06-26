###################### This script constructs the co-variation network of candidate regions based on snATAC-seq data from herring et al, and carries out downstream analyses

###################### 1. Network construction, statistical analyses of module eigengenes and related plots
###################### Set up
library(gplots)
library(data.table)
library(WGCNA)
library(ggplot2)
library(colorRamp2)
library(RColorBrewer)
library(purrr)
library(pheatmap)
library(grid)

setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/2.EnhancerCharact/Astronet_WGCNA/")
rm(list=ls())

###################### Read in and format data
# CRISPRi screen results
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
hits=res[res$HitPermissive, ]
nonhits=res[!res$HitPermissive,]

# Candidate enhancer coverage using data from snATAC-seq from Herring et al.
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Intersections and Coverage.rda")
data<- as.data.frame(log2(herring$CoverageMean_Pooled+0.1)) # data from herring et all only has coverage data for 974 of the 979 enhancers

# Order data in Herring et al by Cell type and then stage
names=colnames(data)
ct=list_transpose(strsplit(names, split="_", fixed=TRUE))[[1]]
stage=list_transpose(strsplit(names, split="_", fixed=TRUE))[[2]]

samples=data.frame(colnames(data), names, ct, stage)
order.ct =c(1:5); names(order.ct)=c("Astro", "Oligo", "Micro",  "Exc", "Inh")
order.stage =c(1:6);names(order.stage)=c("Fetal","Neonatal","Infancy","Childhood","Adolescence","Adult")

samples$orderCT=order.ct[match(samples$ct, names(order.ct))]
samples$orderS=order.stage[match(samples$stage, names(order.stage))]
sorted.ct=samples[order(samples$orderCT,samples$orderS), ]
sorted.stage=samples[order(samples$orderS, samples$orderCT), ]

samples$Astro=samples$ct%in%"Astro"
data=data[, sorted.ct$colnames.data.]
dat=t(data)
elim=which(apply(dat,2,var)==0) # filtering out one enhancer that has 0 coverage across all samples
dat=dat[,-elim]

###################### Construct the network using WGCNA

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, blockSize=40000)

pdf("Soft Thresholding.ATAC.pdf")

par(mfrow = c(1,2));cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

#---- One step network construction using power=10
net=blockwiseModules(dat,
                     power=10,
                     numericLabels=TRUE , 
                     networkType="signed", 
                     #corType = "bicor",
                     minModuleSize=50, 
                     mergeCutHeight=0.15, 
                     deepsplit=1,
                     saveTOMs=FALSE, 
                     verbose=6, 
                     nThreads=24,
                     maxBlockSize=40000, 
                     checkMissingData=FALSE)

#---- Save module information
modules=as.data.frame(table(net$colors)); colnames(modules)=c("Label", "nGenes")
modules$Color=c("grey",labels2colors(modules$Label[-1]))
modules$Label=paste("ME", modules$Label, sep="")
modules=modules[order(modules$nGenes, decreasing = TRUE) , ]

write.csv(modules, "modules.csv")

moduleLabel=paste("ME",net$colors, sep="")
moduleColor=modules$Color[match(moduleLabel, modules$Label)]

#---- Save Module Eigenegens
me<-data.frame(net$MEs)
me=me[,match(modules$Label, colnames(me)) ]
write.csv(me, "ME.csv")

#---- Save KMEs
KMEs<-signedKME(dat, net$MEs,outputColumnName = "M")
kme=data.frame( moduleColor,moduleLabel, KMEs)
rownames(kme)=rownames(KMEs)

write.csv(kme, "kME.csv")

#---- Calculate stats for the association of module eigengenes with cell type and developmental stage

stats=as.data.frame(matrix(NA, nrow=5, ncol=ncol(me)))
colnames(stats)=colnames(me)
rownames(stats)=c("p.anova.ct", "p.anova.stage", "p.w.astro", "rho.astroStage", "p.rho.astroStage")

m=match(rownames(me), sorted.ct$colnames.data.)
for (j in c(1:ncol(me))) 
{
  l=lm(me[,j] ~  sorted.ct$stage[m] + sorted.ct$ct[m])
  a=anova(l)
  stats["p.anova.stage",j]=a$`Pr(>F)`[1]
  stats["p.anova.ct",j]=a$`Pr(>F)`[2]
  w=wilcox.test(me[,j] ~ samples$Astro)
  stats["p.w.astro",j]=w$p.value
  rho=cor.test(me[samples$Astro,j] ,  samples$orderS[samples$Astro], method="s")
  stats["rho.astroStage",j]=rho$estimate
  stats["p.rho.astroStage",j]=rho$p.value
}
write.csv(stats, "moduleStats.csv")


# ###################### Plots
colheat=colorRampPalette(c("#B51A00", "#FFE2D6","white"))
colct <- c("#6F083D", "#F4EAD4","#F3CCDE","#3A5F9A","#A7D2DA")
names(colct) <- unique(samples$ct)
colstage <- brewer.pal(n = 6, name = "BuPu")
names(colstage) <- samples$stage[match(c(1:6), samples$orderS)]
ann_colors <- list(ct = colct, stage = colstage)

#### Exploratory module plots
# 
# pdf("moduleEigengenePlots.pdf", height=10, width=15)
# for (j in c(1:ncol(me)))
# {
#   subtext=paste("p.anova.ct:", round(stats["p.anova.ct",j],5), ";",
#                 "p.anova.stage:",round(stats["p.anova.stage",j],5), ";" ,
#                 "p.w.astro", round(stats["p.w.astro",j],5), ";",
#                 "rho.astroStage", round(stats["rho.astroStage",j],5), ";",
#                 "p.rho.astroStage",round(stats["p.rho.astroStage",j],5)
#   )
#   
#   colors=labels2colors(sorted.ct$ct[m])
#   n=modules$nGenes[match(colnames(me)[j] , modules$Label)]
#   par(las = 2)
#   par(mar = c(10,0,5,0))
#   barplot(me[,j], names.arg=rownames(me),
#           col=colors,
#           main=paste0("Module:", as.character(colnames(me)[j]), "; NEnh=",n, ";",subtext),
#           cex.names=1)
#   e=rownames(kme)[which(kme$moduleLabel%in%colnames(me)[j])]
#   print(pheatmap(as.matrix(data[e,]), color=rev(colheat(20)),
#                  cluster_cols = FALSE,gaps_col=grep("Adult", samples$stage[match(colnames(data), rownames(samples))]),
#                  main=colnames(me)[j]))
#   
#   #grid.newpage()
# }
# dev.off()

#### Figure plots
# Sort the kme data by module assignment
kme_sorted=kme[order(kme$moduleLabel),]
# Within each module, sort data by the kme values of the respective module. 
for (j in c(3:ncol(kme_sorted)))
{
  r=grep(colnames(kme_sorted)[j], gsub("ME", "M", kme_sorted$moduleLabel))
  i=order(kme_sorted[r,j], decreasing = TRUE)
  kme_sorted[r,]=kme_sorted[r[i], ]
}
kme_sorted=kme_sorted[-grep("ME0", kme_sorted$moduleLabel), ]
data_sorted=data[match(rownames(kme_sorted), rownames(data)) ,] 

# SOURCE DATA
write.csv(data_sorted, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/SourceData/SourceData_Fig1F.csv")

ann=kme_sorted[, c(2,1)]

ann_colors$moduleLabel=brewer.pal(n = 6, name = "Dark2")
names(ann_colors$moduleLabel)=unique(kme_sorted$moduleLabel)

rownames(samples)=samples$names
samples=samples[match(colnames(data_sorted), samples$names), ]
pdf("Fig1.ATAC.small.pdf", height=3, width=5) #use this for the main plot
print(pheatmap(as.matrix(data_sorted), 
               scale="row",
               color=rev(colheat(10)),
               cluster_cols = FALSE, cluster_rows = FALSE,
               gaps_col=grep("Adult", samples$stage[match(colnames(data), rownames(samples))]), 
               show_rownames=FALSE, show_colnames=FALSE,
               annotation_row = ann, annotation_col = samples[, c("ct","stage")], annotation_colors = ann_colors
))
grid.newpage() # this is needed for pheatmap to work nicely with pdf()
dev.off()

pdf("Fig1.ATAC.large.pdf", height=7, width=5) #use this for legends
print(pheatmap(as.matrix(data_sorted), 
               scale="row",
               color=rev(colheat(10)),
               cluster_cols = FALSE, cluster_rows = FALSE,
               gaps_col=grep("Adult", samples$stage[match(colnames(data), rownames(samples))]), 
               show_rownames=FALSE, show_colnames=FALSE,
               annotation_row = ann, annotation_col = samples[, c("ct","stage")], annotation_colors = ann_colors
))
grid.newpage() 
dev.off()

# # Subsetting for hits
# h=which(rownames(data_sorted)%in%hits$Enh)
# pdf("Fig3.hitsATAC.pdf", height=3, width=5) 
# print(pheatmap(as.matrix(data_sorted[h,]), 
#                scale="row",
#                color=rev(colheat(10)),
#                cluster_cols = FALSE, cluster_rows = FALSE,
#                gaps_col=grep("Adult", samples$stage[match(colnames(data), rownames(samples))]), 
#                show_rownames=FALSE, show_colnames=FALSE,
#                annotation_row = ann, annotation_col = samples[, c("ct","stage")], annotation_colors = ann_colors
# ))
# grid.newpage() 
# dev.off()

#plot pvalues only (remove rho) and remove the ME0 grey module.
splot=data.frame(t(stats[-which( rownames(stats)%in%"rho.astroStage"), -grep("ME0", colnames(stats))]))
splot=-log10(splot)

# SOURCE DATA
write.csv(splot, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/SourceData/SourceData_Fig1G.csv")

for(j in c(1:ncol(splot))) splot[which(splot[,j]< (-log10(0.05))), j]=0

color_ramp <- colorRamp(c("white", "midnightblue"))
palette <- rgb(color_ramp(seq(0, 1, length.out = 1000)), maxColorValue = 255)
ann2=ann[match(rownames(splot), ann$moduleLabel), ]

colnames(splot)=gsub("p.anova.ct", "Cell Type",colnames(splot))
colnames(splot)=gsub("p.anova.stage", "Stage",colnames(splot))
colnames(splot)=gsub("p.w.astro", "Astrocyte",colnames(splot))
colnames(splot)=gsub("p.rho.astroStage", "Astrocyte Stage",colnames(splot))

rownames(splot)=gsub("ME", "M", rownames(splot))
rownames(ann2)=rownames(splot)

pdf("Fig1.ModuleSig.pdf", height=3, width=4) 
pheatmap(t(splot), cluster_rows = FALSE,cluster_cols = FALSE, color = palette, annotation_col = ann2, annotation_colors = ann_colors)
grid.newpage() 
dev.off()

###################### 2. Comparison of Co-variation Modules vs Hit Enhancers and associated plots
######## Modules vs Hits - Exploratory analysis
# rm(list=ls())
# res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
# power=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/PowerCalculationSuppTable.csv")
# 
# power=power[power$WellPowered015, ]
# keep=c(which(res$HitPermissive==TRUE), which(res$Pair%in%power$Pair))
# 
# res=res[keep,]
# hits=res[res$HitPermissive, ]
# 
# kme=read.csv("kME.csv")
# colnames(kme)=gsub("M", "ME", colnames(kme))
# rownames(kme)=kme[,1]; kme=kme[, -1]
# kme$Hit=rownames(kme)%in% hits$Enh
# 
# kme_powered=kme[which(rownames(kme)%in%res$Enh), ]
# 
# stats=read.csv("moduleStats.csv")
# rownames(stats)=stats[,1]; stats=stats[, -1]
# stats=t(stats)
# 
# names=c(colnames(stats), "p.wilcox.hits", "p.wilcox.hits_vsPowered")
# stats=data.frame(stats, NA, NA); colnames(stats)=names
# 
# for(j in c(1:nrow(stats)))
# {
# stats$p.wilcox.hits[j]=wilcox.test(kme[, rownames(stats)[j]] ~ kme$Hit)$p.value
# stats$p.wilcox.hits_vsPowered[j]=wilcox.test(kme_powered[, rownames(stats)[j]] ~ kme_powered$Hit)$p.value
# }
# write.csv(t(stats), "moduleStats_Hits.csv")

###################### Setup
library(purrr)
library(pheatmap)
library(grid)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/2.EnhancerCharact/Astronet_WGCNA/")
rm(list=ls())

###################### Read in and format data (some steps are repeated from #1 above as the analysis was done subsequently)

# CRISPRi screen results
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv") 
# Subset for hits
hits=res[res$HitPermissive, ] 

# Candidate enhancer coverage using data from snATAC-seq from Herring et al.
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Intersections and Coverage.rda")
data<- as.data.frame(log2(herring$CoverageMean_Pooled+0.1)) # data from herring et all only has coverage data for 974 of the 979 enhancers

# Order data in Herring et al by Cell type and then stage
names=colnames(data)
ct=list_transpose(strsplit(names, split="_", fixed=TRUE))[[1]]
stage=list_transpose(strsplit(names, split="_", fixed=TRUE))[[2]]

samples=data.frame(colnames(data), names, ct, stage)
order.ct =c(1:5); names(order.ct)=c("Astro", "Oligo", "Micro",  "Exc", "Inh")
order.stage =c(1:6);names(order.stage)=c("Fetal","Neonatal","Infancy","Childhood","Adolescence","Adult")

samples$orderCT=order.ct[match(samples$ct, names(order.ct))]
samples$orderS=order.stage[match(samples$stage, names(order.stage))]
sorted.ct=samples[order(samples$orderCT,samples$orderS), ]
sorted.stage=samples[order(samples$orderS, samples$orderCT), ]

samples$Astro=samples$ct%in%"Astro"
data=data[, sorted.ct$colnames.data.]

# Subsetting coverage data for hit enhancers
data=data[which(rownames(data)%in%hits$Enh),]

# Read in ANOVA stats (enhancer level ANOVA for cell type and stage)
stats=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/SuppTable_HerringDevCTreg.csv")
stats=stats[match(rownames(data), stats[,1]) , ]
stats$padj.anova.ct[which(stats$padj.anova.ct >= 0.05)]=1
stats$padj.anova.stage[which(stats$padj.anova.stage >= 0.05)]=1
ann=stats
ann$pval_CT=-log10(stats$padj.anova.ct)
ann$pval_Stage=-log10(stats$padj.anova.stage)
rownames(ann)=ann[,1] ; ann=ann[,-1]

###################### Plots

### Heatmap
colheat=colorRampPalette(c("#B51A00", "#FFE2D6","white"))
colct <- c("#6F083D", "#F4EAD4","#F3CCDE","#3A5F9A","#A7D2DA")
names(colct) <- unique(samples$ct)
colstage <- brewer.pal(n = 6, name = "BuPu")
names(colstage) <- samples$stage[match(c(1:6), samples$orderS)]
ann_colors <- list(ct = colct, stage = colstage)

rownames(samples)=samples$names
samples=samples[match(colnames(data), samples$names), ]

# SOURCE DATA:
write.csv(data, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/SourceData/SourceData_Fig3F.csv")

pdf("Fig3.hitsATAC_v2_small.pdf", height=4, width=6) 
print(pheatmap(t(scale(t(data))),
               scale="none",
               color=rev(colheat(10)),
               cluster_cols = FALSE, 
               gaps_col=grep("Adult", samples$stage[match(colnames(data), rownames(samples))]), 
               show_rownames=FALSE, show_colnames=FALSE,
               annotation_row = ann[, c("pval_CT","pval_Stage")], annotation_col = samples[, c("ct","stage")], annotation_colors = ann_colors
               
))
grid.newpage()
dev.off()

pdf("Fig3.hitsATAC_v2_large.pdf", height=8, width=8) 
print(pheatmap(t(scale(t(data))),
               scale="none",
               color=rev(colheat(10)),
               cluster_cols = FALSE, 
               gaps_col=grep("Adult", samples$stage[match(colnames(data), rownames(samples))]), 
               show_rownames=FALSE, show_colnames=FALSE,
               annotation_row = ann[, c("pval_CT","pval_Stage")], annotation_col = samples[, c("ct","stage")], annotation_colors = ann_colors
               
))
grid.newpage()
dev.off()

### Barplot
barplot_data=as.data.frame(table(stats$padj.anova.ct< 0.05, stats$padj.anova.stage<0.05))
colnames(barplot_data)=c("CT", "Stage", "Count")

# SOURCE DATA:
write.csv(barplot_data, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/SourceData/SourceData_Fig3G.csv")

