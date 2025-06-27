# This script generates the astrocyte regulatory network:
#1. selects the genes, enhancers and TFs included in the network
#2. plots the heatmaps based on pseudobulked snRNA-seq and snATAC-seq data
#3. generates the input files for Cytoscape plots
# It also outputs statistical tests mentioned in the text as a txt file

# NOTE: ignore the following error: Error in dev.off() : cannot shut down device 1 (the null device)
# oddly the pheatmap function needs dev.off() before AND after pdf() to print to a pdf. 

# Load required libraries
library(data.table)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Clear workspace
rm(list=ls())

# Set working directory and sink
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/2.EnhancerCharact/")
#sink("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/Res4Fig5.info_v2.txt")

################# Load data
# snRNA and snATAC, pseudobulked
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Intersections and Coverage.rda")
atac <- as.data.frame(herring$CoverageMean_Pooled)

load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/Pseudobulk2_ATACesquePooling.rda")
exp <- pb2$Final
exp <- exp[, which(colnames(exp) %in% colnames(atac))]
atac <- atac[, match(colnames(exp), colnames(atac))]

# Astrocyte-specific genes based on snRNA-seq data from Herring et al.
herring.genes <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Genes/Astrocyte Markers - Herring 2022.csv")

# Enhancer annotation including astrocyte-specific enhancers based on snATAC-seq data from Herring et al.
ann.enh <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Final - Annotation Logical.csv")
diff.atac=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Herring - Cell-type Specificity Models (Repooled).csv")
rownames(diff.atac)=diff.atac[,1]; diff.atac=diff.atac[,-1]

# CRISPRi hits
res <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
hits <- res[res$HitPermissive, ]

# TF footprinting data
bound <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Summaries/Bound_Matrix_Expressed.csv")
rownames(bound) <- bound$Enh
bound <- bound[, grep("Bound", colnames(bound))]
colnames(bound) <- gsub("Tobias.Bound_", "", colnames(bound))

# TF enrichment based on footprinting
ft <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Summaries/FisherTestsExpressedTFMatrix_oldformat.csv")
ft$TF <- data.table::transpose(strsplit(ft$Variable, split="_", fixed=TRUE))[[2]]
ft <- ft[-grep("Exp", ft$Variable), ]

# Core regulatory circuitry
crc <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/SaintAndre_CRC/CRC_Astro.csv")

# Jaspar cluster annotation
jasp <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/JASPAR_TF_Cluster_Assignments.csv")

############################################################## Core regulatory circuitry
ftb <- ft[grep("Bound", ft$Variable), ]
ftb$crc <- ftb$TF %in% crc$Astrocytes
cat("Fisher test for TFs enriched in hit enhancers vs CRC TFs from Saint Andre et al" , "\n")
print(fisher.test(ftb$crc, ftb$FDR < 0.05))

############################################################## Which TFs show the strongest enrichments?
cat("TFs with FDR<0.05 and OR>3 for enrichment among hit enhancers", "\n")
print(sort(ftb$TF[(ftb$FDR < 0.05 & ftb$T.Stat_or_Odds_ratios > 3)]))

############################################################## Select data for Astrocyte Specific Regulatory network
# Select astro-enriched hit genes based on the Herring data
genes <- intersect(hits$Gene, herring.genes$Gene[(herring.genes$FDR < 0.05) & (herring.genes$log2fc > log2(1.5))])

# Select astro enhancers as enhancers that regulate astro-enriched genes; save Supplementary table
EGnet <- hits[hits$Gene %in% genes, c("Pair", "Enh", "Gene","logfc.vst", "Z","P.SCEPTRE","FDR.N50")]
EGnet <- EGnet[ EGnet$Enh %in% ann.enh$Enh[ann.enh$AstSpecific_Herring_Pooled] , ]
write.csv(EGnet, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/SuppTable.EGdata_v2.csv",row.names=FALSE)

# Select TFs that bind to astro enhancers and have ast-specific expression; save Supplementary table
tfbound.ast <- bound[unique(EGnet$Enh), -grep("Tobias.Exp_Bound_TF_counts", colnames(bound))] #remove the column that counts the total number of binding sites of expressed TFs per enhancer 
k <- colSums(tfbound.ast[, ] > 0)
tfbound.ast <- tfbound.ast[, k > 0]

TEnet <- melt(data.frame(rownames(tfbound.ast), tfbound.ast))
colnames(TEnet) <- c("Enh", "TF", "Bound")
TEnet <- TEnet[TEnet$Bound > 0, ]

TEnet$Enh_AstSpecif <- TEnet$Enh %in% ann.enh$Enh[ann.enh$AstSpecific_Herring_Pooled]
TEnet$TF_AstSpecific <- TEnet$TF %in% herring.genes$Gene[(herring.genes$FDR < 0.05) & (herring.genes$log2fc > log2(2))]

TEnet <- TEnet[TEnet$Enh_AstSpecif & TEnet$TF_AstSpecific,]
write.csv(TEnet, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/SuppTable.TEdata_v2.csv",row.names=FALSE)

############################################################## Cytoscape network
EGnet.cyto <- EGnet[, c("Enh", "Gene")]
colnames(EGnet.cyto) <- c("Node1", "Node2")
TEnet.cyto <- TEnet[, c("TF", "Enh")]
colnames(TEnet.cyto) <- c("Node1", "Node2")

### Plot heatmaps of TFs, Enhancers, and Genes included in the network (Sfig 7B., Fig 4D)
# TF expression
dat.tf <- log2(exp[which(rownames(exp) %in% TEnet.cyto$Node1), ] + 0.1)

# Paired enhancer snATAC and target gene snRNA
h <- EGnet[, c("Enh", "Gene")]
hg <- sort(unique(h$Gene))
b <- rep(NA, length(hg))
for (j in c(1:length(hg))) {
  g <- hg[j]
  e <- grep(g, h$Gene)
  d <- rbind(exp[which(rownames(exp) %in% g), ],
             atac[which(rownames(atac) %in% h$Enh[e]), ])
  d<- data.frame(rownames(d), d)
  if (j == 1) dat.eg <- d else  dat.eg <- rbind(dat.eg, d)
  if (j == 1)  b[j] = nrow(d) else b[j] = b[j - 1] + nrow(d)
}
colnames(dat.eg)[1]<- "Id"
# Heatmap formatting
coldat <- data.frame(data.table::transpose(strsplit(colnames(dat.tf), split = "_", fixed = TRUE))[[1]],
                     data.table::transpose(strsplit(colnames(dat.tf), split = "_", fixed = TRUE))[[2]])
colnames(coldat) <- c("CT", "Stage")
rownames(coldat) <- colnames(dat.tf)

colheat=colorRampPalette(c("#B51A00","#FFE2D6", "#003F01"))
colct <- c("#6F083D", "#F4EAD4","#F3CCDE","#3A5F9A","#A7D2DA")
names(colct) <- unique(coldat$CT)
colstage <- brewer.pal(n = 6, name = "BuPu")
names(colstage) <- unique(coldat$Stage)
ann_colors <- list(CT = colct, Stage = colstage)

dev.off() 
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/AstroNet.HeatmapTF_v2.pdf", height = 7, width = 6)
pheatmap(dat.tf,
         border_color = NA, show_colnames = FALSE,
         fontsize_row = 10, scale = "row", 
         color = rev(colheat(100)), 
         annotation_col = coldat, annotation_colors = ann_colors)
dev.off()
write.csv(dat.tf, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/EData.TFheatmap_v2.csv")

dev.off()
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/AstroNet.HeatmapEG_v2.pdf", height = 8, width = 6)
pheatmap(log2(dat.eg[,-1] + 0.1), fontsize_row = 10, scale = "row",
         border_color = NA, show_colnames = FALSE, cluster_rows = FALSE, gaps_row = b,
         color = rev(colheat(100)), annotation_col = coldat, annotation_colors = ann_colors,
         labels_row=dat.eg$Id)
dev.off()
write.csv(dat.eg, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/EData.EGheatmap_v2.csv", row.names=FALSE)

### Collapse TFs that belong to the same JASPAR cluster based on their binding motifs
cluster <- jasp[match(unique(TEnet.cyto$Node1), jasp$TF), c("TF", "JASPAR_clust")]
f <- function(x) paste(x, collapse = "/")
names <- aggregate(TF ~ JASPAR_clust, data = cluster, FUN = f)
colnames(names)[2] <- "TFcluster"
cluster <- merge(cluster, names, by = "JASPAR_clust")

TEnet.cyto$TFcluster <- cluster$TFcluster[match(TEnet.cyto$Node1, cluster$TF)]
TEnet.cyto$Pair <- paste(TEnet.cyto$TFcluster, TEnet.cyto$Node2, sep = "_")
m <- match(unique(TEnet.cyto$Pair), TEnet.cyto$Pair)
TEnet.cyto <- TEnet.cyto[m, c(3, 2)]
colnames(TEnet.cyto)[1] <- "Node1"

### Assign an arbitrary edge value
EGnet.cyto$Edge <- 1.2
TEnet.cyto$Edge <- 1
TEGnet <- rbind(EGnet.cyto,  TEnet.cyto)

### Make node attribute file
Enodes <- data.frame(unique(EGnet.cyto$Node1), "Enh")
Gnodes <- data.frame(unique(EGnet.cyto$Node2), "Gene")
TFnodes <- data.frame(unique(TEnet.cyto$Node1), "TF")

colnames(Enodes) <- colnames(Gnodes) <- colnames(TFnodes) <- c("Node", "Categ")
TEGnodes <- rbind(Enodes, Gnodes, TFnodes)

### Abbreviate long cluster names
a <- read.csv("AstroNet.TEG.Abv.csv")
for (j in c(1:nrow(a))) {
  TEGnet$Node1 <- gsub(a$TF.cluster[j], a$Abbreviation[j], TEGnet$Node1)
  TEGnodes$Node <- gsub(a$TF.cluster[j], a$Abbreviation[j], TEGnodes$Node)
}
### Save
write.csv(TEGnet, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/AstroNet.TEG_v2.csv", row.names = FALSE, quote = FALSE)
write.csv(TEGnodes, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/AstroNodes.TEG_v2.csv", row.names = FALSE, quote = FALSE)
sink()

#SOURCE DATA:
write.csv(dat.eg, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/SourceData/Fig4D_Heatmap.csv")

herrring.genes.net=herring.genes[which(herring.genes$Gene%in%rownames(dat.eg)) , ]
write.csv(herrring.genes.net, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/SourceData/Fig4D_DiffExp.csv")
          
diff.atac.net=diff.atac[which(rownames(diff.atac)%in%rownames(dat.eg)), ]
write.csv(diff.atac.net, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/SourceData/Fig4D_DiffATAC.csv")

#old=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/Publication/AstroNodes.TEG.csv")
#setdiff(old$Node, TEGnodes$Node)
