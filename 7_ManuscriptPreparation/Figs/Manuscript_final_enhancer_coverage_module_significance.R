# This script generates a heatmap of enhancer coverage (ATAC-seq, H3K27ac, H3K4me3) across cell types
#and assesses module significance using Kruskal-Wallis and Wilcoxon tests;  H3K27ac and H3K4me3 data from  Nott et al. (2019).
#Generates Sup Fig3C and Sup Fig3D

# Load necessary libraries
library(gplots)
library(data.table)
library(WGCNA)
library(ggplot2)
library(colorRamp2)
library(RColorBrewer)
library(purrr)
library(pheatmap)
library(grid)
library(FSA)
library(readxl)
library(tidyverse)
library(valr)

rm(list=ls())

####################################################################################################
##Plot SupFig 3D: Co-accessibility network module significance based on data from Nott et al. (2019)


# Load kME data from WGCNA analysis of candidate enhancer ATAC-seq peaks
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/Astronet_WGCNA")
kme=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/Astronet_WGCNA/kME.csv")

# Sort the kME data by module assignment for better visualization
rownames(kme)=kme[,1]; kme=kme[,-1]
kme_sorted=kme[order(kme$moduleLabel),]
for (j in c(3:ncol(kme_sorted)))
{
  r=grep(colnames(kme_sorted)[j], gsub("ME", "M", kme_sorted$moduleLabel))
  i=order(kme_sorted[r,j], decreasing = TRUE)
  kme_sorted[r,]=kme_sorted[r[i], ]
}
# Remove "ME0" module 
kme_sorted=kme_sorted[-grep("ME0", kme_sorted$moduleLabel), ]


# Load coordinates of candidate enhancers in hg17 format (converted using UCSC liftover tool)
candidate_enh <- read_excel("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable1_ATAC/Supplementary Table 1 - Annotation of candidate peaks.xlsx", sheet = "1A_AllPeakAnnotation")
candidate_enh=as.data.frame(candidate_enh)
candidate_enh=candidate_enh[candidate_enh$EnhancerID!="NA", ]

# Load enhancer coverage data from Nott et al. (2019) 
coverage_nott= read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Coverage All Peaks (Nott).csv")
colnames(coverage_nott)

# Extract genomic coordinates from enhancer and Nott datasets
coordinates_nott=tibble(chrom=sub(":.*","", coverage_nott$X),
start=as.numeric(gsub("chr[0-9XY]+:([0-9]+)-[0-9]+", "\\1", coverage_nott$X)),
end=as.numeric(gsub("chr[0-9XY]+:[0-9]+-([0-9]+)","\\1", coverage_nott$X)))

coordinates_enh=tibble(chrom=sub(":.*","", candidate_enh$`Peak ID`),
                        start=as.numeric(gsub("chr[0-9XY]+:([0-9]+)-[0-9]+", "\\1",candidate_enh$`Peak ID`)),
                        end=as.numeric(gsub("chr[0-9XY]+:[0-9]+-([0-9]+)","\\1", candidate_enh$`Peak ID`)),
                       id=candidate_enh$EnhancerID)

# Identify overlaps between candidate enhancers and Nott et al. dataset (hg19 format)
intersect=bed_intersect(coordinates_enh,coordinates_nott) 
intersect$enh_lenght=intersect$end.x-intersect$start.x
intersect$prop=intersect$enh_lenght/intersect$.overlap
intersect$coordinates=paste0(intersect$chrom,":",intersect$start.y,"-",intersect$end.y)

# Verify that all enhancers fully overlap with Nott peaks
table(intersect$prop) 

# Filter Nott dataset by enhancer regions
coverage_nott_filtered=coverage_nott[coverage_nott$X %in% paste0(intersect$chrom,":",intersect$start.y,"-",intersect$end.y),]
length(unique(coverage_nott_filtered$X)) 
coverage_nott_filtered$enh=intersect[match(coverage_nott_filtered$X, intersect$coordinates),]$id.x
coverage_nott_sorted=coverage_nott_filtered[match(rownames(kme_sorted), coverage_nott_filtered$enh) ,]

# Verify correct ordering
table(rownames(kme_sorted)==coverage_nott_sorted$enh)
rownames(coverage_nott_sorted)=coverage_nott_sorted$enh
sample_order=c("LHX2", "OLIG2", "PU1", "NEUN")

# Compute statistical tests:
# - Kruskal-Wallis test: Compare coverage values across all cell types
# - Wilcoxon rank-sum test: Compare astrocytes vs. all other cell types

stats=list(kruskal=list(atac=list(), H3K27ac=list(), H3K4me3=list()),
           wilcox.astro=list(atac=list(), H3K27ac=list(), H3K4me3=list())
           )
          
for (set in c("atac", "H3K27ac", "H3K4me3"))
  
{ 
  dataset=coverage_nott_sorted[, paste(sample_order, set, sep="_")]
  for (m in unique(kme_sorted$moduleLabel))
  {
    x=dataset[which(kme_sorted$moduleLabel%in%m),]
    x$enh=rownames(x)
    x=gather(x, key="CT", value="value",-"enh")
    
    # Kruskal-Wallis test across cell types
    stats[["kruskal"]][[set]][[m]]=kruskal.test(value ~ as.factor(CT), data=x) 

    # Wilcoxon test comparing astrocytes vs. all other cell types
    stats[["wilcox.astro"]][[set]][[m]]=wilcox.test(x[grepl("LHX2", x$CT),]$value,
                                                   x[!grepl("LHX2", x$CT),]$value,
                                                   alternative = "greater") 
    }}

# Convert lists of statistics to data frames
for (set in c("atac", "H3K27ac", "H3K4me3"))
{stats[["kruskal"]][[set]]=as.data.frame(do.call("rbind",stats[["kruskal"]][[set]]))
stats[["wilcox.astro"]][[set]]=as.data.frame(do.call("rbind",stats[["wilcox.astro"]][[set]]))
}


# Extract and format p-values for plotting
#Extract p values for ATAC
plot.atac=data.frame(unlist(stats$kruskal$atac$p.value), 
                     unlist(stats$wilcox.astro$atac$p.value)
                     )
colnames(plot.atac)=c("Cell Type", "Astrocyte")
for(j in c(1:ncol(plot.atac))) plot.atac[which(plot.atac[,j]>0.05), j]=NA

#Extract p values for H3K27ac
plot.H3K27ac=data.frame( unlist(stats$kruskal$H3K27ac$p.value),
                     unlist(stats$wilcox.astro$H3K27ac$p.value)
                     )
colnames(plot.H3K27ac)=c("Cell Type", "Astrocyte")


for(j in c(1:ncol(plot.H3K27ac))) plot.H3K27ac[which(plot.H3K27ac[,j]>0.05), j]=NA

##Extract p values for H3K4me3
plot.H3K4me3=data.frame(unlist(stats$kruskal$H3K4me3$p.value),
                        unlist(stats$wilcox.astro$H3K4me3$p.value)
                        )
colnames(plot.H3K4me3)=c("Cell Type", "Astrocyte")
for(j in c(1:ncol(plot.H3K4me3))) plot.H3K4me3[which(plot.H3K4me3[,j]>0.05), j]=NA

# Format column names for heatmap
colnames(plot.atac)=paste(colnames(plot.atac), "(ATAC)",sep=" ")
colnames(plot.H3K27ac)=paste(colnames(plot.H3K27ac), "(H3K27ac)",sep=" ")
colnames(plot.H3K4me3)=paste(colnames(plot.H3K4me3), "(H3K4me3)",sep=" ")

# Combine p-value matrices
pvalue_combined=cbind(plot.atac, plot.H3K27ac,plot.H3K4me3)
rownames(pvalue_combined)=sub(".Error: Within.*","",rownames(pvalue_combined))

## Define colors for modules and chromatin marks
ann=data.frame(Module=rownames(pvalue_combined))
rownames(ann)=ann$Module

module_colors <- setNames(
  brewer.pal(length(unique(ann$Module)), "Dark2")[seq_along(unique(ann$Module))],  
  unique(ann$Module)
)

####set marker types format and color 
ann_mark=data.frame(samples=colnames(pvalue_combined))
ann_mark$Marker=gsub(".*\\(([^)]+)\\).*", "\\1", ann_mark$samples)
rownames(ann_mark)=ann_mark$samples
colnames(ann_mark)[2]="Chromatin mark"
Mark_colors <- setNames(
  colct <- c("#6F083D","#F3CCDE","#F4EAD4"),  
unique(ann_mark$`Chromatin mark`)
)

ann_colors<- list(Module = module_colors, `Chromatin mark` = Mark_colors)

#define pallete
color_ramp <- colorRamp(c("lightblue", "midnightblue"))
palette_stats <- rgb(color_ramp(seq(0, 1, length.out = 1000)), maxColorValue = 255)

###########Heatmap ploting stats: It generates Sup Fig 3D
fontsize=18
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/Astronet_WGCNA/Nott_StatsCoverage_heatmap.pdf", width = 11, height = 9)
pheatmap(t(-log10(pvalue_combined)), 
         color=palette_stats,
         cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_col = ann[1], annotation_colors = ann_colors,
         annotation_row = ann_mark[2],
         fontsize = fontsize, fontsize_row = fontsize, fontsize_col = fontsize,
         gaps_row = c(2,4),
         scale = "none")
dev.off()

###################################################################################################################
#Plot SupFig 3C: Heatmap of enhancer coverage  (ATAC-seq, H3K27ac, H3K4me3) for candidate enhancers across cell types

##Format column names
colnames(coverage_nott_sorted)=gsub("LHX2","Astrocytes",colnames(coverage_nott_sorted))
colnames(coverage_nott_sorted)=gsub("NEUN","Neurons",colnames(coverage_nott_sorted))
colnames(coverage_nott_sorted)=gsub("OLIG2","Oligodendrocytes",colnames(coverage_nott_sorted))
colnames(coverage_nott_sorted)=gsub("PU1","Microglia",colnames(coverage_nott_sorted))

#Provide cell type order displayed in heatmap
sample_order=c("Astrocytes", "Oligodendrocytes", "Microglia", "Neurons")

#Define columns/rows format and colors for heatmap
rownames(coverage_nott_sorted)=coverage_nott_sorted$enh
colheat=colorRampPalette(c("#B51A00", "#FFE2D6","white"))
ann=as.data.frame(kme_sorted[, c(2,1)])

## Define colors for modules and cell type 
module_colors <- setNames(
  brewer.pal(length(unique(ann$moduleLabel)), "Dark2")[seq_along(unique(ann$moduleLabel))],  
  unique(ann$moduleLabel)
)

ann_cell=data.frame(samples=colnames(coverage_nott_sorted)[2:13])
ann_cell$CellType=sub("_.*","",ann_cell$samples)
rownames(ann_cell)=ann_cell$samples

# Order  data frame based on factor levels
ann_cell$CellType <- factor(ann_cell$CellType, levels = sample_order)
ann_cell <- ann_cell[order(ann_cell$CellType), ]

cell_colors <- setNames(
  colct <- c("#6F083D", "#F3CCDE", "#F4EAD4", "#3A5F9A"),  
  unique(ann_cell$CellType)
)
ann_colors<- list(moduleLabel = module_colors, CellType = cell_colors)


###########Generates Sup Fig 3C
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/Astronet_WGCNA/Nott_Coverage_heatmap.pdf", height = 7, width = 11)
pheatmap(as.matrix(coverage_nott_sorted[, c(paste0(sample_order, "_atac"), 
                                            paste0(sample_order, "_H3K27ac"),
                                            paste0(sample_order, "_H3K4me3"))]), 
         scale="row",
         color=rev(colheat(10)),
         cluster_cols = FALSE, cluster_rows = FALSE,
         show_rownames=FALSE, show_colnames=FALSE,
         annotation_row = ann[1], annotation_col = ann_cell[2], 
         annotation_colors = ann_colors,
         fontsize = fontsize, fontsize_row = fontsize, fontsize_col = fontsize,
         gaps_col = c(4,8))
dev.off()


###Enhancers from modules 1,2 and 3 accounted for % out of the total putative enhancers
round(nrow(kme_sorted[kme_sorted$moduleLabel %in% c("ME1","ME2", "ME3", "ME6"),])/nrow(kme_sorted),digits=2)
