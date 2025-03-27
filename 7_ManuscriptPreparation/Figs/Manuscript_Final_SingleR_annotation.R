#This script generates figure 1C
# ===============================
# @date: 25-03-27
# This script generates one barplot and three UMAP plots.
# The barplot is a summary of the annotation of NHA CRISPRi cells at the levels of developmental stage and cell type.
# The three UMAPs present the annotation of the same cells at the levels of developmental stage, cell type, and cell cycle stage.
# ===============================

# =================== Setup ===================
setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_LibraryDesign/")
source("../../../FullScale/Scripts/Functions.R")
source("../FinalFigureFunctions.R")

load("../../../FullScale/Data/Preprocessed/NHA Pooled (Final).rda")

outputdir<-"/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/JW_plots/Fig1_250327/"

# Traceable PDF export function
pdf_LibDesign <- function(figNo, title, h, w) {
  pdf(file = paste0(outputdir, figNo, "_Script LibraryDesign_", title, ".pdf"), height = h, width = w) # JW version
}


# =================== Load pseudobulks from Herring 2022 ===================
load("../../../PublicData/snRNAseq/Herring2022_MaturationBrain/RNA/Processed/Pseudobulk_byGJS.rda", verbose = TRUE)

# Pool samples within each stage/cell-type combination
pb$Meta$Category <- paste0(pb$Meta$Celltype, "_", pb$Meta$Stage)
tc <- list()
for (j in unique(pb$Meta$Category)) { tc[[j]] <- rowSums(pb$Exp[,which(pb$Meta$Category == j)])  }
tc <- do.call("cbind", tc)
tc <- apply(tc, 2, function(x) { x / (sum(x) / 10^6) })
tc <- as.data.frame(tc)

# =================== Generate/load singleR related subjects ===================
common <- intersect(rownames(nha), rownames(tc))
tc <- tc[common,]  
load("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/Fig_LibraryDesign/SingleR.rda")                           

# =================== Fig. 1C. Top Barplot: A summary of the singleR statistics ===================
p <- singleR_annots$pruned.labels
p[is.na(p)] <- "NA_NA"
p <- factor(p, levels = c(colnames(tc), "NA_NA"))

q <- as.data.frame(table(p))
q$Stage <- splitter(q$p, "_", 2)
q$Stage <- factor(q$Stage, levels = unique(q$Stage)[c(5,1,3,2,4,6,7)])
levels(q$Stage)[7] <- "Not assigned"
q$Celltype <- splitter(q$p, "_", 1)
q$Celltype <- factor(q$Celltype, levels = unique(q$Celltype)[c(3,1,2,4,5,6,7)])
levels(q$Celltype) <- gsub("Astro", "Astrocyte", levels(q$Celltype)) %>%
  gsub("Exc", "Excitatory\nNeurons", .) %>%
  gsub("Inh", "Inhibitory\nNeurons", .) %>%
  gsub("Micro", "Microglia", .) %>%
  gsub("Oligo", "Oligodendrocytes", .) 
q$Label <- q$Freq
q$Label[which(q$Label == 0)] <- NA

pal <- c("#522569", "#443782", "#3D6B93", "#2F988C", "#98CC42", "#FAE51A", "grey50")

q$Freq.ori<-q$Freq
q$Freq<-log((q$Freq + 1),10)


## Barplot version 1: Acquire the barplot in ideal proportion.
pdf_LibDesign(figNo = "Fig1_Top", title = "SingleR_Barplot", h = 2, w = 4)
ggplot(q, aes(x = Celltype, fill = Stage, y = Freq, colour = Stage)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = Label), position = position_dodge(width = 1), vjust = -0.25,
            hjust = 0.4, size = 2.5, show.legend = FALSE) +
  theme_bw() +
  facet_grid(.~Celltype, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
        # legend.position = "right",
        legend.position = "none",
        # legend.position = c(0.6, 0.7), 
        strip.background = invis,
        axis.text.x = invis, axis.ticks.x = invis) +
  scale_x_discrete(expand = c(0.15,0)) +
  labs(y = "log10(Number of cells + 1)", x = "")
dev.off()


## Barplot version 2: Acquire legend in the optimized format
pdf_LibDesign(figNo = "Fig1_Top", title = "SingleR_Barplot_legend", h = 2, w = 4)
ggplot(q, aes(x = Celltype, fill = Stage, y = Freq, colour = Stage)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = Label), position = position_dodge(width = 1), vjust = -0.25,
            hjust = 0.4, size = 2.5, show.legend = FALSE) +
  theme_bw() +
  facet_grid(.~Celltype, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
        legend.position = "right",
        # legend.position = "none",
        # legend.position = c(0.6, 0.7), 
        strip.background = invis,
        axis.text.x = invis, axis.ticks.x = invis) +
  scale_x_discrete(expand = c(0.15,0)) +
  labs(y = "log10(Number of cells + 1)", x = "")
dev.off()




# =================== Fig. 1C. Bottom Middle UMAP: Annotation by Developmental stage ===================
#### RunPCA using a smaller dimension and lower resolution
nha.15 <- RunPCA(nha, features = VariableFeatures(object = nha), npcs = 15)
nha.15 <- FindNeighbors(nha.15, dims = 1:15)
nha.15 <- FindClusters(nha.15, resolution = 0.05)
nha.15 <- RunUMAP(nha.15, dims = 1:15)

#### Incorporate singleR annotation
nha.15$SingleR<-NA
nha.15$SingleR<- singleR_annots$pruned.labels[match(rownames(nha.15@meta.data),rownames(singleR_annots))]
saveRDS(nha.15, file = paste0(outputdir,"nha.15.rds"))

#### Prepare dataframe for customized DimPlots
p.ori <- data.frame(cbind(nha.15@reductions$umap@cell.embeddings[,c(1:2)]))
length(which((rownames(p.ori) == rownames(nha.15@meta.data))%in%"TRUE")) # [1] 47577
p.ori$SingleR <- nha.15$SingleR
na.rows<-which(is.na(p.ori$SingleR)%in%"TRUE")
p.ori$SingleR[na.rows] <- "Not assigned"
p.ori$Dev_stage<- gsub("^.*_","", p.ori$SingleR)
p.ori$Cell_Type<- gsub("_.*$","", p.ori$SingleR)
p.ori$alpha <- p.ori$SingleR=="Astro_Fetal"
p.ori$Astro <- p.ori$Cell_Type=="Astro"
p.ori$Seurat_cluster<-nha.15$seurat_clusters
p.ori$Cell_Cycle<-nha.15$Cycle_Seurat_Phase

pdf(paste0(outputdir,"Fig1_Bottom_Middle_UMAP_DevStage.pdf"), height = 3, width = 4)
p<-p.ori
Astro_Fetal_row<-which(p$SingleR%in%"Astro_Fetal")
p<-rbind(p[Astro_Fetal_row,], p[na.rows,],p[-c(Astro_Fetal_row,na.rows),])

ggplot(p, aes(x=umap_1, y = umap_2, fill = Dev_stage, colour = Dev_stage, alpha=0.4)) +
  geom_point(aes(fill = Dev_stage), size=1) +
  geom_point(aes(colour=Dev_stage), pch = 21, size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank(), legend.text = (element_text(size = 12)), legend.title = (element_text(size = 14))) +
  guides(colour = guide_legend(override.aes = list(pch=21, colour = c("#98CC42","#522569","#3D6B93","#443782","black"), size = 5, fill=c("#98CC42","#522569","#3D6B93","#443782","white")))) + # This line re-specify the size of the legend.
  scale_colour_manual(values= c("#98CC42","#522569","#3D6B93","#443782","black"))+
  scale_fill_manual(values= c("#98CC42","#522569","#3D6B93","#443782","white"))+
  labs(title = "Developmental stage distribution") +
  guides(alpha= FALSE)
dev.off()


# =================== Fig. 1C. Bottom Left UMAP: Annotation by Cell Type ===================
pdf(paste0(outputdir,"Fig1_Bottom_Left_UMAP_CellType.pdf"), height = 3, width = 4)
p<-p.ori
Astro_row<-which(p$Cell_Type%in%"Astro")
Micro_row <-which(p$Cell_Type%in%"Micro")
p<-rbind(p[Astro_row,],p[-c(Astro_row,Micro_row),],p[Micro_row,])

ggplot(p, aes(x=umap_1, y = umap_2, fill = Cell_Type, colour = Cell_Type, alpha=0.4)) +
  geom_point(aes(fill = Cell_Type), size=1) +
  geom_point(aes(colour=Cell_Type), pch = 21, size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank(), legend.text = (element_text(size = 12)), legend.title = (element_text(size = 14))) +
  guides(colour = guide_legend(override.aes = list(pch=21, colour = c("#FAC858","#484848","black"), size = 5, fill=c("#FAC858","#484848","white")))) + # This line re-specify the size of the legend.
  scale_colour_manual(values= c("#FAC858","#484848","black"))+
  scale_fill_manual(values= c("#FAC858","#484848","white"))+
  labs(title = "Cell Type distribution") +
  guides(alpha= FALSE)
dev.off()


# =================== Fig. 1C. Bottom Left UMAP: Annotation by Cell Cycle genes ===================
pdf(paste0(outputdir,"Fig1_Bottom_Right_UMAP_CellCycle.pdf"), height = 3, width = 4)
p<-p.ori

# ggplot(p, aes(x=umap_1, y = umap_2, fill = Cell_Cycle, colour = Cell_Cycle, alpha=0.05)) +
ggplot(p, aes(x=umap_1, y = umap_2, fill = Cell_Cycle, colour = Cell_Cycle)) +
  geom_point(aes(fill = Cell_Cycle), size=1) +
  geom_point(aes(colour=Cell_Cycle), pch = 21, size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank(), legend.text = (element_text(size = 12)), legend.title = (element_text(size = 14))) +
  guides(colour = guide_legend(override.aes = list(pch=21, size = 5))) + # This line re-specify the size of the legend.
  scale_fill_brewer(type="qual", palette="Set2") +
  scale_colour_brewer(type="qual", palette="Set2") +
  labs(title = "Cell Cycle Phase distribution") +
  guides(alpha= FALSE)
dev.off()

