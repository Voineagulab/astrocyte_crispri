# PCA plot comparing ATAC-seq data from:
# - NHAs (current study)
# - Pseudobulked snATAC-seq data (Herring et al., 2022)
# - Bulk ATAC-seq of immunopanned cell types (Nott et al., 2019)

# Load necessary libraries
library(UpSetR)
library(purrr)
library(ggplot2)

# Load Nott et al. 2019 ATAC-seq data
nott <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Coverage All Peaks (Nott).csv")
rownames(nott) <- nott$X
nott <- nott[,-1]
nott <- nott[, grep("atac", colnames(nott))]

# Load Herring et al. 2022 ATAC-seq data
her <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Coverage All Peaks (Herring).csv")
rownames(her) <- her$X
her <- her[,-1]

# Load NHA ATAC-seq data 
nha <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Chromatin/Coverage/Inhouse.csv")
rownames(nha) <- nha$Id

# Combine Nott and Herring data, matching by peaks
atac <- data.frame(
  her, 
  nott[match(rownames(her), rownames(nott)),],
  NHA_atac = nha$CoverageMean[match(rownames(her), rownames(nha))]
)

# Filter peaks with low coverage across samples
coverage_threshold <- 2
atac_filtered <- atac[rowSums(atac > coverage_threshold) > 2, ]

# Perform PCA on inter-sample Spearman correlation matrix
pca_result <- princomp(cor(atac_filtered, method = "spearman"))
pca_data <- as.data.frame(pca_result$loadings[, 1:5])

# Extract sample information from row names
sample_info <- list_transpose(strsplit(rownames(pca_data), "\\.", fixed = TRUE))[[1]]
sample_info <- gsub("_3", "3", sample_info)
sample_info <- gsub("_6", "6", sample_info)

# Define Stage and Cell type (CT)
pca_data$Stage <- list_transpose(strsplit(sample_info, "_", fixed = TRUE))[[2]]
pca_data$CT <- list_transpose(strsplit(sample_info, "_", fixed = TRUE))[[1]]

# Merge L23","L4","L56" neuronal subtypes
pca_data$CT[pca_data$CT %in% c("L23","L4","L56")] <- "Neurons"

# Define NHA samples explicitly
pca_data$CT[grep("NHA", rownames(pca_data))] <- "NHA"
pca_data$Stage[grep("NHA", rownames(pca_data))] <- "Fetal"

# Remove non-informative CGE/MGE-derived samples
pca_data <- pca_data[-grep("der", pca_data$Stage), ]

# label NHA stage as "Fetal (this study)"
pca_data$Stage <- ifelse(
  pca_data$CT == "NHA", 
  "Fetal (this study)",
  pca_data$Stage
)

# Merge developmental stages to simplify labels
pca_data$Stage <- ifelse(
  pca_data$Stage %in% c("Childhood", "Infancy"), "Infancy/Childhood (Herring et al., 2022)",
  ifelse(pca_data$Stage %in% c("Adult", "Adolescence"), "Adolescence/Adult (Herring et al., 2022)",
         ifelse(pca_data$Stage == "atac", "Postnatal (Nott et al., 2019)",
                ifelse(pca_data$Stage == "Neonatal", "Neonatal (Herring et al., 2022)",
                       ifelse(pca_data$Stage == "Fetal", "Fetal (Herring et al., 2022)", pca_data$Stage)
                )
         )
  )
)

# Merge closely related cell-type labels for consistency
pca_data$CT <- ifelse(
  pca_data$CT %in% c("Astro", "LHX2"), "Astro",
  ifelse(pca_data$CT %in% c("Micro", "PU1"), "Micro",
         ifelse(pca_data$CT %in% c("NEUN", "Neurons"), "Neurons",
                ifelse(pca_data$CT %in% c("OLIG2", "Oligo"), "Oligo", pca_data$CT)
         )
  )
)

# Calculate percentage variance explained by PCs
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
var_PC1 <- round(variance_explained[1] * 100, 0)
var_PC2 <- round(variance_explained[2] * 100, 0)

# Generate PCA plot corresponding to SFig 1b
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/SCRIPTS/PCA_NHA_ATAC validation.pdf", width = 7, height = 5)
ggplot(pca_data, aes(x = Comp.1, y = Comp.2, color = CT, shape = Stage)) +
  geom_point(size = 2.5) +
  theme_bw() +
  xlab(paste0("PC1: ", var_PC1, "% variance explained")) +
  ylab(paste0("PC2: ", var_PC2, "% variance explained")) +
  labs(color = "Cell Type", shape = "Stage - Study")
dev.off()


############
# 
# for (j in c(2:10))
# {th=j
# if (j==2) plot(x=c(1:4), y=apply(atac[, c(1:4)], 2, nth)/1000, pch=20, 
#      col=c("red", "seagreen", "steelblue", "midnightblue"),
#      ylim=c(1,60))
# else points(x=c(1:4), y=apply(atac[, c(1:4)], 2, nth)/1000, pch=20, 
#                col=c("red", "seagreen", "steelblue", "midnightblue"),
#                ylim=c(1,60))
# 
# }
# th=3
# atac_binary=atac
# for (j in c(1:ncol(atac))) atac_binary[,j]=ifelse(atac[,j]>th, 1,0)
# upset(atac_binary, sets = colnames(atac_binary), order.by = "freq")
# 
# th=10
# atac_th=atac
# for (j in c(1:ncol(atac))) atac_th[which(atac[,j]<th),j]=NA
# 
# cormx=cor(atac_th, method = "s", use="pairwise.complete")
# cormx
# pheatmap(cormx)
# 
# atac_plot=atac[rowSums(atac[,]>th)>1 ,]
# pheatmap(cor(log2(atac_plot+1), method="s"), scale="row")
# 
# 
# pca_result <-princomp( atac_th)
# #pca_result <- princomp( log2(atac_plot+1))
# 
# # Prepare data for plotting
# pca_data <- as.data.frame(pca_result$loadings[, 1:5])
# pca_data$CT=gsub("_atac", "", rownames(pca_data))
# 
# pca_plot <- ggplot(pca_data, aes(x = Comp.1, y = Comp.2, color = CT)) +
#   geom_point(size = 3) +
#   theme_bw() 
# print(pca_plot)
# 
# 
# atac_plot=atac[rowSums(atac[,]>th)>1 ,]
# pca_result <- prcomp(cor(atac_plot[,], method="s"))
# pca_data <- as.data.frame(pca_result$rotation[, 1:3])
# pca_data$CT=gsub("_atac", "", rownames(pca_data))
# pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = CT)) +
#   geom_point(size = 3) +
#   theme_bw() 
# print(pca_plot)
# 
