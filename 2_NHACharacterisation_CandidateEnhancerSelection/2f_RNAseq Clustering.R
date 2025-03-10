## Comparison of RNAseq from NHAs with RNA-seq data from immuno-panned astrocytes from Zhang et al. 2016
library(readxl)
library(dplyr)
library(purrr)
library(ggplot2)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq")
rm(list=ls())

## Read in and process metadata from Zhang et al. 2016
zhang_meta <- read_xlsx("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/Zhang2016_ImmunpanningHumanAst/1-s2.0-S0896627315010193-mmc3.xlsx", sheet = "Human data only", range = "A1:AP3", col_names = FALSE)
zhang_meta <- as.data.frame(zhang_meta)

zhang_meta <- t(zhang_meta[,-1]) %>% as.data.frame()
colnames(zhang_meta) <- c("Celltype", "Age", "Gender") # the cell-type column is an attempt to read in a merged cell. here, the first cell of the merge is given the value, and the rest are left NA
zhang_meta$Celltype_Fixed <- "."
for (j in 1:nrow(zhang_meta)) {
  if (is.na(zhang_meta$Celltype[j])) {
    x <- which(!(is.na(zhang_meta$Celltype[1:j]))) %>% max()
    zhang_meta$Celltype_Fixed[j] <- zhang_meta$Celltype[x]
  } else {
    zhang_meta$Celltype_Fixed[j] <- zhang_meta$Celltype[j]
  }
}

zhang_meta$Stage <- "Adult"
zhang_meta$Stage[grep("fetal", zhang_meta$Celltype_Fixed)] <- "Fetal"

zhang_meta$Celltype_Fixed <- gsub("Human ", "", zhang_meta$Celltype_Fixed) %>%
  gsub("/Macrophage", "", .) %>%
  gsub("fetal ", "", .) %>%
  gsub("mature ", "", .) %>%
  gsub("-", "", .) %>%
  gsub("whole ", "", .) %>%
  gsub("hippocampi ", "", .) %>%
  gsub("GBM / ", "", .) %>%
  gsub(" ", "", .)
zhang_meta$Sample <- paste0(zhang_meta$Stage, "_", zhang_meta$Celltype_Fixed)

## Read in and process expression data from Zhang et al. 2016
zhang_exp <- read_xlsx("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq//PublicData/Zhang2016_ImmunpanningHumanAst/1-s2.0-S0896627315010193-mmc3.xlsx", sheet = "Human data only", skip = 3, col_names = FALSE)
dup <- table(zhang_exp$...1) %>% .[. > 1] %>% names() # for duplicated symbols, remove all instances
zhang_exp <- zhang_exp[-which(zhang_exp$...1 %in% dup),]
zhang_exp <- as.data.frame(zhang_exp)
rownames(zhang_exp) <- zhang_exp[,1]
zhang_exp <- zhang_exp[,-1]
colnames(zhang_exp) <- zhang_meta$Sample

#zhang_exp <- zhang_exp[-which(apply(zhang_exp, 1, var) == 0),]
zhang_exp_df <- do.call(cbind, zhang_exp[, ])
rownames(zhang_exp_df)=rownames(zhang_exp)

# Read in RNA-seq data from NHAs
nha=read.delim("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/NHAcharacterisation_resubmission/RNA_seq/star_salmon/salmon.merged.gene_counts.tsv", sep="\t", header=TRUE)

# Select genes present in both datasets
common=intersect(nha$gene_name, rownames(zhang_exp))
nha=nha[match(common, nha$gene_name) ,]

#Format rownames and colnames
rownames(nha)=nha$gene_name
nha_df <- do.call(cbind, nha[, -c(1,2)]) # just dealing with the annoying list type of nha
rownames(nha_df)=rownames(nha)
nha_df<-nha_df[, grep("2019", colnames(nha_df))]

colnames(nha_df)=gsub("_GRE13789A15_22KTVVLT3", "", colnames(nha_df))
colnames(nha_df)=gsub("_GRE13789A13_22KTVVLT3", "", colnames(nha_df))

# merge the two datasets
counts=cbind(zhang_exp_df[match(common, rownames(zhang_exp)) , ],
                  nha_df[match(common, rownames(nha_df)) , ] )

# PCA plot shown in Fig1.
pca_result <- princomp( cor(counts, method="s"))
pca_data <- as.data.frame(pca_result$loadings[, 1:5])

pca_data$sample=list_transpose(strsplit(rownames(pca_data), ".", fixed=TRUE))[[1]]
pca_data$Stage<- list_transpose(strsplit(pca_data$sample, "_", fixed=TRUE))[[1]]
pca_data$CT<- list_transpose(strsplit(pca_data$sample, "_", fixed=TRUE))[[2]]

pca_data$CT[grep("NHA", rownames(pca_data))]<- "NHA"
pca_data$Stage[grep("NHA", rownames(pca_data))]<- "Fetal"

pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/IV_RNAseqClustering.pdf",
    height=3.2, width=4.5)

pca_plot <- ggplot(pca_data, aes(x = Comp.1, y = Comp.2, color = CT, shape=Stage)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
       x = "PC1: 60%",
       y = "PC2: 13%")
print(pca_plot)

dev.off()

# % variance explained as per summary(pca_result). 
summary(pca_result)
# Importance of components:
#                         Comp.1    Comp.2     Comp.3     Comp.4     Comp.5
# Standard deviation     0.4463605 0.2127772 0.17729709 0.13110278 0.10169488
# Proportion of Variance 0.6053845 0.1375656 0.09551313 0.05222563 0.03142375
# Cumulative Proportion  0.6053845 0.7429502 0.83846331 0.89068894 0.92211269
