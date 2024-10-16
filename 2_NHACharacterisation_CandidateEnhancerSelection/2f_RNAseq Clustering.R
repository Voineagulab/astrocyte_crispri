library(readxl)
library(dplyr)
library(purrr)
library(ggplot2)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq")
rm(list=ls())
## Comparison to the Zhang 2016 data
zhang_meta <- read_xlsx("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/Zhang2016_ImmunpanningHumanAst/1-s2.0-S0896627315010193-mmc3.xlsx", sheet = "Human data only", range = "A1:AP3", col_names = FALSE)
zhang_meta <- as.data.frame(zhang_meta)
# rownames(zhang_meta) <- zhang_meta$...1
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

nha=read.delim("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/NHAcharacterisation_resubmission/RNA_seq/star_salmon/salmon.merged.gene_counts.tsv", sep="\t", header=TRUE)
common=intersect(nha$gene_name, rownames(zhang_exp))
nha=nha[match(common, nha$gene_name) ,]
rownames(nha)=nha$gene_name
nha_df <- do.call(cbind, nha[, -c(1,2)]) # just dealing with the annoying list type of nha
rownames(nha_df)=rownames(nha)
nha_df<-nha_df[, grep("2019", colnames(nha_df))]
#nha_df <- nha_df[-which(apply(nha_df, 1, var) == 0),]
colnames(nha_df)=gsub("_GRE13789A15_22KTVVLT3", "", colnames(nha_df))
colnames(nha_df)=gsub("_GRE13789A13_22KTVVLT3", "", colnames(nha_df))

counts=cbind(zhang_exp_df[match(common, rownames(zhang_exp)) , ],
                  nha_df[match(common, rownames(nha_df)) , ] )
# batch=c(rep(1, ncol(zhang_exp)) , rep(2, ncol(nha_df)))  
# adjusted_counts <- ComBat_seq(counts = counts, batch=batch, full_mod = FALSE)
# pheatmap(cor(counts, method="s"))
# pheatmap(cor(adjusted_counts, method="s"))
# var<- apply(adjusted_counts, 1, var)
# mean<- apply(adjusted_counts, 1, mean)
# cv=var/mean
# top_indices <- order(cv, decreasing = TRUE)[1:1000]

pca_result <- princomp( cor(counts, method="s"))

# Prepare data for plotting
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
# % variance explained as per summary(pca_result). Should be possible to code it.

summary(pca_result)
# Importance of components:
#   Comp.1    Comp.2     Comp.3     Comp.4     Comp.5
# Standard deviation     0.4463605 0.2127772 0.17729709 0.13110278 0.10169488
# Proportion of Variance 0.6053845 0.1375656 0.09551313 0.05222563 0.03142375
# Cumulative Proportion  0.6053845 0.7429502 0.83846331 0.89068894 0.92211269
# Comp.6     Comp.7      Comp.8      Comp.9     Comp.10
# Standard deviation     0.07526673 0.06450538 0.056880202 0.048233724 0.040516164
# Proportion of Variance 0.01721337 0.01264305 0.009830648 0.007069057 0.004987886
# Cumulative Proportion  0.93932606 0.95196911 0.961799758 0.968868815 0.973856701
# Comp.11     Comp.12     Comp.13     Comp.14
# Standard deviation     0.036984590 0.033765879 0.030634831 0.028961960
# Proportion of Variance 0.004156248 0.003464304 0.002851615 0.002548683
# Cumulative Proportion  0.978012949 0.981477253 0.984328868 0.986877551
# Comp.15     Comp.16     Comp.17     Comp.18
# Standard deviation     0.025046235 0.020745407 0.019098948 0.018648697
# Proportion of Variance 0.001906096 0.001307687 0.001108355 0.001056713
# Cumulative Proportion  0.988783647 0.990091333 0.991199688 0.992256401
# Comp.19      Comp.20      Comp.21     Comp.22
# Standard deviation     0.018299072 0.0175020103 0.0143884759 0.014123204
# Proportion of Variance 0.001017462 0.0009307558 0.0006290562 0.000606075
# Cumulative Proportion  0.993273862 0.9942046180 0.9948336742 0.995439749
# Comp.23      Comp.24      Comp.25      Comp.26
# Standard deviation     0.0133521877 0.0130356034 0.0122799666 0.0109271578
# Proportion of Variance 0.0005417073 0.0005163238 0.0004581991 0.0003628059
# Cumulative Proportion  0.9959814566 0.9964977804 0.9969559795 0.9973187853
# Comp.27      Comp.28      Comp.29      Comp.30
# Standard deviation     0.010413825 0.0100578692 0.0095743235 0.0091796473
# Proportion of Variance 0.000329519 0.0003073773 0.0002785326 0.0002560424
# Cumulative Proportion  0.997648304 0.9979556817 0.9982342143 0.9984902567
# Comp.31      Comp.32      Comp.33      Comp.34
# Standard deviation     0.0084193592 0.0079855728 0.0072155074 0.0068046327
# Proportion of Variance 0.0002153862 0.0001937635 0.0001581953 0.0001406919
# Cumulative Proportion  0.9987056429 0.9988994064 0.9990576017 0.9991982937
# Comp.35      Comp.36      Comp.37      Comp.38
# Standard deviation     0.0066240184 0.0062763295 0.0060414977 0.0057990669
# Proportion of Variance 0.0001333223 0.0001196937 0.0001109045 0.0001021824
# Cumulative Proportion  0.9993316160 0.9994513097 0.9995622142 0.9996643966
# Comp.39      Comp.40      Comp.41      Comp.42
# Standard deviation     5.661842e-03 5.364165e-03 0.0052242578 4.725103e-03
# Proportion of Variance 9.740368e-05 8.743075e-05 0.0000829295 6.783946e-05
# Cumulative Proportion  9.997618e-01 9.998492e-01 0.9999321605 1.000000e+00
# Comp.43
# Standard deviation           0
# Proportion of Variance       0
# Cumulative Proportion        1


