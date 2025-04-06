# This script generates figure XXX
# ===============================
# @date: 25-03-27
# This script clusters single cells into two clsuters, 
# detects the marker of the minor cluser (referred to as the island cluster hereafter), 
# ID GO terms enriched in the island cluster,
# test for the association between MOI and the cluster identity of cells,
# and test the enrichment of different gRNA treatments in the island clsuter.
# ===============================

# ========== Load required packages ==========
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(purrr)

# ========== Set working directory and output folder ==========
dir.create("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/JW_plots/Rebuttal/Final_data")
setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/JW_plots/Rebuttal/Final_data/")

# ========== Load Seurat object ==========
nha.15 <- readRDS("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/JW_plots/Fig1/JW_nha.15.rds")

# ========== Re-cluster at lower resolution ==========
nha.15.res.new <- FindClusters(nha.15, resolution = 0.01)
nha.15.res.new <- RunUMAP(nha.15.res.new, dims = 1:15)
Idents(nha.15.res.new) <- "RNA_snn_res.0.01"

# ========== Save UMAP plot ==========
pdf("NHA.dim15.resPointZero1.pdf", width = 4, height = 4)
DimPlot(nha.15.res.new)
dev.off()

save(nha.15.res.new, file = "nha.15.res.01.rda")

# ========== Identify island-specific markers ==========
mk <- FindAllMarkers(nha.15.res.new)
write.csv(mk, "Markers_small_island.csv")

# ========== GO enrichment analysis ==========
mk <- mk[sign(mk$avg_log2FC) > 0 , ]
mk <- mk[mk$cluster == 1, ]
mk$diff <- mk$pct.1 - mk$pct.2
mk$rankBy <- mk$diff + mk$pct.1

set <- mk$gene[which((mk$pct.1 > 0.5) & (mk$avg_log2FC > 1))]
bkg <- rownames(nha.15.res.new@assays$RNA)

ont.list <- c("BP", "MF", "CC")
ont.results <- list()
for (i in c(1:length(ont.list))) {
  print(i)
  ont <- ont.list[i]
  
  ont.results[[i]] <- enrichGO(set, OrgDb = "org.Hs.eg.db", 
                               keyType = "SYMBOL", ont = ont, 
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                               minGSSize = 10, maxGSSize = 500, 
                               universe = bkg)
  
  names(ont.results)[i] <- ont
}

save(ont.results, file = "Markers_small_island_GOpct_IV.rda")
bp <- ont.results$BP@result[which(ont.results$BP@result$p.adjust < 0.05), ]
mf <- ont.results$MF@result[which(ont.results$MF@result$p.adjust < 0.05), ]
cc <- ont.results$CC@result[which(ont.results$CC@result$p.adjust < 0.05), ]
write.csv(bp, "Markers_small_island_GOpct_BP_IV.csv")

# ========== Test for MOI enrichment in the island cluster ==========
meta.moi <- nha.15.res.new@meta.data[, c("MOI", "seurat_clusters")]

pval <- t.test(meta.moi$MOI ~ meta.moi$seurat_clusters)$p.value
pval <- formatC(pval, digits = 2)

pdf("MOI_Violin_IV_JW.pdf")
ggplot(meta.moi, aes(y = MOI, x = seurat_clusters, fill = seurat_clusters)) +
  geom_violin(trim = FALSE) +  
  geom_boxplot(width = 0.15, fill = "white") + 
  theme_minimal() +
  labs(
    subtitle = paste0("p = ", pval, " (Welch Two Sample t-test)")
  ) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0, face = "italic", size = 10),
    plot.margin = margin(t = 15, r = 20, b = 10, l = 10)
  )
dev.off()

# ========== Fisher's test: gRNA enrichment in island ==========
nha.15.res.new$Island <- FALSE
nha.15.res.new$Island[which((nha.15.res.new$RNA_snn_res.0.01 == 1) %in% "TRUE")] <- TRUE
Enh_gRNA_col <- seq(24, 4155, 1)

fisher.table <- as.data.frame(cbind(colnames(nha.15.res.new@meta.data)[Enh_gRNA_col],
                                    colnames(nha.15.res.new@meta.data)[Enh_gRNA_col],
                                    colnames(nha.15.res.new@meta.data)[Enh_gRNA_col]))
colnames(fisher.table) <- c("gRNA_Enh", "Fisher_pval", "Fisher_Odds_Ratio")

for (i in c(1:length(Enh_gRNA_col))) {
  print(i)
  col_ID <- Enh_gRNA_col[i]
  
  if (length(unique(nha.15.res.new@meta.data[, col_ID])) > 1) {
    pval <- as.numeric(fisher.test(nha.15.res.new$Island , nha.15.res.new@meta.data[, col_ID])["p.value"])
    OR <- as.numeric(fisher.test(nha.15.res.new$Island , nha.15.res.new@meta.data[, col_ID])[["estimate"]][["odds ratio"]])
    
    fisher.table$Fisher_pval[i] <- pval
    fisher.table$Fisher_Odds_Ratio[i] <- OR
  }
}

write.csv(fisher.table, file = "fisher.table_whole.csv")

# ========== Fisher's test: gRNA category (enhancer-targetting, promoter-targetting, and negative control) enrichment in the island cluster ==========
ft <- fisher.table
res <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")

ft$Categ <- list_transpose(strsplit(ft$gRNA_Enh, split = "_"))[[1]]
ft$Enh <- list_transpose(strsplit(ft$gRNA_Enh, split = "_"))[[1]]

ft$Categ[grep("Enh", ft$Categ)] <- "Enh"
ft <- ft[-grep("_g", ft$gRNA_Enh), ]
ft$Fisher_pval_adj <- p.adjust(ft$Fisher_pval, method = "BH")

sig <- which(ft$Fisher_pval_adj < 0.05)

# Run Fisher tests for each category
Enh.fish<-fisher.test(ft$Categ == "Enh", ft$Fisher_pval_adj < 0.05)
Neg.fish<-fisher.test(ft$Categ == "Neg", ft$Fisher_pval_adj < 0.05)
Pos.fish<-fisher.test(ft$Categ == "Pos", ft$Fisher_pval_adj < 0.05)

ft$neg_pVal_adj <- -(ft$Fisher_pval_adj)

# ========== Visualizing promoter-targeting gRNA enrichment via violin plots ==========
ft$Fisher_Odds_Ratio <- as.numeric(ft$Fisher_Odds_Ratio)

# Prepare label dataframe for violin plot
Enh.p <- formatC(Enh.fish$p.value, digits = 2)
Enh.or <- formatC(Enh.fish$estimate, digits = 2)
Neg.p <- formatC(Neg.fish$p.value, digits = 2)
Neg.or <- formatC(Neg.fish$estimate, digits = 2)
Pos.p <- formatC(Pos.fish$p.value, digits = 2)
Pos.or <- formatC(Pos.fish$estimate, digits = 2)

odds_ratios <- c(paste(Enh.or, "\np = ", Enh.p), paste(Neg.or, "\np = ", Neg.p), paste(Pos.or, "\np = ", Pos.p))
names(odds_ratios) <- c("Enh", "Neg", "Pos")

labels_df <- data.frame(
  Categ = names(odds_ratios),
  neg_pVal_adj = max(ft$neg_pVal_adj, na.rm = TRUE) + 0.5,
  label = paste0("OR = ", odds_ratios)
)

# Final plot
pdf("island_enrichment_in_gRNA_treatments_Fishers.pdf")
ggplot(ft, aes(y = neg_pVal_adj, x = Categ, fill = Categ)) +
  geom_violin(trim = FALSE) +  
  geom_boxplot(width = 0.15, fill = "white") +
  geom_text(data = labels_df, aes(label = label), size = 3, fontface = "italic") +
  theme_minimal() +
  labs(
    x = NULL,
    y = "-pVal.adjusted \n(Fisher's test)",
    caption = paste0("Association with the small island for each individual gRNA \ncategories is indicated at the top of each violin bar.")
  ) +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 15, r = 20, b = 10, l = 10),
    plot.caption = element_text(hjust = 0, face = "italic", size = 10)
  ) +
  ggtitle("Association between individual gRNAs and the small island.")
dev.off()
