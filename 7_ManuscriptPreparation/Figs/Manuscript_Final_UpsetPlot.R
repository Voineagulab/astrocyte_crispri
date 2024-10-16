rm(list=ls())
data=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Final_List/NHA_Peaks_Annotated.csv")

library(UpSetR)

data$CRISPRi=TRUE
data=data[, c("Full.Brain", "Rizz.Brain","Song.Cult","Trev.HCS", "CRISPRi")]
colnames(data)=c("Fullard_Brain NeuN-", "Rizzardi_Brain NeuN-","Song_iPSC Astrocytes","Trevino_hCS Astrocytes", "CRISPRi")
 # Convert your data to a list of sets
list_input <- lapply(names(data), function(col) {
  rownames(data)[data[[col]]]
})
names(list_input) <- names(data)
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/IV_Fig1_UpsetPlot_long.pdf",
    height=3, width=5)
x=upset(fromList(list_input), order.by = "freq", matrix.color="slategray3",
      sets.bar.color=c("indianred",  "midnightblue", "steelblue","chocolate", "sandybrown"))
dev.off()
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/IV_Fig1_UpsetPlot_short.pdf",
    height=2.5, width=2)
upset(fromList(list_input), order.by = "freq", matrix.color="slategray3",
      sets.bar.color=c("indianred",  "midnightblue", "steelblue","chocolate", "sandybrown"))
dev.off()

# ##
# jaccard <- function(set1, set2) {
#   intersection <- sum(set1 & set2)
#   union <- sum(set1 | set2)
#   return(intersection / union)
# }
# 
# # Get the number of categories
# n_categories <- ncol(data)
# 
# # Initialize the Jaccard index matrix
# jaccard_matrix <- matrix(0, nrow = n_categories, ncol = n_categories)
# colnames(jaccard_matrix) <- colnames(data)
# rownames(jaccard_matrix) <- colnames(data)
# 
# # Calculate pairwise Jaccard indices
# for (i in 1:n_categories) {
#   for (j in i:n_categories) {
#     jaccard_matrix[i, j] <- jaccard(data[, i], data[, j])
#     jaccard_matrix[j, i] <- jaccard_matrix[i, j]  # Matrix is symmetric
#   }
# }
