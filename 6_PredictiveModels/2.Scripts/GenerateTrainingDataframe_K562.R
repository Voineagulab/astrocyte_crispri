##############
##
##Generate Data.frame for ENCODE-E2G EGP 
## NOTE: HG38
## Sam Bagot
##Edited by IV and JP 2_ENCODE_rE2G_K562_JP_IV
#############
rm(list=ls())
# Set the working directory to the project directory
path <-"/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised"
setwd(paste0(path, "/1.Data/InputData"))

# Source necessary R scripts for functions
source(paste0(path,"/2.Scripts/Functions.R"))
# Load required libraries
library("readxl")
library("biomaRt")
library("tidyr")
library("PRROC")
library("ggplot2")

# Create a connection to the Ensembl biomart
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

# Read gene information from a file and preprocess it
geneinfo <- read.table("GeneInfo.txt")
geneinfo$TSS <- geneinfo$Start
geneinfo$TSS[geneinfo$Strand == "-"] <- geneinfo$End[geneinfo$Strand == "-"]
colnames(geneinfo)[colnames(geneinfo) == "Symbol"] <- "Gene"

# Read the main ENCODE benchmarking dataset K562_EGPs from a TSV file
K562_EGPs <- read.table("K562/ENCFF968BZL.tsv", sep = "\t", header = TRUE)

# Rename columns for better readability
colnames(K562_EGPs)[1:15] <- c("Enh.chr", "Enh.start", "Enh.end", "name", "EffectSize", "strandPerturbationTarget", 
                               "Enh",  "Gene.chr", "TSS", "Gene.TSS.end", 
                               "strandGene", "EffectSize95ConfidenceIntervalLow", 
                               "EffectSize95ConfidenceIntervalHigh", "Gene", 
                               "EnsID")
K562_EGPs$Original.Gene <- K562_EGPs$Gene

# Replace old HGNC symbols in the dataset
old_HGNC_symbols <- unique(K562_EGPs[!K562_EGPs$Gene %in% geneinfo$Gene, "Gene"])

# Query Ensembl biomart to fetch new HGNC symbols
unique(K562_EGPs[!K562_EGPs$EnsID %in% geneinfo$EnsID, "Gene"])
syns <- getBM(attributes = c("external_synonym", "hgnc_symbol"), filters = c("external_synonym"),
              values = old_HGNC_symbols, mart = ensembl) 
ENSids <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), filters = c("external_gene_name"),
                values = old_HGNC_symbols, mart = ensembl) 
syns <- syns[syns$hgnc_symbol %in% geneinfo$Gene,] # Filter relevant new genes
old_HGNC_symbols <- data.frame(Gene = old_HGNC_symbols, external_synonym = toupper(old_HGNC_symbols))
old_HGNC_symbols <- merge(old_HGNC_symbols, syns, all.x = TRUE)
old_HGNC_symbols <- old_HGNC_symbols[,c("Gene", "hgnc_symbol")]

# Manually fix some symbol mappings
old_HGNC_symbols[old_HGNC_symbols$Gene == "C14orf2", "hgnc_symbol"] <- "ATP5MPL" 
old_HGNC_symbols[old_HGNC_symbols$Gene == "C21orf33", "hgnc_symbol"] <- "GATD3A"
old_HGNC_symbols[old_HGNC_symbols$Gene == "USMG5", "hgnc_symbol"] <- "ATP5MD"
old_HGNC_symbols <- old_HGNC_symbols[match(unique(old_HGNC_symbols$Gene), old_HGNC_symbols$Gene),]

# Merge new HGNC symbols into the main dataset
K562_EGPs <- merge(K562_EGPs, old_HGNC_symbols, all.x = TRUE)
K562_EGPs[!is.na(K562_EGPs$hgnc_symbol), "Gene"] <- K562_EGPs[!is.na(K562_EGPs$hgnc_symbol), "hgnc_symbol"]
K562_EGPs <- K562_EGPs[,colnames(K562_EGPs) != "hgnc_symbol"]
unique(K562_EGPs[!K562_EGPs$Gene %in% geneinfo$Gene, "Gene"])

#####  Perform various data manipulations on K562_EGPs dataset

# Create a new column for the genomic position of enhancers
K562_EGPs$Enh.Pos <- paste0(K562_EGPs$Enh.chr, ":", K562_EGPs$Enh.start, "-", K562_EGPs$Enh.end)
K562_EGPs$Enh <- K562_EGPs$Enh.Pos

# Create a unique identifier for each enhancer-gene pair
K562_EGPs$Pair <- paste0(K562_EGPs$Original.Gene, "_", K562_EGPs$Enh.Pos)

# Calculate and store the size of each enhancer
K562_EGPs$Enh.size <- K562_EGPs$Enh.end - K562_EGPs$Enh.start

# Calculate and store the midpoint of each enhancer
K562_EGPs$Enh.midpoint <- (K562_EGPs$Enh.start + as.integer(K562_EGPs$Enh.size / 2))

# Handle a bug where TRIR has no TSS in their data.frame (assume "-" strand)
K562_EGPs[is.na(K562_EGPs$TSS), c("TSS", "Gene.TSS.end")] <- unlist(geneinfo[geneinfo$Gene == "TRIR",c("End", "End")])

# Calculate the distance between enhancer midpoint and gene TSS
K562_EGPs$Gene.Distance <- abs(K562_EGPs$Enh.midpoint - K562_EGPs$TSS)

# Store gene TSS information
K562_EGPs$Gene.TSS <- K562_EGPs$TSS

# Check for Gasperini genes not present in our hg38 gene annotation file
unique(K562_EGPs[!K562_EGPs$Gene %in% geneinfo$Gene,"Gene"])

# Check for genes not present in either Gene or EnsID columns
unique(K562_EGPs[(!K562_EGPs$Gene %in% geneinfo$Gene) & (!K562_EGPs$EnsID %in% geneinfo$EnsID),"Gene"])

# Figure out which study those genes come from
unique(K562_EGPs[!K562_EGPs$Gene %in% geneinfo$Gene,"Reference"])
# Replace missing values from the ENCODE dataset with 0 
K562_EGPs[is.na(K562_EGPs)] <- 0

# Save a bed file of TSS coordinates
# write.table(unique(K562_EGPs[,c("Gene.chr", "TSS", "TSS", "Gene")]),
#             file = "Data/K562Data/ENCODE_rE2G/Gene_TSS.bed", 
#             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Save a bed file of Enhancer coordinates with and without "chr" in the chromosome names
# write.table(unique(K562_EGPs[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")]), 
#             file = "Data/K562Data/ENCODE_rE2G/ENCODE_Candidates_chr.bed", quote = FALSE,
#             sep = "\t", row.names = FALSE, col.names = FALSE)
K562_EGPs$Enh.chr <- sub("chr", "", K562_EGPs$Enh.chr)

# write.table(unique(K562_EGPs[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")]), 
#             file = "Data/K562Data/ENCODE_rE2G/ENCODE_Candidates.bed", quote = FALSE,
#             sep = "\t", row.names = FALSE, col.names = FALSE)

# Create and save a bed file with a 500bp window around the enhancer midpoint
windowBed <- makeWindowBed(K562_EGPs, window = 500)
# write.table(windowBed, file = "Data/K562Data/ENCODE_rE2G/window.bed", quote = FALSE,
#             sep = "\t", row.names = FALSE, col.names = FALSE)

##### Add nearest gene, and gene stability index information
K562_EGPs <- getGeneNearest(K562_EGPs, geneinfo)
K562_EGPs <- addGeneStability(K562_EGPs, geneinfo, ensembl = ensembl)

##### Add ChIP-seq and open chromatin state data; Files generated in Process_ChIP_DNase.sh
folder <- "K562/Processed_ChIP_ENCODE/"
filenames <- list.files(folder, pattern = ".tab")
# accessions numbers provided by TAP-seq for K562 cells
# hg38
# ENCFF301TVL - H3K27ac
# ENCFF290LQY - H3K4me1
# ENCFF611YPB - H3K4me3
# ENCFF915XIL - H3K27me3
# ENCFF438GBD - POLR2A
# ENCFF538GKX - DNase-seq

# Iterate through each file in the folder
for (file in filenames) {
  averages <- read.table(paste0(folder, file))
  name <- paste0("Chip.", sub("Encode_(.*).tab", "\\1", file))
  averages <- averages[,c(1,6)]  # Extract columns: ID and mean over bigwig
  colnames(averages) <- c("Enh", name)  # Rename columns
  K562_EGPs <- merge(K562_EGPs, averages, by = "Enh")  # Merge data into K562_EGPs
}

# DNAse pileup 
peaks.pileup <- read.table("K562/Enhancer_DNAse_Encode.regions.bed")
peaks.pileup <- peaks.pileup[,4:5]
colnames(peaks.pileup) <- c("Enh", "DNAse.Pileup")
K562_EGPs <- merge(K562_EGPs, peaks.pileup)

# Read rE2G data and perform merges
myrun_E2G <- read.delim("K562/ENCODE-E2G_Predictions.tsv", sep="\t", header=TRUE)
myrun_E2G$Pair <- paste0(myrun_E2G$measuredGeneSymbol, "_", 
                         myrun_E2G$chrom, ":", myrun_E2G$chromStart, "-",
                         myrun_E2G$chromEnd)
myrun_E2G <- myrun_E2G[,c("Pair", "ENCODE.E2G.Score")]
colnames(myrun_E2G)[2] <- "rE2G.DNAseOnly"

# Merge data with K562_EGPs, handling missing pairs
K562_EGPs <- merge(K562_EGPs, myrun_E2G, all.x = TRUE)

# Read rE2G-extended data files and perform merges
myrun_E2G_EXT <- read.delim("K562/ENCODE-E2G_Extended_Predictions.tsv", sep="\t", header=TRUE)
myrun_E2G_EXT$Pair <- paste0(myrun_E2G_EXT$measuredGeneSymbol, "_", 
                             myrun_E2G_EXT$chrom, ":", myrun_E2G_EXT$chromStart, "-",
                             myrun_E2G_EXT$chromEnd)
myrun_E2G_EXT <- myrun_E2G_EXT[,-c(1:21)] # remove CRISPRi benchmarking data and retain training and prediction data
K562_EGPs <- merge(K562_EGPs, myrun_E2G_EXT, all.x = TRUE)

# Rename specific columns
colnames(K562_EGPs)[colnames(K562_EGPs) == "ENCODE.E2G_Extended.Score"] <- "rE2G.Extended"
colnames(K562_EGPs)[colnames(K562_EGPs) == "X3DContact"] <- "hic_contact"
colnames(K562_EGPs)[colnames(K562_EGPs) == "ABCScore"] <- "ABC_Score"

# Extract ENCODE predictrors and write to a file
ENCODE_Ext_predictors <- colnames(K562_EGPs)[which(colnames(K562_EGPs) =="EpiMapScore"): which(colnames(K562_EGPs) == "normalizedEP300_enhActivity")]
#write.table(ENCODE_Ext_predictors, "Results/ENCODE_rE2GPredictions/ENCODE_Predictors.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)

# Count the number of missing values in each column
apply(K562_EGPs, 2, function(x) {sum(is.na(x))})

# Calculate additional columns
K562_EGPs$ENCODE_overlap.Exists <- ! (is.na(K562_EGPs$rE2G.Extended) | is.na(K562_EGPs$rE2G.DNAseOnly))
K562_EGPs$Missing.Values <- K562_EGPs$normalizedDNase_enh == 0 & K562_EGPs$numTSSEnhGene == 0 & signif(K562_EGPs$phyloPMax,6) == 0.596624

# Define Negative and Positive EG pairs
K562_EGPs$Significant_NegEffectSize <- K562_EGPs$Significant & K562_EGPs$EffectSize < 0

# Sort the data by the "Pair" column
K562_EGPs <- K562_EGPs[order(K562_EGPs$Pair),]

# Write the final dataset to a CSV file
write.csv(K562_EGPs, file = "../TrainingData/TrainingDataframe_K562.csv", row.names = FALSE)
