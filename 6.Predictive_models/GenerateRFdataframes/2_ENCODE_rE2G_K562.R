##############
##
##Generate Data.frame for ENCODE-E2G EGP 
## NOTE: HG38
##
#############
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPseq/EnhancerPredictionModels/")
source("Scripts/Annotation_Functions.R")
source("Scripts/Header_functions.R")
library(readxl)
library("biomaRt")
#dir.create("Data/K562Data/ENCODE_rE2G/")
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)


#Want to add ENSID but some symbols are missing 
geneinfo <- read.table("../FullScale/Data/Whitelists/GeneInfo.txt")
geneinfo$TSS <- geneinfo$Start; geneinfo$TSS[geneinfo$Strand == "-"] <- geneinfo$End[geneinfo$Strand == "-"]
colnames(geneinfo)[colnames(geneinfo)  == "Symbol"] <- "Gene"

K562_EGPs <- read.table("../PublicData/ENCODE_rE2G/ENCFF968BZL.tsv", 
                        sep = "\t", header = T)
# K562_EGPs_old <- read.table("../PublicData/ENCODE_rE2G/EPCrisprBenchmark_ensemble_data_GRCh38.tsv", 
#                              sep = "\t", header = T)
colnames(K562_EGPs)[1:15] <- c("Enh.chr", "Enh.start", "Enh.end", "name", "EffectSize", "strandPerturbationTarget", 
                               "Enh",  "Gene.chr", "TSS", "Gene.TSS.end", 
                               "strandGene", "EffectSize95ConfidenceIntervalLow", 
                               "EffectSize95ConfidenceIntervalHigh", "Gene", 
                              "EnsID")
K562_EGPs$Original.Gene <- K562_EGPs$Gene

###########
##Replace old HGNC symbols in their dataset
###########
old_HGNC_symbols<- unique(K562_EGPs[! K562_EGPs$Gene %in% geneinfo$Gene, "Gene"])
#tsvs[tsvs$measuredGeneSymbol %in% old_HGNC_symbols[old_HGNC_symbols != "C19orf43"],]

unique(K562_EGPs[! K562_EGPs$EnsID %in% geneinfo$EnsID, "Gene"])
syns <- getBM(attributes = c("external_synonym", "hgnc_symbol") , filters = c("external_synonym"),
              values = old_HGNC_symbols, mart = ensembl) 
ENSids <- getBM(attributes = c("external_gene_name", "ensembl_gene_id") , filters = c("external_gene_name"),
              values = old_HGNC_symbols, mart = ensembl) 
syns <- syns[syns$hgnc_symbol %in% geneinfo$Gene,]#52/74 new genes picked up
old_HGNC_symbols <- data.frame(Gene = old_HGNC_symbols, external_synonym = toupper(old_HGNC_symbols))
old_HGNC_symbols <- merge(old_HGNC_symbols, syns, all.x = T)
old_HGNC_symbols <- old_HGNC_symbols[,c("Gene", "hgnc_symbol")]
#Fixing manually  we use a different old version not the current one
old_HGNC_symbols[old_HGNC_symbols$Gene == "C14orf2", "hgnc_symbol"] <- "ATP5MPL" 
old_HGNC_symbols[old_HGNC_symbols$Gene == "C21orf33", "hgnc_symbol"] <- "GATD3A"
old_HGNC_symbols[old_HGNC_symbols$Gene == "USMG5", "hgnc_symbol"] <- "ATP5MD"
old_HGNC_symbols <- old_HGNC_symbols[match(unique(old_HGNC_symbols$Gene),old_HGNC_symbols$Gene),]
K562_EGPs <- merge(K562_EGPs, old_HGNC_symbols, all.x = T)
K562_EGPs[! is.na(K562_EGPs$hgnc_symbol), "Gene"]<- K562_EGPs[! is.na(K562_EGPs$hgnc_symbol),"hgnc_symbol"]
K562_EGPs <- K562_EGPs[,colnames(K562_EGPs) != "hgnc_symbol"]
unique(K562_EGPs[! K562_EGPs$Gene %in% geneinfo$Gene, "Gene"])


#write.csv(tsvs[tsvs$measuredGeneSymbol %in% old_HGNC_symbols[is.na(old_HGNC_symbols$hgnc_symbol),"Gene"],], file = "missingEGPs.csv", row.names = F)


#K562_EGPs[! K562_EGPs$Gene %in% geneinfo$Gene, c("Significant", "Pair", "Gene", "Reference")]


K562_EGPs$Enh.Pos <- paste0(K562_EGPs$Enh.chr,":", K562_EGPs$Enh.start, "-",K562_EGPs$Enh.end)
K562_EGPs$Enh <- K562_EGPs$Enh.Pos
K562_EGPs$Pair <- paste0(K562_EGPs$Original.Gene, "_",K562_EGPs$Enh.Pos)
K562_EGPs$Enh.size <- K562_EGPs$Enh.end - K562_EGPs$Enh.start
K562_EGPs$Enh.midpoint <- (K562_EGPs$Enh.start + as.integer(K562_EGPs$Enh.size / 2))
#Bug TRIR has no TSS in their data.frame so doing it like this. "-" strand so end is TSS
K562_EGPs[is.na(K562_EGPs$TSS),c("TSS", "Gene.TSS.end")] <- unlist(geneinfo[geneinfo$Gene == "TRIR",c("End", "End")])
K562_EGPs$Gene.Distance <- abs(K562_EGPs$Enh.midpoint - K562_EGPs$TSS)
K562_EGPs$Gene.TSS <- K562_EGPs$TSS


gene_promoters <- K562_EGPs[match(unique(K562_EGPs$Gene),K562_EGPs$Gene),c('Gene', "Gene.chr", "Gene.TSS", "Gene.TSS")]
gene_promoters$Gene.TSS <-  as.numeric(gene_promoters$Gene.TSS) - 500
gene_promoters$Gene.TSS.1 <-  as.numeric(gene_promoters$Gene.TSS.1) + 500
write.bed(gene_promoters[,c(2,3,4,1)], file = "Data/K562Data/Gene_Promoters.bed")



#Fixed all but one Gasperini gene that doesn't seem to be in HG38
unique(K562_EGPs[! K562_EGPs$Gene %in% geneinfo$Gene,"Gene"])
unique(K562_EGPs[(! K562_EGPs$Gene %in% geneinfo$Gene) & (!K562_EGPs$EnsID %in% geneinfo$EnsID),"Gene"]) 
unique(K562_EGPs[! K562_EGPs$Gene %in% geneinfo$Gene,"Reference"])

write.table(unique(K562_EGPs[,c("Gene.chr", "TSS", "TSS", "Gene")]) ,
            file = "Data/K562Data/ENCODE_rE2G/Gene_TSS.bed", 
            quote =F, row.names = F, col.names = F, sep ="\t" )
write.table(unique(K562_EGPs[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")]), 
            file = "Data/K562Data/ENCODE_rE2G/ENCODE_Candidates_chr.bed", quote = F,
            sep = "\t", row.names = F, col.names = F)
windowBed <- makeWindowBed(K562_EGPs, window = 500)
write.table(windowBed, file = "Data/K562Data/ENCODE_rE2G/window.bed", quote = F,
            sep = "\t", row.names = F, col.names = F)
K562_EGPs$Enh.chr <- sub("chr", "",K562_EGPs$Enh.chr)
write.table(unique(K562_EGPs[,c("Enh.chr", "Enh.start", "Enh.end", "Enh")]), 
            file = "Data/K562Data/ENCODE_rE2G/ENCODE_Candidates.bed", quote = F,
            sep = "\t", row.names = F, col.names = F)


#Add tobias
K562_EGPs <- addExpTobiasBound(K562_EGPs, 
                  exp_bound_overlaps = "Results/Tobias/ENCODEFootprint/BINDetect/exp_bound_overlaps.bed", 
                  exp_unbound_overlaps = "Results/Tobias/ENCODEFootprint/BINDetect/exp_unbound_overlaps.bed")
#Previously this analysis had a slight bias towards Gasperini peaks (since I used their peaks) 
#so have replaced with new ENCODE dataset... 
#turned out to become more biased but a better predictor so whatever
fisher.test(K562_EGPs$Tobias.Exp_Bound_TF_counts != 0, K562_EGPs$Reference == "Gasperini et al., 2019")
aucPlots(K562_EGPs$Tobias.Exp_Bound_TF_counts, K562_EGPs$Significant)

#SuperEnhancer annotation
#bedtools intersect -a /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/Candidates.bed -b /mnt/Data0/PROJECTS/CROPSeq/PublicData/Hnisz_Cell2013_SuperEnhancers/PROCESSED/liftOver_hg38/K562.hg38.bed | awk '{ print $4 }' > /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/K562/ENCODE_E2G_OverlappingSuperEnhancers.txt
superEnhancers <- read.table("Results/K562/ENCODE_E2G_OverlappingSuperEnhancers.txt")
K562_EGPs$Hnisz_SuperEnhancers <- K562_EGPs$Enh %in% superEnhancers$V1

#Gene expression, nearest gene, gene stability
K562_EGPs <- addFeatureCounts(K562_EGPs, "Data/K562Data/RNAseq/Ref_RNAseq_FeatCounts.csv",try.genecol = T)
colnames(K562_EGPs)[colnames(K562_EGPs)== "Gene.Ref_RNAseq_RPKM"] <- "Gene.RNAseq_RPKM"
K562_EGPs <- getGeneNearest(K562_EGPs, geneinfo)
K562_EGPs <- addGeneStability(K562_EGPs, geneinfo, ensembl = ensembl)

###### Add TTseq gene counts
TTseq_counts <-read.table("../PublicData/K562TTseq/salmon.merged.gene_counts.tsv", header = T)
TTseq_counts<- TTseq_counts[,c("gene_id", "gene_name", "TTseq_R1", "RNAseq_R1")]
colnames(TTseq_counts) <- c("EnsID","Gene",  "Gene.TTseq_Counts", "Gene.RNAseq_Counts2")
TTseq_counts <- TTseq_counts[match(unique(TTseq_counts$Gene), TTseq_counts$Gene),]
K562_EGPs <- merge(K562_EGPs, TTseq_counts[,c("Gene", "Gene.TTseq_Counts", "Gene.RNAseq_Counts2")], all.x = T)
tmp <-  merge(K562_EGPs[is.na(K562_EGPs$Gene.TTseq_Counts),colnames(K562_EGPs) != "Gene.TTseq_Counts"], TTseq_counts[,c("EnsID", "Gene.TTseq_Counts")], all.x = T)
K562_EGPs<- rbind(K562_EGPs[!is.na(K562_EGPs$Gene.TTseq_Counts),], tmp)

#Add BigWigs
#K562_EGPs <- addChipBigWigs(K562_EGPs, filenames, folder )
folder <- "Results/Process_Chip/ENCODE/"
filenames <- list.files(folder, pattern = ".tab")
for (file in filenames) {
  averages <- read.table(paste0(folder, file))
  name <-  paste0("Chip.",sub("Encode_(.*).tab", "\\1", file))
  averages <- averages[,c(1,6)] #column 1 is the ID and column 6 is the mean over bigwig
  colnames(averages) <- c("Enh", name)
  K562_EGPs <- merge(K562_EGPs, averages, by = "Enh")
}
K562_EGPs$ABC_Approx <- K562_EGPs$Chip.H3K27ac / K562_EGPs$Gene.Distance

#Use DNAse pileup instead of ATAC.pileup
peaks.pileup <- read.table("Results/Process_Chip/ENCODE/Enhancer_DNAse_Encode.regions.bed")
peaks.pileup <- peaks.pileup[,4:5]
colnames(peaks.pileup) <- c("Enh", "DNAse.Pileup")
K562_EGPs <- merge(K562_EGPs, peaks.pileup)
#Add ttseq enhancer data 
ttseq <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/4_EnhancerTranscription/K562TTseq/Results_Table_K562_ENCODE.csv")
colnames(ttseq) <- paste0("TTseq.", colnames(ttseq))
colnames(ttseq)[1] <- "Enh"
K562_EGPs <- merge(K562_EGPs, ttseq[,c("Enh", "TTseq.TTseq_Total", "TTseq.TPM_TTseq_Total", "TTseq.Ratio_TTversusRNA")])
#"TTseq.TT_Enrich", "TTseq.RNAseq_Total",


myrun_E2G <- read.csv("ENCODE_E2G_clone/ENCODE-E2G/FullPredictions_ENCODE.csv")
myrun_E2G$Pair <- paste0(myrun_E2G$measuredGeneSymbol, "_", 
                    myrun_E2G$chrom, ":", myrun_E2G$chromStart, "-",
                    myrun_E2G$chromEnd)
myrun_E2G <- myrun_E2G[,c("Pair","ENCODE.E2G.Score")]
colnames(myrun_E2G)[2] <- c("rE2G.DNAseOnly")

#37 pairs not in their predictions
K562_EGPs <- merge(K562_EGPs, myrun_E2G, all.x = T)

#This file also has all their predictors
myrun_E2G_EXT <- read.csv("ENCODE_E2G_clone/ENCODE-E2G/FullPredictions_ENCODE_EXT.csv")
their_E2G <- read.table("ENCODE_E2G_clone/ENCODE-E2G/data/crispri/ENCODE-E2G_Extended_Predictions.tsv", header = T,sep = "\t")
myrun_E2G_EXT$Pair <- paste0(myrun_E2G_EXT$measuredGeneSymbol, "_", 
                             myrun_E2G_EXT$chrom, ":", myrun_E2G_EXT$chromStart, "-",
                             myrun_E2G_EXT$chromEnd)
myrun_E2G_EXT <- myrun_E2G_EXT[,- c(1:22)]
K562_EGPs <- merge(K562_EGPs, myrun_E2G_EXT, all.x = T)
colnames(K562_EGPs)[colnames(K562_EGPs) == "ENCODE.E2G_Extended.Score"] <- "rE2G.Extended"
colnames(K562_EGPs)[colnames(K562_EGPs) == "X3DContact"] <- "hic_contact"
colnames(K562_EGPs)[colnames(K562_EGPs) == "ABCScore"] <- "ABC_Score"
ENCODE_Ext_predicitors <-  colnames(K562_EGPs)[which(colnames(K562_EGPs) =="EpiMapScore"):
                                                 which(colnames(K562_EGPs) == "normalizedEP300_enhActivity")]
write.table(ENCODE_Ext_predicitors, "Results/ENCODE_rE2GPredictions/ENCODE_Predictors.txt",col.names = F, quote = F, row.names = F)

apply(K562_EGPs, 2, function(x) {sum(is.na(x))})
#Set them to 0 in the comparisons function for now 
#They have 932 EGPs where they have no infromation

#They seem to have missing values for this and Hi-C (and I use this in my code so quite annoying) 
#so I'm re-running it and replacing with the same name
#K562_EGPs$numTSSEnhGene <- addnumTSSbetween(K562_EGPs, geneinfo)

#Both of these must be removed
K562_EGPs$ENCODE_overlap.Exists <- ! (is.na(K562_EGPs$rE2G.Extended) | is.na(K562_EGPs$rE2G.DNAseOnly))
K562_EGPs$Missing.Values <- K562_EGPs$normalizedDNase_enh == 0 & K562_EGPs$numTSSEnhGene == 0 & signif(K562_EGPs$phyloPMax,6) == 0.596624

K562_EGPs$Significant_NegEffectSize <- K562_EGPs$Significant & K562_EGPs$EffectSize < 0

K562_EGPs <- K562_EGPs[order(K562_EGPs$Pair),]
write.csv(K562_EGPs, file = "Results/K562/ENCODE_EGP.csv", row.names = F)
#K562_EGPs <- read.csv("Results/K562/ENCODE_EGP.csv")

#They have pre-predicted dataset with NAs filled (and non-NAs filled)
#tsvs <-  read.table("ENCODE_E2G_clone/ENCODE-E2G/data/crispri/EPCrisprBenchmark_ensemble_data_GRCh38.K562_AllFeatures_NAfilled.FullModel.tsv",header = T, sep = "\t")
