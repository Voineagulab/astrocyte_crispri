
library(biomaRt)
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
attr=listAttributes(ensembl)
genes <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position', 'entrezgene_id'), mart = ensembl)

data=read.csv("Enhancers - All Characterisation.csv")
data$EnsGene=genes$ensembl_gene_id[match(data$Gene, genes$hgnc_symbol)]


table(is.na(data$EnsGene))
# FALSE  TRUE 
# 17613  4406 

table(is.na(data$EnsGene[data$Hit]))
# FALSE  TRUE 
# 88     4 

library(disgenet2r)

disgenet_api_key <- get_disgenet_api_key(
  email = "i.voineagu@unsw.edu.au", 
  password = "UNSWlab123" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

gene_list=as.character(data$Gene[data$Hit])

disgenes<- gene2disease( gene = gene_list, database= "GENOMICS_ENGLAND", score =c(0.2, 1), verbose= TRUE)

pdf("DisGenes_Hits_IV.pdf", height=10, width=10)
#for (j in c("network", "Piechart", "DiseaseClass", "Venn", "Heatmap", "ProteinClass", "Barplot",  "Lollipop"))
plot(disgenes)
plot(disgenes, class="Heatmap")
plot(disgenes, class="DiseaseClass")
dev.off()