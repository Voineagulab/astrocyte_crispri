# setup
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")
load("../2_DE/Enhancers - All Results Summary.rda")
pe.res <- read.csv("../3_HitEnrichment/eQTL - PE - Data.csv", row.names = 1)
gtex.res <- read.csv("../3_HitEnrichment/eQTL - GTEx - Data.csv", row.names = 1)
disg.res <- read.csv("../3_HitEnrichment/Disgenet - Peak Variant Overlap.csv", row.names = 1)

library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl)
genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"), mart = ensembl)

# define hit
res <- res.all[which(res.all$High.Expression),]
res$NegBinom.FDR <- p.adjust(res$NegBinom.P, method = "fdr")
res$Hit <- res$NegBinom.FDR < 0.05
res.hit <- res[which(res$Hit),]
hits <- unique(res$Enh.Pos[which(res$Hit)])

# pair characteristics
length(unique(res.hit$Enh)) # number of enhancers
length(unique(res.hit$Gene)) # number of genes
nrow(res.hit) # number of genes

# pe eqtl
pe.res$Hit <- pe.res$Peak.id %in% hits
table(pe.res$Hit)
m <- match(pe.res$eQTL.gene, genes$ensembl_gene_id)
pe.res$eQTL.gene <- genes$hgnc_symbol[m]
m <- match(pe.res$Peak.id, res.hit$Enh.Pos)
pe.res$Peak.gene <- res.hit$Gene[m]
pe.res$SameTarget <- pe.res$eQTL.gene == pe.res$Peak.gene
table(pe.res$SameTarget)


# gtex eqtl
gtex.res$Hit <- gtex.res$Peak.id %in% hits
table(gtex.res$Hit)

gtex.res.b <- gtex.res[grep("Brain", gtex.res$eQTL.tissue),]
length(unique(gtex.res$eQTL.id)) # 661 
length(unique(gtex.res.b$eQTL.id)) # 166
z <- gtex.res.b[which(gtex.res.b$SameTarget),]
length(unique(z$eQTL.id))

z <- gtex.res[which(gtex.res$SameTarget),]
length(unique(z$eQTL.id))

m <- match(gtex.res$Peak.id, res.hit$Enh.Pos)
gtex.res$Peak.gene <- res.hit$Gene[m]
gtex.res$SameTarget <- gtex.res$eQTL.gene == gtex.res$Peak.gene
table(gtex.res$SameTarget)

# disgenet variants
disg.res$Hit <- disg.res$Peak.id %in% hits
table(disg.res$Hit)

# disgenet genes
library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(
  email = "i.voineagu@unsw.edu.au", 
  password = "UNSWlab123" )
Sys.setenv(DISGENET_API_KEY = disgenet_api_key)
disg.genes <- gene2disease(gene = as.character(x$Symbol), 
                      database = "ALL", 
                      score = c(0.2, 1), 
                      verbose= TRUE)
View(disg.genes@qresult)
length(unique(disg.genes@qresult$gene_symbol))


## disgenet genes v2
disg <- read.delim("../../../PublicData/DisGeNet_070422/all_gene_disease_associations.tsv")
x$Disgenet <- x$Symbol %in% disg$geneSymbol
disg.hit <- disg[which(disg$geneSymbol %in% hits$Gene),]
save(disg.hit, file = "../../../NHMRC2022_IV/Disgenet.hits.rda")
