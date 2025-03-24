library(data.table)

rm(list=ls())
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.Results.rda")
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.DEP.rda")
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
vars=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/ROSMAP_snRNAseq_PFC-main/Results/Differential_gene_expression_analysis/Description_of_variables_Categorised.csv")
vars=vars[-grep("##", vars$Variable), ]

genes=data.frame(unique(res$Gene), FALSE)
colnames(genes)=c("Symbol", "Hit")
genes$Hit[which(genes$Symbol%in%res$Gene[res$HitPermissive])]=TRUE

enh=data.frame(unique(res$Enh), FALSE)
colnames(enh)=c("Symbol", "Hit")
enh$Hit[which(enh$Symbol%in%res$Enh[res$HitPermissive])]=TRUE

fisher.tests=rbind(rosmap.results$DE$fisher.tests, 
                rosmap.results$FourGroup$fisher.tests,
                rosmap.results$DisProg$fisher.tests)

sig.genes=rbind(rosmap.results$DE$sig.genes, 
                   rosmap.results$FourGroup$sig.genes,
                   rosmap.results$DisProg$sig.genes)

#Focus on results in the entire Astrocyte population (rather than Astrocyte sub-types)
sig.genes=sig.genes[which(sig.genes$cluster_id%in%"Ast"), ]
fisher.tests=fisher.tests[which(fisher.tests$CT%in%"Ast"), ]
# #adjust pvals across all comparisons
fisher.tests$loc.fisher.padj=p.adjust(fisher.tests$loc.fisher.p, method="BH")
fisher.tests$glb.fisher.padj=p.adjust(fisher.tests$glb.fisher.p, method="BH")

write.csv(fisher.tests,"/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.Ast.FisherTests.csv",
          row.names=FALSE)

#subset
sig.genes.hit=sig.genes[sig.genes$Hit, ]
sig.genes.dep=sig.genes[which(sig.genes$Symbol%in%res.hits.dep$Gene), ]
sig.genes.hit.grp41=sig.genes.hit[grep("group_fourth_vs_first", sig.genes.hit$coef), ]
sig.genes.hit.Ast=sig.genes.hit[which(sig.genes.hit$cluster_id%in%"Ast"), ]
sig.genes.hit.AstGRM3=sig.genes.hit[which(sig.genes.hit$cluster_id%in%"Ast GRM3"), ]

genes$DE_AD=genes$Symbol%in%sig.genes$Symbol
genes$DE_AD_Ast=genes$Symbol%in%sig.genes$Symbol[which(sig.genes$cluster_id%in%"Ast")]
genes$DE_AD_AstGRM3=genes$Symbol%in%sig.genes$Symbol[which(sig.genes$cluster_id%in%"Ast GRM3")]

table(genes$Hit, genes$DE_AD)

#       FALSE TRUE
# FALSE  1327 1682
# TRUE     36   80
table(genes$Hit, genes$DE_AD_Ast)
# FALSE TRUE
# FALSE  1358 1651
# TRUE     38   78

# 80 of the 116 Hit genes are significant in Astrocytes or Ast subtypes in at least one comparison in the ROSMAP data.
# 78 of the 116 Hit genes are significant in Astrocytes overall in at least one comparison in the ROSMAP data.
fisher.test(genes$Hit, genes$DE_AD)
# Fisher's Exact Test for Count Data
# 
# data:  genes$Hit and genes$DE_AD
# p-value = 0.005529
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.160367 2.693394
# sample estimates:
# odds ratio 
#     1.7529 


###### Core regulatory circuitries https://genome.cshlp.org/content/26/3/385.long 
crc=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/SaintAndre_CRC/Extended_RC_Astro.csv")
genes$CRC=genes$Symbol%in%crc$Astrocytes
fisher.test(genes$Hit, genes$CRC)

# Fisher's Exact Test for Count Data
# 
# data:  genes$Hit and genes$CRC
# p-value = 1.867e-09
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   3.693403 11.336100
# sample estimates:
# odds ratio 
#   6.597118 


#Assigning categories to variables was done manually in excel:

sig.genes.hit.Ast$Category=vars$Category[match(sig.genes.hit.Ast$coef, vars$Variable)]
sig.genes.hit.AstGRM3$Category=vars$Category[match(sig.genes.hit.AstGRM3$coef, vars$Variable)]


ROSMAP.net=sig.genes.hit.Ast[-grep("Other", sig.genes.hit.Ast$Category), c("Category", "gene")]
colnames(ROSMAP.net)=c("Node1", "Node2")
ROSMAP.net$Edge=1
#write.csv(ROSMAP.net, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.net.csv", row.names = FALSE)
ROSMAP.nodes=as.data.frame(table(c(ROSMAP.net$Node1, ROSMAP.net$Node2)))
ROSMAP.nodes=ROSMAP.nodes[ROSMAP.nodes$Freq>1, ]
ROSMAP.net=ROSMAP.net[ROSMAP.net$Node2%in%ROSMAP.nodes$Var1 , ]
write.csv(ROSMAP.net, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.net.csv", row.names = FALSE)

sig.genes.hit.Ast=sig.genes.hit.Ast[-grep("Other", sig.genes.hit.Ast$Category), ]
n=as.data.frame(table(paste(sig.genes.hit.Ast$gene, sig.genes.hit.Ast$Category, sep=";")))
n$Var1=as.character(n$Var1)
ROSMAP.net=data.frame( transpose(strsplit(n$Var1, split=";", fixed=TRUE))[[1]],
                       transpose(strsplit(n$Var1, split=";", fixed=TRUE))[[2]],
                       n$Freq
                       )
colnames(ROSMAP.net)=c("Node1", "Node2", "Edge")

write.csv(ROSMAP.net, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.net.csv", row.names = FALSE, quote = FALSE)
ROSMAP.nodes=as.data.frame(table(c(sig.genes.hit.Ast$gene, sig.genes.hit.Ast$Category)))
write.csv(ROSMAP.nodes, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.nodes.csv", row.names = FALSE, quote= FALSE)
