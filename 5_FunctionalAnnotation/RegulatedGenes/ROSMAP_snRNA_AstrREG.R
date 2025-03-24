rm(list=ls())
library(data.table)

###################################################### Functions
makePlaceholders=function(astfiles) 
{
  astfiles$Ngenes=NA
  astfiles$loc.fisher.p=NA
  astfiles$loc.fisher.or=NA
  astfiles$loc.Nintersect=NA
  astfiles$glb.fisher.p=NA
  astfiles$glb.fisher.or=NA
  astfiles$glb.Nintersect=NA
  return(astfiles)
}

extractResults=function(astfiles)
{
  results=list()
  #Loop through ROSMAP results files and do fisher tests for Hit vs Nonhit genes / DE sig vs Nonsig using either the local or global pvalue for significance
  #Combine all significant genes across files 
  #Note: Fisher test gives an error if there is a single category (for example none of our tested genes has p<0.05), that is dealt with by tryCatch and the output will have an NA value.
  for (j in c(1:nrow(astfiles)))
  {
    tryCatch({
      data=read.csv(astfiles$Path[j])
      data=cbind(genes, data[match(genes$Symbol, data$gene) , ])
      data=data[!is.na(data$gene), ]
      # replace coef with the specific comparison
      data$coef=astfiles$Trait[j]
      
      astfiles$Ngenes[j]=nrow(data)
      
      astfiles$loc.fisher.p[j]=fisher.test(data$Hit, data$p_adj.loc < 0.05)$p.value
      astfiles$loc.fisher.or[j]=fisher.test(data$Hit, data$p_adj.loc < 0.05)$estimate
      astfiles$loc.Nintersect[j]=table(data$Hit, data$p_adj.loc < 0.05)["TRUE", "TRUE"]
      
      astfiles$glb.fisher.p[j]=fisher.test(data$Hit, data$p_adj.glb < 0.05)$p.value
      astfiles$glb.fisher.or[j]=fisher.test(data$Hit, data$p_adj.glb < 0.05)$estimate
      astfiles$glb.Nintersect[j]=table(data$Hit, data$p_adj.glb < 0.05)["TRUE", "TRUE"]
      
      s=data[((data$p_adj.loc < 0.05)|(data$p_adj.glb <0.05)) , ]
      if (j==1) sig.genes=s else sig.genes=rbind(sig.genes , s)
      
      print(j)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  # pvalue adjustment
  astfiles$loc.fisher.padj=p.adjust(astfiles$loc.fisher.p, method="BH")
  astfiles$glb.fisher.padj=p.adjust(astfiles$glb.fisher.p, method="BH")
  #results
  results$fisher.tests=astfiles
  results$sig.genes=sig.genes
  return(results)
}

###################################################### Load crispri results
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
genes=data.frame(unique(res$Gene), FALSE)
colnames(genes)=c("Symbol", "Hit")
genes$Hit[which(genes$Symbol%in%res$Gene[res$HitPermissive])]=TRUE

rosmap.results=list()

###################################################### 1. ROSMAP Differential_gene_expression_analysis results
######################################################
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/ROSMAP_snRNAseq_PFC-main/Results/Differential_gene_expression_analysis")
#Load ROSMAP DE data
x=list.files(, recursive = TRUE, pattern=".csv") # sometimes this function does not list all files!! does it time out?? It should be a total 3268 of files and 196 Ast files
#select astrocyte files
astfiles=data.frame(x[grep("Ast", x)], NA)
colnames(astfiles)=c("Path", "File")
names=transpose(strsplit(astfiles[,1], split="/", fixed=TRUE))
length(names)
#most csv files are in a subdirectory; a few subdirs have one more level -> extract file name
y=which(is.na(names[[3]])==FALSE)
astfiles$File[y]=names[[3]][y]
astfiles$File[-y]=names[[2]][-y]
#extract ct and trait
astfiles$CT=transpose(strsplit(astfiles$File, split="_", fixed=TRUE))[[1]]
astfiles$Trait=transpose(strsplit(astfiles$Path, split="/", fixed=TRUE))[[1]]

#placeholders
astfiles=makePlaceholders(astfiles)
#extract results
rosmap.results$DE=extractResults(astfiles)


###################################################### 2. ROSMAP Four Group comparisons
######################################################
rm(list=setdiff(ls(), c("genes", "rosmap.results", "makePlaceholders", "extractResults")))
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/ROSMAP_snRNAseq_PFC-main/Results/Four-group/")
#Load ROSMAP Four Group comparison data
x=list.files(, recursive = TRUE, pattern=".csv") 
#select astrocyte files
astfiles=data.frame(x[grep("Ast", x)], NA)
colnames(astfiles)=c("Path", "File")
names=transpose(strsplit(astfiles[,1], split="/", fixed=TRUE))
length(names)
astfiles$File=names[[2]]
#extract ct and trait
astfiles$CT=names[[1]]
astfiles$Trait=gsub("Ast_", "", astfiles$File)
astfiles$Trait=gsub(".csv", "", astfiles$Trait)

#placeholders
astfiles=makePlaceholders(astfiles)
#extract results
rosmap.results$FourGroup=extractResults(astfiles)

###################################################### 3. ROSMAP Disease progression
######################################################
rm(list=setdiff(ls(), c("genes", "rosmap.results", "makePlaceholders", "extractResults")))
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/ROSMAP_snRNAseq_PFC-main/Results/Disease_progression/")

#Load ROSMAP Disease progression data
x=list.files(, recursive = TRUE, pattern=".csv") 
#select astrocyte files
astfiles=data.frame(x[grep("Ast", x)], NA)
colnames(astfiles)=c("Path", "File")

names=transpose(strsplit(astfiles[,1], split="/", fixed=TRUE))
length(names)
astfiles$File=names[[3]]
#extract ct and trait
astfiles$CT=names[[2]]
astfiles$Trait=gsub("Ast_", "", astfiles$File)
astfiles$Trait=gsub(".csv", "", astfiles$Trait)

#placeholders
astfiles=makePlaceholders(astfiles)
#extract results
rosmap.results$DisProg=extractResults(astfiles)

save(rosmap.results, file="/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.Results.AstroREG.rda")
