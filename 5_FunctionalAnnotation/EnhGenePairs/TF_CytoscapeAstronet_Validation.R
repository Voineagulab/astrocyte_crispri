rm(list=ls())
setwd("/Users/irina/Downloads/TF_Validation/")
library(tidyr)
library(dplyr)
library(pheatmap)
library(readxl)
library(data.table)
library(stringr)

############## Read in network data 
  nodes=read.csv("STables/Nodes.csv") # node information from Astronet
  edges=read.csv("STables/Edges.csv") # edge information from Astronet
  coord=read.csv("STables/ST1.csv") # Supplementary Table 1A
  
  edges$Node1=gsub("**", "", edges$Node1, fixed=TRUE)
  edges$Node1=gsub("TEAD1/2/3", "TEAD1/TEAD2/TEAD3", edges$Node1, fixed=TRUE)
  edges_expanded <- edges %>% separate_rows(Node1, sep = "/")
  
  edges_tf=edges_expanded [which(edges_expanded$Edge==1),]
  colnames(edges_tf)[1:2]=c("TF", "Enh")
  edges_genes=edges_expanded [which(edges_expanded$Edge==1.2),]
  colnames(edges_genes)[1:2]=c("Enh", "Gene")
  edges_tf$Gene=edges_genes$Gene[match(edges_tf$Enh, edges_genes$Enh)]
  edges_tf$pair=paste(edges_tf$TF, edges_tf$Enh, sep="_")

############## Generate a bed file of enhancer coordinates for enhancers included in the network
  enh=coord[which(coord$EnhancerID%in%edges_tf$Enh), ]
  enh$Peak.ID=enh$EnhancerID
  enh=enh[which(enh$EnhancerID%in%edges_tf$Enh) , ]
  enh=enh[,1:6]
  enh[,5]=0
  enh[,6]="."
  #write.table(enh, "enh.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

############## ############## ############## Process ChIP data
  ############## ChIP Atlas
    atlas=read.csv("ChIP_Atlas/experimentList.csv")
    #########subset for:
    #human hg38
    atlas=atlas[grep("hg38", atlas$Genome.assembly, ignore.case=TRUE),]
    #TF ChIP
    atlas=atlas[grep("TF", atlas$Track.type.class, ignore.case=TRUE),]
    #TFs in our network
    atlas=atlas[which(atlas$Track.type%in%c(edges_tf$TF)),]
    #Neural cell types OR astrocytes (iPSC derived astrocytes are under the Pluripotent stem cell track class)
    atlas=atlas[c(which(atlas$Cell.type.class%in%"Neural"),
                  grep("astrocytes", atlas$Cell.type)),]
    #Subset to glial cells
    gli=c(grep("iPSC derived astrocytes", atlas$Cell.type),
          grep("U-87 MG", atlas$Cell.type),
          grep("SF268", atlas$Cell.type),
          grep("hTERT RPE-1", atlas$Cell.type))
    atlas=atlas[gli,]
    
    ############### Download the following files from https://chip-atlas.org/peak_browser 
    # 99687 Downloaded/Oth.Neu.05.AllAg.Astrocytes.bed
    # 3225 Downloaded/Oth.Neu.05.IRF3.AllCell.bed
    # 368204 Downloaded/Oth.Neu.05.JUN.AllCell.bed
    # 60396 Downloaded/Oth.Neu.05.MXI1.AllCell.bed
    # 101039 Downloaded/Oth.Neu.05.NFIB.AllCell.bed
    # 166606 Downloaded/Oth.Neu.05.REST.AllCell.bed
    # 54175 Downloaded/Oth.Neu.05.TEAD1.AllCell.bed
    # 4306 Downloaded/Oth.Neu.05.TP53.AllCell.bed
    
    ###bash
    # cd TF_Validation/ChIP_Atlas
    # cat Downloaded/Oth.* >> TF.ChIP.bed
    # bedtools intersect -wa -wb -a TF.ChIP.bed -b ../enh.bed > Overlaps.Atlas.bed
    # wc -l Overlaps.Atlas.bed
    # 109 Overlaps.Atlas.bed

  ############## Loupe et al.
    ab=read.csv("Loupe/TFs_41593_2024_1658_MOESM5_ESM.csv")
    sheet_names <- excel_sheets("Loupe/Edited_41593_2024_1658_MOESM6_ESM.xlsx")
    
    for (s in c(1:length(sheet_names)))
    {
      print(sheet_names[s])
      data=read_excel("Loupe/Edited_41593_2024_1658_MOESM6_ESM.xlsx", sheet = sheet_names[s])
      data$Categ=sheet_names[s]
      if (s==1) tfdata=data else tfdata=rbind(tfdata, data)
    }
    
    tfdata$start <- format(tfdata$start, scientific = FALSE, trim = TRUE)
    tfdata$end <- format(tfdata$end, scientific = FALSE, trim = TRUE)
    tfdata$Id=paste(tfdata$seqnames, tfdata$start, tfdata$end, sep="_")
    tfbed=tfdata[, c("seqnames", "start", "end", "Id")]
    tfbed$score=0; tfbed$strand="."
    
    #write.table(tfbed, "Loupe/LoupeTFdata.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    #bash
    # cd Loupe
    # bedtools intersect -wa -wb -a LoupeTFdata.bed -b ../enh.bed > OverlapsLoupeTF.bed
    # wc -l OverlapsLoupeTF.bed
    # 139 OverlapsLoupeTF.bed

############## ############## ############## Overlap analysis
  ############## ChIP Atlas
    # Read in chip peaks overlapping network enhancers
    overlaps.atlas=read.table("ChIP_Atlas/Overlaps.Atlas.bed")
    colnames(overlaps.atlas)=c(paste0("ChIP_", c("chr", "start", "end", "desc", "score", "strand" ,"V7", "V8", "V9") ),
                               paste0("Enh_", c("chr", "start", "end", "id", "score", "strand") ))
                               
    overlaps.atlas$Name = str_extract(overlaps.atlas$ChIP_desc, "(?<=Name=)[^;]*")
    
    overlaps.atlas$Name = gsub("%20", "", overlaps.atlas$Name, fixed=TRUE)
    overlaps.atlas$Name = gsub("(", "", overlaps.atlas$Name, fixed=TRUE)
    overlaps.atlas$Name = gsub(")", "", overlaps.atlas$Name, fixed=TRUE)
    names=transpose(strsplit(overlaps.atlas$Name, split="@", fixed=TRUE))
    overlaps.atlas$TF=names[[1]]
    overlaps.atlas$CellType=names[[2]]
    
    #Keep ChIP for glial cells only 
    overlaps.atlas=overlaps.atlas[which(overlaps.atlas$CellType%in%c("iPSCderivedastrocytes",  
                                                                     "SF268", 
                                                                     "U-87MG", 
                                                                     "hTERTRPE-1")),]
    
    #Add TRUE/FALSE columns to the edges data frame 
    edges_tf$Atlas=edges_tf$TF%in%atlas$Track.type # TF is present in the ChIP Atlas
    edges_tf$AtlasOverlap=FALSE # The linked enhancer overlaps the TF chip peak (filled in in the for loop below)
    edges_tf$AtlasOverlapData=NA # Store overlap data decriptors such as cell type etc. (filled in in the for loop below)
    
    agg_atlas = aggregate(Cell.type ~ Track.type, data = atlas, FUN = function(x) paste(x, collapse = ";"))
    m=match(edges_tf$TF, agg_atlas$Track.type)
    edges_tf$Atlas=agg_atlas$Cell.type[m]
    
    Nedges.atlas=0
    Nvalidated.atlas=0
    t=unique(intersect(edges_tf$TF, atlas$Track.type))
    for(k in  c(1:length(t)))
    {
      tf=t[k]
      e=edges_tf[grep(tf, edges_tf$TF) , ] 
      o=overlaps.atlas[grep(tf, overlaps.atlas$ChIP_desc) , ]  
      
      if (nrow(o)>0)  
      {
      o$pair=paste(tf, o$Enh_id, sep="_"); 
      o=o[which(o$pair%in%e$pair), ]
      }
      if (nrow(o)>0) 
      {
      edges_tf$AtlasOverlap[which(edges_tf$pair%in%o$pair)]=TRUE
      agg_o= aggregate(Name ~ pair, data = o, FUN = function(x) paste(x, collapse = ";"))
      m=match(agg_o$pair, edges_tf$pair)
      edges_tf$AtlasOverlapData[m]=agg_o$Name
      }
      Nedges.atlas=Nedges.atlas+nrow(e)
      Nvalidated.atlas=Nvalidated.atlas+length(unique(o$pair))
    }
    
    t.edges_in_atlas=as.data.frame(table(edges_tf$TF[which(is.na(edges_tf$Atlas)==FALSE)]))
    t.valid_in_atlas=as.data.frame(table(edges_tf$TF[which(edges_tf$AtlasOverlap==TRUE)]))
    
    colnames(t.edges_in_atlas)=c("TF", "Nedges")
    t.edges_in_atlas$NvalidatedAtlas=t.valid_in_atlas$Freq[match(t.edges_in_atlas$TF,
                                                                 t.valid_in_atlas$Var1)]
    
    summary_atlas=t.edges_in_atlas

  ################## Loupe et al.
    overlaps.loupe=read.table("Loupe/OverlapsLoupeTF.bed")
    colnames(overlaps.loupe)=c(paste0("ChIP_", c("chr", "start", "end", "id", "score", "strand") ),
                               paste0("Enh_", c("chr", "start", "end", "id", "score", "strand") )
    )
    m=match(overlaps.loupe$ChIP_id, tfdata$Id)
    overlaps.loupe=cbind(overlaps.loupe, tfdata[m,])
    
    edges_tf$Loupe=edges_tf$TF%in%ab$Target
    edges_tf$LoupeOverlap=FALSE
    edges_tf$LoupeOverlapData=NA
    
    Nedges.loupe=0
    Nvalidated.loupe=0
    t=unique(intersect(edges_tf$TF, ab$Target))
    for(k in  c(1:length(t)))
    {
      tf=t[k]
      e=edges_tf[grep(tf, edges_tf$TF) , ]  
      o=overlaps.loupe[grep(tf, overlaps.loupe$TF) , ]  
      if (nrow(o)>0)  
      {
        o$pair=paste(tf, o$Enh_id, sep="_"); 
        o=o[which(o$pair%in%e$pair), ]
      }
      
      if (nrow(o)>0) 
      {
        edges_tf$LoupeOverlap[which(edges_tf$pair%in%o$pair)]=TRUE
        agg_o= aggregate(Categ ~ pair, data = o, FUN = function(x) paste(x, collapse = ";"))
        m=match(agg_o$pair, edges_tf$pair)
        edges_tf$LoupeOverlapData[m]=agg_o$Categ
      }
      Nedges.loupe=Nedges.loupe+nrow(e)
      Nvalidated.loupe=Nvalidated.loupe+length(unique(o$pair))
      #
     }
    
    t.edges_in_loupe=as.data.frame(table(edges_tf$TF[which(edges_tf$Loupe==TRUE)]))
    t.valid_in_loupe=as.data.frame(table(edges_tf$TF[which(edges_tf$LoupeOverlap==TRUE)]))
    
    colnames(t.edges_in_loupe)=c("TF", "Nedges")
    t.edges_in_loupe$Nvalidatedloupe=t.valid_in_loupe$Freq[match(t.edges_in_loupe$TF,
                                                                 t.valid_in_loupe$Var1)]
    
    summary_loupe=t.edges_in_loupe

############## ############## ############## Final output of comparison to Astronet
  edges_tf$EitherData=(!is.na(edges_tf$Atlas)) + edges_tf$Loupe
  edges_tf$EitherDataValidated=edges_tf$AtlasOverlap + edges_tf$LoupeOverlap
  edges_tf.chip=edges_tf[which(edges_tf$EitherData >0 ) , ]
  
  #write.csv(edges_tf.chip, "edges_tf.chip.csv")
  
  complexes=unique(edges$Node1[grep("/", edges$Node1)])
  # "FOS/JUN"           "GLI2/ZBTB7B"       "IRF2/IRF3"         "KLF15/KLF9"       
  # "MXI1/TGIF1"        "RORA/RORB"         "SOX21/SOX9"        "TEAD1/TEAD2/TEAD3"
  
  unique(edges_tf.chip$TF)
  # "FOS"    "JUN"    "ZBTB7B" "MEIS2"  "NFIB"   "REST"   "SOX9"   "TEAD1"  "TP53" 
  # For FOS/JUN where both TFs are present in the ChIP validation data, aggregate the data for the two TFs
  
  edges_tf.chip$pair=gsub("JUN_", "FOS/JUN_", edges_tf.chip$pair)
  edges_tf.chip$pair=gsub("FOS_", "FOS/JUN_", edges_tf.chip$pair)
  
  validation=edges_tf.chip[, c("pair" , "Atlas"  , "AtlasOverlap",   "Loupe" ,"LoupeOverlap"  ,"EitherData", "EitherDataValidated")]
  validation$Atlas=1-as.numeric(is.na(validation$Atlas))
  validation_agg= aggregate(. ~ pair, data = validation, FUN = sum, na.rm = TRUE)


