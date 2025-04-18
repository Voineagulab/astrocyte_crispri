####### This script compares Disease impact score (DIS) values between SNPs within 1kb of functional hits and powered non-hit enhancers
####### and generates the figure displaying variant effect predictions for all variants within Enh427, which regulates CCL2
## ISM data generated using https://hb.flatironinstitute.org/deepsea/?analysis=insilico
## DIS scores generated using https://hb.flatironinstitute.org/sei/
library(pheatmap)
####### Fisher test
####### SNP DIS scores: hits vs powered non-hit
setwd("~/Library/CloudStorage/GoogleDrive-ivlabunsw@gmail.com/My Drive/MANUSCRIPTS_IN PROGRESS/CROPseq_MS/Manuscript/Resubmission_NatNeuro/SuppFigs_tables")
dis=readxl::read_xlsx("Supplementary Table 8 - Deep learning-derived disease scores for variants in functional enhancers.xlsx",
                  sheet="8A_SNPs")
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
hits=res[res$HitPermissive, ]
power=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/PowerCalculationSuppTable.csv")
keep=(power$Hit)|(power$WellPowered015)
power=power[keep, ]

dis.powered=dis[which(dis$Enh%in%power$Enhancer), ]

fisher.test(dis$DIS_max_score > -log10(0.05), dis$Enh%in%hits$Enh)
# Fisher's Exact Test for Count Data
# 
# data:  dis$DIS_max_score > -log10(0.05) and dis$Enh %in% hits$Enh
# p-value = 0.0009354
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.078732 1.353302
# sample estimates:
# odds ratio 
#   1.208186 

# fisher.test(dis.powered$DIS_max_score > -log10(0.05), dis.powered$Enh%in%hits$Enh)
# data:  dis.powered$DIS_max_score > -log10(0.05) and dis.powered$Enh %in% hits$Enh
# p-value = 0.009237
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.037987 1.309945
# sample estimates:
#   odds ratio 
# 1.166016 

####### Enh427
setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/")
rm(list=ls())
library(grid)

########################
enh="Enh427"

#read in data
res=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv")
hits=res[res$HitPermissive, ]
  
# AD variants (p<0.05 in ADSP) that overlap Hit enhancers
# Note: in this bed file the variant position is stored in the "Start" column of the bed file. See : mnt/Data0/PROJECTS/CROPSeq/PublicData/ADSP/DataProcessing.R
overlap=read.table("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/PublicData/ADSP/Processed/adsp.hits.bed", sep="\t")
colnames(overlap)=c("Var_chr", "Var_start", "Var_End", "Var_Id", "Var_pval", "Var_strand","Var_source", "Enh_chr", "Enh_start", "Enh_end", "Enh_Id")
overlap=cbind(overlap, hits[match(overlap$Enh_Id, hits$Enh.Pos), c(1:9)])
write.csv(overlap, "ADvars_HitEnh.csv") 
# Beluga disease impact scores for all nucleotide positions in hit enhancers (DIS ISM)
b=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Beluga/BelugaVariants/BelugaVariants.csv")
# Note: the SNP position needs to be converted from 0-base to 1-base (+1)
b$pos=b$pos+1
b.enh=b[which(b$Enh%in%enh), ]

#ISM data (effects on chromatin marks) for the enhancer
ism=read.table(paste0("/Volumes/share/mnt/Data0/PROJECTS/GWP/crispri/", enh, "_logits.tsv"), header=TRUE)
ism=ism[, c(1:3,grep("Astro", colnames(ism)))]

#### Select variants in Enh 
vars=overlap[which(overlap$Enh%in%enh) , ]
vars=vars[order(vars$Var_pval), ]

#### Convert Variant positions to the relative coordinates in the ISM data 
e=which(overlap$Enh%in%enh)[1]

end=overlap$Enh_end[e]
start=overlap$Enh_start[e]
offset=(1999-(end-start+1))/2

vars$ismRelativePos=ceiling(offset + (vars$Var_start - start))
b.enh$ismRelativePos=ceiling(offset + (b.enh$pos - start))

#### Aggregate the ISM data to get the maximum value at each position for consistency with the DIS score
max.ism=aggregate(ism[, c(4:15)], by=list(ism[,1]), FUN="max")

# Fix colnames(for plotting)
colnames(max.ism)[1]="relative_pos"
colnames(max.ism)=gsub("NH_A_Astrocytes.", "", colnames(max.ism))
colnames(max.ism)=gsub(".None", "", colnames(max.ism))

# Keep only the enhancer coordinates for plotting (i.e remove the ISM flanks)
max.ism=max.ism[which(max.ism$relative_pos%in%b.enh$ismRelativePos), ]

# Add DIS and AD variant pvals to the ISM data
max.ism$DIS=b.enh$DIS_max_score[match(max.ism$relative_pos, b.enh$ismRelativePos)]
max.ism$ADpval=vars$Var_pval[match(max.ism$relative_pos,vars$ismRelativePos)]
max.ism$ADvar=as.numeric(!(is.na(max.ism$ADpval)))

#Plot heatmap
colheat=colorRampPalette(c("red", "white", "blue"))
ann=max.ism[, c("DIS",  "ADvar")]
pdf(paste0(enh, "_ISM.pdf"), height=2.5, width=8)
pheatmap(t(max.ism[,2:13]), color=rev(colheat(200)), cluster_cols = FALSE,
         scale="row", annotation_col = ann, show_colnames = FALSE)
grid.newpage() 
dev.off()

#Export ism data
max.ism$genomic_pos=paste0("chr",unique(b.enh$chr),"_", b.enh$pos[match(max.ism$relative_pos, b.enh$ismRelativePos)])
write.csv(max.ism, paste0(enh, "_maxISM.csv"))

# Data for manuscript text
# AD variants (p<0.05 in ADSP) that overlap hit enhancers
advars=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Resubmission/IV/ADvars_HitEnh.csv")

length(unique(advars$Enh))
#[1] 84
min(advars$Var_pval)
#[1] 2.05754e-06
length(unique(advars$Var_Id))
#[1] 200

t=advars[grep("TSC22D1", advars$Gene) , ]
length(unique(t$Var_Id))
#[1] 12

