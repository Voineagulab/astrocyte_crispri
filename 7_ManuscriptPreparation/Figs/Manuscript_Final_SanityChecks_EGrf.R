########################################################################################################################################################################
#Script for ploting Fig 7c Using EGrf predictions
########################################################################################################################################################################

remove(list=ls());gc()

setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub")
library(tidyverse)
library(rcartocolor)
library(ggsci)
library(ggplot2)
library(ggbeeswarm) 
library(readxl)

source("Scripts/6.Predictive_models/Finalised/2.Scripts/Functions.R")
source("Scripts/6.Predictive_models/Finalised/2.Scripts/FigureFunctions.R")

#Load TAD overlapping data
tad <- read.csv("Scripts/6.Predictive_models/Sanity_checks/Pair-TAD Annotation_allintergenic.csv", row.names = 1)

#Load Superenhancer overlapping data
superenhancers=read.csv("Scripts/6.Predictive_models/Sanity_checks/SuperEnh_Vs_UniqueIntergenicPeaks.csv")

#Load ovlerapping known data
Yao=read.csv("Scripts/6.Predictive_models/Sanity_checks/UniqueIntergenicPeaks_vs_YaoEnhancers.csv")

#Load CRISPRi screen data
res.final <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv") #Load Crispri hits

##Load EGrf data
All_EGPs <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")

#Define significant EGRF predictions as those peaks successfully predicted by EGRF with a RPKM >0.5, and not overlapping with any CRISPRi screen peaks
All_EGPs$EGRF.sign <- ifelse(
  All_EGPs$pass_rf == TRUE & 
    All_EGPs$Gene.RNAseq_RPKM > 0.5 & 
    !All_EGPs$Enh.Pos %in% res.final$Enh.Pos, 
  TRUE, 
  FALSE
)

#Add supernehancer and known enhancer data
All_EGPs$superenhancer=ifelse(All_EGPs$Enh %in% superenhancers$V4, TRUE, FALSE)
All_EGPs$knownenh=ifelse(All_EGPs$Enh %in% Yao$V4, TRUE, FALSE)

#Extract enh level data for each EGP
enh_level_data=unique(All_EGPs[c("Enh","superenhancer","knownenh", "EGRF.sign")])
table(enh_level_data$EGRF.sign)

##add rpkm values to TAD data and define pass filter as predictions for which genes had RPKM >= 0.5 and peaks were not tested in CRISPRi screening
tad$geneRPKM=All_EGPs$Gene.RNAseq_RPKM[match(All_EGPs$Pair, tad$Pair)]
tad$EGRF.sign=ifelse(tad$pass_rf==TRUE & tad$geneRPKM >= 0.5 & !tad$Enh.Pos %in% res.final$Enh.Pos, TRUE, FALSE)

table(tad[tad$EGRF.sign==TRUE,]$Enh.Pos %in% All_EGPs[All_EGPs$EGRF.sign==TRUE,]$Enh.Pos) #sanity check

#define variables for storing fisher test statistics
x_ttl_final=character()
y_ttl_final=character()
z_ttl_final=character()


x <- table(tad$EGRF.sign, tad$CrossTAD)
fish <- fisher.test(x)
x <- x / rowSums(x) * 100
x <- as.data.frame(x)
x <- x[which(as.logical(x$Var2)),]


x$Resource <- "TAD\nBoundary (EGP-level)"

x_ttl <- paste0("p=", signif(fish$p.value, 1), "\n", 
                "OR=", signif(fish$estimate, 2))
x_ttl_final=c(x_ttl_final,x_ttl)


# Prepare superenhancers data
y <- table(enh_level_data$EGRF.sign, enh_level_data$superenhancer)
fish=fisher.test(y) ###Check that provide same statistic than previous computed fisher test (see candidate.enrich)
y <- y / rowSums(y) * 100
y <- as.data.frame(y)
y <- y[which(as.logical(y$Var2)),]
# levels(y$Var1) <- c("Inactive\ncandidates", "Functional\nenhancers")
y$Resource <- "Super-\nenhancer"

y_ttl <- paste0("p=", signif(fish$p.value, 2), "\n", 
                "OR=", signif(fish$estimate, 2))
y_ttl_final=c(y_ttl_final,y_ttl)


# Prepare known enhancers data
#z <- candidate.annot[which(candidate.annot$Tested),]
z <- table(enh_level_data$EGRF.sign, enh_level_data$knownenh)
fish=fisher.test(z)###Check that provide same statistic than previous computed fisher test (see candidate.enrich)
z <- z / rowSums(z) * 100
z <- as.data.frame(z)
z <- z[which(as.logical(z$Var2)),]


z$Resource <- "K562 enhancer"
z_ttl <- paste0("p=", signif(fish$p.value, 2), "\n", 
                "OR=", signif(fish$estimate, 2))
z_ttl_final=c(z_ttl_final,z_ttl)

## Plot
po<- rbind(x,y,z)
po$Var1 <- factor(po$Var1)
levels(po$Var1) <- c("NHA intergenic peaks", "EGrf predictions")
po$Resource <- gsub(" ", "\n", po$Resource)

#This plot Fig 6C)
pdf("Scripts/6.Predictive_models/Sanity_checks/NHA_EGRF_OverlapBarplot_0.5RPKM.pdf",width = 3.5, height = 4)
ggplot(po, aes(x = Resource, y = Freq, colour = Var1, fill = Var1)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_colour_manual(values = pals$Hits_Darker) +
  scale_fill_manual(values = pals$Hits) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y=element_text(size =14), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(size = 10),
        axis.text.x = element_text(size=12), axis.text.y=element_text(size=12),
        legend.position = c(0.35,0.7), legend.text = element_text(size=11)) + 
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0), breaks = seq(0, 60, by = 10)) +
  labs(y = "Percentage", title = "", color = "", fill = "") +
  
  geom_text(data = subset(po, Resource == "K562\nenhancer"), 
            aes(x = 1, y = 14, label = z_ttl_final[1]),
            color = "gray30", # Change color as desired
            size = 4) + # Adjust text size as needed
  geom_text(data = subset(po,Resource == "Super-\nenhancer"), 
            aes(x = 2, y = 22, label = y_ttl_final[1]),
            color = "gray30", # Change color as desired
            size = 4) +
  geom_text(data = subset(po, Resource == "TAD\nBoundary\n(EGP-level)"), 
            aes(x = 3, y = 54, label = x_ttl_final[1]),
            color = "gray30", # Change color as desired
            size = 4) 
  # geom_segment(aes(x = 0.8, xend = 1.2, y = 10.5, yend = 10.5), color = "black") +
  # geom_segment(aes(x = 0.8, xend = 0.8, y = 10, yend = 10.5), color = "black") +
  # geom_segment(aes(x = 1.2, xend = 1.2, y = 10, yend = 10.5), color = "black") +
  # 
  # geom_segment(aes(x = 1.8, xend = 2.2, y = 14.5, yend = 14.5), color = "black") +
  # geom_segment(aes(x = 1.8, xend = 1.8, y = 14, yend = 14.5), color = "black") +
  # geom_segment(aes(x = 2.2, xend = 2.2, y = 14, yend = 14.5), color = "black") +
  # 
  # geom_segment(aes(x = 2.8, xend = 3.2, y = 52.5, yend = 52.5), color = "black") +
  # geom_segment(aes(x = 2.8, xend = 2.8, y = 52, yend = 52.5), color = "black") +
  # geom_segment(aes(x = 3.2, xend = 3.2, y = 52, yend = 52.5), color = "black") 


dev.off()
