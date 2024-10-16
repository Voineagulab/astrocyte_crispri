###########################################################################################
##Sanity checks: Distance density plot between Egrf Predicted enhancer midpoint and the genes TSS
##########################################################################################


## Setup
remove(list=ls());gc()
setwd("/Users/javierperez/Documents/PostDoc_UNSW/UNSW/CROPseq/revision/all_intergenic_ATACpeaks/Fig3A")
library(tidyverse)
library(rcartocolor)
library(ggplot2)
invis <- element_blank()
maxh <- 29.7 / 2.54
maxw <- 21.0 / 2.54
source("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Scripts/Functions.R")
source("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/FinalFigureFunctions.R")

## Read EGRF data
All_EGPs <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")

#read Crispri hits
res.final <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv") 
length(unique(res.final$Enh.Pos))
# #Remove CRISPRi screening peaks from EGRF data
# All_EGPs=All_EGPs[!All_EGPs$Enh.Pos %in% res.final$Enh.Pos,] 

#Define significant EGRF predictions as those peaks successfully predicted by EGRF with a RPKM >0.5, and not in CRISPRi screen peaks
All_EGPs$EGRF.sign <- ifelse(
  All_EGPs$pass_rf == TRUE & 
    All_EGPs$Gene.RNAseq_RPKM > 0.5 & 
    !All_EGPs$Enh.Pos %in% res.final$Enh.Pos, 
  TRUE, 
  FALSE
)

#sanity checks
table(All_EGPs$EGRF.sign)

table(All_EGPs[All_EGPs$EGRF.sign==TRUE,]$Enh %in% All_EGPs[All_EGPs$pass_rf==TRUE & All_EGPs$Gene.RNAseq_RPKM>0.5,]$Enh)

table(All_EGPs[All_EGPs$EGRF.sign==TRUE,]$Enh %in% c(All_EGPs[All_EGPs$pass_rf==TRUE & All_EGPs$Gene.RNAseq_RPKM>0.5,]$Enh, 
                                                     All_EGPs[!All_EGPs$Enh.Pos %in% res.final$Enh.Pos,]$Enh))
table(res.final$Enh.Pos %in% All_EGPs$EGRF.sign)

#Create object for plotting and define levels
p <- All_EGPs
p$EGRF.sign <- factor(p$EGRF.sign)
levels(p$EGRF.sign) <-  c("NHA intergenic peaks", "EGrf predictions")
table(p$EGRF.sign)


# Provide median distance for comparing significant EGRF predictions from all other ATAC peaks
m <- aggregate(p$Distance, list(p$EGRF.sign), median)[,2]
m <- m / 1000
names(m) <- levels(p$EGRF.sign)

# pdf(file = "EGPs - Distance Density.pdf", height = 2, width = 2.8)
bandwidth.adjust <- c(1, 0.5, 0.25, 0.1)
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Sanity_checks/EGRF_distance_18Sept.pdf", width = 4, height = 4)
ggplot(p, aes(x = Distance / 1000, colour = EGRF.sign, fill = EGRF.sign)) + 
  geom_density(alpha = 0.25, adjust = bandwidth.adjust[1]) +
  theme_bw() +
  scale_fill_manual(values = pals$Hits, name=NULL) +
  scale_colour_manual(values = pals$Hits_Darker, name=NULL) +
  guides(fill = guide_legend(title = ""), colour = guide_legend(title = "")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.002, 0.004, 0.006, 0.008, 0.010, 0.012)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 510), breaks = c(0, m[2], 100, 200, m[1],300,400, 500),
                     labels = scales::label_number(accuracy = 1)) +
  geom_vline(xintercept = m, linetype = 2, colour = pals$Hits_Darker) +
  labs(x = "Distance between EGP (kb)", y = "Density", title = "") +
  theme(panel.border = invis, panel.grid = invis, legend.position = c(0.7, 0.8),# legend.title = invis,
         axis.text.x = element_text(size=10), axis.text.y=element_text(size=10),
        legend.text = element_text(size=12))


dev.off()
