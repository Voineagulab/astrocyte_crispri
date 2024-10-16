#######################################################################################
####Sanity checks on TTseq using EGrf predictions 
#######################################################################################

################################################################################################################################################################################################
# TTseq read comparing EGRF predicted vs non predicted peaks
###############################################################################################################################################################################################

###Data preparation
remove(list=ls())

load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/4_EnhancerTranscription/TTseq/FeatureCounts_intergenic.rda")

#create dataframe with positive and negative peaks counts
intergenic_ttseq=cbind(counts_intergenic$Pos$counts, counts_intergenic$Neg$counts)
colnames(intergenic_ttseq)=c("Pos_ttseq","Pos_rnaseq","Neg_ttseq","Neg_rnaseq")
intergenic_ttseq=as.data.frame(intergenic_ttseq)


#prepare data
intergenic_ttseq$TTseq_Total <- intergenic_ttseq$Pos_ttseq + intergenic_ttseq$Neg_ttseq
intergenic_ttseq$RNAseq_Total <- intergenic_ttseq$Pos_rnaseq + intergenic_ttseq$Neg_rnaseq
intergenic_ttseq$NascentEnriched <- intergenic_ttseq$TTseq_Total > intergenic_ttseq$RNAseq_Total


################################################################################################################################ 
## Binary classification of peaks to be expressing eRNA 
################################################################################################################################
threshes <- c(2, 3, 5, 10)

for (use.thresh in threshes) {
  u <- paste0("Unidirectional_", use.thresh) 
  b <- paste0("Bidirectional_", use.thresh)
  c <- paste0("Category_", use.thresh)
  
  intergenic_ttseq[,u] <- (rowSums(intergenic_ttseq[,c("Pos_ttseq", "Neg_ttseq")] >= use.thresh) >= 1) & intergenic_ttseq$NascentEnriched
  intergenic_ttseq[,b] <- (rowSums(intergenic_ttseq[,c("Pos_ttseq", "Neg_ttseq")] >= use.thresh) == 2) & intergenic_ttseq$NascentEnriched
  
  intergenic_ttseq[,c] <- "Not transcribed"
  intergenic_ttseq[which(intergenic_ttseq[,u]),c] <- "Unidirectional"
  intergenic_ttseq[which(intergenic_ttseq[,b]),c] <- "Bidirectional"
  intergenic_ttseq[,c] <- factor(intergenic_ttseq[,c], levels = c("Not transcribed", "Unidirectional", "Bidirectional"))
}


## Save
write.csv(intergenic_ttseq, file = "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Sanity_checks/IntergenicPeaks_Transcriptional_classification.csv", row.names = TRUE)


#######################
TTseq comparison
######################


remove(list=ls());gc()


#load requiered functions
source("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Scripts/Functions.R")
source("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Figs/FinalFigureFunctions.R")


#load data
intergenic_ttseq <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Sanity_checks/IntergenicPeaks_Transcriptional_classification.csv", row.names = 1)

###create Hit column with Predicred and non predicted Peaks
library(readxl)
All_EGPs <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")

#read Crispri hits
res.final <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv") 

#Define significant EGRF predictions as those peaks successfully predicted by EGRF with a RPKM >0.5, and not in CRISPRi screen peaks
predictions=unique(All_EGPs[All_EGPs$pass_rf==TRUE & All_EGPs$Gene.RNAseq_RPKM>0.5 & !All_EGPs$Enh.Pos %in% res.final$Enh.Pos,]$Enh)
length(predictions)

#Provide EGRF predictions into TTseq data
intergenic_ttseq$peak=rownames(intergenic_ttseq)
intergenic_ttseq$pass_rf=ifelse(intergenic_ttseq$peak %in% predictions, TRUE, FALSE)
table(intergenic_ttseq$pass_rf)

#check that intergenic pass rf are the one predicted
table(intergenic_ttseq[intergenic_ttseq$pass_rf==TRUE,]$peak %in% All_EGPs[All_EGPs$pass_rf==TRUE & All_EGPs$Gene.RNAseq_RPKM>0.5 & !All_EGPs$Enh.Pos %in% res.final$Enh.Pos,]$Enh)


## Select threshold
thresh <- 3


fish_trns <- table(intergenic_ttseq$pass_rf, intergenic_ttseq[,paste0("Category_", thresh)] != "Not transcribed") %>% fisher.test()
lab_trns <- paste0("Transcribed: p=",
                   signif(fish_trns$p.value, 2),
                   ", OR=", signif(fish_trns$estimate, 2))

fish_uni <- table(intergenic_ttseq$pass_rf, intergenic_ttseq[,paste0("Category_", thresh)] == "Unidirectional") %>% fisher.test()
lab_uni <- paste0("Unidirectional: p=",
                  signif(fish_uni$p.value, 2),
                  ", OR=", signif(fish_uni$estimate, 2))

fish_bi <- table(intergenic_ttseq$pass_rf, intergenic_ttseq[,paste0("Category_", thresh)] == "Bidirectional") %>% fisher.test()
lab_bi <- paste0("Bidirectional: p=",
                 signif(fish_bi$p.value, 2),
                 ", OR=", signif(fish_bi$estimate, 2))

# title
ttl <- paste0("\n", lab_trns,
              "\n", lab_uni,
              "\n", lab_bi)


#tab <- table(intergenic_ttseq$pass_rf, intergenic_ttseq[,paste0("Category_", thresh)])
tab <- table(intergenic_ttseq$pass_rf, intergenic_ttseq[,paste0("Category_", thresh)])
tab <- tab / rowSums(tab)
tab <- t(tab)
tab <- as.data.frame.matrix(tab)
tab$Transcribed <- rownames(tab)

p <- melt(tab)
colnames(p)[2] <- "pass_rf"
p$pass_rf <- factor(p$pass_rf, levels = c("FALSE", "TRUE"))  
levels(p$pass_rf) <- (c("NHA intergenic peaks", "EGrf Predictions")  )
p$Transcribed <- factor(p$Transcribed, levels = c("Not transcribed", "Unidirectional", "Bidirectional"))

#Barplot based on TTseq information of EGrf predicted peaks
pdf("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Sanity_checks/TTseq_barplot_August_0.5RPKM.pdf", width = 5.5, height = 4.5)  
pal <- c("grey90", pals$Primary[2], pals$Primary_Darker[2])
ggplot(p, aes(x = pass_rf, fill = Transcribed, y = value*100, linetype = Transcribed, alpha = Transcribed)) +
  geom_col(colour = "black", width = 0.7) +
  theme_bw() +
  scale_alpha_manual(values = c(1, 0.9, 1)) +
  theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
        axis.title.x = invis, axis.text.x = element_text(size=12), plot.title = element_text(size = 8),
        axis.text.y = element_text(size=10), axis.title.y = element_text(size=12),  legend.position = "right") +
  scale_fill_manual(values = pal) +
  scale_linetype_manual(values = c("blank", "dashed", "solid")) +
  scale_y_continuous(expand = c(0,1)) +
  labs(y = "Percent of peaks", title = ttl) 

dev.off()

