####Prepare MAGMA GENELOC_FILE for distances 
remove(list=ls());gc()
source("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/2.Scripts/Functions.R")

###########################################################################
#Script for preparing Loc and Set files used for EGrf MAGMA analyses 
############################################################################

#Prepare background data including Predicted enhancers corresponding to EGPs with RPKM > 0.5 (loc. file)
All_EGPs <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")
All_EGPs=All_EGPs[All_EGPs$Gene.RNAseq_RPKM>= 0.5,]
All_EGPs=All_EGPs[c("Pair", "Enh.Pos", "TSS", "Enh.start", "Enh.end")]
All_EGPs$chr= splitter(All_EGPs$Enh.Pos, ":", 1)
All_EGPs <- All_EGPs[,c("Pair", "chr","Enh.start", "Enh.end")]
All_EGPs$chr=gsub("chr","",All_EGPs$chr)
All_EGPs$Pair=sub("^(.*?_[0-9]+)_.*", "\\1", All_EGPs$Pair)
All_EGPs=unique(All_EGPs)

dim(unique(All_EGPs))

#write loc file
write.bed(All_EGPs, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/Allintergenic_peaks_background_0.5RPKM.loc")

#####Prepare set files: Grouping Predicted Enhancers according to distance with closest gene:  0 to 50, 50 to 200 and 200 to 500KB

#Load Crispri hits
res.final <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv") 
#Load EGRF
Predictions_RF <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")

##Check how many CRISPRi screening peaks are among our EGRF data
table(res.final[res.final$HitPermissive==TRUE,]$Enh.Pos %in% Predictions_RF[Predictions_RF$pass_rf==TRUE & Predictions_RF$Gene.RNAseq_RPKM>0.5,]$Enh.Pos)
table(Predictions_RF[Predictions_RF$pass_rf==TRUE & Predictions_RF$Gene.RNAseq_RPKM>0.5,]$Enh.Pos %in% res.final$Enh.Pos)
table(Predictions_RF[Predictions_RF$pass_rf==TRUE,]$Enh.Pos %in% res.final$Enh.Pos)

#Provide ATAC id for consistency with background set
res.final$Enh=Predictions_RF[match(res.final$Enh.Pos,Predictions_RF$Enh.Pos),]$Enh 

#filter EGRF by those that pass RF and have RPKM >0.5
Predictions_RF=Predictions_RF[Predictions_RF$pass_rf==TRUE & Predictions_RF$Gene.RNAseq_RPKM>0.5,] 

#Remove ATAC Peaks corresponding to ATAC Peaks already tested in the CRISPRi screening
Predictions_RF=Predictions_RF[!Predictions_RF$Enh %in% res.final$Enh,] 

#Create category containing all EGrf predictions (RPKM > 0.5)
Predictions_RF_set=Predictions_RF
Predictions_RF_set$category="All_Egrf_predictions" 


#Create set containing EGRF within 2-50kb
Predictions_RF_0_50=Predictions_RF[Predictions_RF$Distance < 50000,] #note that peaks were previously filtered by > 2 KB, so there are no peaks within 0-2KB range
Predictions_RF_0_50$category="0-50kb"
dim(unique(Predictions_RF_0_50[c("category","Enh")]))

#Create set containing EGRF within 50-200kb
Predictions_RF_50_200=Predictions_RF[Predictions_RF$Distance >= 50000 & Predictions_RF$Distance < 200000,]
Predictions_RF_50_200$category="50_200KB"
dim(unique(Predictions_RF_50_200[c("category","Enh")]))

#Create set containing EGRF within 200-500kb
Predictions_RF_200_500=Predictions_RF[Predictions_RF$Distance >= 200000 & Predictions_RF$Distance <= 500000,]
Predictions_RF_200_500$category="200_500KB"
dim(unique(Predictions_RF_200_500[c("category","Enh")]))

#Gather all the groups
set_data=rbind(unique(Predictions_RF_0_50[c("category","Enh")]), unique(Predictions_RF_50_200[c("category","Enh")]),
               unique(Predictions_RF_200_500[c("category","Enh")]), unique(Predictions_RF_set[c("category","Enh")]))


#Perform some sanity checks
table(set_data$category)

table(set_data[set_data$category=="0-50kb",]$Enh %in% 
        Predictions_RF[Predictions_RF$pass_rf==TRUE & Predictions_RF$Gene.RNAseq_RPKM>0.5 & Predictions_RF$Distance <50000,]$Enh)

table(set_data[set_data$category=="50_200KB",]$Enh %in% 
        Predictions_RF[Predictions_RF$pass_rf==TRUE & Predictions_RF$Gene.RNAseq_RPKM>0.5 & Predictions_RF$Distance >= 50000 & Predictions_RF$Distance < 200000,]$Enh)


table(set_data[set_data$category=="200_500KB",]$Enh %in% 
        Predictions_RF[Predictions_RF$pass_rf==TRUE & Predictions_RF$Gene.RNAseq_RPKM>0.5 & Predictions_RF$Distance >= 200000,]$Enh)

#All Looks good!


#Save file as a set  file
write.bed(set_data, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/EGrf_distance_05RPKM.set")

