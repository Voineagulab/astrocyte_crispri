remove(list=ls());gc()
library(ggplot2)

#####################################################################
##STable 12b: Magma gene-set level analyses 
#####################################################################

#Load magma outputs
data_path="/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/AstroREG_Extended/Outputs/MAGMA_queries/EGrf_Distance_05RPKM_queries"
files=file.path(data_path, list.files(data_path, pattern = "gsa.out"))

data_list= lapply(files, function(file){
  
  #read data
  data_0=read.table(file, header=TRUE, quote="\"")
  #extract dataset info from filename
  dataset_name <- strsplit(basename(file), "_")[[1]][4]
  #add dataset column
  data_0$dataset=dataset_name
  return(data_0)
})


# Combine all data frames into one
data_snpmean <- do.call(rbind, data_list)
data_snpmean$model="mean"
data=data_snpmean



#Prepare data for plotting
data_crispri.egrf=data[data$model=="mean",]
data_crispri.egrf$dataset=gsub(".gsa.out","", data_crispri.egrf$dataset)

#Extract full GWAS dataset names
library(readxl)
GWAS <- read_excel("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/GWAS.xlsx")
data_crispri.egrf$dataset_names=GWAS[match(data_crispri.egrf$dataset, GWAS$`Study pheno-type`),]$Description_JP

#renaming dataset 
data_crispri.egrf$dataset_names[data_crispri.egrf$dataset_names=="Attention deficit/hyperactivity disorder"]="Attention deficit\nhyperactivity disorder"

#sort data to plot
data_crispri.egrf$order=GWAS[match(data_crispri.egrf$dataset, GWAS$`Study pheno-type`),]$Order_plot
data_crispri.egrf=data_crispri.egrf[order(data_crispri.egrf$order, decreasing = TRUE),]

#set order
data_crispri.egrf$dataset_names <- as.character(data_crispri.egrf$dataset_names)
data_crispri.egrf$dataset_names <- factor(data_crispri.egrf$dataset_names, levels=unique(data_crispri.egrf$dataset_names))

#filter by set groups to be plotted
table(data_crispri.egrf$VARIABLE)

data_crispri.egrf=data_crispri.egrf[data_crispri.egrf$VARIABLE %in% c("All_Egrf_predictions","0-50kb","200_500KB",
                                                                      "50_200KB"),]
table(data_crispri.egrf$VARIABLE)
# data_crispri.egrf=data_crispri.egrf[data_crispri.egrf$VARIABLE %in% c("AstroREG_0-50Kb","AstroREG_200-500Kb",
#                                                                       "AstroREG_50-200Kb","All_AstroReg"),]

#Leave only 1 insomnia dataset
data_crispri.egrf=data_crispri.egrf[data_crispri.egrf$dataset != "INS",]

#renaming SANG insomnia DATASET 
data_crispri.egrf$dataset_names[data_crispri.egrf$dataset_names=="INSang"]="Insomnia"

#re-order and re-name columns
data_crispri.egrf=data_crispri.egrf[c(1,3,4,5,6,7,9,8,10)]
colnames(data_crispri.egrf)=c("Gene set", "NEnhancers", "BETA", "BETA_STD", "SE", "P", "Model", "Dataset_code", "Dataset_fullname")
data_crispri.egrf$Model="SNP-wise Mean"
data_crispri.egrf$`Gene set`= ifelse(data_crispri.egrf$`Gene set`=="0-50kb", "EGrf (2−50kb)", ifelse(data_crispri.egrf$`Gene set`=="50_200KB", "EGrf (50−200kb)", 
                                                                       ifelse(data_crispri.egrf$`Gene set`=="200_500KB", "EGrf (200−500kb)", "EGrf")))

#Export Table
write.csv(data_crispri.egrf, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/Stable12/STable12b_MAGMA.csv")




#####################################################################
##STable 12c: Magma gene level analyses 
#####################################################################


#Load magma outputs
data_path="/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/Outputs/MAGMA_bkg/General_Background_bkg/"
files=file.path(data_path, list.files(data_path, pattern = "genes.out"))

data_list= lapply(files, function(file){
  
  #read data
  data_0=read.table(file, header=TRUE, quote="\"")
  #extract dataset info from filename
  dataset_name <- strsplit(basename(file), "_")[[1]][5]
  #add dataset column
  data_0$dataset=dataset_name
  return(data_0)
})


# Combine all data frames into one
data_genelevel <- do.call(rbind, data_list)



#Load sete data
set_data=read.delim("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/EGrf_distance_05RPKM.set", header = FALSE)

#Filter MAGMA gene level analyses to keep only EGrf predictions
data_genelevel=data_genelevel[data_genelevel$GENE %in% set_data[set_data$V1 =="All_Egrf_predictions",]$V2,]

#sanity check: verify Peaks correspond only to EGrf predictions
length(unique(data_genelevel$GENE))
table(data_genelevel$GENE %in% set_data[set_data$V1 =="All_Egrf_predictions",]$V2)


#change dataset names

#Extract full GWAS dataset names
library(readxl)
GWAS <- read_excel("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/GWAS.xlsx")

#renaming dataset 
data_genelevel$dataset=gsub(".genes.out","", data_genelevel$dataset)
data_genelevel$dataset=GWAS[match(data_genelevel$dataset, GWAS$`Study pheno-type`),]$Description_JP
data_genelevel$dataset[data_genelevel$dataset=="Attention deficit/hyperactivity disorder"]="Attention deficit-hyperactivity disorder"

#Leave only 1 insomnia dataset
data_genelevel=data_genelevel[data_genelevel$dataset != "Insomnia",]

#renaming SANG insomnia DATASET 
data_genelevel$dataset[data_genelevel$dataset=="INSang"]="Insomnia"
table(data_genelevel$dataset)


#change column names and save data
colnames(data_genelevel)[1] =("EGRF Enhancer")

#Load EGRF
Predictions_RF <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/3.Predictions/RF_Results/Astrocytes_All_Intergenic_Predictions/All_EGPs.csv")

#Load Crispri hits
res.final <- read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Final.csv") 

#filter EGRF by those that pass RF and have RPKM >0.5
Predictions_RF=Predictions_RF[Predictions_RF$pass_rf==TRUE & Predictions_RF$Gene.RNAseq_RPKM>0.5,] 

#Remove ATAC Peaks corresponding to ATAC Peaks already tested in the CRISPRi screening
Predictions_RF=Predictions_RF[!Predictions_RF$Enh %in% res.final$Enh,] 

#FiLlter to just EGRF
Predictions_RF=Predictions_RF[Predictions_RF$Enh %in% data_genelevel$`EGRF Enhancer`,]
Predictions_RF=unique(Predictions_RF[c(1,2)])


#merge information
Predictions_RF=merge(Predictions_RF, data_genelevel, by.x = "Enh", by.y = "EGRF Enhancer", all.x = TRUE)
colnames(Predictions_RF)[1] =("EGRF Enhancer")

#Filter to just AD
Predictions_RF=Predictions_RF[Predictions_RF$dataset=="Alzheimer’s disease",]

write.csv(Predictions_RF, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/Stable12/STable12c_MAGMA.csv")






