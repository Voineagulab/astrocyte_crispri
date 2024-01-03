##########
##
##Aggregate Signal plots in R
##
## @author: Gavin Sutton & Sam Bagot
## @23-11-01
#########

setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels")
library(data.table)
library("ggplot2")
source("Scripts/Header_functions.R")
filenames <- list.files("Results/Tobias/Footprint/Aggregate_txt/", full.names = T)

res <- lapply(filenames, function(x) { 
  data <- as.data.table(read.table(x))
  rows <- nrow(data)
  colnames(data) <- c("signal", "TFBS", "type", "regions", "values")
  data <- data[ , list( values = unlist(strsplit(values, ","))) , by = c("signal","TFBS","type", "regions") ]
  data$regions <- factor(data$regions, levels = c("Hits", "Non-Hits", "All_Peaks"))
  data$TF <- sub(".*\\.(.*?)_.*","\\1",data$TFBS)
  data$bound_status <- sub(".*_(.*)","\\1",data$TFBS)
  data$values <- as.numeric(data$values)
  #MANually encoding 60bp window because that is default
  data$pos <- rep(-59:60, rows)
  #data[ , list( pep = unlist( strsplit( pep , ";" ) ) ) , by = pro ]
  #data$values <- lapply(strsplit(data$values, ","), as.numeric)
  return(data)
})
names(res) <- sub(".*../(.*).txt","\\1",filenames)

saveRDS(res, file = "Results/Tobias/Summaries/Aggregated_Signal.rds")

