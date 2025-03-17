#This script define the functions sourced in RunRFmodels.R

###Function for Integrating Tobias Bound and Unbound Transcription Factor Data
addTobiasBound <- function (
    df,            # df: A data frame to which Tobias binding data will be added.
    bound_file,    # bound_file: File path to the file containing data on bound transcription factors.
    unbound_file,  # unbound_file: File path to the file containing data on unbound transcription factors.
    merge_col = "Enh.Pos",  # merge_col: Column name in 'df' used for merging bound and unbound data. Default is "Enh.Pos".
    expressed = T, # expressed: Logical flag to indicate if only expressed transcription factors should be considered. Default is TRUE.
    exp.file = "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/1.Data/InputData/Astrocytes/Expressed_TFs.csv" # exp.file: File path to the file containing data on expressed transcription factors. Used only if 'expressed' is TRUE.
) {
  # Read the bound and unbound overlaps data
  bound_overlaps <- read.table(bound_file)
  unbound_overlaps <- read.table(unbound_file)
  
  # Convert the overlaps data into matrices
  unbound_matrix <- TFmatrix(unbound_overlaps, outnumeric = T)
  bound_matrix <- TFmatrix(bound_overlaps, outnumeric = T)
  
  # Filter matrices for expressed transcription factors if expressed == TRUE
  if (expressed == T) {
    homo_TFs <- read.csv(exp.file)
    bound_matrix <- bound_matrix[,colnames(bound_matrix) %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names | colnames(bound_matrix) == "Group.1"]
    unbound_matrix <- unbound_matrix[,colnames(unbound_matrix) %in% homo_TFs[homo_TFs$TF.1_Expressed & homo_TFs$TF.2_Expressed,]$concat_names | colnames(unbound_matrix) == "Group.1"]
  }
  
  # Modify column names of matrices for clarity
  colnames(unbound_matrix) <- paste0("Tobias.Unbound_", colnames(unbound_matrix))
  colnames(bound_matrix) <- paste0("Tobias.Bound_", colnames(bound_matrix))
  
  # Merge the bound and unbound data with the original data frame
  df <- merge(df, bound_matrix,by.x = merge_col, by.y = "Tobias.Bound_Group.1", all.x = T)
  df <- merge(df, unbound_matrix, by.x = merge_col, by.y = "Tobias.Unbound_Group.1", all.x = T)
  
  # Replace NA values in Tobias columns with 0
  df[is.na(df)] <- 0 
  
  # Add columns to count the number of bound and unbound transcription factors
  if (expressed == T) {
    df$Tobias.Exp_Bound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Bound")])
    df$Tobias.Exp_Unbound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Unbound")])
  } else {
    df$Tobias.Bound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Bound")])
    df$Tobias.Unbound_TF_counts <- rowSums(df[,str_detect(colnames(df), "Tobias.Unbound")])
  }
  
  # Return the modified data frame
  return(df)
}


# TFmatrix function: A function for processing overlaps data and generating a matrix
TFmatrix <- function (overlaps, # The input data frame containing overlap information.
                      outnumeric = FALSE, #If TRUE, the output matrix will contain numeric values instead of logical.
                      expressed_genes = NULL, # A vector of gene names to filter the output matrix columns.
                      TF_col = "V8", #The name of the column containing transcription factor information in 'overlaps'.
                      target_col = "V4") { #The name of the column containing target gene information in 'overlaps'.
  # Display a warning message about potential bugs in the function
  warning("This function was recently changed and results may be bugged")
  
  # Extract the transcription factor names from TF_col using regular expressions
  overlaps[,TF_col] <- sub(".*\\.(.*)_.*","\\1",overlaps[,TF_col])
  
  # Aggregate data based on target gene and transcription factor columns
  agg <- aggregate(by = list(overlaps[,target_col], overlaps[,TF_col]), overlaps[,TF_col], length)
  
  # Spread the aggregated data to create a matrix
  matrix<- spread(agg, Group.2, x)
  
  # Store the names from the first column of the matrix for later use
  names <- matrix[,"Group.1"]
  
  # Replace missing values in the matrix with 0
  matrix[is.na(matrix)] <- 0
  
  # If outnumeric is FALSE, convert the matrix to a logical format
  if (outnumeric == FALSE) {
    matrix[matrix > 0] <- TRUE
    matrix[,2:ncol(matrix)] <- apply(matrix[,2:ncol(matrix)], MARGIN = 2, as.logical)
    matrix[,"Group.1"] <- names
  } 
  
  # If expressed_genes is provided, filter the matrix columns accordingly
  if (! is.null(expressed_genes)) {
    expressed_genes <- c(expressed_genes, "Group.1")
    matrix <- matrix[,colnames(matrix) %in% expressed_genes]
  }
  
  # Return the resulting matrix
  return(matrix)
}

##### Function to extract chromosome, start, end, and size from position data.
extract_pos <- function(data, # data: The dataframe containing the data.
                        pos_col, # pos_col: The name of the column in 'data' that contains the position information. 
                        #The position data is expected to be in the format 'chr<chromosome>:<start>-<end>'.
                        out_name = "")# out_name: It is used to prefix the names of the new columns. 
{
  data[,paste0(out_name, "chr")] <-  sub("chr(.*):.*","\\1",data[,pos_col])
  data[,paste0(out_name, "start")] <- as.numeric(sub(".*:(.*)-.*","\\1",data[,pos_col]))
  data[,paste0(out_name, "end")] <- as.numeric(sub(".*:.*-(.*)","\\1",data[,pos_col]))
  data[,paste0(out_name, "size")] <- data[,paste0(out_name, "end")] - data[,paste0(out_name, "start")]
  return(data)
}
###Function for Integrating Chip data
addChipBigWigs <- function (
    df,         # df: The data frame to which ChIP-seq data will be added.
    filenames,  # filenames: A vector of file names containing ChIP-seq data to be added to the data frame.
    folder      # folder: The directory path where the ChIP-seq files are located.
) {
  # Loop through each file in the filenames list
  for (file in filenames) {
    # Read the data from each file, assuming it's in tabular format
    averages <- read.table(paste0(folder, file))
    
    # Extract the accession number from the file name and convert it to a readable name
    name=sub(".*(ENCFF.*).tab", "\\1", file)
    name <- accession2Name(name)
    
    # Select specific columns from the data and rename them
    averages <- averages[,c(1,6)]  # Assuming columns 1 and 6 are of interest
    colnames(averages) <- c("Enh", name)  # Renaming columns for clarity
    
    # Merge the current file's data with the main data frame
    df <- merge(df, averages, by = "Enh")  
  }
  
  # Return the updated data frame with all ChIP-seq data added
  return(df)
}

# Function to add gene stability and housekeeping gene information to EGPs (Enhancer-Gene Pairs)
addGeneStability <- function(EGPs, # EGPs: DataFrame containing Enhancer-Gene Pair data
                             geneinfo, # geneinfo: DataFrame containing gene information 
                             ensembl# ensembl: Ensembl dataset used for gene data retrieval 
){
  
  # Read housekeeping genes data
  housekeeping_genes <- read.table("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/1.Data/InputData/Eisenberg2013_Genes.txt")
  housekeeping_genes <- data.frame(Gene = housekeeping_genes$V1)
  
  # Fetch synonyms for housekeeping genes from Ensembl
  syns <- getBM(attributes = c("external_synonym", "hgnc_symbol"), filters = c("external_synonym"),
                values = housekeeping_genes$Gene, mart = ensembl)
  colnames(syns) <- c("Gene", "hgnc_symbol")
  
  # Merge housekeeping genes data with synonyms and gene information
  housekeeping_genes <- merge(housekeeping_genes, syns, all.x = TRUE)
  housekeeping_genes <- merge(housekeeping_genes, geneinfo[,c("Gene", "EnsID")], all.x = TRUE)
  
  # Add a column indicating if a gene is a housekeeping gene
  EGPs$Gene.Housekeeping <- EGPs$Gene %in% housekeeping_genes$Gene | 
    toupper(EGPs$Gene) %in% housekeeping_genes$hgnc_symbol | 
    EGPs$EnsID %in% housekeeping_genes$EnsID
  
  # Read Lin2019 Single-Cell gene stability data
  lin_hk <- read_excel("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/Finalised/1.Data/InputData/Lin2019_SingleCell.xlsx")
  lin_hk <- as.data.frame(lin_hk[c("GeneSymbol","Stability index")])
  
  # Fetch synonyms for Lin2019 genes from Ensembl
  syns <- getBM(attributes = c("external_synonym", "hgnc_symbol"), filters = c("external_synonym"),
                values = lin_hk$GeneSymbol, mart = ensembl)
  colnames(syns) <- c("Gene", "hgnc_symbol")
  colnames(lin_hk) <- c("Gene", "Gene.StabilityIndex")
  
  # Merge Lin2019 gene stability data with synonyms and gene information
  lin_hk <- merge(lin_hk, syns, all.x = TRUE)
  lin_hk <- merge(lin_hk, geneinfo[,c("Gene", "EnsID")], all.x = TRUE)
  
  # Merge Lin2019 gene stability information with EGPs based on Gene or hgnc if gene stability is not there
  EGPs <- merge(EGPs, unique(lin_hk[,c("Gene", "Gene.StabilityIndex")]), all.x = TRUE)
  
  # Handle cases where gene stability information is missing in EGPs
  lin_hk_hgnc <- lin_hk[match(unique(lin_hk$hgnc_symbol), lin_hk$hgnc_symbol), c("hgnc_symbol", "Gene.StabilityIndex")]
  NAcols <- merge(EGPs[is.na(EGPs$Gene.StabilityIndex), colnames(EGPs) != "Gene.StabilityIndex"], 
                  unique(lin_hk_hgnc[,c("hgnc_symbol", "Gene.StabilityIndex")]),
                  by.x = "Gene", by.y = "hgnc_symbol", all.x = TRUE)
  EGPs <- rbind(EGPs[!is.na(EGPs$Gene.StabilityIndex),], NAcols)
  
  # Handle cases where EnsID match is not found (kept for reference)
  lin_hk <- lin_hk[match(unique(lin_hk$EnsID), lin_hk$EnsID), c("EnsID", "Gene.StabilityIndex")]
  NAcols <- merge(EGPs[is.na(EGPs$Gene.StabilityIndex), colnames(EGPs) != "Gene.StabilityIndex"], 
                  unique(lin_hk[,c("EnsID", "Gene.StabilityIndex")]),
                  all.x = TRUE)
  EGPs <- rbind(EGPs[!is.na(EGPs$Gene.StabilityIndex),], NAcols)
  
  # Return EGPs with added gene stability information
  return(EGPs)
}

# Function to calculate the number of Transcription Start Sites (TSS) between enhancers and genes
addnumTSSbetween <- function (
    all_rf_factors, # all_rf_factors: A data frame containing information about enhancer-gene pairs. Expected to have columns for gene chromosome ('Gene.chr') and enhancer midpoint ('Enh.Midpoint').
    geneinfo        # geneinfo: A data frame containing information about genes, including their chromosome ('Chr') and TSS location ('TSS').
) {
  # Extract the numeric part of TSS coordinates from geneinfo
  geneinfo$TSS <- sub(".*:", "", geneinfo$TSS)
  
  # Apply a function to each row of all_rf_factors to calculate the number of TSS between enhancers and genes
  numTSSEnhGene <- apply(all_rf_factors, 1, function(pair) {
    # Filter  to include only genes on the same chromosome as the current pair-- Gene.chr changed to Enh.chr
    chr_genes <- geneinfo[geneinfo$Chr == pair["Enh.chr"],]
    
    # Convert enhancer midpoint and TSS to numeric values for comparison
    midpoint <- as.numeric(pair["Enh.Midpoint"])
    TSS <- as.numeric(pair["TSS"])
    
    # Calculate the number of TSS between the enhancer midpoint and gene TSS
    if (TSS > midpoint) {
      # Count TSS that fall between the enhancer midpoint and the gene TSS
      out <- sum(midpoint < chr_genes[,"TSS"] & chr_genes[,"TSS"] < TSS)
    } else {
      # Count TSS that fall outside of this range
      out <- sum(midpoint > chr_genes[,"TSS"] & chr_genes[,"TSS"] > TSS)
    }
    
    # Return the count for each pair
    return(out)
  })
  
  # Return the vector of TSS counts corresponding to each row in all_rf_factors
  return(numTSSEnhGene)
}

# Function to create a BED format data frame with windowed regions around enhancer midpoints
makeWindowBed <- function (
    EGPs,           # EGPs: Data frame containing Enhancer-Gene pair information.
    midpoint.col = "Enh.midpoint",  # midpoint.col: Name of the column in EGPs representing the midpoint of the enhancer. Default is "Enh.midpoint".
    window = 500    # window: Numeric value specifying the size of the window to create around each enhancer midpoint. Default is 500 bases.
) {
  # Create a new data frame 'bed' with BED format columns
  bed <- data.frame(
    Enh.chr = EGPs$Enh.chr,  # Chromosome information for each enhancer
    Enh.start = EGPs[,midpoint.col] - window,  # Start position of the window around the enhancer midpoint
    Enh.end = EGPs[,midpoint.col] + window,    # End position of the window around the enhancer midpoint
    Enh = EGPs$Enh  # Enhancer identifier
  )
  
  # Return the 'bed' data frame
  return(bed)
}

# Function to identify the nearest gene to each enhancer in the EGPs data frame
getGeneNearest <- function (EGPs, geneInfo) { 
  # Initialize lists to store nearest gene names and Ensembl IDs
  nearest_Gene <- list()
  nearest_EnsID <- list()
  
  # Loop through each unique enhancer in the EGPs data frame
  for (enh in unique(EGPs$Enh)) {
    # Extract the midpoint and chromosome of the current enhancer
    enh.mid <- unique(EGPs[EGPs$Enh == enh, "Enh.midpoint"])
    enh.chr <- paste0("chr", unique(EGPs[EGPs$Enh == enh, "Enh.chr"]))
    
    # Filter gene information for genes on the same chromosome as the enhancer
    chr_genes <- geneInfo[geneInfo$Chr == enh.chr,]
    
    # Find the nearest gene based on the minimum distance to the enhancer midpoint
    nearest_Gene[[enh]] <- chr_genes[which.min(abs(enh.mid - chr_genes$TSS)), "Gene"]
    nearest_EnsID[[enh]] <- chr_genes[which.min(abs(enh.mid - chr_genes$TSS)), "EnsID"]
  }
  
  # Convert the lists of nearest genes and Ensembl IDs to vectors
  nearest_EnsID <- do.call("c", nearest_EnsID)
  nearest_Gene <- do.call("c", nearest_Gene)
  
  # Initialize a logical column in EGPs to indicate if a gene is the nearest to an enhancer
  EGPs$Gene.Nearest <- FALSE
  
  # Loop through each row in EGPs to update the 'Gene.Nearest' column
  for (i in 1:nrow(EGPs)) {
    EGPs$Gene.Nearest[i] <- EGPs[,"EnsID"][i] == nearest_EnsID[names(nearest_EnsID) == EGPs$Enh[i]] | 
      EGPs[,"Gene"][i] %in% nearest_Gene[names(nearest_Gene) == EGPs$Enh[i]]
  }
  
  # Return the updated EGPs data frame
  return(EGPs)
}

##### Function to convert ENCODE Astrocyte sample Ids  to a name specifying assay and cell type
accession2Name <- function (mylist) {
  accessionNames <- list(
    "ENCFF499UDS" = "DNASE.astrocyte_hippocampus", 
    "ENCFF901UBX" = "DNASE.astrocyte_spinal_cord",
    "ENCFF382FZE" = "DNASE.astrocyte_cerebellum",
    "ENCFF320CHE" = "DNAse",
    "ENCFF033TTE" = "CHIP.CTCF.astrocyte",
    "ENCFF058GEM" = "Chip.H3K27ac",
    "ENCFF184NZS" =  "Chip.H3K4me3", 
    "ENCFF656QVF" = "Chip.H3K9ac",
    "ENCFF736AZD" = "Chip.H4K20me1",
    "ENCFF424JNY" = "CHIP.CTCF.astrocyte_spinal_cord",
    "ENCFF510YMH" = "CHIP.H3K4me3.astrocyte_spinal_cord",
    "ENCFF072YMW" = "CHIP.H3K4me3.astrocyte_cerebellum",
    "ENCFF757YRI" =  "CHIP.CTCF.astrocyte_cerebellum",
    "ENCFF476HBN" = "Chip.H3K27me3",
    "ENCFF242JDN" = "Chip.H3K4me1",
    "ENCFF702AYE" = "Chip.POLR2A",
    "ENCFF640EWX" = "CHIP.H3K4me3.K562"
  )
  outList <- mylist
  names(outList) <- mylist
  outList[outList[outList %in% names(accessionNames)]] <- accessionNames[outList[outList %in% names(accessionNames)]]
  return(outList)
}


# # This is for running TT.tests & Fisher tests on the entire data.frame
# # remove variables that should not be used for calculating stats
# removed_variables <- c("X", "Enh_Gene","ENSG.targetgene","HitCore", "HitCategory", "Fasta", "Enh",  "Enh.Pos", "Pair", "Gene", "Gene.Pos", "Gene.TSS",
#                        "FDR.N50", "FDR.SCEPTRE", "FDR.N250", "P.N250", "P.N50", "P.NB", "P.SCEPTRE", "X", "X.1",# 
#                        "nCells", "logfc", "logfc.vst", "Enh.chr", "Enh.start", "Enh.end", "Gene.chr", "Gene.start", "Tobias.Intersecting_TFs",
#                        "Gene.end", "Enh.Pos_Gene", "PeakId", "Gene.Distance.Bin", "Exp_Intersecting_TFs", "Intersecting_TFs", "start", "end", "chr", 'TSS', 'EnsID', "Hnisz_SuperEnhancers_OLD")


#Performs a Fisher's Exact Test on a specified column of a data frame 
run.EnhancerFisher <- function(data, data.column, Hit = "HitPermissive") {
  # Count the total number of rows in the data
  total <- nrow(data)
  
  # Convert the specified column to logical values (TRUE/FALSE)
  logical_data <- data[,data.column] %>% as.logical()
  
  # Extract the hit column values
  Hits <- data[,Hit]
  
  # Perform Fisher's Exact Test between the logical data and hits
  ft <- table(logical_data, Hits) %>% fisher.test()
  
  # Compile the output data frame with statistics and test results
  out <- data.frame(
    Variable = data.column,  # Name of the tested variable
    Total_TRUE = sum(Hits),  # Total number of hits
    fraction_total_TRUE = sum(Hits) / total,  # Fraction of hits in total data
    Total_Hit_TRUE = sum(Hits & logical_data),  # Number of true hits for the logical data
    Fraction_Hit_TRUE = sum(Hits & logical_data) / sum(logical_data),  # Fraction of true hits in the logical data
    p = ft$p.value,  # P-value from the Fisher's test
    OR = ft$estimate,  # Odds ratio
    Lower = ft$conf.int[1],  # Lower bound of the confidence interval for the odds ratio
    Upper = ft$conf.int[2],  # Upper bound of the confidence interval
    row.names = data.column  # Set row names to the variable name
  )
  
  # Return the results
  return(out)
}

###########################################
#####Functions used in TrainEGPmodels.R########################
############################################
#Header File 
library(ggplot2)
library(ranger)
library(caret)
library(data.table)
#library(pROC)
library(PRROC)
library(stringr)
library("biomaRt")
library("xgboost")
library("tidyr")
library("rapportools")
library("parallel")
library(readxl)
#install.packages("xgboost")

seed = 135643
set.seed(seed)


####### Functions for executing random forest #######

# Random Forest Functions

getTrainControl <- function(method = "cv",   # method: The resampling method (e.g., cross-validation)
                            number = 10,     # number: Number of resampling iterations
                            modelsTested = 6, # modelsTested: Number of models to test (for seed generation)
                            rf.method) {      # rf.method: Random forest method (not directly used here)
  
  # Setting the seed for reproducibility
  set.seed(seed)
  
  # Creating a list of seeds for each iteration of the training process
  seeds <- vector(mode = "list", length = (number + 1))
  for(i in 1:number) seeds[[i]] <-  sample.int(n=1000, modelsTested )
  seeds[[number+1]] <-  sample.int(n=1000, 1)
  
  # Configuring training control parameters
  tc <-  trainControl(method= method, number= number, savePredictions = TRUE, seeds = seeds)
  return(tc)
}

# Function to train a Random Forest model using 'caret' and 'ranger'
getRFmodel <- function (all_rf_factors, # all_rf_factors: Dataset including features and target variable
                        subset,         # subset: Subset of variables for model training
                        train_var = "Z", # train_var: Target variable name (default "Z")
                        importance = NULL, # importance: File path for saving variable importance
                        hyperparameters=NULL, # hyperparameters: File path for saving best hyperparameters
                        train.method = "ranger", # train.method: Method for training the model (e.g., "ranger")
                        weights = NULL) { # weights: Optional weights for training
  
  # Preparing data for training
  random_forest_data <- all_rf_factors[,c(train_var, subset)]
  colnames(random_forest_data)[colnames(random_forest_data) == train_var] <- "train_var"
  
  # Converting target variable to factor if it contains specific strings
  if(str_detect(train_var, "Hit|Significant")) {
    random_forest_data$train_var <- as.factor(random_forest_data$train_var > 0)
  }
  
  # Setting seed and getting training control parameters
  set.seed(seed)
  train_control <- getTrainControl(rf.method = train.method)
  
  # Training the model using caret package
  model <- train(train_var ~., data=random_forest_data, trControl=train_control, 
                 method= train.method, weights = weights, seed = seed)
  
  # If importance file is provided, calculate variable importance using ranger
  if (!is.null(importance)) {
    set.seed(seed)
    rf <- ranger(train_var ~ ., random_forest_data, seed = seed, mtry = model$bestTune$mtry, 
                 splitrule = model$bestTune$splitrule, num.trees = 500, 
                 min.node.size = model$bestTune$min.node.size, importance = "permutation")
    write.csv(rf$variable.importance[order(rf$variable.importance, decreasing = T)], file = importance)
    write.csv(data.frame(mtry = model$bestTune$mtry, splitrule = model$bestTune$splitrule, min.node.size = model$bestTune$min.node.size), file = hyperparameters)
  }
  return(model)
}

# Function to get predictions from the trained Random Forest model


getRFmodelPreds <- function (all_rf_factors, # all_rf_factors: Dataset for making predictions
                             model,          # model: Trained Random Forest model
                             colname = "rf_pred", # colname: Column name for storing predictions
                             pred.method = "oob", # pred.method: Method for prediction ("oob" or "cv")
                             train_var = NULL, # train_var: Target variable name (required for some prediction methods)
                             probability = F, # probability: Flag for returning class probabilities
                             save_models= NULL) { # save_models: File path to save the model
  
  # Initializing the prediction column
  all_rf_factors[,colname] <- 0 
  vars <- colnames(model$trainingData)[colnames(model$trainingData) != ".outcome"]
  
  # Ensuring the target variable is specified
  if (is.null(train_var)) {
    stop("Must set target var train_var")
  }
  
  # Setting seed and predicting using ranger
  set.seed(seed)
  rf <- ranger(all_rf_factors[,train_var] ~ ., data = all_rf_factors[,vars],num.trees = 1000,
               seed = seed,mtry = model$bestTune$mtry, 
               splitrule = model$bestTune$splitrule, 
               min.node.size = model$bestTune$min.node.size, 
               probability = probability)
  
  # Applying cross-validation or out-of-bag prediction method
  if (pred.method == "cv") {
    all_rf_factors <- crossFoldGeneLevel(all_rf_factors, model, vars = vars,
                                         train_var = train_var, colname  = colname, 
                                         probability = probability)
  } else if (pred.method == "oob") {
    if (probability == F) {
      all_rf_factors[,colname] <- rf$predictions
    } else {
      all_rf_factors[,colname] <- rf$predictions[,"TRUE"]
    }
  }
  
  # Saving the model if required
  if (!is.null(save_models)) {
    saveRDS(rf, file = save_models)
  }
  
  return(all_rf_factors)
}


# Function to perform cross-validation
getFolds <- function(data, #data: Dataset to be divided into folds
                     nfolds = 10) { # nfolds: Number of folds for cross-validation
  # Randomizing data and dividing it into folds
  data <- sample(data, size = length(data))
  foldsize <- length(data) / nfolds
  folds <- list()
  for (i in 1:nfolds) {
    lower <- ceiling((foldsize * (i-1)) + 1)
    upper <- ceiling(foldsize * i)
    folds[[i]] <- data[lower:upper]
  }
  return(folds)
}

# Function to apply cross-validation at gene level
crossFoldGeneLevel <- function (df,         # df: DataFrame containing the data
                                model,      # model: Trained model for making predictions
                                vars,       # vars: Variables used in the model
                                train_var = "logfc.vst", # train_var: Target variable name for training
                                nfolds = 10,             # nfolds: Number of folds for cross-validation
                                colname = "rf_pred",     # colname: Column name for storing predictions
                                probability = F) {      # probability: Flag for returning class probabilities
  # Preparing data and folds
  df_order <- df[,"Pair"]
  set.seed(seed)
  if (probability == T) {
    df[,train_var] <- as.factor(df[,train_var])
  }
  genefolds <- getFolds(unique(df[,"Gene"]), nfolds = 10)
  enhfolds <- getFolds(unique(df[,"Enh"]), nfolds = 10)
  out <- data.frame()
  
  # Running cross-validation
  message(paste0("Running CV performance\nRunning:", nfolds, "x", nfolds, "folds"))
  for (genefold in 1:nfolds) {
    for (enhfold in 1:nfolds) {
      train_set <- df[! df[,"Gene"] %in% genefolds[[genefold]] & ! df[,"Enh"] %in% enhfolds[[enhfold]],]
      test_set <- df[df[,"Gene"] %in% genefolds[[genefold]] & df[,"Enh"] %in% enhfolds[[enhfold]],]
      if (nrow(test_set) > 0) {
        set.seed(seed)
        rf <- ranger(train_set[,train_var] ~ ., data = train_set[,vars], num.trees = 1000,
                     seed = seed, mtry = model$bestTune$mtry, 
                     splitrule = model$bestTune$splitrule, 
                     min.node.size = model$bestTune$min.node.size, 
                     probability = probability)
        pred <- predict(rf, test_set)
        if (probability == F) {
          test_set[,colname] <- pred$predictions
        } else {
          test_set[,colname] <- pred$predictions[,"TRUE"]
        }
        out <- rbind(out, test_set)
      }
    }
  }
  # Adjusting the output
  if (probability == T) {
    out[,train_var] <- as.logical(out[,train_var])
  }
  out <- out[match(df_order, out[,"Pair"]),]
  return(out)
}


#####Comparison functions####
##These are primarily focused on handling specific data structures and formats, and preparing them for analysis with the Random Forest models##


library("biomaRt")
library("cowplot")
library("stringr")
library("seqinr")

###################
#Define variables for each model
###################
#########Astrocytes

# # EGrf variables
# selected_vars_exp <- c("Chip.H3K4me3", "Chip.H3K27ac","Gene.Nearest", 
#                        "Distance","ATACseq.Pileup", "numTSSEnhGene", "Gene.StabilityIndex")
# # Variables for variants of the EGrf model, adding or replacing variables
# selected_vars_exp_replace_ATAC <- c(selected_vars_exp[selected_vars_exp != "ATACseq.Pileup"], "DNAse") #Note this is not quite the same I use bigwig here
# plus_hic_vars <- c(selected_vars_exp, "HiC_interaction_NeuN_neg") #May need to replace this with actual contacts not from ABC pipeline #HiC_interaction_NeuN_neg
# 
# # TAPseq.rf variables
# tap_seq_vars <- c("Distance","HiC_interaction_NeuN_neg", 
#                   "Chip.H3K27ac", "Chip.H3K4me3","DNAse", 
#                   "Chip.POLR2A", "Chip.H3K4me1", "Chip.H3K27me3")
# 
# #########K562
# # EGrf variables (using DNase-seq rather than ATAC-seq for open chromatin state)
# selected_vars_exp_encode <- c(selected_vars_exp[selected_vars_exp != c( "ATACseq.Pileup") ], 
#                               "DNAse.Pileup")
# # TAPseq.rf variables (using the variable names form the K562 dataframe)
# tap_seq_vars_encode<- c(tap_seq_vars[!tap_seq_vars %in% c("DNAse","HiC_interaction_NeuN_neg")],
#                         "DNAse.Pileup", "hic_contact")
# 
# # EGrf-extended variables
# combined_model_encode <- c(selected_vars_exp_encode,"hic_contact",
#                            "normalizedEP300_enhActivity", "activity_prom")  
# # rE2G-extended.rf variables
# encode_ext_model <- read.table("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/ENCODE_rE2GPredictions/ENCODE_Predictors.txt")$V1

#This function trains all models in a list
#Variable list is a named list of variables
# dir.create("Results")
# dir.create("SavedModels")
trainAllModels <- function(all_rf_factors, variables.list, noexp_name = "noexp_rf",probability = F,
                           suffix = "", train_var = "logfc.vst", save_models = F, pred_method = "cv", folder = " ") {
  #OOB with CV splitting unique enhancers and genes for training so test is unique from training (10x10)
  final_suffix <- paste0(suffix, "_",train_var,"_",pred_method )
  
  for (name in names(variables.list)) {
    message(name)
    message(paste("Vars:", paste0(variables.list[[name]], collapse = ",")))
    model<- getRFmodel(all_rf_factors, subset = variables.list[[name]], train_var = train_var, 
                       importance = paste0(folder,"/EG_Variable_Importance_",name,final_suffix ,".csv"), hyperparameters= paste0(folder,"/hyperparameters_",name,".csv")) #Hyperparamters is added by JP to retrieve best hyperparameters deduced by the caret
    all_rf_factors <- getRFmodelPreds(all_rf_factors, model , colname = name, pred.method = pred_method, 
                                      train_var = train_var, probability = probability, 
                                      save_models =  paste0(folder,"/EGPModel",name,final_suffix,".rds"))
  }
  return(all_rf_factors)
}


updateVarNames <- function(names) {
  names <- sub("Chip.","",names)
  names[names == "activity_prom"] <- "Promoter.Activity"
  names[names == "hic_contact"] <- "HiC.Contact"
  names[names =="normalizedEP300_enhActivity"] <- "EP300"
  names[names =="DNAse.Pileup"] <- "DNAse"
  names[names =="Tobias.Exp_Bound_TF_counts"] <- "Bound.TFfootprints"
  names[names =="Tobias.Exp_Unbound_TF_counts"] <- "Unbound.TFfootprints"
  names[names =="TTseq.TTseq_Total"] <- "TTseq.Counts"
  names[names =="TTseq.Ratio_TTversusRNA"] <- "TTseq.RNAseqRatio"
  return(names)
}
getVarsExcel<- function (varslist) {
  
  variables <- data.frame(variables = sort(unique(unlist(varslist))))
  for (model  in names(varslist)) {
    variables[,model] <- variables$variables %in% varslist[[model]]
  }
  variables$variables<- updateVarNames(variables$variables)
  return(variables)
}





