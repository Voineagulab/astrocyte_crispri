## In this script, we analyse Naostring validation data

################################################################################################################################ #
## Setup ----

## Setup
  rm(list = ls())
  gc()
  setwd("/mnt/Data0/PROJECTS/CROPSeq/Validation_RNAseq/Nanostring/Results/Final/")
  
## Functions and libraries
  library(NanoTube)
  library(readxl)
  library(tidyverse)

  # to process character vectors and split at desired values
  splitter <- function(
    x, # x is a character vector
    split, # split is the character at which the split shall occur (and a split occurs at EVERY match)
    side # side is which side of the split to be returned, any value from 1 - X where X is the number of splits + 1
  ) { 
    sapply(strsplit(as.character(x), split = split), "[", side)
  }

  # nanostring qc  
  runQC <- function(l) { # l is a list object from readRCC
    
    ## Run positive control QC
      QC_pos <- positiveQC(l)$tab
      colnames(QC_pos) <- c("Sample", "Pos_ScaleFactor", "Pos_R2")
      QC_pos$Pos_Flag <- (QC_pos$Pos_ScaleFactor < 0.3 | QC_pos$Pos_ScaleFactor > 3 | QC_pos$Pos_R2 < 0.95)
      l$samples <- data.frame(l$samples, QC_pos[,-1])
  
    ## QC negative control probes. these should be visualised on a distribution
      QC_neg <- negativeQC(l)$tab
      colnames(QC_neg) <- c("Neg_Mean", "Neg_Max", "Neg_sd", "Neg_BgThresh", "Neg_GenesBelowBg")
      l$samples <- data.frame(l$samples, QC_neg)
  
    ## QC the negative to positive control ratio
      posE <- l$dict$Accession[which(l$dict$Name == "POS_E(0.5)")] # index of the second lowest positive control probe
      posE <- l$exprs.raw[posE,]
      l$samples$LoD_Flag <- posE < l$samples$Neg_BgThresh
      
    ## QC housekeeping genes
      l$samples$Housekeeping_ScaleFactor <- l$hk.scalefactors
      l$samples$Housekeeping_Flag <- (l$hk.scalefactors < 0.1 | l$hk.scalefactors > 10 )
  
    ## Cartridge-based QC metrics
      QC_cart <- as.data.frame(l$qc)
      QC_cart <- data.frame(FOV_Count = as.numeric(QC_cart$FovCount),
                            FOV_Counted = as.numeric(QC_cart$FovCounted),
                            FOV_Ratio = NaN,
                            FOV_Flag = NA,
                            Binding_Dens = as.numeric(QC_cart$BindingDensity),
                            Binding_Flag = NA,
                            Cart_StagePosition = factor(QC_cart$StagePosition),
                            Cart_ID = factor(QC_cart$CartridgeID),
                            Cart_Barcode = factor(QC_cart$CartridgeBarcode))
    
      QC_cart$FOV_Ratio <- QC_cart$FOV_Counted / QC_cart$FOV_Count
      QC_cart$FOV_Flag <- QC_cart$FOV_Ratio < 0.75
      QC_cart$Binding_Flag <- (QC_cart$Binding_Dens < 0.1 | QC_cart$Binding_Dens > 2)
      
      l$samples <- data.frame(l$samples, QC_cart)
      
      l$samples$AnyFlag <- rowSums(l$samples[,grep("Flag", l$samples)]) > 0
      
      ## Output
        return(l)
  }
  
  # to visualise positive controls
  visPosQC <- function(l) {
    ctlProbes <- which(l$dict$CodeClass %in% c("Positive", "Negative")) # index of the second lowest positive control probe
    ctlExp <- l$exprs.raw[ctlProbes,]
    ctlExp <- data.frame(CodeClass = l$dict$CodeClass[ctlProbes],
                         ID = l$dict$Name[ctlProbes],
                         Concentration = NA,
                         ctlExp)
    ctlExp$Concentration <- splitter(ctlExp$ID, "\\(", 2) %>%
      splitter(")", 1) %>%
      as.numeric()
    
    ctlExp <- melt(ctlExp, id.vars = c("CodeClass", "ID", "Concentration"))
    
    ggplot(ctlExp, aes(x = Concentration, y = value, colour = CodeClass)) +
      geom_point() +
      facet_wrap(~variable)
  }
  
  # for ggplot2
  invis <- element_blank()
  
## List to store data
  nano <- list()
  meta <- list()
  
# Files
  # id2enh <- read_xlsx("../../../SelectionOfGuides/Oligos to order for CRISPRi validation.xlsx", sheet = "Sheet1", range = "G1:J18")
  # id2guide <- read.csv("../../../SelectionOfGuides/res.sig.combined.topguides_20bp.csv")
  
  id2enh <- read_xlsx("../../Resources/new_oligos to order for validation.xlsx", sheet = 1, range = "A1:G21")
  id2enh <- as.data.frame(id2enh)
  colnames(id2enh) <- c("Pair", "Guide", "PS", "median.unsuppressed", "suppression", "Gene", "ID")
  id2enh$ID <- gsub("New_|_Fwd_CSO_sgRNA", "", id2enh$ID)
  id2enh$Enh <- splitter(id2enh$Pair, "_", 1)
  id2enh <- relocate(id2enh, c("ID", "Enh", "Gene", "Pair"))
  
## Quick load
  # load("ProcessedData.rda")
  
################################################################################################################################ #
## Read in batch1 ----
  
## Read in and normalise expression
  # paths1 <- list.files("../../Ramaciotti/20230726_210582350125_RCC_GOK12461/", pattern = "300ng", full.names = TRUE) # keep only samples at 300ng concentration, not 150
  paths1 <- "../../Ramaciotti/20230726_210582350125_RCC_GOK12461/"
  
  nano$b1 <- processNanostringData(nsFiles = paths1, # path to directory or file
                                    normalization = "nSolver", # choice of nSolver, RUV, and none
                                    bgType = "threshold", bgThreshold = 2, bgProportion = 0.5, # parameters for nSolver normalisation, remain as default. 
                                    skip.housekeeping = FALSE, # housekeeping normalisation works across samples
                                    includeQC = TRUE, # QC from nanostring rcc file. seems to be coerced to character...
                                    output.format = "list") # or an expdata data type. i prefer to work with lists
  
## Read in metadata
  meta$b1 <- read_xlsx("../../Ramaciotti/SampleOrder_V2.xlsx", sheet = "Run1")
  meta$b1$ChipBatch <- colnames(meta$b1)[1] # batch is encoded by this column name, strangely.
  colnames(meta$b1) <- c("Sample", "Concentration", "LIMS_ID", "Batch")
  meta$b1 <- as.data.frame(meta$b1)

  # trim sample names
  meta$b1$Chip_ID <- gsub("_", "", meta$b1$Sample)
  meta$b1$Sample <- gsub("HD_|NG_", "", meta$b1$Sample)
  meta$b1$Input_ng <- stringr::str_sub(meta$b1$Sample, -5, -3) %>% as.numeric()
  
  meta$b1$Sample <- gsub("_300ng|_150ng", "", meta$b1$Sample)
  meta$b1$Group <- gsub("_A$|_B$", "", meta$b1$Sample)
  
  # extract target and batch information
  meta$b1$Target <- splitter(meta$b1$Sample, "_", 1)
  meta$b1$Batch <- sub(" ", "_", meta$b1$Batch)
  meta$b1$Pool <- grepl("Pool", meta$b1$Sample)
  
  meta$b1$Experimenter <- meta$b1$TargetNo <- NA
  
  for (j in 1:nrow(meta$b1)) {
    nSampleFields <- str_count(meta$b1$Sample[j], "_")
    
    if (nSampleFields == 1) {
      meta$b1$Experimenter[j] <- splitter(meta$b1$Sample[j], "_", 2)
      meta$b1$TargetNo[j] <- 1
    } else {
      meta$b1$Experimenter[j] <- splitter(meta$b1$Sample[j], "_", 3)
      meta$b1$TargetNo[j] <- splitter(meta$b1$Sample[j], "_", 2)
    }

  }
  
  # add enhancer and guide name
  m <- match(meta$b1$Group, id2enh$ID)
  meta$b1$Enh <- id2enh$Pair[m]
  
  # m <- match(meta$b1$Enh, id2guide$X)
  meta$b1$Guide <- id2enh$Guide[m]
  meta$b1$ScreenSuppression <- id2enh$suppression[m]
  
  
  # match order
  m <- match(splitter(colnames(nano$b1$exprs), "_", 3), meta$b1$Chip_ID)
  meta$b1 <- meta$b1[m,]
  meta$b1$ID <- colnames(nano$b1$exprs)
  meta$b1 <- dplyr::relocate(meta$b1, c("ID", "Sample", "Chip_ID", "LIMS_ID", "Enh", "Guide"))
  rownames(meta$b1) <- meta$b1$ID

  nano$b1$samples <- meta$b1

## Run QC
  nano$b1 <- runQC(nano$b1)
  
  # noting that runQC augments the metadata stored in $samples, extract that!
  meta$b1 <- nano$b1$samples
  
################################################################################################################################ #
## Read in batches 2 and 3, as these contain biological replicates for the same sets of samples ----

  
## Read in and normalise expression
  paths23 <- c("../../Ramaciotti/GOK12461 batch 2rpt 20230911_210611690225_RCC/",
             "../../Ramaciotti/GOK12461 Batch3 20230830_210611660225_RCC/")
  
  nano$b23 <- processNanostringData(nsFiles = paths23, # path to directory or file
                                    normalization = "nSolver", # choice of nSolver, RUV, and none
                                    bgType = "threshold", bgThreshold = 2, bgProportion = 0.5, # parameters for nSolver normalisation, remain as default. 
                                    skip.housekeeping = FALSE, # housekeeping normalisation works across samples
                                    includeQC = TRUE, # QC from nanostring rcc file. seems to be coerced to character...
                                    output.format = "list") # or an expdata data type. i prefer to work with lists
  
## Read in metadata
  meta$b23 <- list(read_xlsx("../../Ramaciotti/SampleOrder_V2.xlsx", sheet = "Run2"),
                 read_xlsx("../../Ramaciotti/SampleOrder_V2.xlsx", sheet = "Run3"))
                 
  meta$b23 <- lapply(meta$b23, function(x) {
    x$ChipBatch <- colnames(x)[1] # batch is encoded by this column name, strangely.
    colnames(x) <- c("Sample", "Concentration", "LIMS_ID", "Batch")
    x <- as.data.frame(x)
    return(x)
  })
  
  meta$b23 <- do.call("rbind", meta$b23)
                 
  # trim sample names
  meta$b23$Chip_ID <- meta$b23$LIMS_ID
  meta$b23$Input_ng <- 150
  meta$b23$Sample <- gsub("HD_|NG_", "", meta$b23$Sample)
  meta$b23$Group <- gsub("_A$|_B$", "", meta$b23$Sample)
  
  # extract target and batch information
  meta$b23$Target <- splitter(meta$b23$Sample, "_", 1)
  meta$b23$Batch <- sub(" ", "_", meta$b23$Batch)
  meta$b23$Pool <- grepl("Pool", meta$b23$Sample)
  
  meta$b23$Experimenter <- meta$b23$TargetNo <- NA
  
  for (j in 1:nrow(meta$b23)) {
    nSampleFields <- str_count(meta$b23$Sample[j], "_")
    
    if (nSampleFields == 1) {
      meta$b23$Experimenter[j] <- splitter(meta$b23$Sample[j], "_", 2)
      meta$b23$TargetNo[j] <- 1
    } else {
      meta$b23$Experimenter[j] <- splitter(meta$b23$Sample[j], "_", 3)
      meta$b23$TargetNo[j] <- splitter(meta$b23$Sample[j], "_", 2)
    }

  }
  
  # add enhancer and guide name
  m <- match(meta$b23$Group, id2enh$ID)
  meta$b23$Enh <- id2enh$Pair[m]
  
  # m <- match(meta$b1$Enh, id2guide$X)
  meta$b23$Guide <- id2enh$Guide[m]
  meta$b23$ScreenSuppression <- id2enh$suppression[m]
  
  
  # match order
  m <- match(splitter(colnames(nano$b23$exprs), "_", 3), meta$b23$LIMS_ID)
  meta$b23 <- meta$b23[m,]
  meta$b23$ID <- colnames(nano$b23$exprs)
  meta$b23 <- dplyr::relocate(meta$b23, c("ID", "Sample", "Chip_ID", "LIMS_ID", "Enh", "Guide"))
  rownames(meta$b23) <- meta$b23$ID

  nano$b23$samples <- meta$b23

## Run QC
  nano$b23 <- runQC(nano$b23)
  
  # noting that runQC augments the metadata stored in $samples, extract that!
  meta$b23 <- nano$b23$samples


################################################################################################################################ #
## As per the previous sections, but for batches 4 and 5 ----
  
## Read in and normalise expression
  paths45 <- c("../../Ramaciotti/GOK12461 Batch 4 20230911_210611680225_RCC/",
             "../../Ramaciotti/GOK12461 cartridge 5 20230915_210611700225_RCC/")
  
  nano$b45 <- processNanostringData(nsFiles = paths45, # path to directory or file
                                    normalization = "nSolver", # choice of nSolver, RUV, and none
                                    bgType = "threshold", bgThreshold = 2, bgProportion = 0.5, # parameters for nSolver normalisation, remain as default. 
                                    skip.housekeeping = FALSE, # housekeeping normalisation works across samples
                                    includeQC = TRUE, # QC from nanostring rcc file. seems to be coerced to character...
                                    output.format = "list") # or an expdata data type. i prefer to work with lists
  
## Read in metadata
  meta$b45 <- list(read_xlsx("../../Ramaciotti/SampleOrder_V2.xlsx", sheet = "Run4"),
                 read_xlsx("../../Ramaciotti/SampleOrder_V2.xlsx", sheet = "Run5"))
                 
  meta$b45 <- lapply(meta$b45, function(x) {
    x$ChipBatch <- colnames(x)[1] # batch is encoded by this column name, strangely.
    colnames(x) <- c("Sample", "Concentration", "LIMS_ID", "Batch")
    x <- as.data.frame(x)
    return(x)
  })
  
  meta$b45 <- do.call("rbind", meta$b45)
                 
  # trim sample names
  meta$b45$Sample <- gsub("HD_|NG_", "", meta$b45$Sample)
  meta$b45$Chip_ID <- meta$b45$LIMS_ID
  meta$b45$Input_ng <- 150
  meta$b45$Group <- gsub("_A$|_B$", "", meta$b45$Sample)
  
  # extract target and batch information
  meta$b45$Target <- splitter(meta$b45$Sample, "_", 1)
  meta$b45$Batch <- sub(" ", "_", meta$b45$Batch)
  meta$b45$Pool <- grepl("Pool", meta$b45$Sample)
  
  meta$b45$Experimenter <- meta$b45$TargetNo <- NA
  
  for (j in 1:nrow(meta$b45)) {
    nSampleFields <- str_count(meta$b45$Sample[j], "_")
    
    if (nSampleFields == 1) {
      meta$b45$Experimenter[j] <- splitter(meta$b45$Sample[j], "_", 2)
      meta$b45$TargetNo[j] <- 1
    } else {
      meta$b45$Experimenter[j] <- splitter(meta$b45$Sample[j], "_", 3)
      meta$b45$TargetNo[j] <- splitter(meta$b45$Sample[j], "_", 2)
    }

  }
  
  # add enhancer and guide name
  m <- match(meta$b45$Group, id2enh$ID)
  meta$b45$Enh <- id2enh$Pair[m]
  
  # m <- match(meta$b1$Enh, id2guide$X)
  meta$b45$Guide <- id2enh$Guide[m]
  meta$b45$ScreenSuppression <- id2enh$suppression[m]
  
  
  # match order
  m <- match(splitter(colnames(nano$b45$exprs), "_", 3), meta$b45$LIMS_ID)
  meta$b45 <- meta$b45[m,]
  meta$b45$ID <- colnames(nano$b45$exprs)
  meta$b45 <- dplyr::relocate(meta$b45, c("ID", "Sample", "Chip_ID", "LIMS_ID", "Enh", "Guide"))
  rownames(meta$b45) <- meta$b45$ID

  nano$b45$samples <- meta$b45

## Run QC
  nano$b45 <- runQC(nano$b45)
   
  # noting that runQC augments the metadata stored in $samples, extract that!
  meta$b45 <- nano$b45$samples
  
  
################################################################################################################################ #
## First: visualise sample QC ----


## QC metrics
  ## Batch 1
    p <- nano$b1$samples
    p <- melt(p[,c("Batch", "Neg_BgThresh", "Housekeeping_ScaleFactor", "Pos_ScaleFactor", "Pos_R2", "FOV_Count", "FOV_Ratio", "Binding_Dens")], id.vars = "Batch")
    
    pdf(file = "QC - Batch 1.pdf", height = 6, width = 6.5)
    ggplot(p, aes(x = Batch, y = value, colour = Batch)) +
      geom_violin(draw_quantiles = 0.5, scale = "width") +
      geom_beeswarm() +
      facet_wrap(~variable, scales = "free_y") +
      theme_bw() +
      theme(axis.title.x = invis, axis.text.x = invis, panel.grid = invis, legend.position = c(0.8, 0.2)) +
      labs(y = "Score on QC Metric")
    dev.off()

  ## Batches 2 and 3
    p <- nano$b23$samples
    p <- melt(p[,c("Batch", "Neg_BgThresh", "Housekeeping_ScaleFactor", "Pos_ScaleFactor", "Pos_R2", "FOV_Count", "FOV_Ratio", "Binding_Dens")], id.vars = "Batch")
    
    pdf(file = "QC - Batches 2 and 3.pdf", height = 6, width = 6.5)
    ggplot(p, aes(x = Batch, y = value, colour = Batch)) +
      geom_violin(draw_quantiles = 0.5, scale = "width") +
      geom_beeswarm() +
      facet_wrap(~variable, scales = "free_y") +
      theme_bw() +
      theme(axis.title.x = invis, axis.text.x = invis, panel.grid = invis, legend.position = c(0.8, 0.2)) +
      labs(y = "Score on QC Metric")
    dev.off()
    
  ## Batches 4 and 5
    p <- nano$b45$samples
    p <- melt(p[,c("Batch", "Neg_BgThresh", "Housekeeping_ScaleFactor", "Pos_ScaleFactor", "Pos_R2", "FOV_Count", "FOV_Ratio", "Binding_Dens")], id.vars = "Batch")
    
    pdf(file = "QC - Batches 4 and 5.pdf", height = 6, width = 6.5)
    ggplot(p, aes(x = Batch, y = value, colour = Batch)) +
      geom_violin(draw_quantiles = 0.5, scale = "width") +
      geom_beeswarm() +
      facet_wrap(~variable, scales = "free_y") +
      theme_bw() +
      theme(axis.title.x = invis, axis.text.x = invis, panel.grid = invis, legend.position = c(0.8, 0.2)) +
      labs(y = "Score on QC Metric")
    dev.off()

  

################################################################################################################################ #
## Run statistics ----
  
## Function
  runDE <- function(batches, gene, group, plot = TRUE) {
    
    x <- nano[[batches]]
    
    # get samples to compare
    sampsTest <- which(x$samples$Group == group & x$samples$Input_ng == 150) # the test group, based on the function's input
    targ <- unique(x$samples$Target[sampsTest]) # the main target of the enhancer. this is distinct from the gene variable, which can be set to anything
    sampsBg <- which((x$samples$Group != group) & !(x$samples$Target %in% gene) & x$samples$Input_ng == 150) # the control group, which is also filtered to remove samples with a different group but same expected target gene
    
    # get expression
    exp <- x$exprs[which(x$dict$Name == gene),]
    exp <- log2(exp + 1)
    
    # get design matrix information 
    design <- x$samples[,c("Batch", "Group","FOV_Ratio", "Binding_Dens")]
    design$GroupBinary <- (design$Group == group) %>% as.numeric()
    
    # filter
    exp <- exp[c(sampsTest, sampsBg)]
    design <- design[c(sampsTest, sampsBg),]
    
    # linear model
    formula <- "exp ~ Batch+FOV_Ratio+Binding_Dens+GroupBinary"
    if (batches == "b1") formula <- sub("Batch+", "", formula) # as this level has but a single batch
    
    mod <- lm(as.formula(formula), data = design)
    sum <- summary(mod)
    res <- data.frame(Group = group,
                      Batches = batches,
                      Gene = gene,
                      nGroup = sum(design$GroupBinary),
                      nBg = sum(design$GroupBinary == 0),
                      log2fc = sum$coefficients["GroupBinary", "Estimate"],
                      p = sum$coefficients["GroupBinary", "Pr(>|t|)"])
    
    # output
    if (plot) {
      p <- data.frame(design, Exp = exp)
      p$Group[which(p$GroupBinary == 0)] <- "Other Target"
      
      print(ggplot(p, aes(x = Group, y = Exp, shape = Batch)) +
        geom_beeswarm(colour = "black", size = 3) +
        stat_summary(geom = "point", shape = "-", colour = "red", size = 10, fun = mean) +
        theme_bw() +
        # scale_y_continuous(limits = c(0, NA)) +
        labs (y = "log2(Normalised Expression)") +
        theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis))
      
    }
    
    return(res)
    
  }
  
## Setup comparisons
  comps <- lapply(names(meta), function(x) {
    y <- data.frame(batches = x,
               gene = (nano[[x]]$samples$Target),
                group = (nano[[x]]$samples$Group))
    y <- unique(y)
    return(y)
  })
  comps <- do.call("rbind", comps)
  comps <- comps[-grep("Neg", comps$group),]
  # comps <- comps[-which(comps$batches == "b23" & comps$gene == "FTH1"),] # remove as tested in batch1
  
  # add a custom test
  custom_comp <- data.frame(batches = "b45", gene = "PCGF5", group = "ANKRD1")
  comps <- rbind(comps, custom_comp)
  
  # add enhancer
  m <- match(comps$group, meta_nano$Group)
  comps$Enh <- splitter(meta_nano$Enh[m], "_", 1)

  
## Run
  pdf(file = "Basic Dotplots.pdf", height = 4, width = 4)
  # nanostringResults <- apply(comps$b23, 1, function(x) runDE(batches = x[1], group = x[2]))
  nanostringResults <- apply(comps, 1, function(x) runDE(batches = x[1], gene = x[2], group = x[3]))
  dev.off()
  nanostringResults <- do.call("rbind", nanostringResults)
  
## Multiple testing
  nanostringResults$FDR <- p.adjust(nanostringResults$p, method = "fdr")
  nanostringResults$Bonferroni <- p.adjust(nanostringResults$p, method = "bonf")
  
## Save
  write.csv(nanostringResults, "LinearModels.csv")
  
################################################################################################################################ #
## Compare to scRNA-seq screen ----

## Prepare data
  p <- nanostringResults
  p$Replicated <- p$p < 0.05

  m <- match(p$Group, meta_nano$Group)
  p$ScreenFC <- meta_nano$ScreenSuppression[m] %>% exp()
  p <- p[-which(is.na(p$ScreenFC)),]

  p$NanoFC <- 2 ^ p$log2fc
  source("../../../../Manuscript/Figs/FinalFigureFunctions.R")
  p$Replicated <- factor(p$Replicated, levels = c("TRUE", "FALSE"))
  levels(p$Replicated) <- c("FDR < 0.05", "ns")
  
  r <- cor(p$ScreenFC, p$NanoFC, method = "p") %>% signif(2)
  txt <- paste0("r = ", r)
  
## Plot
  pdf(file = "scRNAseq vs Nanostring.pdf", height = 3.2, width = 3.5)
  ggplot(p, aes(x = ScreenFC, y = NanoFC, fill = Replicated, label = Group)) +
    geom_point(shape = 21, colour = "black", size = 2.5) +
    theme_bw() +
    annotate("text", label = txt, x = 1.15, y = 0.05) +
    # geom_text(size = 3) +
    scale_fill_manual(values = pal_iv_discrete2[c(2,1)]) +
    labs(x = "scRNA-seq fold-change", y = "Nanostring replication fold-change") +
    theme(panel.grid = invis, legend.position = c(0.25, 0.85), legend.background = element_rect(colour = "black")) +
    theme(panel.border = invis, axis.line = element_line()) +
    guides(fill = guide_legend(title = "Nanostring")) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.2, linetype = 2) +
    geom_vline(xintercept = 1, alpha = 0.2, linetype = 2) +
    geom_hline(yintercept = 1, alpha = 0.2, linetype = 2) +
    scale_y_continuous(limits = c(0,1.3), expand = c(0,0), breaks = c(0.5, 1)) +
    scale_x_continuous(limits = c(0,1.3), expand = c(0,0), breaks = c(0, 0.5, 1))
  dev.off()
  
################################################################################################################################ #
## Save ----
  
exp_nano <- lapply(nano, function(x) {
  y <- as.data.frame(x$exprs) # expression matrix
  # colnames(y) <- x$samples$Group # rename columns to a better sample id
  rownames(y) <- x$dict$Name # rename rows to gene name
  y <- y[unique(comps$gene),] # filter genes
  
}) 
exp_nano <- do.call("cbind", exp_nano)
meta_nano <- do.call("rbind", meta)

save(nano, meta_nano, exp_nano, file = "ProcessedData.rda")
  

