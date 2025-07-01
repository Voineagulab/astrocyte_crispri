## This script examines the power to detect DE

################################################################################################################################ #
## Setup ----


## Generic
  rm(list = ls()); gc()
  setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Power/")
  options(stringsAsFactors = FALSE)

## Packages
  library(sceptre)
  
## Loading
  load("../Sceptre Input Files.rda", verbose = TRUE)
  res.final <- read.csv("../2_DE/Enh/Results Final.csv")
  source("../../../Scripts/Functions.")
  invis <- element_blank()
  source("../../../../Manuscript/Figs/FinalFigureFunctions.R")
  

################################################################################################################################ #
## Run simulation ----
  
  
## Functions
  # bin the cell count
    toBins <- function(x) {
      
      # define bin boundaries: 100-500 with steps of 50, then up to 1000 by steps of 100
      bins <- c(seq(100, 500, 50),
                seq(600, 1000, 100))
      
      # assign
      assigned <- list()
      for (j in 1:length(bins)) {
        diffs <- abs(x[j]-bins)
        b <- bins[which.min(diffs)]
        assigned[[j]] <- b
      }
      
      assigned <- do.call("c", assigned)  
      assigned <- unique(assigned)
      assigned <- sort(assigned)
      
      return(assigned)
      
    }
  
  # to downsample counts in a vector
    ds_fun <- function(vec, end_sum) { # from https://rdrr.io/github/mirvie/mirmisc/src/R/downsample.R
      checkmate::assert_integerish(vec, lower = 0, min.len = 1, any.missing = FALSE)
      checkmate::assert_count(end_sum)
      vec_sum <- sum(vec)
      n_to_draw <- vec_sum - end_sum
      draws <- 0
      if (n_to_draw > 0) {
        draws <- detrendr::rfromboxes(n_to_draw, vec, vec)
      }
      vec - draws
    }
  
## Parameters
  nRep <- 10 # number of times to repeat the simulation per gene/foldchange/ncell parameter combination
  fcs <- c(0.15, 0.25) # the fold-changes to use. Note that these get converted to downregulations
  simDE <- list() # to store data
  
  # and the covariate matrix
  meta_covar <- sceptre.covar
  # rownames(meta_covar) <- colnames(mx)
 
  
## Loop
  # toRun <- unique(res.final$Gene)
  toRun <- sapply(simDE, function(x) {
    y <- splitter(x$gene_id, "_", 2) %>%
      gsub("Rep", "", .) %>%
      as.numeric() %>% 
      max()
    return(y != 20)
  })
  
  toRun <- names(toRun)[which(toRun)]
  
  # for (gene in unique(res.final$Gene)) {
  for (gene in toRun) {

    # check the gene; if already run, skip 
    alreadyRun <- gene %in% names(simDE)
    # if(alreadyRun) next
    # actually, above was replaced by adding simulations
    print(gene)
      
    
  ## Get the number of cells to use in the test group
    # for this gene, what are the number of cells for the enhancers targeting it?
    e <- res.final$Enh[which(res.final$Gene == gene)] # the enhancers for the gene
    
    if (length(e) > 1) { # if there's more than one...
      cellCounts <- rowSums(sceptre.guide.pooled[e,])  
    } else {
      cellCounts <- sum(sceptre.guide.pooled[e,])  
    }
    
    cellCounts_binned <- toBins(cellCounts) # function defined earlier in script
    
  ## Dataframes and list to store and read data
    mx <- list() # simulated expression
    meta_guide <- list() # guide metadata
    
    gene_exp <- sceptre.exp[gene,] # actual expression for the gene of interest

  ## Your simulation loop
    # check if there is already a simulation for this gene. add to it if so
    if (alreadyRun) {
      currentN <- splitter(simDE[[gene]]$gene_id, "_", 2) %>%
        gsub("Rep", "", .) %>%
        as.numeric() %>%
        max()
      reps <- (currentN + 1) : (currentN + nRep)
    } else {
      reps <- 1:nRep
    }
    
    for (j in reps) { # for each of the replications
      
      for (k in cellCounts_binned) { # for each of the numbers of cells to perturb
        
        for (l in fcs) { # for each fold change
          
          # label containing all the above information
          lab <- paste(gene, 
                       paste0("Rep", j),
                       paste0("n", k), 
                       paste0("fc", l), 
                       sep = "_")
          
          # sample k cells to be perturbed
          s <- sample(1:length(gene_exp), k, replace = FALSE) 
          
          # annotate as being perturbed
          g <- rep(0, length(gene_exp))
          g[s] <- 1
          
          # downsample expression of the k cells
          x <- gene_exp
          ds_count <- sum(x[s]) * (1-l)
          ds_count <- round(ds_count)
          x[s] <- ds_fun(x[s], ds_count)
          
          # save
          mx[[lab]] <- x
          meta_guide[[lab]] <- g
          
        }
        
     
        
      }
      

  }

  ## Reformat simulation data to suit SCEPTRE
    mx <- do.call("rbind", mx) # in mx, each row is a different simulation, and is treated as though it has a different gene id
    mx <- as(mx, "sparseMatrix")
    
    meta_guide <- do.call("rbind", meta_guide)
    colnames(meta_guide) <- colnames(mx)
    meta_guide <- as(meta_guide, "sparseMatrix")
  
    # finally, create a "pairs" matrix, which indicates which guides should be tested against each gene
    meta_pairs <- data.frame(gene_id = rownames(mx),
                             gRNA_id = rownames(meta_guide),
                             pair_type = "positive_control")
    
 
  
  ## Run DE
    a <- Sys.time()
    run <- run_sceptre_high_moi(gene_matrix = mx,
                                 combined_perturbation_matrix = meta_guide,
                                 covariate_matrix = meta_covar,
                                 gene_gRNA_group_pairs = meta_pairs,
                                 side = "both",
                                 B = 1000)
    b <- Sys.time()
    gc()
    
    print(b-a)
  
  ## Save
    if (alreadyRun) {
      simDE[[gene]] <- rbind(simDE[[gene]], run)
    } else {
      simDE[[gene]] <- run
    }
    
    # simDE[[gene]] <- run
  
    # to disk
    if (which(unique(res.final$Gene) == gene) %% 100 == 0) save(simDE, file = paste0("Temp_simDE_", gene, ".rda"))   

  
  }
    
## Save
  save(simDE, file = "Power Simulation - All Simulations, FC 0.15 and 0.25, n=20.rda")
  
  
################################################################################################################################ #
## Process simulation results ----
  
## Make dataframe from list 
  simRes <- do.call("rbind", simDE)
  
## Process columns to add pertinent metadata
  rownames(simRes) <- simRes$gene_id

  simRes$SceptreP <- simRes$p_value
  simRes$Z <- simRes$z_value
  
  # simulation information
  simRes$Gene <- splitter(simRes$gene_id, "_", 1)
  simRes$Rep <- splitter(simRes$gene_id, "_", 2)
  simRes$nCells <- splitter(simRes$gene_id, "_", 3) %>% gsub("n", "", .) %>% as.numeric()
  simRes$nCells <- splitter(simRes$gene_id, "_", 3) %>% gsub("n", "", .) %>% as.numeric()
  simRes$FC <- splitter(simRes$gene_id, "_", 4) %>% gsub("fc", "", .) %>% as.numeric()
  
  simRes <- simRes[,c("Gene", "Rep", "nCells", "FC", "SceptreP", "Z")]
    
 
## Calculate empirical p-value
  ## Load negative control guide-gene pairs
    load("../Neg/SCEPTRE Output (NegE, Guide-level).rda")
    neg <- de.negE
    
  ## Calculate negative binomial p
    neg$NB.p <- 2 * pnorm(abs(neg$z_value), lower.tail = FALSE) # in negative controls
    simRes$NegBinomP <- 2 * pnorm(abs(simRes$Z), lower.tail = FALSE) # in simulation

  ## Use this to calculate an empirical P for negative binomial tests on enh-gene pairs
    total.ntc <- nrow(neg) + 1
      
    # for enh-level analysis
    simRes$EmpiricalP <- NA

    for (j in 1:nrow(simRes)) {
      
      if (j %% 100 == 0) print(j)
      
      lower <- length(which(neg$NB.p < simRes$NegBinomP[j])) # number of neg-gene pairs more significant 
      e <- (lower + 1) / total.ntc # p-value
      simRes$EmpiricalP[j] <- e
      
        
    }
  
  ## Base your p cut-off on the full-scale screen (rather than FDR correct, as they will have different strengths)
    thresh <- res.final$P.N50[which(res.final$HitPermissive)] %>% max() # the maximum empirical p to be called a hit
    simRes$EmpiricalSig <- simRes$EmpiricalP < thresh

  ## Save
    write.csv(simRes, file = "Processed 1 - Simulation results.csv")
    
## Convert to power
  # in other words, take all the reps for a given FC/nCells/gene, and get the fraction significant
  x <- simRes
  x$Sim <- paste(x$Gene, x$FC, x$nCells, sep = "_")
  x <- x[order(x$Gene, x$nCells, x$FC, x$Rep),]

  power <- list()
  for (j in unique(x$Sim)) {
    print(j)
    y <- x[which(x$Sim == j),]
    z <- y[1,1:4]
    z$Power <- sum(y$EmpiricalSig) / nrow(y)
    power[[j]] <- z
  }

  power <- do.call("rbind", power)  
  power$FC <- factor(power$FC)
    
## Add mean expression
  m <- match(power$Gene, names(meanExp_NHA))
  power$MeanExp <- meanExp_NHA[m]
  
  power$ExpressionBin <- cut(power$MeanExp, c(0,0.2,0.5,1,10))
  levels(power$ExpressionBin) <- c("Low (0-0.2)",
                                  "Moderate (0.2-0.5)",
                                  "High (0.5-1)",
                                  "Extreme (1+)")
  levels(power$ExpressionBin) <- paste0(levels(power$ExpressionBin), " (n=", table(power$ExpressionBin), ")")
  
  power$ExpressionQuintile <- cut(power$MeanExp, 
                                breaks = quantile(power$MeanExp, probs = seq(0, 1, by = 0.2)),
                                labels = c("0-20th", "20-40th", "40-60th", "60-80th", "80-100th"),
                                include.lowest = TRUE)
    
  
## Save
  write.csv(power, file = "Processed 2 - Power.csv")
  
## Plot
  ## Effect of gene expression
    # ggplot(power, aes(x = MeanExp, y = Power, colour = (FC))) +
    #   geom_point() +
    #   geom_smooth() +
    #   theme_bw() +
    #   theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) +
    #   geom_hline(yintercept = 0.8, linetype = 2) +
    #   scale_y_continuous(limits = c(0, 1)) +
    #   labs(x = "Gene normalised expression", y = "Sensitivity")
    # 
    # ggplot(power, aes(y = MeanExp, x = as.factor(Power), colour = (FC))) +
    #   geom_violin(scale = "width") +
    #   
    #   theme_bw() +
    #   theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) +
    #   geom_hline(yintercept = 0.8, linetype = 2) +
    #   # scale_y_continuous(limits = c(0, 1)) +
    #   labs(x = "Power", y = "Normalised expression")
    
  pdf(file = "Power Simulation - Effect of expression quintile.pdf", height = 4, width = 8)
  ggplot(power, aes(x = ExpressionQuintile, y = Power, colour = ExpressionQuintile)) +
      # geom_boxplot() +
    geom_quasirandom() +
    stat_summary(fun = mean, geom = "point", shape = "+", colour = "black", size = 10) +
    # geom_violin() +
      facet_wrap(~FC) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
            axis.text.x = invis) +
      geom_hline(yintercept = 0.8, linetype = 2) +
      # scale_y_continuous(limits = c(0, 1)) +
      labs(x = "", y = "Sensitivity")  
  dev.off()
  
    
  ## Effect of fold-change
    # process data
    p <- power[,c(1,3,4,5)]
    p$Sim <- paste0(p$Gene, p$nCells)
    
    p <- pivot_wider(p[,3:5], names_from = Sim, values_from = Power)
    p <- as.data.frame(p)
    rownames(p) <- p[,1]
    p <- p[,-1]
    p <- t(p)
    p <- as.data.frame(p)
    colnames(p) <- c("LowFC", "HighFC")
    
    p$LowFC <- factor(p$LowFC, levels = as.character(seq(0, 1, 0.05)))
    p$HighFC <- factor(p$HighFC, levels = as.character(seq(0, 1, 0.05)))
    p <- table(p$LowFC, p$HighFC)
    p <- as.data.frame(p)
    
    colnames(p) <- c("LowFC", "HighFC", "Freq")
    p$Label <- p$Freq
    p$Label[which(p$Label == 0)] <- "."
    
    # plos
    pdf(file = "Power Simulation - Foldchange comparison.pdf", height = 6, width = 6)
    ggplot(p, aes(x = LowFC, y = HighFC, label = Label, fill = Freq)) +
      geom_tile() +
      geom_text(size = 3) +
      scale_fill_carto_c(palette = "Geyser", trans = "log2") +
      geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
      labs(x = "Sensitivity for fc of -15%", y = "Sensitivity for fc of -25%") +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
            legend.position = "none") 
    dev.off()
    
  ## Effect of nCells
    pdf(file = "Power Simulation - Effect of nCells.pdf", height = 6, width = 8)
    ggplot(power, aes(x = as.factor(nCells), y = Power, colour = ExpressionQuintile)) +
      geom_boxplot() +
      facet_wrap(~FC, ncol = 1) +
      theme_bw() +
      labs(x = "Number of cells", y = "Sensitivity") +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) 
    dev.off()
    
  ## Combined
    pdf(file = "Power Simulation - All variables.pdf", height = 12, width = 8)
    ggplot(power, aes(x = ExpressionQuintile, y = Power, colour = ExpressionQuintile)) +
      # geom_boxplot() +
      geom_quasirandom() +
      stat_summary(fun = mean, geom = "point", shape = "+", colour = "black", size = 10) +
      # geom_violin() +
      facet_grid(nCells~FC) +
      theme_bw() +
      theme(axis.line.y = element_line(), panel.grid = invis,
            axis.text.x = invis) +
      geom_hline(yintercept = 0.8, linetype = 2) +
      # scale_y_continuous(limits = c(0, 1)) +
      labs(x = "", y = "Sensitivity", title = "Columns: Fold-changes\nRows: nCells")  
    
    ggplot(power, aes(x = MeanExp, y = as.numeric(Power))) +
      # geom_boxplot() +
      geom_point() +
      
      # stat_summary(fun = mean, geom = "point", shape = "+", colour = "black", size = 10) +
      # geom_violin() +
      facet_grid(nCells~FC) +
      theme_bw() +
      theme(axis.line.y = element_line(), panel.grid = invis) +
      geom_hline(yintercept = 0.8, linetype = 2) +
      # scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Mean expression", y = "Sensitivity", title = "Columns: Fold-changes\nRows: nCells")  
    dev.off()
    
  ## The above, but nicer
    p <- power
    levels(p$ExpressionBin) <- gsub(") ", ")\n", levels(p$ExpressionBin))
    p$nCells <- factor(p$nCells)
    levels(p$nCells) <- paste0(levels(p$nCells), "\ncells")
    
    p$FC <- factor(p$FC)
    levels(p$FC) <- paste0("Fold-change -", levels(p$FC))
    
    
    levels(p$ExpressionQuintile) <- gsub("th", "", levels(p$ExpressionQuintile))
    #Supplementary Figure 4A
    pdf(file = "Power Simulation - All variables V2.pdf", height = 7.5, width = 7)
    ggplot(p, aes(x = ExpressionQuintile, y = Power*100, colour = ExpressionQuintile)) +
      # geom_boxplot() +
      geom_hline(yintercept = 80, linetype = 2, alpha = 0.5) +
      geom_quasirandom(width = 0.3) +
      stat_summary(fun = median , geom = "point", shape = "+", colour = "black", size = 5) +
      # geom_violin() +
      facet_grid(nCells~FC) +
      theme_bw() +
      scale_colour_manual(values = pals$Primary[c(2,5,1,6,7)]) +
      # scale_colour_manual(values = pals$grn2orng[c(3,3,2,2,1)]) +
      theme(axis.line.y = element_line(), panel.grid = invis,
            legend.position = "none") +
      
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 100)) +
      labs(x = "Expression quintile", y = "Power")  
    dev.off()
    
  ## Summarise as the proportion that are moderately and well-powered
    p <- power
    p$Group <- paste0(p$FC, "_", p$nCells, "_", p$ExpressionQuintile)

    p$Powered50 <- p$Power >= 0.5
    p$Powered80 <- p$Power >= 0.8
    
    prop50 <- table(p$Group, p$Powered50) %>% proportions(1)
    prop80 <- table(p$Group, p$Powered80) %>% proportions(1)
    
    prop <- data.frame(FC = splitter(rownames(prop50), "_", 1),
                       nCells = splitter(rownames(prop50), "_", 2),
                       Quintile = splitter(rownames(prop50), "_", 3),
                       Powered50 = prop50[,"TRUE"],
                       Powered80 = prop80[,"TRUE"])
    
    prop$nCells <- factor(prop$nCells, levels = c(sort(unique(as.numeric(prop$nCells)))))
    
    prop$FC <- factor(prop$FC)
    levels(prop$FC) <- paste0("Fold-change -", levels(prop$FC))
    
    #Supplementary Figure 4B
    pdf(file = "Proportion well-powered (wide).pdf", width = 7, height = 3)
    
    # ggplot(prop, aes(x = nCells, y = Quintile, fill = Powered80*100, label = round(Powered80*100, 0))) +
    #   geom_tile() +
    #   geom_text(size = 3) +
    #   scale_y_discrete(position = "left", expand = c(0,0)) +
    #   theme_bw() +
    #   facet_wrap(~FC, ncol = 1, strip.position = ) +
    #   theme(panel.grid = invis, legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    #         axis.ticks.x = invis) +
    #   scale_x_discrete(expand = c(0,0)) +
    #   scale_fill_gradientn(colours = c("#b51a00","#ee5900","#ff9d68","#feceb8","grey95"),
    #                      values = c(1.0,0.7,0.5,0.4,0.2,0)) +
    #   guides(fill = guide_colourbar(title = "Proportion\nwell-powered")) +
    #   labs(x = "Number of cells", y = "Expression quintile")
    
    ggplot(prop, aes(y = nCells, x = Quintile, fill = Powered80*100, label = round(Powered80*100, 0))) +
      geom_tile() +
      geom_text(size = 3) +
      scale_y_discrete(position = "left", expand = c(0,0)) +
      theme_bw() +
      facet_wrap(~FC, ncol = 2, strip.position = ) +
      theme(panel.grid = invis, legend.position = "right", axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
      scale_x_discrete(expand = c(0,0)) +
      scale_fill_gradientn(colours = c("#b51a00","#ee5900","#ff9d68","#feceb8","grey95"),
                         values = c(1.0,0.7,0.5,0.4,0.2,0)) +
      guides(fill = guide_colourbar(title = "Proportion\nwell-powered")) +
      labs(y = "Number of cells", x = "Expression quintile")
    dev.off()
    
    
    pdf(file = "Proportion well-powered (wide v2).pdf", width = 7, height = 2.5)
    ggplot(prop, aes(x = nCells, y = Quintile, fill = Powered80*100, label = round(Powered80*100, 0))) +
      geom_tile() +
      # geom_text(size = 3) +
      scale_y_discrete(position = "left", expand = c(0,0)) +
      theme_bw() +
      facet_wrap(~FC, ncol = 2, strip.position = ) +
      theme(panel.grid = invis, legend.position = "right", 
            # axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_x_discrete(expand = c(0,0)) +
      scale_fill_gradientn(colours = c("#b51a00","#ee5900","#ff9d68","#feceb8","grey90"),
                         values = c(1.0,0.7,0.5,0.4,0.2,0)) +
      guides(fill = guide_colourbar(title = "Percent\nwell-powered")) +
      labs(x = "Number of cells", y = "Expression quintile")
   
    dev.off() 
    
    # pdf(file = "Proportion well-powered.pdf", height = 8, width = 3)
    # ggplot(prop, aes(y = nCells, x = Quintile, fill = Powered80*100, label = round(Powered80*100, 0))) +
    #   geom_tile() +
    #   geom_text(size = 3) +
    #   scale_y_discrete(position = "right", expand = c(0,0)) +
    #   theme_bw() +
    #   facet_wrap(~FC, ncol = 1, strip.position = ) +
    #   theme(panel.grid = invis, legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    #         axis.ticks.x = invis) +
    #   scale_x_discrete(expand = c(0,0)) +
    #   scale_fill_gradientn(colours = c("#b51a00","#ee5900","#ff9d68","#feceb8","grey95"),
    #                      values = c(1.0,0.7,0.5,0.4,0.2,0)) +
    #   guides(fill = guide_colourbar(title = "Proportion\nwell-powered")) +
    #   labs(y = "Number of cells", x = "Expression quintile")
    # 
    # dev.off() 
    
    
    
    
    
################################################################################################################################ #
## Integrate power analyses with other project analyses ----
    
## How does power relate to the actual screen DE results?
  powerInScreen <- list()
    
  for (j in 1:nrow(res.final)) {
    print(j)
    
    # simplify your screen results mx
    x <- res.final[j,c("Pair", "Enh", "Gene","Gene.Exp","nCells","logfc.vst","HitPermissive")]
    
    # if (!(x$Gene %in% power$Gene)) next # only useful if the simulation is incompletely run
    
    # for the egp, get the simulation that most closely matches it in terms of gene and nCells
    cellBin <- toBins(x$nCells) # round the cell count to be similar to those used in simulation
    y <- power[which(power$Gene == x$Gene & power$nCells == cellBin),]
    
    # extract
    x$Sensitivity_0.15 <- y$Power[which(y$FC == "0.15")]
    x$Sensitivity_0.25 <- y$Power[which(y$FC == "0.25")]
    
    powerInScreen[[j]] <- x
    
  }
      
  powerInScreen <- do.call("rbind", powerInScreen)
  powThresh <- 0.8
  powerInScreen$WellPowered015 <- powerInScreen$Sensitivity_0.15 >= powThresh
  powerInScreen$WellPowered025 <- powerInScreen$Sensitivity_0.25 >= powThresh
  
  table(powerInScreen$WellPowered015, powerInScreen$HitPermissive)
  table(powerInScreen$WellPowered015, powerInScreen$HitPermissive) %>% fisher.test()
  table(powerInScreen$WellPowered015, powerInScreen$HitPermissive) %>% proportions(margin = 2)
  
  table(powerInScreen$WellPowered025, powerInScreen$HitPermissive)
  table(powerInScreen$WellPowered025, powerInScreen$HitPermissive) %>% fisher.test()
  table(powerInScreen$WellPowered025, powerInScreen$HitPermissive) %>% proportions(margin = 2)
  
  p <- powerInScreen
  p$WellPowered015 <- factor(p$WellPowered015)
  p$WellPowered025 <- factor(p$WellPowered025)
  p$HitPermissive <- factor(p$HitPermissive)
  levels(p$WellPowered015) <- levels(p$WellPowered025) <- c("<80%", ">80%")
  levels(p$HitPermissive) <- c("Inactive EGP", "Functional EGP")
  
  p <- melt(p[,7:9])
  
  colnames(p)[1] <- "Hit"
  levels(p$variable) <- c("-15%", "-25%")
  pdf(file = "Comparison to screen - Wellpowered hits versus nonhits.pdf", height = 3, width = 3)
  ggplot(p, aes(x = variable, y = value, fill = Hit, colour = Hit)) +
    # geom_col(position = "dodge") +
    # geom_quasirandom() +
    geom_violin(scale = "width", width = 0.7, draw_quantiles = 0.5, alpha = 0.9) +
    # geom_quasirandom(dodge.width = 0.7, width = 0.15) +
    theme_bw() +
    labs(x = "Simulated fold-change", y = "Power") +
    guides(fill = guide_legend(nrow = 2), colour = "none") +
    theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis, legend.position = "bottom") +
    scale_fill_manual(values = pals$Hits) +
    scale_colour_manual(values = pals$Hits_Darker) 
  dev.off()
  
  # prop50 <- table(p$HitPermissive, p$WellPowered015) %>% proportions(1)
  #   prop80 <- table(p$HitPermissive, p$WellPowered015) %>% proportions(1)
  #   
  #   prop <- data.frame(FC = splitter(rownames(prop50), "_", 1),
  #                      nCells = splitter(rownames(prop50), "_", 2),
  #                      Quintile = splitter(rownames(prop50), "_", 3),
  #                      Powered50 = prop50[,"TRUE"],
  #                      Powered80 = prop80[,"TRUE"])
  # 
  # as.data.frame(table(p$HitPermissive, p$WellPowered015, p$WellPowered025))
  
  ## Save
    write.csv(powerInScreen, file = "Power Simulation - Power per EGP.csv")
  
  ## For non-powered hits, what was the fold-change?
    x <- powerInScreen
    x <- x[which(x$HitPermissive),]
    x$WellPowered015 <- factor(x$WellPowered015)
    x$WellPowered025 <- factor(x$WellPowered025)
    levels(x$WellPowered025) <- levels(x$WellPowered015) <- c("< 80%", ">= 80%")
    
    pdf(file = "Comparison to screen - Wellpowered hits versus foldchange.pdf", height = 3, width = 3)
    ggplot(x, aes(x = WellPowered015, y = log2(exp(logfc.vst)))) +
      geom_hline(yintercept = c(0,log(1-0.15)), linetype = c(1,2), alpha = 0.5) +
      geom_violin() +
      geom_quasirandom(size = 1, width = 0.3) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis ) +
      labs(x = "Power to detect -15% downregulation", y = "Log2 fold-change\nat functional EGPs")
    # dev.off()
    
    ggplot(x, aes(x = WellPowered025, y = log2(exp(logfc.vst)))) +
      geom_hline(yintercept = c(0,log(1-0.15)), linetype = c(1,2), alpha = 0.5) +
      geom_violin() +
      geom_quasirandom(size = 1, width = 0.3) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis ) +
      labs(x = "Power to detect -25% downregulation", y = "Log2 fold-change\nat functional EGPs")
    dev.off()
      
## Condense to the enhancer level
  powerAtEnh <- list(Fifteen = aggregate(powerInScreen$Sensitivity_0.15~powerInScreen$Enh, FUN = max),
                       Twentyfive = aggregate(powerInScreen$Sensitivity_0.25~powerInScreen$Enh, FUN = max))
  powerAtEnh <- do.call("cbind", powerAtEnh)
  powerAtEnh <- powerAtEnh[,c(1,2,4)]
  colnames(powerAtEnh) <- c("Enh", "Sensitivity15", "Sensitivity25")
  powerAtEnh$Hit <- powerAtEnh$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
  powerAtEnh$Hit <- factor(powerAtEnh$Hit)
  levels(powerAtEnh$Hit) <- c("Inactive candidate", "Functional enhancer")
  
  # save
  write.csv(powerAtEnh, file = "Power Simulation - Power per Enh.csv")
  
  # categorise as well power
  wellPoweredThresh <- 0.8
  powerAtEnh$Powered15 <- powerAtEnh$Sensitivity15 >= wellPoweredThresh   
  powerAtEnh$Powered25 <- powerAtEnh$Sensitivity25 >= wellPoweredThresh   
  
  # table
  table(powerAtEnh$Hit, powerAtEnh$Powered15)
  table(powerAtEnh$Hit, powerAtEnh$Powered15) %>% proportions(1)
  table(powerAtEnh$Hit, powerAtEnh$Powered15) %>% fisher.test()
  
  table(powerAtEnh$Hit, powerAtEnh$Powered25)
  table(powerAtEnh$Hit, powerAtEnh$Powered25) %>% proportions(1)
  table(powerAtEnh$Hit, powerAtEnh$Powered25) %>% fisher.test()
   
    
# ## Compare power to the technical score from an enhancer prediction model
#   # read in
#   tech <- read.csv("../../../../EnhancerPredictionModels/Data/Enhancer_factors_pluspredictions.csv")
#   tech$Hit <- factor(tech$HitPermissive > 0)
#   levels(tech$Hit) <- c("Inactive candidate", "Functional enhancer")
#   
#   # for each enhancer, get the maximum power
#   enhMax15  <- aggregate(powerInScreen$Sensitivity_0.15~powerInScreen$Enh, FUN = max)
#   m <- match(tech$Enh, enhMax15$`powerInScreen$Enh`)
#   tech$PowerMax15 <- enhMax15$`powerInScreen$Sensitivity_0.15`[m]
#   
#   enhMax25  <- aggregate(powerInScreen$Sensitivity_0.25~powerInScreen$Enh, FUN = max)
#   m <- match(tech$Enh, enhMax25$`powerInScreen$Enh`)
#   tech$PowerMax25 <- enhMax25$`powerInScreen$Sensitivity_0.25`[m]
#   
#   # scatterplot
#   pdf(file = "Power Simulation - Versus Technical Score.pdf", height = 5, width = 8)
#   ggplot(tech, aes(y = all_genes_tech_rf_ref, x = as.factor(PowerMax15), colour = Hit )) +
#     geom_hline(yintercept = c(0,1,-1,-2,-3), linetype = 2, alpha = 0.5) +
#     # geom_violin(scale = "width", draw_quantiles = 0.5) +
#     geom_boxplot(outlier.shape = NA, width = 0.8, position = position_dodge(width = 0.8)) +
#     
#     geom_quasirandom(dodge.width = 0.8, size = 1, alpha = 0.5, width = 0.1) +
#     # stat_summary(fun = mean, geom = "point", shape = "+", size = 10, position = position_dodge(width = 0.7)) +
#     theme_bw() +
#     theme(panel.border = invis, axis.line = element_line(), panel.grid = invis) +
#     labs(x = "Power", y = "Technical score", title = "FC = -0.15")
#   
#   ggplot(tech, aes(y = all_genes_tech_rf_ref, x = as.factor(PowerMax25), colour = Hit )) +
#     geom_hline(yintercept = c(0,1,-1,-2,-3), linetype = 2, alpha = 0.5) +
#     # geom_violin(scale = "width", draw_quantiles = 0.5) +
#     geom_boxplot(outlier.shape = NA, width = 0.8, position = position_dodge(width = 0.8)) +
#     
#     geom_quasirandom(dodge.width = 0.8, size = 1, alpha = 0.5, width = 0.1) +
#     # stat_summary(fun = mean, geom = "point", shape = "+", size = 10, position = position_dodge(width = 0.7)) +
#     theme_bw() +
#     theme(panel.border = invis, axis.line = element_line(), panel.grid = invis) +
#     labs(x = "Power", y = "Technical score", title = "FC = -0.25")
#   dev.off()
#   
#   
# ## Compare the effect of power to nearest gene 
#   # specifically, their plots!
#   distCat <- read.csv(file = "../../3_HitEnrichment/EnhGenePairs/Intervening Gene Classification Between EGPs.csv", row.names = 1)
#   
#   m <- match(distCat$Pair, powerInScreen$Pair)
#   distCat$Powered <- powerInScreen$WellPowered015[m]
#   
#   
#   
#    ## Original version
#     p <- table(distCat$HitPermissive, distCat$Distance.Category)
#     p <- p / rowSums(p)
#     p <- t(p)
#     p <- as.data.frame.matrix(p)
#     
#     colnames(p) <- c("ns", "Hit")
#     p$Cat <- rownames(p)
#     p <- melt(p)
#     
#     levels(p$variable) <- c("Inactive\nEGPs", "Functional\nEGPs")
#   
#   p$Cat <- factor(p$Cat)
#   levels(p$Cat) <- gsub("Cond", "Type ", levels(p$Cat)) %>%
#     sub("_", "\n", .) %>%
#     sub("_", " \n", .) %>%
#     sub("Distal", "Not nearest", .) %>%
#     sub("NonExp", " off TSSs only", .) %>%
#     sub("SkipExp", "Skip 1+ on TSS", .) %>%
#     sub("NoSkip", "No intervening TSS", .) 
#   
#   legtitle <- guide_legend(title = "EGP type")
#   pal_ord <- c(5,7,4,1)
#   
#   # pdf(file = "Nearest Gene - Stacked Barplot.pdf", height = 1.8, width = 1.8)
#   ggplot(p, aes(x = variable, fill = Cat, colour = Cat, y = value*100)) +
#     geom_col(width = 0.7) +
#     scale_fill_manual(values = pals$Primary[pal_ord]) +
#     scale_colour_manual(values = pals$Primary_Darker[pal_ord]) +
#     labs(y = "Percent of EGPs") +
#     theme_bw() +
#     scale_y_continuous(expand = c(0,0)) +
#     theme(panel.grid = invis, axis.line.y = element_line(), panel.border = invis, 
#           axis.title.x = invis, legend.position = "none") 
#   # dev.off()
# 
#   
#   distCat2 <- distCat[which(distCat$Powered),]
# 
#   ## Collect data
#     p <- table(distCat2$HitPermissive, distCat2$Distance.Category)
#     p <- p / rowSums(p)
#     p <- t(p)
#     p <- as.data.frame.matrix(p)
#     
#     colnames(p) <- c("ns", "Hit")
#     p$Cat <- rownames(p)
#     p <- melt(p)
#     
#     levels(p$variable) <- c("Inactive\nEGPs", "Functional\nEGPs")
#   
#   p$Cat <- factor(p$Cat)
#   levels(p$Cat) <- gsub("Cond", "Type ", levels(p$Cat)) %>%
#     sub("_", "\n", .) %>%
#     sub("_", " \n", .) %>%
#     sub("Distal", "Not nearest", .) %>%
#     sub("NonExp", " off TSSs only", .) %>%
#     sub("SkipExp", "Skip 1+ on TSS", .) %>%
#     sub("NoSkip", "No intervening TSS", .) 
#   
#   legtitle <- guide_legend(title = "EGP type")
#   pal_ord <- c(5,7,4,1)
#   
#   # pdf(file = "Nearest Gene - Stacked Barplot.pdf", height = 1.8, width = 1.8)
#   ggplot(p, aes(x = variable, fill = Cat, colour = Cat, y = value*100)) +
#     geom_col(width = 0.7) +
#     scale_fill_manual(values = pals$Primary[pal_ord]) +
#     scale_colour_manual(values = pals$Primary_Darker[pal_ord]) +
#     labs(y = "Percent of EGPs") +
#     theme_bw() +
#     scale_y_continuous(expand = c(0,0)) +
#     theme(panel.grid = invis, axis.line.y = element_line(), panel.border = invis, 
#           axis.title.x = invis, legend.position = "none") 
# 
#   
#   table(distCat$Distance.Category, distCat$HitPermissive) 
#   table(distCat$Distance.Category, distCat$HitPermissive) %>% proportions(margin = 2)
#   
#   table(distCat2$Distance.Category, distCat2$HitPermissive)
#   table(distCat2$Distance.Category, distCat2$HitPermissive) %>% proportions(margin = 2)
#   
#   table(distCat$Distance.Category[-which(distCat$Powered)], distCat$HitPermissive[-which(distCat$Powered)]) %>% proportions(margin = 2)
#   
#   # tab <- as.data.frame.table(table(distCat$Distance.Category, distCat$HitPermissive, distCat$Powered) )
#   
#   
# ## And now check TTseq plots
#   transcribed <- read.csv("../../4_EnhancerTranscription/TTseq/Transcriptional classification.csv")
#   
#   
#   # for each enhancer, get the maximum power
#   enhMax15  <- aggregate(powerInScreen$Sensitivity_0.15~powerInScreen$Enh, FUN = max)
#   m <- match(transcribed$Enh, enhMax15$`powerInScreen$Enh`)
#   
#   transcribed$PowerMax15 <- enhMax15$`powerInScreen$Sensitivity_0.15`[m]
#   
#   enhMax25  <- aggregate(powerInScreen$Sensitivity_0.25~powerInScreen$Enh, FUN = max)
#   m <- match(transcribed$Enh, enhMax25$`powerInScreen$Enh`)
#   transcribed$PowerMax25 <- enhMax25$`powerInScreen$Sensitivity_0.25`[m]
# 
#   table(transcribed$Category_3, transcribed$Hit)    
#   table(transcribed$Category_3, transcribed$Hit)    %>% proportions(margin = 2)
#   
#   table(transcribed$Category_3[which(transcribed$PowerMax15 >= 0.8)], transcribed$Hit[which(transcribed$PowerMax15 >= 0.8)])    
#   table(transcribed$Category_3[which(transcribed$PowerMax15 >= 0.8)], transcribed$Hit[which(transcribed$PowerMax15 >= 0.8)])    %>% proportions(margin = 2)
# 
#   
#    final.thresh <- 3
#   
#   ## Functions
#     plot.cats <- function(dat = transcribed, thresh, complex = FALSE) {
#       # stats
#       fish_trns <- table(dat$Hit, dat[,paste0("Category_", thresh)] != "Not transcribed") %>% fisher.test()
#       lab_trns <- paste0("Transcribed: p=", 
#                          signif(fish_trns$p.value, 2),
#                          ", OR=", signif(fish_trns$estimate, 2))
#       
#       fish_uni <- table(dat$Hit, dat[,paste0("Category_", thresh)] == "Unidirectional") %>% fisher.test()
#       lab_uni <- paste0("Unidirectional: p=", 
#                          signif(fish_uni$p.value, 2),
#                          ", OR=", signif(fish_uni$estimate, 2))
#       
#       fish_bi <- table(dat$Hit, dat[,paste0("Category_", thresh)] == "Bidirectional") %>% fisher.test()
#       lab_bi <- paste0("Bidirectional: p=", 
#                          signif(fish_bi$p.value, 2),
#                          ", OR=", signif(fish_bi$estimate, 2))
#     
#       # title
#       ttl <- paste0("Threshold: ", thresh,
#                     "\n", lab_trns,
#                     "\n", lab_uni,
#                     "\n", lab_bi)
#       
#       # plot values
#       tab <- table(dat$Hit, dat[,paste0("Category_", thresh)])
#       tab <- tab / rowSums(tab)
#       tab <- t(tab)
#       tab <- as.data.frame.matrix(tab)
#       tab$Transcribed <- rownames(tab)
#       
#       p <- melt(tab)
#       colnames(p)[2] <- "Hit"
#       p$Hit <- factor(p$Hit, levels = c("FALSE", "TRUE"))  
#       levels(p$Hit) <- (c("Inactive\ncandidates", "Functional\nenhancers")  )
#       p$Transcribed <- factor(p$Transcribed, levels = c("Not transcribed", "Unidirectional", "Bidirectional"))
#       # pal <- c("grey90", pals$Primary[1:2])
#       # pal <- c("grey90", pals$grn2orng[c(3,2)])
#       
#       if (complex) {
#         pal <- c("grey90", pals$Hits[1], pals$Hits_Darker[1], # non-hit colours
#                "grey91", pals$Hits[2], pals$Hits_Darker[2]) # hit colours
#         
#          ggplot(p, aes(x = Hit, fill = interaction(Transcribed, Hit), y = value*100, alpha = Transcribed)) +
#           geom_col(width = 0.7) +
#           theme_bw() +
#            scale_alpha_manual(values = c(1, 0.8, 1)) +
#           theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
#                 axis.title.x = invis, plot.title = element_text(size = 6),
#                 axis.text.y = text90) +
#           scale_fill_manual(values = pal) +
#           scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100)) +
#           labs(y = "Percent of peaks", title = ttl)
#         
#         
#       } else {
#         # pal <- c("grey90", pals$Primary[1], pals$Primary[8])
#         pal <- c("grey90", pals$Primary[2], pals$Primary_Darker[2])
#         ggplot(p, aes(x = Hit, fill = Transcribed, y = value*100, linetype = Transcribed, alpha = Transcribed)) +
#           geom_col(colour = "black", width = 0.7) +
#           theme_bw() +
#           scale_alpha_manual(values = c(1, 0.9, 1)) +
#           theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis,
#                 axis.title.x = invis, plot.title = element_text(size = 6),
#                 axis.text.y = text90) +
#           scale_fill_manual(values = pal) +
#           scale_linetype_manual(values = c("blank", "dashed", "solid")) +
#           scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100)) +
#           labs(y = "Percent of peaks", title = ttl)
#       }
#       
#       
#       
#     }
#     
#     plot.cats(dat = transcribed, thresh = 3)
#     plot.cats(dat = transcribed[which(transcribed$PowerMax15 >= 0.8),], thresh = 3)
#   
    
