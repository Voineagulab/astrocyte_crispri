setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/")


################################################################################################################################ #
## On nGenes per enhancer, as a function of distance ----

res.final <- read.csv("../2_DE/Enh/Results Final.csv")

## Number of EGPs at each threshold
  table(res.final$Gene.Distance < 200000) # 3246 true of 7759
  res.final$Gene[which(res.final$Gene.Distance < 200000)] %>% unique() %>% length() # 1598
  res.final$Enh[which(res.final$Gene.Distance < 200000)] %>% unique() %>% length() # 855
  sum(res.final$HitPermissive[which(res.final$Gene.Distance < 200000)]) # 120
  
  table(res.final$Gene.Distance < 50000) # 739 true of 7759, interestingly OR of 10 for being a hit!
  res.final$Gene[which(res.final$Gene.Distance < 50000)] %>% unique() %>% length() # 461
  res.final$Enh[which(res.final$Gene.Distance < 50000)] %>% unique() %>% length() # 465
  sum(res.final$HitPermissive[which(res.final$Gene.Distance < 50000)]) # 79

  
## Distribution for each enhancer
   dist <- 50
  pdf(file = paste0("nTested Genes ", dist, "kb.pdf"), height = 3.5, width = 5)
  x <- res.final
  tab <- table(x$Enh, x$Gene.Distance < (dist*1000)) %>% as.data.frame()
  tab$isHit <- tab$Var1 %in% x$Enh[which(x$HitPermissive)]
  tab <- tab[which(tab$Var2 == "TRUE"),]
  # y <- table(tab$Freq, tab$isHit) %>% as.data.frame()
  # tab <- table(tab$Freq) %>% as.data.frame()
  # tab <- tab[,c("isHit", "Freq")]
  tab$isHit <- factor(tab$isHit)
  tab$Scaled <- scale(tab$Freq)
  levels(tab$isHit) <- c("ns", "Hit Enhancer")
  
  t.test(tab$Freq ~ tab$isHit)
  
  # boxplot(tab$Freq ~ tab$isHit)
  pA <- ggplot(tab, aes(x = isHit, y = Freq)) +
    geom_violin(draw_quantiles = 0.5) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line()) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = paste0("Number of Genes < ", dist, "kb per Enhancer"), x = "Enhancer is Hit")
  
  pB <- ggplot(tab, aes(x = isHit, y = Scaled)) +
    geom_violin(draw_quantiles = 0.5) +
    theme_bw() +
    theme(panel.border = invis, axis.line.y = element_line()) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = paste0("Number of Genes < ", dist, "kb per Enhancer"), x = "Enhancer is Hit")
  
  plot_grid(pA, pB)
  
  ggplot(tab, aes(x = Freq)) +
    geom_bar() +
    facet_wrap(~isHit, scales = "free_y", ncol = 1) +
    theme_bw() +
    # theme(panel.border = invis, axis.line.y = element_line()) +
    scale_y_continuous(expand = c(0,0)) +
    # labs(y = paste0("Number of Genes < ", dist, "kb per Enhancer"), x = "Enhancer is Hit")
    labs(x = paste0("Number of Genes < ", dist, "kb per Enhancer"), y = "Count of Enhancers")
  
  dev.off()
  
  
################################################################################################################################ #
## On sgRNA binding efficiency ----
  
  
## Load
  guides <- read.csv("../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  gLong <- read.csv("../../../FullLibrary_Selection/Results/FromSefi/Final Guides, Long Regions.csv")
  gShort <- read.csv("../../../FullLibrary_Selection/Results/FromSefi/Final Guides, Short Regions.csv")
  g <- rbind(gShort, gLong)
  
  
## Process
  # fix names
  g$Input <- splitter(g$Input, "-exp-", 1)
  
  # add enhancer and guide id
  m <- match(g$sgRNA.Sequence, guides$GuideSequence)
  g$Enh <- guides$TargetID[m]
  g$GuideID <- guides$GuideID[m]
  
  # add n
  m <- match(g$GuideID, colnames(nha@meta.data))
  g$nCells <- colSums(nha@meta.data[,m])
  
  # add hit
  g$Hit <- g$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
  
  # filter columns
  g <- g[,c("GuideID", "Enh", "Hit", "On.Target.Efficacy.Score", "nCells")]
  

  
## Plot overall
  
  
## Summarise per enhancer
  x <- data.frame(Enh = unique(res.final$Enh),
                  Max = NaN,
                  WeightedMean = NaN)
  
  # get summary
  for (j in x$Enh) {
    i <- which(x$Enh == j)
    y <- g[which(g$Enh == j),]
    x$Max[i] <- max(y$On.Target.Efficacy.Score)
    x$WeightedMean[i] <- sum(y$On.Target.Efficacy.Score * (y$nCells / sum(y$nCells)))
  }
  
  x$Hit <- x$Enh %in% res.final$Enh[which(res.final$HitPermissive)]
  
  # x$Max <- scale(x$Max)
  # x$WeightedMean <- scale(x$WeightedMean)
  
  # ttest
  t.test(x$Max ~ x$Hit)
  t.test(x$WeightedMean ~ x$Hit)
  
  # plot
  y <- melt(x[,-1])
  
  pdf(file = "sgRNA On Target Score.pdf", height = 3, width = 4)
  ggplot(y, aes(x = variable, y = value, colour = Hit, fill = Hit)) +
    geom_violin(draw_quantiles = 0.5, scale = "width", width = 0.75, colour = "black") +
    # geom_quasirandom(dodge.width = 0.8, alpha = 0.2, colour = "black") +
    labs(x = "Summary method", y = "On Target Score per Enhancer (Scaled)") +
    theme_bw() +
    theme(panel.border = element_blank(), axis.line.y = element_line())
  dev.off()
    
################################################################################################################################ #
## Compare gene expression and gene variance within negative and enhancer libraries ----
  
## Compare gene expression and gene variance within negative and enhancer libraries
  index.enh <- which(nha$TransductionPool == "Primary")
  index.neg <- which(nha$TransductionPool == "Negative")
  index.enh <- sample(index.enh, length(index.neg))
  
  mean.enh <- rowMeans(nha@assays$VST@data[,index.enh])
  mean.neg <- rowMeans(nha@assays$VST@data[,index.neg])
  
  var.enh <- apply(nha@assays$VST@data[,index.enh], 1, var)
  var.neg <- apply(nha@assays$VST@data[,index.neg], 1, var)
  
  
  used.genes <- which(rownames(nha) %in% res.final$Gene)
  
  
  pdf(file = "Gene Expression Comparison Neg vs Enh Pool.pdf", height = 3, width = 9)
  x <- data.frame(Enh = mean.enh,
                  Neg = mean.neg)
  x$GeneType <- "Not Tested"
  x$GeneType[used.genes] <- "Tested ns"
  x$GeneType[(rownames(x) %in% res.final$Gene[which(res.final$HitPermissive)])] <- "Hit"
  
  
  
  # qplot(var.enh[used.genes], var.neg[used.genes]) +
  ggplot(x, aes(x = Enh, y = Neg)) +
    geom_point() +
    facet_wrap(~GeneType) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
    theme() +
    labs(x = "Enhancer Pool Mean Expression (VST)", y = "Negative Control Pool Mean Expression (VST)")
  
  x <- data.frame(Enh = var.enh,
                  Neg = var.neg)
  x$GeneType <- "Not Tested"
  x$GeneType[used.genes] <- "Tested ns"
  x$GeneType[(rownames(x) %in% res.final$Gene[which(res.final$HitPermissive)])] <- "Hit"
  
  
  
  # qplot(var.enh[used.genes], var.neg[used.genes]) +
  ggplot(x, aes(x = Enh, y = Neg)) +
    geom_point() +
    facet_wrap(~GeneType) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
    theme() +
    labs(x = "Enhancer Pool Variance (VST)", y = "Negative Control Pool Variance (VST)")
  dev.off()  
  
  
  