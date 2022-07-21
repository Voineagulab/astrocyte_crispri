## Setup
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/1a_GuideAssignment/")
load("../2_DE/Sceptre NHA Input.rda")

## Dataframes
  dat.guide <- as.data.frame(sceptre.guide)
  dat.guide <- dat.guide[grep("Enh", rownames(dat.guide)),]
  
  # pooled at enhancer level
  dat.enh <- as.data.frame(sceptre.perturb)
  dat.enh <- dat.enh[grep("Enh", rownames(dat.enh)),]
  
## Stats
  
  
  
  
  
## Plot 1: Cells per guide and enhancer
  # add level for gasperini
  x <- data.frame(Type = "Guide", Guide = rownames(dat.guide), N = apply(dat.guide, 1, function(x) length(which(x == 1))))
  y <- data.frame(Type = "Enhancer\nTarget", Guide = rownames(dat.enh), N = apply(dat.enh, 1, function(x) length(which(x == 1))))
  x <- rbind(x, y)
  x$Bin <- cut(x$N, c (0, 30, 100, 200, 500, 750, 1000, 2000))
  levels(x$Bin) <- c("<30", "31-100", "101-200", "201-500", "501-1000", "1001-2000")
  
  pdf(file = "Cells Per Guide or Enhancer.pdf", height = 2.5, width = 3.5)
  ggplot(x, aes(x = Type, y = N)) +
    geom_violin(scale = "width", draw_quantiles = c(0.5), fill = "darkorange1", alpha = 0.7) +
    theme_bw() +
    coord_flip() +
    theme(panel.border = invis, axis.line.x = element_line(), axis.title.y = invis) +
    geom_hline(yintercept = 915, linetype = 2) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 2100), breaks = c(0, 500, 1000, 1500, 2000)) +
    labs(y = "Number of Cells")
  dev.off()  
  
  
  ## Plot 3 has been copied to a formal script
  
# ## Plot 3: UMIs per cell
#   x <- data.frame(UMI = nha$nCount_RNA, Dummy = ".")
#   pdf(file = "UMIs Per Cell.pdf", height = 2, width = 3.5)
#   m <- median(x$UMI)
#   ggplot(x, aes(x = UMI)) +
#     geom_density(fill = "darkorange1", alpha = 0.7) +
#     theme_bw() +
#     # coord_flip() +
#     theme(panel.border = invis, axis.line.x = element_line(), axis.text.y = invis, axis.ticks.y = invis) +
#     scale_x_continuous(expand = c(0,0), limits = c(0, 250000)) +
#     geom_vline(xintercept = m, linetype = 2) +
#     scale_y_continuous(expand = c(0,0)) +
#     labs(x = "UMIs Per Cell", y = "Density")
#   dev.off()    
#   
# ## Plot 4: UMIs per cell, compariing to Gasperini
#   a <- readRDS("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/GSE120861_at_scale_screen.cds.rds")
#   x <- data.frame(UMI = nha$nCount_RNA, Study = "Voineagu")
#   y <- data.frame(UMI = colSums(a@assayData$exprs), Study = "Gasperini")
#   x <- rbind(x, y)
#   pdf(file = "UMIs Per Cell (vs. Gasperini).pdf", height = 2, width = 3.5)
#   ggplot(x, aes(x = Study, y = UMI)) +
# 
#     geom_violin(scale = "width", draw_quantiles = c(0.5), fill = "darkorange1", alpha = 0.7) +
#     theme_bw() +
#     coord_flip() +
#     theme(panel.border = invis, axis.line.x = element_line(), axis.title.y = invis) +
#     geom_hline(yintercept = 915, linetype = 2) +
#     scale_y_continuous(expand = c(0,0), limits = c(0, 250000)) +
#   labs(y = "UMIs Per Cell")
#   
#   dev.off()    
  
## Test: 
  s <- data.frame(gene_id = rep(rownames(nha), each = 2),
                  gRNA_id = rep(c("Enh53", "Enh54"), times = 2),
                  pair_type = "candidate")
  
  ID3 <- list()
  sceptre.direction <- "both"

  nPerRun <- 200
  chunks <- ceiling(nrow(s) / nPerRun)
    for (j in 264:chunks) { # this chunks sceptre in runs with nPerRun test pairs
      range <- ((nPerRun*(j-1)) + 1) : (nPerRun*j) # get nPerRun test pairs
      range <- range[which(range < nrow(s))]
      
      a <- Sys.time()
      ID3[[paste0("Genes", min(range), "_", max(range))]] <- run_sceptre_high_moi(gene_matrix = sceptre.gene,
                                                                                  combined_perturbation_matrix = sceptre.perturb,
                                                                                  covariate_matrix = sceptre.covar,
                                                                                  gene_gRNA_group_pairs = s[range,],
                                                                                  side = sceptre.direction,
                                                                                  B = 500)
      b <- Sys.time()
      
      print(paste0(j, ": ", b-a))
      save(ID3, file = "Temp.rda")
      gc()
    }
    
  ## UP TO 262 of 274
   mean <- rowMeans(nha@assays$RNA@data)
   exp.thresh <- 2^-6
   mean <- names(mean)[which(mean > exp.thresh)]
  
   x <- do.call("rbind", ID3)
   x <- x[which(x$gene_id %in% mean),]
  y <- x[grep("4", x$gRNA_id),]
  x <- x[grep("3", x$gRNA_id),]
  
  y$FDR <- p.adjust(y$p_value, method = "fdr")
  x$FDR <- p.adjust(x$p_value, method = "fdr")
  
  plot(-log10(x$p_value), -log10(y$p_value))
  
  
  
## QQPLOT
  load("../2_DE/Sceptre.rda")
  x <- do.call("rbind", sceptre$neg)
  y <- sample(1:nrow(x), 100000, FALSE)   
  make_qq_plot(x$p_value[y])
  
  z <- 2 * pnorm(abs(x$z_value[y]), lower.tail = FALSE)
  make_qq_plot(z)

  
  my.pvalues=runif(10000)

#Calculate expectations
  p <- data.frame(SCEPTRE = x$p_value[y], NB = z)
  p$SCEPTRE <- rev(sort(p$SCEPTRE))
  p$NB <- rev(sort(p$NB))
  p$Expected <- (rank(p$SCEPTRE, ties.method = "first") + 0.5) / (length(p$SCEPTRE) + 1) # as $Expected is rank based, it is the same as when calculated from $NB 
  # p$Expected2 <- (rank(p$NB, ties.method = "first") + 0.5) / (length(p$NB) + 1)
  
  
# exp.pvalues <- (rank(my.pvalues, ties.method = "first") + 0.5) / (length(my.pvalues) + 1)
  p <- melt(p, id.vars = "Expected")
  
  pdf(file = "../2_DE/Negative Control - qqplot, 10000 sampled tests (V2).pdf", height = 3, width = 4)
  ggplot(p, aes(x = -log10(Expected), y = -log10(value), colour = variable)) +
    geom_point() +
    theme_bw() +
    theme(panel.border = invis, axis.line = element_line(), legend.title = invis) +
    scale_y_continuous(limits = c(0, 50), expand = c(0,0)) +
    scale_x_continuous(limits = c(0, 5), expand = c(0,0)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = "Observed -log10(P)", x = "Expected -log10(P)")
  dev.off()
 

## Empirical p
  r2 <- r[grep("Neg", names(r))]
  
  load("../2_DE/Enhancers - All Results Summary.rda")
  

  hits <- res[which(res$Hit),]
  
  # test the (equal) top hit: Enh53 vs ID3
  neg <- do.call("rbind", sceptre$neg)
  neg.id3 <- neg[which(neg$gene_id == "ID3"),]
    