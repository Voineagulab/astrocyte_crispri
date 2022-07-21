## Load and save
  # save(s, file = "Temp.rda") # s contains statistics for enhancer-gene pairs, for all highly-expressed genes
  load("Temp.rda")

## Function
  # to apply, choose an n, as the number of hits you want. the function will plot a summary comparing the enhancer-gene pairs of the 
  # top n hits, ranked by p-value, against all other genes
  TopHitQC <- function(n) {
    p <- s[order(s$LocalP),]
    p$Hit <- FALSE  
    p$Hit[1:n] <- TRUE  
    
    # p-value (FDR) threshold equivalence
    equiv <- signif(max(p$LocalFDR[1:n]),3)
    
    # report numbers of hits
    q <- data.frame(Local = n,
                    Sceptre = length(which(p$ExpressionFilteredFDR < 0.1 & p$Hit)))
    q <- melt(q)
    levels(q$variable) <- c(paste0("Local Perturbation\nFDR equivalent of ", equiv), "Overlap With Sceptre\n(92 Hits at FDR < 0.1)")
    
    pA <- ggplot(q, aes(x = variable, y = value)) +
      geom_col() +
      labs(y = "Number of Genes", x = "Category") +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), axis.title.x = invis) +
      scale_y_continuous(expand = c(0,0))
    
    # compare fold-change
    pB <- ggplot(p[which(p$Hit),], aes(x = Localfc, y = logfc_VST)) +
      geom_point() +
      theme_bw() +
      theme(panel.border = invis) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs(x = "Local Perturbation Fold-change\nHits Only", y = "VST Fold-change") +
      geom_abline(slope = -1, intercept = 0, linetype = 2, colour = "red")
    
    # volcano
    pC <- ggplot(p, aes(x = Localfc, y = -log10(LocalP), colour = Hit))  +
      geom_point() +
      theme_bw() +
      theme(panel.border = invis, panel.grid = invis, axis.line.y = element_line()) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      # scale_y_continuous(limits = c(0,17), expand = c(0,0)) +
      # scale_x_continuous(limits = c(-27,27)) +
      labs(x = "log2 fold-change", y = "-log10 Unadjusted P") +
      scale_colour_manual(values = sig.colours)
    
    # fold-change bias
    q <- table(p$Hit, sign(p$Localfc))
    q <- q / rowSums(q)
    q <- as.data.frame(q)
    colnames(q) <- c("Hit", "Sign", "Freq")
    levels(q$Hit) <- c("ns", paste0("Top ", n, " Genes\n(FDR < ", equiv, ")"))
    levels(q$Sign) <- c("Downregulated", "Upregulated")
    
    # pdf(file = "Enhancers - Sign of Fold-Change vs Hit Rate.pdf", height = 3, width = 3)
    pD <- ggplot(q, aes(x = Hit, y = Freq, fill = Sign)) +
      geom_col(colour = "black", width = 0.7) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = c("dodgerblue1", "firebrick1")) +
      labs(y = "Fraction of Enhancer-Gene Pairs", x = "Statistical association")
    
    # distance
    m <- aggregate(p$Distance, list(p$Hit), median)[,2]
    m <- m / 1000
    
    # density plot, comparing hit and non-hit pairings
    pE <- ggplot(p, aes(x = Distance / 1000, colour = Hit, fill = Hit)) +
      geom_density(alpha = 0.1) +
      theme_bw() +
      scale_colour_manual(values = sig.colours) +
      scale_fill_manual(values = sig.colours) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0), limits = c(0, 1010), breaks = c(50, 100, 200, 500, 1000)) +
      geom_vline(xintercept = m, colour = sig.colours, linetype = 2) +
      labs(x = "Distance Between Enhancer-Gene Pair (kilobases)", y = "Density") +
      theme(panel.border = invis, panel.grid = invis, legend.position = c(0.8, 0.8)) 
    
    # nearest gene
    q <- table(p$Hit, p$Nearest)
    q <- q / rowSums(q)
    q <- as.data.frame(q)
    colnames(q) <- c("Hit", "NearestGene", "Freq")
    levels(q$Hit) <- c("ns", paste0("Top ", n, " Genes\n(FDR < ", equiv, ")"))
    
    pF <- ggplot(q, aes(x = Hit, y = Freq, fill = NearestGene)) +
      geom_col(colour = "black", width = 0.7) +
      theme_bw() +
      theme(panel.border = invis, axis.line.y = element_line(), panel.grid = invis) +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = c("white", "black")) +
      labs(y = "Fraction of Enhancer-Gene Pairs", x = "Statistical association")
    
    
    # output plot
    plot_grid(pA, pB, pC, pD, pE, pF, ncol = 2)
  }

## Apply function, runs in ~ 3s
  pdf(file = "Enhancers - Top N Hit Sanity Check (Local Perturbation).pdf", height = 7, width = 8)
  TopHitQC(20)
  TopHitQC(50)
  TopHitQC(92)
  TopHitQC(100)
  TopHitQC(150)
  TopHitQC(250)
  TopHitQC(500)
  TopHitQC(1000)
  TopHitQC(2000)
  dev.off()

