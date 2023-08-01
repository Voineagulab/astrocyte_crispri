## Load

test <- read.csv("../Enh/Results Summary.csv", row.names = 1)
load("SCEPTRE Output (Neg, Guide-level).rda")
load("SCEPTRE Output (NegE, Guide-level).rda")
neg <- list()
neg$Neg <- de.neg.guidelevel
neg$Enh <- de.negE

## Correct
total.ntc <- lapply(neg, function(x) nrow(x) + 1)
    
    test$EmpiricalP.Neg <- NA
    test$EmpiricalP.Enh <- NA
  
    for (j in 1:nrow(test)) {
      print(j)
      
      # I note an ambiguity in Gasperini et al.'s phrasing, hence I do two versions
      
      # ## Version 1: a futile endeavour, using only the neg-gene pairs for each given gene
      #   g <- bg$gene_id[j] # get gene
      #   lower <- length(which(bg.p[[g]] < s$p[j])) # number of neg-gene pairs more significant 
      #   e <- (lower + 1) / total.ntc1 # p-value
      #   s$EmpiricalP[j] <- e
      
      ## Version 2: more reasonably,using the neg-gene pairs for all genes as a background
        # for the smaller enh pool of 50
        lower <- length(which(neg$Enh$p_value < test$P[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Enh # p-value
        test$EmpiricalP.Enh[j] <- e
        
        # for the larger neg pool of 250
        lower <- length(which(neg$Neg$p_value < test$P[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc$Neg # p-value
        test$EmpiricalP.Neg[j] <- e
        
    }

    
    test$EmpiricalFDR.Enh <- p.adjust(test$EmpiricalP.Enh, method = "fdr")
    test$EmpiricalFDR.Neg <- p.adjust(test$EmpiricalP.Neg, method = "fdr")
    
    
    
## Plot
  thresh <- -log10(max(test$P[which(test$Hit)]))
  pdf(file = "../../Scratchspace/EmpiricalP SCEPTRE.pdf", height = 4, width = 4)
  ggplot(test, aes(x = -log10(P), y = -log10(EmpiricalP.Neg))) +
    geom_point() +
    theme_bw() +
    labs(x = "Raw SCEPTRE -log10P", y = "EmpiricalN SCEPTRE -log10P") +
    scale_y_continuous(limits = c(0,17)) +
    geom_vline(xintercept = thresh, linetype = 2, colour = "red") +
    geom_hline(yintercept = thresh, linetype = 2, colour = "red")
  
  
  ggplot(test, aes(x = -log10(P), y = -log10(EmpiricalP.Enh))) +
    geom_point() +
    theme_bw() +
    labs(x = "Raw SCEPTRE -log10P", y = "EmpiricalE SCEPTRE -log10P") +
    scale_y_continuous(limits = c(0,17)) +
    geom_vline(xintercept = thresh, linetype = 2, colour = "red") +
    geom_hline(yintercept = thresh, linetype = 2, colour = "red")
  dev.off()
  