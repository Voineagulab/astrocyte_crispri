################################################################################################################################ #
## Calculate empirical p-value ----

## Load    
  bg <- do.call("rbind", sceptre$neg) # bg for "background"
    
  # qq plot as a sanity check  
  pdf(file = "Negative Control - qqplot, 10000 sampled tests.pdf", height = 4, width = 4)
  make_qq_plot(bg$p_value[sample(1:nrow(bg), 10000)])
  dev.off()  
  
## For the enh-gene pairs, compare observed p-value to the background distribution of neg-gene pairs
  
  ## This is how Gasperini et al. (2019) calculated it:
  
    # [(the number of NTCs with a smaller P-value than that test’s raw P-value) + 1] divided by [the total number of NTCs tests + 1].
    # These empirical P-values were Benjamini-Hochberg corrected, and those < 0.1 were kept for 10% empirical FDR sets.
  
  ## Loop
    total.ntc <- nrow(bg) + 1
    
    res$EmpiricalP <- "."
  
    for (j in 1:nrow(res)) {
      print(j)
      
      # I note an ambiguity in Gasperini et al.'s phrasing, hence I do two versions
      
      # ## Version 1: a futile endeavour, using only the neg-gene pairs for each given gene
      #   g <- bg$gene_id[j] # get gene
      #   lower <- length(which(bg.p[[g]] < s$p[j])) # number of neg-gene pairs more significant 
      #   e <- (lower + 1) / total.ntc1 # p-value
      #   s$EmpiricalP[j] <- e
      
      ## Version 2: more reasonably,using the neg-gene pairs for all genes as a background
        lower <- length(which(bg$p_value < res$p[j])) # number of neg-gene pairs more significant 
        e <- (lower + 1) / total.ntc2 # p-value
        res$EmpiricalP[j] <- e
        
    }

    
    res$EmpiricalFDR <- p.adjust(res$EmpiricalP, method = "fdr")