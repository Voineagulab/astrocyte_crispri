
load("../2_DE/Enhancers - All Results Summary.rda")


## Transcription factor binding: binary
  jasp <- read.csv("../3_HitEnrichment/Jaspar - Overlaps.csv")
  jasp <- dcast(jasp, Peak~tf) 
  
  # rename peaks
  jasp$Peak <- sub("_", ":", jasp$Peak)
  jasp$Peak <- sub("_", "-", jasp$Peak)
  
  rownames(jasp) <- jasp$Peak
  
  # hit peaks
  jasp$Hit <- jasp$Peak %in% res$Enh.Pos[which(res$Hit)]
  # input <- as.data.frame(t(jasp[,-1]))
  

  
  # logistic
  mod <- glm(Hit ~., family = binomial(link = "logit"), data = jasp[,-1])
  sum <- summary(mod)
  
  # 1: glm.fit: algorithm did not converge
  # 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
  
## Overenrichment
  jasp[,-c(1,331)] <- jasp[,-c(1,331)] > 0
  jasp <- jasp[,-1]
  
  ## Analyse
  # x <- apply(vld[,annot.cols], 2, table) %>% t() %>% as.data.frame
  # colnames(x) <- c("NoOverlap", "Overlap")

  jasp.fish <- apply(jasp[,-330], 2, function(y) { 
    if (length(table(y)) == 1) return(NA) 
    f <- fisher.test(table(y, jasp$Hit))
    data.frame(OR = f$estimate, P = f$p.value)
  })
  
  jasp.fish <- do.call("rbind", jasp.fish)
  jasp.fish$Bonf <- p.adjust(jasp.fish$P, method = "bonferroni")  
  jasp.fish$FDR <- p.adjust(jasp.fish$P, method = "fdr")  
  
## One for technical factors
  
## One for biological factors