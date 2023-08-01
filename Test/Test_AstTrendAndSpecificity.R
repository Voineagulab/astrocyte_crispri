tf.biol <- data.frame(TF = unique(tf.robust$TF),
                       EnsID = ".",
                       AstTrend = ".",
                       AstSpecificN = ".",
                       FunBiology = ".")

m <- match(tf.biol$TF, geneInfo$Symbol)
tf.biol$EnsID <- geneInfo$EnsID[m] 

## Trend in Lister snRNA-seq
  load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/5.DevTimecourseAnalysis/ProcessedData/trends.rda")
  m <- match(tf.biol$TF, trends$Symbol)
  tf.biol$AstTrend <- trends$Astro[m]
  
  tf.biol$AstTrendBroad <- "."
  tf.biol$AstTrendBroad[which(as.numeric(tf.biol$AstTrend) <= 7 & tf.biol$AstTrend != "0")] <- "Up"
  tf.biol$AstTrendBroad[which(as.numeric(tf.biol$AstTrend) >= 8)] <- "Dn"

## Astrocyte specificity
  load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/Annotation/sigsBrain.rda")
  sigsBrain <- sigsBrain[c("IP", "DM", "VL", "NG", "CA", "LK")]
  
  
  spec <- lapply(sigsBrain, function(a) { # gives list of all tfs in the dataset which are most highly expressed in Ast
    a <- a[intersect(rownames(a), tf.biol$EnsID),]
    b <- apply(a, 1, which.max)
    names(b)[which(b == grep("Astrocytes", colnames(a)))]
  })
  
  spec <- do.call("c", spec)
  spec <- table(spec) # number of datasets in which the tf is ast-highest (zero excluded)
  
  m <- match(tf.biol$EnsID, names(spec))
  tf.biol$AstSpecificN <- spec[m]
  
  
## Save
  write.csv(tf.biol, file = "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/TF Biological Information.csv")
  
  
  
  
## Hit gene level analysis
  hitGene.biol <- data.frame(Gene = unique(res.final$Gene[which(res.final$HitPermissive)]),
                       EnsID = ".",
                       AstTrend = ".",
                       AstSpecificN = ".")
  
  m <- match(hitGene.biol$Gene, geneInfo$Symbol)
  hitGene.biol$EnsID <- geneInfo$EnsID[m] 
  
  # add trend
  m <- match(hitGene.biol$Gene, trends$Symbol)
  hitGene.biol$AstTrend <- trends$Astro[m]
  
    hitGene.biol$AstTrendBroad <- "."
  hitGene.biol$AstTrendBroad[which(as.numeric(hitGene.biol$AstTrend) <= 7 & hitGene.biol$AstTrend != "0")] <- "Up"
  hitGene.biol$AstTrendBroad[which(as.numeric(hitGene.biol$AstTrend) >= 8)] <- "Dn"
  
  # add specificity
  spec <- lapply(sigsBrain, function(a) { # gives list of all tfs in the dataset which are most highly expressed in Ast
    a <- a[intersect(rownames(a), hitGene.biol$EnsID),]
    b <- apply(a, 1, which.max)
    names(b)[which(b == grep("Astrocytes", colnames(a)))]
  })
  
  spec <- do.call("c", spec)
  spec <- table(spec) # number of datasets in which the tf is ast-highest (zero excluded)
  
  m <- match(hitGene.biol$EnsID, names(spec))
  hitGene.biol$AstSpecificN <- spec[m]
  
  

  
## Extra fun TFs
  fun <- read.csv("~/Downloads/AstFunTFs.csv")
  colnames(fun) <- "Gene"
  
  m <- match(fun$Gene, geneInfo$Symbol)
  fun$EnsID <- geneInfo$EnsID[m] 
  
  spec <- lapply(sigsBrain, function(a) { # gives list of all tfs in the dataset which are most highly expressed in Ast
    a <- a[intersect(rownames(a), fun$EnsID),]
    b <- apply(a, 1, which.max)
    names(b)[which(b == grep("Astrocytes", colnames(a)))]
  })
  
  spec <- do.call("c", spec)
  spec <- table(spec) # number of
  
  m <- match(fun$EnsID, names(spec))
  fun$AstSpecificN <- spec[m]
  
  m <- match(fun$Gene, trends$Symbol)
  fun$AstTrend <- trends$Astro[m]
  