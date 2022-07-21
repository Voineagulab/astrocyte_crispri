# setup
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/")
load("Enhancers - All Results Summary.rda")

# convert sceptre Z-score to p-values
x <- res
x$Sceptre.ZP <- 2 * pnorm(abs(x$Sceptre.Z), lower.tail = FALSE) # this is a two-sided test
x$Sceptre.ZFDR <- p.adjust(x$Sceptre.ZP, method = "fdr") 

table(x$Sceptre.ZFDR < 0.1)
table(x$Sceptre.ZFDR < 0.1, x$Sceptre.FDR < 0.1)

y <- x[which(x$Sceptre.ZFDR < 0.1 & x$Sceptre.FDR > 0.1),]
summary(y$Distance)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 2632   55839  182569  209791  365003  496047

 table(sign(y$logfc.vst))

# -1  1 
# 81 23 

