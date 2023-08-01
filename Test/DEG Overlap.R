setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Neg/")


# I will need:
  # a nice table of overlaps at 0.1 and 0.05
  # sanity checks
  # correspondence to raw NB p?

## Make overlaps
  a <- "NB.FDR" # NB
  b <- "EmpiricalFDR.Enh" # N50
  c <- "EmpiricalFDR.Neg" # N250
  d <- "FDR" # S
  keyCols <- enh[,c(a,b,c,d)]
  colnames(keyCols) <- c("NB", "N50", "N250", "S")
  

  calc.overlap <- function(cols, thresh, nRep) {
    if (length(cols) == 1) {
      e <- length(which(as.numeric(keyCols[,cols]) < thresh))
    } else {
      e <- keyCols[,cols]
      e <- e < thresh
      e <- rowSums(e)
      e <- length(which(e >= nRep))  
    }
    
    return(e)
  }
  
  # calc.overlap(c("EmpiricalFDR.Enh", "EmpiricalFDR.Neg", "FDR"), 0.1, 3)
  # calc.overlap(c("EmpiricalFDR.Enh", "EmpiricalFDR.Neg", "FDR"), 0.1, 2)
  # calc.overlap(c("EmpiricalFDR.Enh", "EmpiricalFDR.Neg", "FDR"), 0.1, 2)
  # calc.overlap(c("FDR", "FDR"), 0.1, 1)
  # 
  # cols <- c("EmpiricalFDR.Enh", "EmpiricalFDR.Neg", "FDR")
  # 
  # calc.overlap(c("EmpiricalFDR.Enh", "FDR"), 0.1, 2)
  # calc.overlap(c("EmpiricalFDR.Enh", "FDR"), 0.1, 1)
  # calc.overlap(c("EmpiricalFDR.Enh", "FDR"), 0.05, 2)
  # calc.overlap(c("EmpiricalFDR.Enh", "FDR"), 0.05, 1)
  # 
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.1, 2)
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.1, 1)
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.05, 2)
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.05, 1)
  # 
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.1, 2)
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.05, 2)
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.1, 1)
  # calc.overlap(c("EmpiricalFDR.Neg", "FDR"), 0.05, 1)
  

z <- data.frame(nAlgs = c(rep(1, 4), rep(2, 6), rep(3, 4), rep(4,1)), 
                NB = ".", N50 = ".", N250 = ".", S = ".", 
                FDR10 = NaN, FDR10_inAny = NaN, FDR10_in2 = NaN,
                FDR05 = NaN, FDR05_inAny = NaN, FDR05_in2 = NaN)
rownames(z) <- c("NB", "N50", "N250", "S",
                 "NB+N50", "NB+N250", "NB+S", "N50+N250", "N50+S", "N250+S",
                 "NB+N50+N250", "NB+N50+S", "NB+N250+S", "N50+N250+S",
                 "NB+N50+N250+S")
# rownames(a) <- paste0(a$nAlgs, "_", rownames(a))  

## Loop
  for (j in rownames(z)) {
    # parse rowname
    
    k <- strsplit(j, "\\+") %>% do.call("c", .)
    
    # get overlap at FDR < 0.1 
    z[j,"FDR10"] <- calc.overlap(k, 0.1, length(k)) # in all algorithms shown
    z[j,"FDR10_inAny"] <- calc.overlap(k, 0.1, 1) # in at least one algorithm
    z[j,"FDR10_in2"] <- calc.overlap(k, 0.1, min(length(k), 2)) # in at least two algorithm, if available
    
    
    
    # get overlap at FDR < 0.05, in all algorithms
    z[j,"FDR05"] <- calc.overlap(k, 0.05, length(k))
    z[j,"FDR05_inAny"] <- calc.overlap(k, 0.05, 1) # in at least one algorithm
    z[j,"FDR05_in2"] <- calc.overlap(k, 0.05, min(length(k), 2)) # in at least two algorithm, if available
  }



## In excel, white text the useless ones!!
write.csv(z, file = "../Algorithm Consistency (Enh-level).csv")
  
################################################################################################################################ #
## Guide-level analyses!! ----
  
## This output table needs to be a number of hits + number novel in bracket
m <- match(gl$Pair.enh, enh$Pair)
gl[,paste0(a, "_pool")] <- enh[m,a]
gl[,paste0(b, "_pool")] <- enh[m,b]
gl[,paste0(c, "_pool")] <- enh[m,c]
gl[,paste0(d, "_pool")] <- enh[m,d]


x <- gl[,c(a,b,c,d, paste0(c(a,b,c,d), "_pool"))]
colnames(x) <- c("NB", "N50", "N250", "S", paste0(c("NB", "N50", "N250", "S"), "_pool"))

## Function
  calc.overlap.gl <- function(cols, thresh, nRep) {
    if (length(cols) == 1) {
      e1 <- length(which(as.numeric(x[,cols]) < thresh))
      e2 <- length(which(as.numeric(x[,cols]) < thresh & as.numeric(x[,paste0(cols, "_pool")]) < thresh))
      e <- paste0(e1, " (+", e1-e2, ")")
    } else {
      e1 <- x[,cols]
      e1 <- e1 < thresh
      e1 <- rowSums(e1)
      e1 <- ((e1 >= nRep))
      
      e2 <- x[,paste0(cols, "_pool")]
      e2 <- e2 < thresh
      e2 <- rowSums(e2)
      e2 <- ((e2 >= nRep))
      
      e <- paste0(length(which(e1)), " (+", length(which(e1 & !(e2))), ")")
    }
    
    return(e)
    }
    
## Calculate overlaps
    w <- data.frame(nAlgs = c(rep(1, 4), rep(2, 6), rep(3, 4), rep(4,1)), 
                NB = ".", N50 = ".", N250 = ".", S = ".", 
                FDR10 = NA, FDR10_inAny = NA, FDR10_in2 = NA,
                FDR05 = NA, FDR05_inAny = NA, FDR05_in2 = NA)
rownames(w) <- c("NB", "N50", "N250", "S",
                 "NB+N50", "NB+N250", "NB+S", "N50+N250", "N50+S", "N250+S",
                 "NB+N50+N250", "NB+N50+S", "NB+N250+S", "N50+N250+S",
                 "NB+N50+N250+S")
# rownames(a) <- paste0(a$nAlgs, "_", rownames(a))  

## Loop
  for (j in rownames(w)) {
    # parse rowname
    
    k <- strsplit(j, "\\+") %>% do.call("c", .)
    
    # get overlap at FDR < 0.1 
    w[j,"FDR10"] <- calc.overlap.gl(k, 0.1, length(k)) # in all algorithms shown
    w[j,"FDR10_inAny"] <- calc.overlap.gl(k, 0.1, 1) # in at least one algorithm
    w[j,"FDR10_in2"] <- calc.overlap.gl(k, 0.1, min(length(k), 2)) # in at least two algorithm, if available
    
    
    
    # get overlap at FDR < 0.05, in all algorithms
    w[j,"FDR05"] <- calc.overlap.gl(k, 0.05, length(k))
    w[j,"FDR05_inAny"] <- calc.overlap.gl(k, 0.05, 1) # in at least one algorithm
    w[j,"FDR05_in2"] <- calc.overlap.gl(k, 0.05, min(length(k), 2)) # in at least two algorithm, if available
  }

## Save
write.csv(w, file = "../Algorithm Consistency (Guide-level).csv")



######################## V2
## This will count the number of hits which are also supported by >2 guides

## Setup
y <- gl[,c("Pair.enh",a,b,c,d, paste0(c(a,b,c,d), "_pool"))]
colnames(y) <- c("Pair","NB", "N50", "N250", "S", paste0(c("NB", "N50", "N250", "S"), "_pool"))

## Function
  calc.overlap.glSupport <- function(cols, thresh, nRep) {
    if (length(cols) == 1) {
      # e1 is the number of pairs with  at least 2 sig guides
      e1 <- (which(as.numeric(y[,cols]) < thresh))
      e1 <- length(which(table(y$Pair[e1]) >= 2))
      
      # e2 is the number of pairs sig at the pool level
      e2 <- (which(as.numeric(y[,paste0(cols, "_pool")]) < thresh))
      e2 <- length(unique(y$Pair[e2]))
      
      # output  
      e <- paste0(e2, " (", e1, ")")
    } else {
      # e1 is the number of pairs with at least 2 sig guides in all algs...
      e1 <- y[,cols]
      e1 <- e1 < thresh
      e1 <- rowSums(e1)
      e1 <- which(e1 >= nRep)
      e1 <- length(which(table(y$Pair[e1]) >= 2))
      
      # e2 is the number of pairs sig at the pool level
      e2 <- y[,paste0(cols, "_pool")]
      e2 <- e2 < thresh
      e2 <- rowSums(e2)
      e2 <- which(e2 >= nRep)
      e2 <- length(unique(y$Pair[e2]))
      
      e <- paste0(e2, " (", e1, ")")
    }
    
    return(e)
  }
  
  
## Calculate overlaps
    q <- data.frame(nAlgs = c(rep(1, 4), rep(2, 6), rep(3, 4), rep(4,1)), 
                NB = ".", N50 = ".", N250 = ".", S = ".", 
                FDR10 = NA, FDR10_inAny = NA, FDR10_in2 = NA,
                FDR05 = NA, FDR05_inAny = NA, FDR05_in2 = NA)
rownames(q) <- c("NB", "N50", "N250", "S",
                 "NB+N50", "NB+N250", "NB+S", "N50+N250", "N50+S", "N250+S",
                 "NB+N50+N250", "NB+N50+S", "NB+N250+S", "N50+N250+S",
                 "NB+N50+N250+S")
# rownames(a) <- paste0(a$nAlgs, "_", rownames(a))  

## Loop
  for (j in rownames(q)) {
    # parse rowname
    
    k <- strsplit(j, "\\+") %>% do.call("c", .)
    
    # get overlap at FDR < 0.1 
    q[j,"FDR10"] <- calc.overlap.glSupport(k, 0.1, length(k)) # in all algorithms shown
    q[j,"FDR10_inAny"] <- calc.overlap.glSupport(k, 0.1, 1) # in at least one algorithm
    q[j,"FDR10_in2"] <- calc.overlap.glSupport(k, 0.1, min(length(k), 2)) # in at least two algorithm, if available
    
    
    
    # get overlap at FDR < 0.05, in all algorithms
    q[j,"FDR05"] <- calc.overlap.glSupport(k, 0.05, length(k))
    q[j,"FDR05_inAny"] <- calc.overlap.glSupport(k, 0.05, 1) # in at least one algorithm
    q[j,"FDR05_in2"] <- calc.overlap.glSupport(k, 0.05, min(length(k), 2)) # in at least two algorithm, if available
  }

## Save
write.csv(q, file = "../Algorithm Consistency (Pool- and Guide-level).csv")

######################## V3
## This will count the number of hits which are also supported by >2 guides at a (bonferroni-adjusted) p < 0.05/0.01/0.001 in the same direction

v3 <- gl[c(6,7,5,
           4,14,21,
           9,12,20,
           10,11,19)]
colnames(v3) <- c("Pair.Pool", "Pair.GL", "Sign", "S.P", "S.FDR", "S.Pool", "N250.P", "N250.FDR", "N250.Pool", "N50.P", "N50.FDR", "N50.Pool")
v3$Sign <- sign(v3$Sign)

v3.fun <- function(alg, thresh = 0.05) {
  # get hits at the pool level
  # subset <- v3[which(v3[,paste0(alg, ".Pool")] < 0.1),]
  subset <- v3
  
  # for each 
  l <- list()
  
  for (j in unique(subset$Pair.Pool)) {
    rep <- subset[which(subset$Pair.Pool == j),]
    nRep <- length(which(rep[,paste0(alg, ".P")] < (thresh / nrow(rep))))
    l[[j]] <- nRep >= 2
  }
  
  l <- do.call("c", l)
  l <- names(l)[which(l)]
  # table(l)
  
}

# of the 90 SCEPTRE hits at FDR 0.1
v3.fun("S", thresh = 0.05) # 60
v3.fun("S", thresh = 0.01) # 44
v3.fun("S", thresh = 0.001) # 27

# of the 158 SCEPTRE hits at N50 0.1
v3.fun("N50", thresh = 0.05) # 81
v3.fun("N50", thresh = 0.01) # 57
v3.fun("N50", thresh = 0.001) # 29

# of the 63 SCEPTRE hits at N250 0.1
v3.fun("N250", thresh = 0.05) # 54
v3.fun("N250", thresh = 0.01) # 37
v3.fun("N250", thresh = 0.001) # 16



n50 <- v3.fun("N50", 0.05)
s <- unique(v3$Pair.Pool[which(v3$S.FDR < 0.1)])
s2 <- unique(res$Pair[which(res$FDR < 0.1)])

union(s, n50)
union(s2, n50)

