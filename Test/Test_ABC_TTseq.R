# # load res
# 
# 
# # rm(list=ls())
# 
# getStats=function(thresholds, var)
# {
#   th.stats=as.data.frame(thresholds); colnames(th.stats)[1]="Threshold"
#   th.stats$TPR=NA; th.stats$TNR=NA; th.stats$FPR=NA; th.stats$FNR=NA; th.stats$PPV=NA; th.statsNPV=NA; th.stats$ACC=NA
#   
#   for (j in c(1:nrow(th.stats)))
#   {
#     t=th.stats[j,1]
#     realpos=which(screen$Hit ==TRUE); realneg= which(screen$Hit ==FALSE)
#     pos=which(var >= t); neg=which(var < t)
#     
#     th.stats$TPR[j]=length(intersect(realpos, pos))/length(realpos)
#     th.stats$TNR[j]=length(intersect(realneg, neg))/length(realneg)
#     
#     th.stats$FPR[j]=length(intersect(realneg, pos))/length(pos)
#     th.stats$FNR[j]=length(intersect(realneg, neg))/length(neg)
#     
#     th.stats$PPV[j]=length(intersect(realpos, pos))/length(pos)
#     th.stats$NPV[j]=length(intersect(realneg, neg))/length(neg)
#     
#     th.stats$ACC[j]=(length(intersect(realpos, pos)) + length(intersect(realneg, neg))) /( length(realneg)+ length(realpos))
#     
#   }
#   
#   return(th.stats)
# }
# 
# plotpairs=function(data, x, y, title)
# {
#   plot(data[,x] ~ data[,1], pch=20, col="red" ,ylim=c(0,1), 
#        main=paste(title, x, y, sep=","),
#        ylab=paste(x, ":red;", y, ":blue"),
#        xlab="Significance Threshold")
#   points(data[,y] ~ data[,1], pch=20, col="steelblue")
#   abline(h=0.5, col="red")
# }
#################

## Load

  # tt-seq data for screen enhancers
  load("/mnt/Data0/PROJECTS/CROPSeq/TTseq/RESULTS/TTseq_data.rda")
  ttseq <- data; rm(data)
  
  m <- match(ttseq$Id, guides$TargetCoord)
  ttseq$Enh <- guides$TargetID[m]
  write.csv(ttseq, file = "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/TTseq.csv")
  
  #CRISPRi data
  # load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Archive_V1/Enhancers - All Results Summary.rda")
  # res.sig=res[which(res$Hit==TRUE) , ]
  # res.sig <- res
  
  #abc data
  abc <- read.delim("/mnt/Data1/PROJECTS/Farbod_ABC/ABC_results/predictions/EnhancerPredictions_NHAPeaks.bed", sep="\t")
  colnames(abc) <- c("NHA_chr", "NHA_start", "NHA_end", "NHA_id", "NHA_score", "NHA_strand",
                  "ABC_chr", "ABC_start", "ABC_end", "ABC_id", "Gene", "Score", "CellType", "ABC_Score")
  
  abc$ABC_id <- gsub("intergenic|", "", abc$ABC_id, fixed = TRUE)
  abc$NHA_id <- paste0(abc$NHA_chr, ":", abc$NHA_start, "-", abc$NHA_end)
  abc <- abc[which(abc$NHA_id %in% x$Enh.Pos), ]
  # abc.screen <- abc[which(abc$NHA_id %in% x$Enh.Pos) , ]
  # abc.screen.sig <- abc[which(abc$NHA_id %in% x.sig$Enh.Pos) , ]
  
  
  # abc.sorted <- abc[order(abc$ABC_Score, decreasing = TRUE) , ]
  # res.sorted <- x[order(x$Z) , ]

## Evaluate ABC pairs
  # get enhancer and pair annotation for ABC
  m <- match(abc$NHA_id, guides$TargetCoord)
  abc$Enh <- guides$TargetID[m]
  abc$Pair <- paste0(abc$Enh, "_", abc$Gene)  
  
  # match
  abc.replicate <- abc[which(abc$Pair %in% x$Pair), c("Enh", "Gene", "Pair", "NHA_id", "ABC_id", "ABC_Score")]
  m <- match(abc.replicate$Pair, x$Pair)
  abc.replicate$Screen_FDR <- x$FDR[m]
  abc.replicate$Screen_Z <- x$Z[m]
  abc.replicate$Screen_FDR <- x$FDR[m]
  abc.replicate$Hit <- x$Hit[m]
  
  table(abc.replicate$Hit)
  
  write.csv(abc.replicate, file = "ABC Replication.csv")
  
  # reverse that match...
  m <- match(x$Pair, abc.replicate$Pair)
  x$ABC.Score <- abc.replicate$ABC_Score[m]
  x$ABC.Hit <- x$Pair %in% abc.replicate$Pair
  table(x$ABC.Hit, x$Hit) %>% fisher.test()
  
  table(x$ABC.Hit, x$Hit) %>% fisher.test()
  h <- x[which(x$Gene.Distance > 10000),]
  table(h$ABC.Hit, h$Hit) %>% fisher.test()
  h <- x[which(x$Gene.Distance < 10000),]
  table(h$ABC.Hit, h$Hit) %>% fisher.test()
  
## Evaluate the general "strength" of the enhancer in ABC 
  e <- data.frame(Enh = unique(guides$TargetID[which(guides$TargetCat == "Enh")]))
  e$Hit <- e$Enh %in% x$Enh[which(x$Hit)]
  e$Hit <- factor(e$Hit)
  levels(e$Hit) <- c("ns", "Hit")
  e$ABC.any <- e$Enh %in% abc.replicate$Enh
  x <- abc.replicate[order(abc.replicate$ABC_Score, decreasing = TRUE),]
  e$ABC.top <- x$ABC_Score[match(e$Enh, x$Enh)]
  
## Evaluate TTseq expression of enhancers...
  m <- match(e$Enh, guides$TargetID)
  e$Coord <- guides$TargetCoord[m]  
  m <- match(e$Coord, ttseq$Id)
  e$TT.count <- ttseq$TTseq_NHA[m]
  e$TT.exp <- ttseq$TTseq_exp[m]
  e$TT.enrich <- ttseq$TTseq_enriched[m]
  
  table(e$Hit, e$TT.exp) %>% fisher.test()
  table(e$Hit, e$TT.count > 2) %>% fisher.test()
  table(e$Hit, e$TT.enrich) %>% fisher.test()
  

