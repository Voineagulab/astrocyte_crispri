# load res from sceptre
x <- read.csv("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Enh/Results Summary.csv")

# load d from local cooper
# load("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Scratchspace/Cooper (Local) DE Method.rda")
#     mean <- rowMeans(nha@assays$RNA@data) # rna assay, not sct
#     mean <- names(mean)[which(mean > 2^-6)]
# localpb <- lapply(local.cooper, function(x) {
#       if (class(x$DE) == "character") {
#         x <- data.frame(Mean = "Error",
#                         log2fc = "Error",
#                         p = "Error")
#       } else {
#         x <- x$DE
#         x <- data.frame(Mean = x$baseMean,
#                         log2fc = x$log2FoldChange,
#                         p = x$pvalue,
#                         Gene = rownames(x))
#       }
#       
#       x <- x[which(x$Gene %in% mean),] # mean was made in the temp sceptre script
#       
#       return(x)
#  })
#  
#  localpb <- do.call("rbind", localpb)
#  localpb$Enh <- splitter(rownames(localpb), "\\.", 1)
#  localpb$Pair <- paste0(localpb$Enh, "_", localpb$Gene)
#  localpb <- localpb[which(localpb$Pair %in% x$Pair),]
  load("../Scratchspace/Cooper (Signed V2, Processed).rda")
  load("../Scratchspace/Cooper (Local V2, New PCA) DE Method.rda")
  localpb <- cooper.signed.v2

## Add new FDRs
 x$NegBinomFDR <- 2 * pnorm(abs(x$Z), lower.tail = FALSE) %>% p.adjust(p = ., method = "fdr")
m <- match(x$Pair, localpb$Pair)
x$LocalPB.FDR <- localpb$FDR[m]
x$LocalPB.FC <- localpb$Mean[m]

## Categorise hits
p.thresh <- 0.05
  x$Cat <- "."
  x$Cat[which(x$Hit & x$NegBinomFDR < p.thresh & x$LocalPB.FDR < p.thresh)] <- "3: All"
  x$Cat[which((x$Hit) & !(x$NegBinomFDR) < p.thresh & !(x$LocalPB.FDR) < p.thresh)] <- "1: SCEPTRE"
  x$Cat[which(!(x$Hit) & (x$NegBinomFDR) < p.thresh & !(x$LocalPB.FDR) < p.thresh)] <- "1: NB"
  x$Cat[which(!(x$Hit) & !(x$NegBinomFDR) < p.thresh & (x$LocalPB.FDR) < p.thresh)] <- "1: LocalPB"
  x$Cat[which(!(x$Hit) & (x$NegBinomFDR) < p.thresh & (x$LocalPB.FDR) < p.thresh)] <- "2: NB+LocalPB"
  x$Cat[which((x$Hit) & !(x$NegBinomFDR) < p.thresh & (x$LocalPB.FDR) < p.thresh)] <- "2: Sceptre+LocalPB"
  x$Cat[which((x$Hit) & (x$NegBinomFDR) < p.thresh & !(x$LocalPB.FDR) < p.thresh)] <- "2: Sceptre+NB"
  x$Cat[which(!(x$Hit) & !(x$NegBinomFDR) < p.thresh & !(x$LocalPB.FDR) < p.thresh)] <- "0: None"
  
  x$Cat <- factor(x$Cat)
  levels(x$Cat) <- paste0(levels(x$Cat), " (", table(x$Cat), ")")
  
  
## Compare on basic sanity check of direction and distance
  pdf(file = "Core Hit Definition - Violin Sanity Check (FDR05).pdf", height = 4, width = 7)
  ggplot(x, aes(x = Cat, y = Z)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
    labs(y = "Z (from NB/SCEPTRE)") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
  
  ggplot(x, aes(x = Cat, y = LocalPB.FC)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
    labs(y = "Fold Change (from LocalPB)") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
  
   ggplot(x, aes(x = Cat, y = Gene.Distance)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
     labs(y = "Enh-gene Distance") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
   
   ggplot(x, aes(x = Cat, y = Gene.Exp)) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
     labs(y = "Gene Exp") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis)
  dev.off()  
  
## Compare on P-values 
   pdf(file = "Core Hit Definition - Violin P-values.pdf", height = 4, width = 7)
  ggplot(x, aes(x = Cat, y = -log10(P))) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      # geom_hline(yintercept = 0) +
     labs(y = "SCEPTRE P") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis) 
  
  ggplot(x, aes(x = Cat, y = -log10(NegBinomFDR))) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      # geom_hline(yintercept = 0) +
     labs(y = "NB FDR") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis) 
  
  ggplot(x, aes(x = Cat, y = -log10(NegBinomFDR))) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      # geom_hline(yintercept = 0) +
     labs(y = "NB FDR") +
    scale_y_continuous(limits = c(0,10)) +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis) 
  
  ggplot(x, aes(x = Cat, y = -log10(LocalPB.FDR))) +
      geom_violin(scale = "width", draw_quantiles = 0.5, colour = "red") +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.2) +
      theme_bw() +
      # geom_hline(yintercept = 0) +
     labs(y = "LocalPB FDR") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis) 
  dev.off()
  
## Look at TTSeq and ABC
  # abc
  abc.replicate <- read.csv("ABC Replication.csv")
  m <- match(x$Pair, abc.replicate$Pair)
  x$ABC.Score <- abc.replicate$ABC_Score[m]
  x$ABC.Hit <- x$Pair %in% abc.replicate$Pair
  
  table(x$Cat, x$ABC.Hit)
  
  p <- table(x$Cat, x$ABC.Hit)
  p <- p / rowSums(p)
  p <- data.frame(p[,2])
  
  p$Cat <- rownames(p)
  
  ggplot(p, aes)
  
  # ttseq
  ttseq <- read.csv("TTseq.csv")
  m <- match(x$Enh, ttseq$Enh)
  x$TT.count <- ttseq$TTseq_NHA[m]
  x$TT.exp <- ttseq$TTseq_exp[m]
  x$TT.enrich <- ttseq$TTseq_enriched[m]
  
  table(x$Cat, x$TT.exp)
  table(x$Cat, x$TT.count > 9)
  
  table(x$Cat, x$TT.enrich)
  
  
  
#   
# ## Quick look at TTseq and ABC
#   source("mnt/Data0/PROJECTS/CROPSeq/IV/SCRIPTS/utils/Functions.R")
# #tt-seq data for screen enhancers
# load("/mnt/Data0/PROJECTS/CROPSeq/TTseq/RESULTS/TTseq_data.rda")
# ttseq=data; rm(data)
# #CRISPRi data
# # load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/2_DE/Archive_V1/Enhancers - All Results Summary.rda")
# # res.sig=res[which(res$Hit==TRUE) , ]
# #abc data
# abc=read.delim("/mnt/Data1/PROJECTS/Farbod_ABC/ABC_results/predictions/EnhancerPredictions_NHAPeaks.bed", sep="\t")
# colnames(abc)=c("NHA_chr", "NHA_start", "NHA_end", "NHA_id", "NHA_score", "NHA_strand","ABC_chr", "ABC_start", "ABC_end", "ABC_id", "Gene", "Score", "CellType", "ABC_Score")
# abc$ABC_id=gsub("intergenic|", "", abc$ABC_id, fixed=TRUE)
# abc$NHA_id=paste0(abc$NHA_chr, ":", abc$NHA_start, "-", abc$NHA_end)
# abc.screen=abc[which(abc$NHA_id%in%res$Enh.Pos) , ]
# # abc.screen.sig=abc[which(abc$NHA_id%in%res.sig$Enh.Pos) , ]
# 
# 
# abc.sorted=abc[order(abc$ABC_Score, decreasing = TRUE) , ]
# res.sorted=res[order(res$Z) , ]
# 
# # Combine data at enhancer level
# screen=res.sorted[match(unique(res.sorted$Enh), res.sorted$Enh) , ]
# 
# screen$Sig=screen$Enh.Pos%in%res$Enh.Pos
# screen$ABCSig=screen$Enh.Pos%in%abc$NHA_id
# screen$ABCHighSig=screen$Enh.Pos%in%abc$NHA_id[which(abc$ABC_Score > 0.1)]
# screen$MaxABCScore=abc.sorted$ABC_Score[match(screen$Enh.Pos, abc.sorted$NHA_id)]
# screen$MaxABCScore[which(is.na(screen$MaxABCScore)==TRUE)]=0
# 
# screen=cbind(screen, ttseq[match(screen$Enh.Pos, ttseq$Id) , ])
# 
# th.abc=getStats(seq(from=0.02, to=0.1, by=0.01) ,screen$MaxABCScore )
# th.ttseq.counts=getStats(seq(from=2, to=10, by=1) ,screen$TTseq_NHA )
# th.ttseq.ratio=getStats(seq(from=2, to=10, by=1) ,screen$cpm_Ratio )

  