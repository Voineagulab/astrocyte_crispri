## This script preprocesses all counts data output from CellRanger using R

################################################################################################################################ #
## Setup ----


## Generic
rm(list = ls()); gc()
setwd("/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/EnhGenePairs/")
options(stringsAsFactors = FALSE)

## Packages, functions, and libraries
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(reshape2)
  library(readxl)
  library(tidyr)
  library(liftOver)
  library(rtracklayer)
  library(rcartocolor)

  source("../../../Scripts/Functions.R")

## Load
  # expression
  load("../../../Data/Preprocessed/NHA Pooled (Final).rda")
  guides <- read.csv("../../../Data/Whitelists/Protospacer Whitelist (NHA).csv", row.names = 1)
  guides <- guides[which(guides$Celltype == "NHA"),]
  
  # results from earlier scripts
  res.final <- read.csv("../../2_DE/Enh/Results Final.csv")  
  annot.enh <- read.csv("../Chromatin/Final - Annotation Logical.csv")
  annot.gene <- read.csv("../Genes/Final - Annotation Logical.csv")
  
  
  
## Plotting
  maxh <- 24 / 2.54
  maxw <- 14.9/ 2.54
  invis <- element_blank()
  
  sig.colours <- c("black", "firebrick1")

## Setup dataframe
  ## Generic dataframe to store annotations pertinent at the pair level
    x <- res.final[,c("Pair", "Enh", "Gene", "HitPermissive")]
    rownames(x) <- x$Pair
    colnames(x)[4] <- "Hit"
  
    # make a list, with a level for coenrichment analysis, and replication of the pair link
    # egp <- list()
    # egp$Coenrich <- egp$Replicate <- x
    egp_replicate <- x
    egp_coannot <- x
    
  ## Functions to get annotation
    isGene <- function(col) {
      annot.gene$Gene[which(annot.gene[,col])]
    }
    
    isEnh <- function(col) {
      annot.enh$Enh[which(annot.enh[,col])]
    } 
    
  ## Functions to calculate enrichment
    calc.coenrich <- function(enhCol, geneCol, label, rbind = FALSE, rbind.to = coannot_stats) {
      # setup data
      x <- egp_coannot[,c("Hit", geneCol, enhCol)]
      colnames(x) <- c("Hit", "Gene", "Enh")
      x$Concord <- x$Gene & x$Enh
      total <- nrow(x)
      
      # run stats
      f <- table(x$Hit, x$Concord) %>% fisher.test()
      
      
      # calculate stats
      out <- data.frame(EnhSet = enhCol,
                        GeneSet = geneCol,
                        Total_Enh_TRUE = sum(x$Enh),
                        Total_HitEnh_TRUE = sum(x$Enh & x$Hit),
                        Total_Gene_TRUE = sum(x$Gene),
                        Total_HitGene_TRUE = sum(x$Gene & x$Hit),
                        Total_Both_TRUE = sum(x$Concord),
                        Total_BothHit_TRUE = sum(x$Concord & x$Hit),
                        p = f$p.value,
                        OR = f$estimate,
                        Lower = f$conf.int[1],
                        Upper = f$conf.int[2],
                        row.names = label)
      
      # return
      if (rbind) {
        return(rbind(rbind.to, out))
      } else {
        return(out)
      }
      
      
    }
     
  
  


################################################################################################################################ #
## Comparison to ABC ----
  
# ## Read inABC
#   abc <- read.delim("/mnt/Data1/PROJECTS/Farbod_ABC/ABC_results/predictions/EnhancerPredictions_NHAPeaks.bed", sep="\t")
#   colnames(abc) <- c("NHA_chr", "NHA_start", "NHA_end", "NHA_id", "NHA_score", "NHA_strand",
#                   "ABC_chr", "ABC_start", "ABC_end", "ABC_id", "Gene", "Score", "CellType", "ABC_Score")
#   
#   abc$ABC_id <- gsub("intergenic|", "", abc$ABC_id, fixed = TRUE)
#   abc$NHA_id <- paste0(abc$NHA_chr, ":", abc$NHA_start, "-", abc$NHA_end)
#   abc <- abc[which(abc$NHA_id %in% x$Enh.Pos), ]
# 
# ## Evaluate ABC pairs
#   # get enhancer and pair annotation for ABC
#   m <- match(abc$NHA_id, guides$TargetCoord)
#   abc$Enh <- guides$TargetID[m]
#   abc$Pair <- paste0(abc$Enh, "_", abc$Gene)  
#   
#   # match
#   abc <- abc[which(abc$Pair %in% x$Pair), c("Enh", "Gene", "Pair", "NHA_id", "ABC_id", "ABC_Score")]
# 
#   # save
#   write.csv(abc, file = "ABC Raw.csv")
#   
# ## Plot
#   m <- match(y$Pair, abc$Pair)
#   y$ABC.Score <- abc$ABC_Score[m]
#   y$ABC.Hit <- y$Pair %in% abc$Pair
#   
#   p <- table(y$HitCategory, y$ABC.Hit)
#   p <- p / rowSums(p)
#   p <- data.frame(ABC.Rate = p[,2])
#   
#   p$Cat <- factor(rownames(p), levels = levels(y$HitCategory))
#   sig.colours2 <- c(pal_lancet()(2), "grey80")
#   
#   pdf(file = "ABC Replication Rate.pdf", height = 3, width = 5)
#   ggplot(p, aes(x = Cat, y = ABC.Rate, fill = Cat)) +
#     geom_col(colour = "black") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = invis, legend.position = "none") +
#     scale_fill_manual(values = sig.colours2) +
#     labs(y = "Fraction of Pairs Called by ABC") +
#     scale_y_continuous(limits = c(0,1), expand = c(0,0))
#   dev.off()  
#   
    
################################################################################################################################ #
## Comparison to K562 screens ----
  
  
## Preprocess input files
  ## Gasperini
    gasp <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx", sheet = 2) # data from at-scale 
    gasp <- gasp[grep("candidate_enhancer", gasp$Category),] # filter to enhancers
    gasp <- gasp[,c(3,4,5,2)] # coordinates (hg19) and id
    gasp <- unique(gasp) # duplicate rows due to different guides are removed
    
    write.bed(gasp, "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/Candidates_AtScale_hg19.bed")
    
    gasp.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/Candidates_AtScale_hg38.bed"
    lift.bed(bed.in.19 = "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/Candidates_AtScale_hg19.bed", 
             bed.out.38 = gasp.in,
             isList = FALSE) # lift to hg38
    
  ## Xie
    xie <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Xie2019/mmc2.xlsx", sheet = 1) 
    x <- sub(":", "-", xie$`region pos (hg38)`)
    xie <- data.frame(chr = splitter(x, "-", 1),
                      start = splitter(x, "-", 2),
                      end = splitter(x, "-", 3), 
                      id = splitter(xie$`Position (hg19)`, "\\|", 1) %>% gsub(">", "", .))    
    xie <- unique(xie)
    xie.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Xie2019/Candidates_hg38.bed"
    write.bed(xie, xie.in)    

## Run intersection
  candidates <- "/mnt/Data0/PROJECTS/CROPSeq/TransferFromRNA/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  k562.out <- "/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/BED_Annotations_K562.bed"
  k562.in <- paste(gasp.in, xie.in, collapse = " ")
  
  # run  
  call <- paste("intersectBed",
                "-a", candidates,
                "-b", k562.in,
                "-wa",
                "-wb", # write the original entry in B for each overlap. 
                "-names gasp xie",
                ">", k562.out)

  system(call, intern = FALSE, wait = TRUE) 
  
## Read in and wrangle
  k562 <- read.delim(k562.out, sep = "\t", header = FALSE)
  colnames(k562) <- c("chr", "start", "end", "id", "study", "k562.chr", "k562.start", "k562.end", "k562.id", "V10", "V11")
  
## Annotate with the target gene in our Astrocytes
  k562$id <- paste0(k562$chr, ":", k562$start, "-", k562$end)
  k562$HitVoineagu <- k562$id %in% enh.hits  
  
## Compare hits
  x <- res[which(res$Hit),]  
  y <- k562[which(k562$HitVoineagu),]
  
  x <- x[which(x$Enh.Pos %in% y$id),]
  
  m <- match(x$Enh.Pos, y$id)
  x$GaspEnh <- y$k562.id[m]
  
  # all are from Gasperini, so find those hits!
  gasp.hits <- read_xlsx("/mnt/Data0/PROJECTS/CROPSeq/PublicData/K562Screens/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx", sheet = 3) # data from at-scale 
  gasp.hits <- gasp.hits[which(gasp.hits$Target_Site %in% x$GaspEnh),] # all are unique, can use a simpler match
  m <- match(x$GaspEnh, gasp.hits$Target_Site)
  x$GaspGene <- gasp.hits$target_gene_short[m]

  # remove rows and columns for clarity
  x <- x[-which(is.na(x$GaspGene)),c(1:4, 24:25)]  
  
  # save
  write.csv(x, file = "K562 - Overlapping Peaks With Target Genes.csv")
  
################################################################################################################################ #
## Comparison to Dong et al Nat Genet 2022 ----
  
## Used lasso regression in >1000 samples to link eRNA to mRNA, 500kb window
  
## Process input
  dong <- readxl::read_xlsx("../../../PublicData/Dong2022_LassoEnhPredictions/41588_2022_1170_MOESM3_ESM.xlsx", sheet = 5)
  dong.bed <- gsub("glia_", "", dong$enhancer) %>% strsplit(., "_")
  dong.bed <- data.frame(chr = sapply(dong.bed, "[", 1),
                         start = sapply(dong.bed, "[", 2),
                         end = sapply(dong.bed, "[", 3),
                         id = dong$enhancer_name)  
  
  dong.in <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Dong2022_LassoEnhPredictions/ST5.bed"
  dong.out <- "BED_Annotations_Dong.bed"
  write.bed(dong.bed, dong.in)  

## Intersect
  call <- paste("intersectBed",
                "-a", candidates,
                "-b", dong.in,
                "-wa",
                "-wb", # write the original entry in B for each overlap.
                ">", dong.out)

  system(call, intern = FALSE, wait = TRUE) 
  
## Read
  x <- read.delim(dong.out, header = FALSE)
  x <- unique(x)
  colnames(x) <- c("chr", "start", "end", "id", "dong.chr", "dong.start", "dong.end", "dong.id")
  x$id <- sub("_", ":", x$id)
  x$id <- sub("_", "-", x$id)

## Match
  y <- list()
  
  # for every enhancer
  for (j in 1:nrow(x)) {
    if (!(x$dong.id %in% res$Enh.Pos)) next
    
    # screen results
    z <- res[which(res$Enh.Pos == x$id[j]), c("Enh", "Gene", "Z", "FDR", "Hit")]
    
    # Dong results
    d <- dong[which(dong$enhancer_name == x$dong.id[j]),]
    
    # append
    m <- match(z$Gene, d$gene_name)
    z$Dong.ID <- unique(d$enhancer)
    z$Dong.Hit <- z$Gene %in% d$gene_name
    z$Dong.Coef <- d$coef[m]
   
    # account for cases where Dong identifies a hit but your screen failed to test it
    if (any(!(d$gene_name %in% z$Gene))) {
      w <- setdiff(d$gene_name, z$Gene)
      m <- match(w, d$gene_name)
      m <- data.frame(Enh = unique(z$Enh),
                 Gene = w,
                 Z = NA,
                 FDR = NA,
                 Hit = FALSE,
                 Dong.Hit = TRUE,
                 Dong.Coef = d$coef[m],
                 Dong.ID = unique(d$enhancer_name))
      z <- rbind(z, m)
    }
    
    # fin
    y[[j]] <- z
  }
  
  y <- do.call("rbind", y)
  
  # get all genes within 500kb
  
  # tested in screen
  
  # hit by screen?
  
  # hit by Dong?
  
  
################################################################################################################################ #
## Comparison to Yao 2022 ----
  
  
## Intersect
  dir.yao.known <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Yao_NatBiotech2022/ST2_KnownEnhancers.bed"
  dir.yao.non <- "/mnt/Data0/PROJECTS/CROPSeq/PublicData/Yao_NatBiotech2022/ST3_NonEnhancers.bed"
  dir.nha <- "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"
  # dir.nha.hg19 <- "/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Final_List/NHA_Peaks"

  # known enhancers
   call <- paste("intersectBed",
                "-a", dir.nha,
                "-b", dir.yao.known,
                "-wb", 
                ">", "Yao_Intersect_Known.bed")

  system(call, intern = FALSE, wait = TRUE)
  
  # non-enhancers
  call <- paste("intersectBed",
                "-a", dir.nha,
                "-b", dir.yao.non,
                "-wb", 
                ">", "Yao_Intersect_Non.bed")

  system(call, intern = FALSE, wait = TRUE)
  
  # note that the above has no window to either side; this consistent with other analyses in this study
  
## Read in
  yao_known <- read.delim("Yao_Intersect_Known.bed", header = FALSE)
  yao_non <- read.delim("Yao_Intersect_Non.bed", header = FALSE)
  
## Annotate
  # results from our study
  yao_annotated <- data.frame(Enh = unique(guides$TargetID[which(guides$TargetCat == "Enh")]),
                        Coord = unique(guides$TargetCoord[which(guides$TargetCat == "Enh")]))
  
  yao_annotated$TestPerformed <- yao_annotated$Enh %in% res.final$Enh
  yao_annotated$HitCore <- yao_annotated$Enh %in% res.final$Enh[which(res.final$HitCore)]
  yao_annotated$HitPerm <- yao_annotated$Enh %in% res.final$Enh[which(res.final$HitPermissive)]

  yao_annotated$HitCat <- "ns"
  yao_annotated$HitCat[yao_annotated$HitCore] <- "Core Hit"
  yao_annotated$HitCat[!(yao_annotated$HitCore) & yao_annotated$HitPerm] <- "Permissive Hit"
  yao_annotated$HitCat <- factor(yao_annotated$HitCat, levels = c("Core Hit", "Permissive Hit", "ns"))
  
  # results from Yao, using matched enh
    # I note that some of our enh (20) have multiple intersects to the known enhancer set, but the code will not note this in the output
  yao_annotated$Known <- yao_annotated$Coord %in% yao_known$V4
  yao_annotated$Non <- yao_annotated$Coord %in% yao_non$V4
  yao_annotated$Yao <- "No Intersection"
  yao_annotated$Yao[which(yao_annotated$Known)] <- "Yao Known Enh"
  yao_annotated$Yao[which(yao_annotated$Non)] <- "Yao Non Enh"
  
## Save 
  write.csv(yao_annotated, file = "Yao 2022 - Annotation of Screen Enh.csv")
  
## Stats
    yao_annotated_filt <- yao_annotated[which(yao_annotated$TestPerformed),] 
  
    x <- as.data.frame(table(yao_annotated_filt$HitCat, yao_annotated_filt$Known)) %>% dcast(Var1~Var2)
    colnames(x) <- c("Hit", "Yao_Known_False", "Yao_Known_True")
    
    y <- as.data.frame(table(yao_annotated_filt$HitCat, yao_annotated_filt$Non)) %>% dcast(Var1~Var2)
    colnames(y) <- c("Hit", "Yao_Non_False", "Yao_Non_True")
    
    z <- cbind(x, y[,-1])
    
    write.csv(z, file = "Yao 2022 - Overlap Counts.csv", row.names = FALSE)
    
    # output fisher test statistics
    write.fisher <- function(tab, name) {
    f <- fisher.test(tab)
    data.frame(OR = round(f$estimate, 2), p = f$p.value, row.names = name)
    }
    
    x <- rbind(table(yao_annotated_filt$HitPerm, yao_annotated_filt$Known) %>% write.fisher(name = "Permissive_With_Known"),
               table(yao_annotated_filt$HitCore, yao_annotated_filt$Known) %>% write.fisher(name = "Core_With_Known"),
               table(yao_annotated_filt$HitPerm, yao_annotated_filt$Non) %>% write.fisher(name = "Permissive_With_Non"),
               table(yao_annotated_filt$HitCore, yao_annotated_filt$Non) %>% write.fisher(name = "Core_With_Non"))
    write.csv(x, file = "Yao 2022 - Fisher.csv")
    
## Plot
  x <- table(yao_annotated_filt$HitCat, yao_annotated_filt$Yao)
  x <- x / rowSums(x) * 100
  x <- as.data.frame(x)
  
  pdf(file = "Yao 2022 - Barplot 1.pdf", height = 3, width = 3)
  ggplot(x, aes(x = Var2, y = Freq, colour = Var1, fill = Var1)) +
    geom_col(position = "dodge") +
    scale_fill_lancet() +
    scale_colour_lancet() +
    NoExpand +
    theme_bw() +
    theme_gjs +
    theme(axis.title.x = invis) +
    labs(y = "Fraction Within Enhancer Set", x = "Intersection With Enhancer Annotation in Yao 2022")
  dev.off()
  
    # version 2, with hit permissive only
    x <- table(yao_annotated_filt$HitPerm, yao_annotated_filt$Yao)
    x <- x / rowSums(x) * 100
    x <- as.data.frame(x)
    
    pdf(file = "Yao 2022 - Barplot 2.pdf", height = 3, width = 2.5)
    ggplot(x, aes(x = Var2, y = Freq, fill = Var1)) +
      geom_col(position = "dodge") +
      scale_fill_lancet() +
      NoExpand +
      theme_bw() +
      theme_gjs +
      guides(fill = guide_legend(title = "Screen Hit\n(Permissive)")) +
      theme(axis.title.x = invis) +
      labs(y = "Fraction Within Enhancer Set", x = "Intersection With Enhancer Annotation in Yao 2022")
    dev.off()
  
    
################################################################################################################################ #
## TADs ----
  
  
## Read in
  tad <- read.delim("../../../../FullLibrary_Selection/PublicData_forLibrarySelectionOnly/CulturedCells/Rajarajan2018/RajarajanScience2018_Synapse/Glia.100000_hg38.bed", header = FALSE)
  
    
## Some characteristics
  ## Size distribution
    p <- data.frame(Size = tad$V3 - tad$V2)
    p$Size <- p$Size / 10^6
    
    pdf(file = "TAD/Size Distribution.pdf", height = 3, width = 3)
    ggplot(p, aes(x = ".", y = Size)) +
      geom_quasirandom(alpha = 0.5) +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
      labs(y = "Size Distribution of 1131 TADs (MB)", x = "") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.text.x = invis, panel.grid.major.x = invis)
    dev.off()
  
## Do TADs overlap?
  ## Setup
    x <- tad
    x$ID <- paste0("TAD_", 1:nrow(x))
    write.bed(x, "TAD/SelfIntersection_In.bed")
    
  ## Intersect with self
    call <- paste("intersectBed",
                  "-a", "TAD/SelfIntersection_In.bed",
                  "-b", "TAD/SelfIntersection_In.bed",
                  "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                  ">", "TAD/SelfIntersection_Out.bed")
    system(call, intern = FALSE, wait = TRUE) 
    
  ## Read in
    self <- read.delim("TAD/SelfIntersection_Out.bed", header = FALSE) # the logic here is to intersect bed with itself
    colnames(self) <- c("A_chr", "A_start", "A_end", "A_ID", "B_chr", "B_start", "B_end", "B_ID", "Overlap")
  
    # remove duplicates
    self <- self[-which(self$A_ID == self$B_ID),] # where the a and b peaks are the same
    self <- self[-which(as.numeric(splitter(self$A_ID, "_", 2)) > as.numeric(splitter(self$B_ID, "_", 2))),] # where the a and b peak pair is the same as a b and a peak pair
  
  ## Calculate the fraction of overlap, where 1 means that a tad is entirely within another
    self$OverlapFraction <- apply(self[,-c(1,4,5,8)], 1, function (x) {
      sizeA <- x[2] - x[1]
      sizeB <- x[4] - x[3]
      frac <- x[5] / min(sizeA, sizeB)
      
    })
  
  ## Plot
    self$OverlapComplete <- factor(self$OverlapFraction == 1)
    levels(self$OverlapComplete) <- c("2 TAD Windows\nPartially Overlap", "Smaller TAD Wholly\nInside Larger TAD")
    pdf(file = "TAD/Self Overlap Distribution.pdf", height = 3, width = 3)
    ggplot(self, aes(x = OverlapComplete)) +
      geom_bar() +
      theme_bw() +
      scale_y_continuous(expand = c(0,0)) +
      labs(y = "Count of TAD Intersections") +
      theme(panel.border = invis, axis.line.y = element_line(), axis.title.x = invis)
    dev.off()
  
 
    
    
## Intersect EGPs with TADs
  ## Create a bed file whose coordinates are the intervening section between E and P
    x <- res.final[,c("Pair", "Enh.Pos", "Gene.TSS")]
    
    # get the chromosome
    x$chr <- splitter(x$Enh.Pos, ":", 1)
    
    # get the start, in this case the gene TSS
    x$Coord1 <- splitter(x$Gene.TSS, ":", 2) %>% as.numeric()
    
    # get the end, which is Enh's end or start (whichever is closer!)
    x$Coord2 <- apply(x, 1, function(x) {
      e.coord <- splitter(x[2], ":", 2)
      e.left <- splitter(e.coord, "-", 1) %>% as.numeric() 
      e.right <- splitter(e.coord, "-", 2) %>% as.numeric() 
      t.start <- as.numeric(x[5])
      use.left <- (abs(t.start - e.left)) < (abs(t.start - e.right))
      
      if (use.left) {
        return(e.left)
      } else {
        return(e.right)
      }
    })
    
    # the start coordinate is whichever is smaller of Coord1 and Coord2
    x$start <- apply(x[,c("Coord1", "Coord2")], 1, min)
    x$end <- apply(x[,c("Coord1", "Coord2")], 1, max) # and the opposite for the end
    
    # write
    x <- x[,c("chr", "start", "end", "Pair")]
    write.bed(x, "TAD/Pair_Intervening_Coord.bed")
    
    
  # here, reuse the file from the CTCF section, which contains the intervening section between E and P
  call <- paste("intersectBed",
                "-a", "TAD/Pair_Intervening_Coord.bed",
                "-b", "TAD/SelfIntersection_In.bed", # despite the name, this is simply TAD coordinates with an annotation
                "-wao", # Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
                ">", "TAD/Tad_Vs_Pair_Overlap.bed")
  
  system(call, intern = FALSE, wait = TRUE) 
  
  ## Read in
    tad_overlap <- read.delim("TAD/Tad_Vs_Pair_Overlap.bed", header = FALSE)
    colnames(tad_overlap) <- c("EGP_Chr", "EGP_Start", "EGP_End", "EGP", "TAD_Chr", "TAD_Start", "TAD_End", "TAD_ID", "Overlap_bp")
    
  ## How much of the EGP is contained within a tad?
    tad_overlap$EGP_Distance <- (tad_overlap$EGP_End - tad_overlap$EGP_Start)
    tad_overlap$Overlap_Frac <- tad_overlap$Overlap_bp / tad_overlap$EGP_Distance
    
  
  ## Now categorise each EGP
    # first, get all EGPs
    pair2tad <- res.final[,c("Pair", "Gene.Distance.Bin", "HitPermissive")]
    pair2tad$Gene.Distance.Bin <- factor(pair2tad$Gene.Distance.Bin, levels = c("2-10kb", "10-50kb", "50-100kb", "100-500kb"))
    
    # define basic true/false statements to determine categories
    within <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac == 1)] 
    cross <- tad_overlap$EGP[which(tad_overlap$Overlap_Frac < 1)] 
    
    # is the EGP fully contained within a single TAD?
    pair2tad$WithinTAD <- pair2tad$Pair %in% within
  
    # does the EGP cross a tad boundary?
    pair2tad$CrossTAD <- pair2tad$Pair %in% cross
    
    # combinations
    pair2tad$WithinAndCrossTAD <- pair2tad$WithinTAD & pair2tad$CrossTAD
    pair2tad$WithinNotCrossTAD <- pair2tad$Pair %in% setdiff(within, cross)
    pair2tad$CrossNotWithinTAD <- pair2tad$Pair %in% setdiff(cross, within)
    pair2tad$NoTAD <- rowSums(pair2tad[,4:8]) == 0 # fortunately 0, as was defined this way...
  
  ## Save
    write.csv(pair2tad, "TAD/Pair-TAD Annotation.csv")
    
    
    
################################################################################################################################ #
## Concordance: Are both members of the pair astrocyte markers? ----
          
  
    
################################################################################################################################ #
## Concordance: Are both members of the pair astrocyte markers? ----
  
  
## Add to annotation
  ## Genes
    egp_coannot$AstMarker_Gene_Herring2022 <- egp_coannot$Gene %in% isGene("AstMarkers_Herring2022")
    egp_coannot$AstMarker_Gene_Zhong2018 <- egp_coannot$Gene %in% isGene("AstMarkers_Foetal_Zhong2018")
    
  ## Enhancers
    egp_coannot$AstMarker_Enh_NottAtac <- egp_coannot$Enh %in% isEnh("AstSpecific_NottAtac")
    egp_coannot$AstMarker_Enh_NottH3K27ac <- egp_coannot$Enh %in% isEnh("AstSpecific_NottH3K27ac")
    egp_coannot$AstMarker_Enh_Herring2022 <- egp_coannot$Enh %in% isEnh("AstSpecific_Herring")
    
    
## Compare
  coannot_stats <- calc.coenrich("AstMarker_Enh_Herring2022", "AstMarker_Gene_Herring2022", label = "AstMarker_HerringVsHerring", rbind = FALSE)
  coannot_stats <- calc.coenrich("AstMarker_Enh_NottAtac", "AstMarker_Gene_Herring2022", label = "AstMarker_NottAtacVsHerring", rbind = TRUE)
  coannot_stats <- calc.coenrich("AstMarker_Enh_NottH3K27ac", "AstMarker_Gene_Herring2022", label = "AstMarker_NottH3K27acVsHerring", rbind = TRUE)
  
  coannot_stats <- calc.coenrich("AstMarker_Enh_Herring2022", "AstMarker_Gene_Zhong2018", label = "AstMarker_HerringVsZhong", rbind = TRUE)
  coannot_stats <- calc.coenrich("AstMarker_Enh_NottAtac", "AstMarker_Gene_Zhong2018", label = "AstMarker_NottAtacVsZhong", rbind = TRUE)
  coannot_stats <- calc.coenrich("AstMarker_Enh_NottH3K27ac", "AstMarker_Gene_Zhong2018", label = "AstMarker_NottH3K27acVsZhong", rbind = TRUE)
  

################################################################################################################################ #
## Concordance: Are both members of the involved in AD? ----
  
## Add to annotation
  ## Genes
    egp_coannot$AD_Gene_Saddick2022 <- egp_coannot$Gene %in% isGene("AstDisease_AD_Sadick2022")
    
  ## Enhancers
    egp_coannot$AD_Enh_Morabito2021 <- egp_coannot$Enh %in% isEnh("AD_Morabito2021")
    
    x <- apply(annot.enh[,grep("Bendl2022", colnames(annot.enh))], 1, any)
    x <- annot.enh$Enh[x]
    egp_coannot$AD_Enh_Bendl2022 <- egp_coannot$Enh %in% x
    
    
## Compare
  coannot_stats <- calc.coenrich("AD_Enh_Morabito2021", "AD_Gene_Saddick2022", label = "AD_MorabitoVsSaddick", rbind = TRUE)
  coannot_stats <- calc.coenrich("AD_Enh_Bendl2022", "AD_Gene_Saddick2022", label = "AD_BendlVsSaddick", rbind = TRUE)
  
  
################################################################################################################################ #
## Concordance: Are both members of the involved in GBM? ----
  
## Add to annotation
  ## Genes
    egp_coannot$GBM_Gene_Zhang2016 <- egp_coannot$Gene %in% isGene("AstDisease_GBM_Zhang2016")
    
  ## Enhancers
    x <- apply(annot.enh[,grep("Xu2021", colnames(annot.enh))], 1, any)
    x <- annot.enh$Enh[x]
    egp_coannot$GBM_Enh_Xu2021 <- egp_coannot$Enh %in% x
    
    
## Compare
  coannot_stats <- calc.coenrich("GBM_Enh_Xu2021", "GBM_Gene_Zhang2016", label = "GBM_XuVsZhang", rbind = TRUE)
  
  
  
################################################################################################################################ #
## Concordance: With Disgenet ----
  
  
## Once-off authorisation
  library(disgenet2r)
  disgenet_api_key <- get_disgenet_api_key(
    email = "i.voineagu@unsw.edu.au", 
    password = "UNSWlab123" )
  Sys.setenv(DISGENET_API_KEY = disgenet_api_key)
  
  
## Get pre-calculated enrichments
  x <- read.csv("../Genes/Disgenet - Enrichment.csv")
  
## Get genes!
  
  dg_genes <- list()
  dg_genes$Glioma <- disease2gene("C0017638", score = c(0, 1))@qresult$gene_symbol
  dg_genes$AD <- disease2gene("C0002395", score = c(0, 1))@qresult$gene_symbol 
  annot.gene$Disg_Glioma <- 
  
  
  
################################################################################################################################ #
## Save ----
  
write.csv(egp_replicate, file = "Final - Pair Replication.csv")    
write.csv(egp_coannot, file = "Final - Pair Coannotation.csv")
write.csv(coannot_stats, file = "Final - Pair Coannotation Stats.csv")
  