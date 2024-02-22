setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable1_ATAC/")
options(stringsAsFactors = FALSE)
source("../../../FullScale/Scripts/Functions.R")
res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
guides <- read.csv("../../../FullScale/Data/Whitelists/NHA Enh Library.csv")


update_enr_fractions <- function(dat) {
  # the function to output enrichment statistics outputs a misleading total
    # the "Total_TRUE" and "fraction_Bg_TRUE" both give refer to the complete set of bg + hit, rather than just bg
    # this does not affect the calculation of statistics
    # this script adjusts these numbers
  
  # infer n
  total <- dat$Total_TRUE / dat$Fraction_Bg_TRUE 
  nHit <- dat$Total_Hit_TRUE / dat$Fraction_Hit_TRUE 
 
  # adjust 
  dat$Total_TRUE <- dat$Total_TRUE - dat$Total_Hit_TRUE
  dat$Fraction_Bg_TRUE <- dat$Total_TRUE  / (total - nHit)
  
  # colnames!
  colnames(dat)[2:5] <- c("Nonhit_True", "Nonhit_Fraction", "Hit_True", "Hit_Fraction")
  
  # round
  dat$Nonhit_Fraction <- round(dat$Nonhit_Fraction, 3)
  dat$Hit_Fraction <- round(dat$Hit_Fraction, 3)
  dat[,c(6:9)] <- apply(dat[,c(6:9)], 2, signif, 3)
  
  
  # output
  return(dat)
  
}

## 1A: all ATAC information
  ## Read in
    # atac <- read.csv("../../../TransferFromRNA/FullLibrary_Selection/Results/NHA/AllAnnotatedPeaks.csv")
    atac <- read.csv("../../../FullLibrary_Selection/Results/Peaks_Annotated.csv")

  ## Sort
    o <- gsub("NHA_ATAC_S3.filtered.BAM_peak_", "", atac$name)
    o <- as.numeric(o)
    o <- order(o)
    atac <- atac[o,]
  
    atac$name <- gsub("NHA_ATAC_S3.filtered.BAM_peak_", "Peak_", atac$name)
    
  ## Add whether it is screened as a candidate
    atac$Candidate <- (atac$GENCODE32 == FALSE & 
                         atac$TAD &
                         atac$nRep >= 2 & 
                         atac$PE & 
                         atac$high.cpm)
    
    # further reordering: put the screened enhancers at the top
    atac <- atac[order(atac$Candidate, decreasing = TRUE),]
    m <- match(atac$id, splitter(guides$TargetCoord, "-exp", 1))
    atac$name <- guides$TargetID[m]
  
  ## Finalise column names
    colnames(atac) <- c("Chromosome",
                      "Start",
                      "End",
                      "Peak ID",
                      "Absolute Summit",
                      "Pileup",
                      "P-value (-log10)",
                      "Fold Enrichment",
                      "Q-value (-log10)",
                      "EnhancerID",
                      "Genic",
                      "Fullard2018_Brain_NeuN-_ATAC",
                      "Rizzardi2019_Brain_NeuN-_ATAC",
                      "Song2019_Culture_Ast_ATAC",
                      "Trevino2020_Organoid_Glial_ATAC",
                      "Rajarajan2018_Culture_Ast_HiC",
                      "PsychENCODE_Enhancer",
                      "CPM",
                      "CPM > 2",
                      "ATAC_Replication_N",
                      "Candidate")
    
    atac <- atac[,c(1,2,3,4,10,
                6,9, # removes p-value and fold-enrichment and absolute summit
                18, 11, 17, 16, 12:15, 19, 20, 21)]
    
  ## Write
    write.csv(atac, "1A_AllPeakAnnotation.csv", row.names = FALSE)
  
## 1B: candidate annotation using public resources
  ## Read in
    load("../../../FullScale/Results/3_HitEnrichment/Chromatin/Final.rda", verbose = TRUE)
    cand <- candidate.annot
    
  ## Add genes
    cand$Gene <- ""
    h <- res.final[which(res.final$HitPermissive),]
    for (j in which(cand$Hit)) {
      h1 <- h[which(h$Enh == cand$Enh[j]),]
      cand$Gene[j] <- paste(h1$Gene, collapse = "/")
    }
    
  ## Filter columns
    cn <- colnames(cand)
    cand <- cand[,c("Enh", "Coord", "Tested", "Hit", "Gene", 
                    cn[grep("Altius", cn)], 
                    "ValidatedEnh_K562_Yao2022", 
                    "AstSpecific_NottAtac", "AstSpecific_NottH3K27ac", "AstSpecific_Morabito", 
                    "AstSpecific_Herring_Pooled", "AstSpecific_Glia_Pooled",
                    # cn[grep("GBM", cn)],
                    "Superenhancer")]
    
    colnames(cand) <- gsub("Superenhancer", "NHA_Superenhancer_Hnisz2013", colnames(cand)) %>%
      gsub("_Xu", "_Supenhancer_Xu", .) %>%
      gsub("Altius", "ENCODEv3_Meuleman2020", .) %>%
      gsub("Gene", "LinkedGene", .) %>%
      gsub("AstSpecific_Glia_Pooled", "GliaSpecific_Herring_Pooled", .) %>%
      gsub("Herring", "Herring2022", .) %>%
      gsub("Morabito", "Morabito2021", .) %>%
      gsub("Nott", "Nott2019_", .)
      
                    
  ## Save
    write.csv(cand, "1B_CandidatePeakAnnotation.csv", row.names = FALSE)
                    
## 1C: TTseq annotation
  ## Load 
    tt <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/TTseq/Results Table.csv", row.names = 1)
    tt <- tt[which(rownames(tt) %in% res.final$Enh),]
    
    transcribed <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/TTseq/Transcriptional classification.csv")
    
  ## Process columns
    tt$Enh <- rownames(tt)
    colnames(tt) <- gsub("RNAseq", "RiboDepRNAseq", colnames(tt))
    tt <- tt[,c("Enh", "Hit", "TTseq_Pos", "TTseq_Neg", "TTseq_Total", "TPM_TTseq_Total", "RiboDepRNAseq_Pos", "RiboDepRNAseq_Neg", "RiboDepRNAseq_Total", "TPM_RiboDepRNAseq_Total", "TT_Enrich")]
    
    colnames(tt) <- gsub("TT_Enrich", "NascentEnriched", colnames(tt)) %>%
      gsub("^TT", "Count_TT", .) %>%
      gsub("^Ribo", "Count_Ribo", .) %>%
      gsub("Pos", "PosStrand", .) %>%
      gsub("Neg", "NegStrand", .) %>%
      gsub("TPM", "CPM", .)
    
  ## Is it transcribed?
    m <- match(tt$Enh, transcribed$Enh)
    tt$TTseq_eRNA_Unidirectional <- transcribed$Unidirectional_3[m]
    tt$TTseq_eRNA_Bidirectional <- transcribed$Bidirectional_3[m]
    tt$TTseq_eRNA <- transcribed$Category_3[m]
   
  ## Add gene
    tt$LinkedGene <- ""
    h <- res.final[which(res.final$HitPermissive),]
    for (j in which(tt$Hit)) {
      h1 <- h[which(h$Enh == tt$Enh[j]),]
      tt$LinkedGene[j] <- paste(h1$Gene, collapse = "/")
    }

  ## Reorder columns
    tt <- relocate(tt, c("Enh", "Hit", "LinkedGene", "TTseq_eRNA", "TTseq_eRNA_Unidirectional", "TTseq_eRNA_Bidirectional", "NascentEnriched"))
    
  ## Reorder rows
    ord <- order(as.numeric(gsub("Enh", "", tt$Enh)))
    tt <- tt[ord,]
     
  ## Add fantom5
    f5 <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/FANTOM5/Peak Annotation.csv")
    m <- match(tt$Enh, f5$Enh)
    tt$FANTOM5_Ast_CAGE <- f5$UsedAst[m] > 0
    # tt$FANTOM5_GBM_CAGE <- f5$UsedGBM[m] > 0
    tt$FANTOM5_AnySample_CAGE <- f5$FANTOM5[m]
    
  ## Update the names of eRNA calls
    tt$TTseq_eRNA <- gsub("Not transcribed", "eRNA-", tt$TTseq_eRNA) %>%
      gsub("directional", "directional eRNA+", .) 
    
  ## Filter columns
    tt <- tt[,-c(5:6)]
    
  ## Save
    write.csv(tt, "1C_eRNA.csv", row.names = FALSE)
    
## 1D: Enrichment for annotations and transcription
  ## Load for annotations
    load("../../../FullScale/Results/3_HitEnrichment/Chromatin/Final.rda", verbose = TRUE)
    enr <- candidate.enrich
    enr$Resource <- rownames(enr)
    enr <- relocate(enr, "Resource")
    
    enr$Resource <- gsub("Superenhancer", "NHA_Superenhancer_Hnisz2013", enr$Resource) %>%
      gsub("_Xu", "_Supenhancer_Xu", .) %>%
      gsub("Altius", "ENCODEv3_Meuleman2020", .) %>%
      gsub("Gene", "LinkedGene", .) %>%
      gsub("AstSpecific_Glia_Pooled", "GliaSpecific_Herring_Pooled", .) %>%
      gsub("Herring", "Herring2022", .) %>%
      gsub("Morabito", "Morabito2021", .) %>%
      gsub("Nott", "Nott2019_", .)
    
    enr <- enr[which(enr$Resource %in% colnames(cand)),]
    
  ## Add transcription
    eRNA.enrich <- read.csv("../../../FullScale/Results/4_EnhancerTranscription/Hit Enrichments.csv")
    colnames(eRNA.enrich)[1] <- "Resource"
    
    # filter rows
    eRNA.enrich <- eRNA.enrich[c(1,2,4,5),]
    eRNA.enrich$Resource <- gsub("FANTOM5_Ast", "FANTOM5_Ast_CAGE", eRNA.enrich$Resource) %>%
      gsub("FANTOM5_Any", "FANTOM5_AnySample_CAGE", .) %>%
      gsub("TTseq", "TTseq_eRNA+", .) %>%
      gsub("Unidirectional", "Uni_or_Bi", .)
    
      
  ## Combined
    enr <- rbind(enr, eRNA.enrich)
    
  ## Fix fractions
    enr <- update_enr_fractions(enr)
    
  ## Save
    write.csv(enr, "1D_CandidatePeakEnrichments.csv", row.names = FALSE)
    