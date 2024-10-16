setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable4_GeneFunctionalAnnot//")
options(stringsAsFactors = FALSE)
source("../../../FullScale/Scripts/Functions.R")
res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")
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
## 4A: GO
  ## Read in results from clusterprofiler
    go <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Wellpowered/GO - Clusterprofiler.csv", row.names = 1)
  
  ## Columns
    colnames(go) <- c("Ontology", "GO_ID", "GO_Description", "HitRatio", "NonhitRatio", "P", "PAdj", "Q", "HitGeneSymbols", "Count")
    go <- go[,c("Ontology", "GO_ID", "GO_Description", "HitRatio", "NonhitRatio", "P", "PAdj", "HitGeneSymbols")]
    
  ## Prevent date formatting in excel
    go$HitRatio <- gsub("/", " of ", go$HitRatio)
    go$NonhitRatio <- gsub("/", " of ", go$NonhitRatio)
    
  ## Save (no filter)
    write.csv(go, file = "4A_GO_wellpowered.csv", row.names = FALSE)


## 4B: Gene functional annotations
    load("../../../FullScale/Results/3_HitEnrichment/Genes/Wellpowered/Final.rda", verbose = TRUE)
    
  ## Signed gene annotations
    # get downregulations
    dn <- annot.logi.down[,-c(1,2)]  
    dn <- -apply(dn, 2, as.numeric) # the - before apply converts TRUE to -1
  
    # get upregulations
    up <- annot.logi.up[,-c(1,2)]  
    up <- apply(up, 2, as.numeric) 
    
    # combined the above
    signed <- up + dn
    
    # rename columns
    colnames(signed) <- gsub("_TIC_hiPSC_", "_Leng2022_Protocol", colnames(signed))
  
    # get unsigned (but, in essence, these are all up)
    unsigned <- annot.logi[,-c(1,2)]
    unsigned <- unsigned[,-which(colnames(unsigned) %in% colnames(dn))] 
    unsigned <- apply(unsigned, 2, as.numeric) # code a TRUE as +1
    
    # get the "combined" annotations, in which like resources are summed
    comb <- annot.combined[,c("AstMarker", "AstAgeing", "AstMaturation", "Housekeeping", "AstActivation_HumanOnly")]
    colnames(comb) <- gsub("_HumanOnly", "", colnames(comb)) %>%
      paste0(., "_CombinedStudies") %>%
      gsub("Marker", "Markers", .)
      
      
    
    # combine all
    out <- data.frame(Gene = annot.logi$Gene,
                      Hit = annot.logi$Hit,
                      signed, 
                      unsigned,
                      comb)
    
    # add enhancer
    out$LinkedEnh <- NA
    h <- res.final[which(res.final$HitPermissive),]
    for (j in which(out$Hit)) {
      h1 <- h[which(h$Gene == out$Gene[j]),]
      out$LinkedEnh[j] <- paste(h1$Enh, collapse = "/")
    }
    # out <- relocate(out, c("Gene", "Hit", "Enh"))
    
    # filter columns
    keepCols <- sort(colnames(out))
    keepCols <- c("Gene", "LinkedEnh", "Hit", # annotation
                  keepCols[grep("Maturation|Ageing|Marker|Housekeeping", keepCols)], # from the four functional categories
                  keepCols[grep("Leng$|TCW$|Krencik$|Li$|Activation_CombinedStudies", keepCols)]) # subset of activation
    
    out <- out[,keepCols]
    
    # save as table
    write.csv(out, "4B_SignedAnnotations_Wellpowered.csv", row.names = FALSE)
  
## 4C: Enrichments
  load("../../../FullScale/Results/3_HitEnrichment/Genes/Wellpowered/Final.rda", verbose = TRUE)

  ## Combine the two enrichment dataframes
    # relabel resources
    enrichments.combined$Resource <- paste0(enrichments.combined$Resource, "_CombinedStudies") %>% 
      gsub("_HumanOnly", "", .)  %>% 
      gsub("Marker", "Markers", .)
    
    # combine
    enr <- rbind(enrichments, enrichments.combined)
  
  ## Filter rows
    enr$Resource <- gsub("_TIC_hiPSC_", "_Leng2022_Protocol", enr$Resource)
    enr <- enr[which(enr$Resource %in% colnames(out)),]
    
  ## Fix fractions
    enr <- update_enr_fractions(enr)

  # Save
    write.csv(enr, "4C_SignedEnrichments_wellpowered.csv", row.names = FALSE)
    