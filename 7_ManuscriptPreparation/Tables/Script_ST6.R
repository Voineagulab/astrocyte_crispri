setwd("/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable6_Variants/")
options(stringsAsFactors = FALSE)
source("../../../FullScale/Scripts/Functions.R")
res.final <- read.csv("../../../FullScale/Results/2_DE/Enh/Results Final.csv")


## 6A: all functional enhancer to SNP annotations
  ## Read in
    snp2peak <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/Final - SNP-to-Peak List.csv")
    
  ## Just hits
    snp2peak <- snp2peak[which(snp2peak$Hit),]
    
  ## Column filtering
    snp2peak <- snp2peak[,-c(2,11)] # "Hit" and "in153"
    
  ## Clean the synonym column
    snp2peak$Synonym <- gsub("\\|$", "", snp2peak$Synonym) # removes a trailing "|"
    
  ## Clarify the coordinate column
    colnames(snp2peak)[4] <- "EnhCoord"
    snp2peak <- relocate(snp2peak, c("Enh", "EnhCoord"))
    
  ## Save
    write.csv(snp2peak, "6A_SNPlist.csv", row.names = FALSE)

## 6B: eQTL annotation
  ## Read in
    eSig <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs pooled.csv", row.names = 1)
    eFm  <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs pooled, Fine-mapped.csv", row.names = 1)
    eAst  <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/eQTL Final - EGPs pooled, Ast-specific.csv", row.names = 1)
    
  ## Update colnames
    # function
    resources <- c("Metabrain", "Bryois", "GTEx", "mmQTL")
    renameCols <- function(dat) {
      s <- which(colnames(dat) %in% resources)
      colnames(dat)[s] <- paste0(colnames(dat)[s], "_eQTLs")
      
      dat[,s] <- apply(dat[,s], 2, function(x) {
        x[which(is.na(x))] <- "No eQTL"
        x[which(x == "TRUE")] <- "Same gene"
        x[which(x == "FALSE")] <- "Different gene(s)"
        return(x)
      })  
      
      colnames(dat)[which(colnames(dat) == "eQTL_Category")] <- "Overall_category"
      return(dat)
    }
    
    # apply
    eSig <- renameCols(eSig)
    eFm <- renameCols(eFm)
    eAst <- renameCols(eAst)
      
  ## Save
    write.csv(eSig, file = "6B_eQTLs_significant.csv", row.names = FALSE)
    write.csv(eFm, file = "6B_eQTLs_finemapped.csv", row.names = FALSE)
    write.csv(eAst, file = "6B_eQTLs_astro.csv", row.names = FALSE)
    
    
## 6C: GWAS catalogue
  ## Read in
    gwascat <- read.csv("../../../FullScale/Results/3_HitEnrichment/Variants/Final - SNPs of Interest.csv")
    colnames(gwascat) <- c("Enhancer", "Hit", "SNP", "Resource", "SNP.Category", "Phenotype", "P")
    
  ## Filter rows
    gwascat <- gwascat[which(gwascat$Hit),] # hits only
    gwascat <- gwascat[which(gwascat$Resource == "GWAS_Catalogue"),] # GWAS_catalogue only
    
  ## Add relevant genes
    gwascat$LinkedGenes <- NA
    h <- res.final[which(res.final$HitPermissive),]
    for (j in 1:nrow(gwascat)) {
      h1 <- h[which(h$Enh == gwascat$Enhancer[j]),]
      gwascat$LinkedGenes[j] <- paste(h1$Gene, collapse = "/")
    }
    
  ## Filter columns
    gwascat <- gwascat[,c("Enhancer", "LinkedGenes", "SNP", "SNP.Category", "Phenotype", "P"),]
    
  ## Save
    write.csv(gwascat, "6C_GWASCatalogue.csv", row.names = FALSE)
    