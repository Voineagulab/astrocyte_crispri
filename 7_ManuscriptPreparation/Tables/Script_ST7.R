setwd("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/Manuscript/Tables/STable7_Disease/")
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
  if (any(is.na(nHit))) {
    nHit[which(is.na(nHit))] <- mean(nHit, na.rm = TRUE) # this triggers when no genes are TRUE in the annotation (just GBM_Zhang2016, Downregulated). Note that all values of nHit are the same, anyway
  }
  
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



## 7A: Disgenet
# read in  
disg <- read.csv("../../../FullScale/Results/3_HitEnrichment/Genes/Disgenet - Enrichment.csv", row.names = 1)

# filter
disg <- disg[which(disg$FDR < 0.05),]

# colnames
colnames(disg)[7:10] <- c("Nonhit_count", "Nonhit_fraction", "Hit_count", "Hit_fraction")

# adjust the non-hit count columns to be total - hits
disg$Nonhit_count <- disg$Nonhit_count - disg$Hit_count
disg$Nonhit_fraction <- disg$Nonhit_count / (length(unique(res.final$Gene)) - length(unique(res.final$Gene[which(res.final$HitPermissive)])))

# clean
disg$Disease_Type <- str_to_sentence(disg$Disease_Type)
disg$Disease <- str_to_sentence(disg$Disease)

# remove columns
disg <- disg[,-c(11, 12)] # upper and lower bounds of the confidence interval

# save
write.csv(disg, file = "7A_Disgenet.csv", row.names = FALSE)

## 7B: Disease-associated DEGs in astrocytes
## Load
load("../../../FullScale/Results/3_HitEnrichment/Genes/Final.rda", verbose = TRUE)

## Signed gene annotations
# get downregulations
dn <- annot.logi.down[,-c(1,2)]  
dn[,"AstDisease_MS_Jakel2019"] <- 0 # here, I note that it is NA prior to conversion. this is because the study only supplied DEGs that were upregulated in MS, not downregulated
dn <- -apply(dn, 2, as.numeric) # the - before apply converts TRUE to -1

# get upregulations
up <- annot.logi.up[,-c(1,2)]  
up <- apply(up, 2, as.numeric) 

# combined the above
signed <- up + dn
signed[,"AstDisease_MS_Jakel2019"] <- up[,"AstDisease_MS_Jakel2019"] 


# combine all
out <- data.frame(Gene = annot.logi$Gene,
                  Hit = annot.logi$Hit,
                  signed)

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
              keepCols[grep("Disease", keepCols)])
keepCols <- keepCols[-grep("FocalCortical|DAA|Velmeshev", keepCols)]

out <- out[,keepCols]

colnames(out) <- gsub("AstDisease_", "", colnames(out)) 

# save
write.csv(out, file = "7B_AstDEGs.csv", row.names = FALSE)

## 7C: enrichments for the above
## Load
load("../../../FullScale/Results/3_HitEnrichment/Genes/Final.rda", verbose = TRUE)
enr <- enrichments

## Filter rows
enr$Resource <- gsub("AstDisease_", "", enr$Resource) 
enr <- enr[which(enr$Resource %in% colnames(out)),]

## Update the fraction and nIntersect
enr <- update_enr_fractions(enr)

# Save
write.csv(enr, "7C_AstDEGsEnrichments.csv", row.names = FALSE)

## 7D: Rosmap annotations
## Load
load("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.Results.rda", verbose = TRUE)
rosmap <- lapply(rosmap.results, function(x) x$sig.genes)
rosmap <- do.call("rbind", rosmap)

## Clean rows
rosmap <- rosmap[which(rosmap$cluster_id == "Ast"),] # analyses performed across all astrocytes, rather than the individual sub-types

## Save
write.csv(rosmap[, -c(1:3])], file = "7D_RosmapAnnotation_v2.csv", row.names = FALSE)

## 7E: Rosmap enrichments
## Load
load("/mnt/Data0/PROJECTS/CROPSeq/IV/RESULTS/3.SigGeneCharact/ROSMAP.Results.rda", verbose = TRUE)
rosmap <- lapply(rosmap.results, function(x) x$fisher.tests)
rosmap <- do.call("rbind", rosmap)

## Clean rows
rosmap <- rosmap[which(rosmap$CT == "Ast"),] # analyses performed across all astrocytes, rather than the individual sub-types
rosmap <- rosmap[-which(is.na(rosmap$loc.Nintersect)),] # nas in this column indicate no test

# and sort
rosmap <- rosmap[order(rosmap$glb.fisher.padj),]

## Clean columns
rosmap$ROSMAP_Class <- splitter(rownames(rosmap), "\\.", 1) # add class, based on the list's level in rosmap.results

# filter columns
glb_or_loc <- "loc."
rosmap <- rosmap[,c("ROSMAP_Class", "Trait",  paste0(glb_or_loc, c("Nintersect", "fisher.p", "fisher.or", "fisher.padj")))]
colnames(rosmap) <- c("ROSMAP_Class", "Trait", "N_intersection", "P", "OR", "FDR")
# remove l=

## Save
write.csv(rosmap, file = "7E_RosmapEnrich.csv", row.names = FALSE)


