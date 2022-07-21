## Read in NHA_1 for with either NovaSeq1 (old), or NovaSeq1+NovaSeq2+hnPCR2
new <- as.Seurat("../Sequencing/Cat/NHA_1_ENRICHED/outs/filtered_feature_bc_matrix/", h5 = FALSE, sample.id = "New", min.cells = 3, min.features = 200)
old <- as.Seurat("../Sequencing/GIMR_GWCCG_200802_IRIVOI_10X/211215_A00152_0509_BHWJ5MDSX2_summary/NHA_1/summary/sample_feature_bc_matrix", h5 = FALSE, sample.id = "Old", min.cells = 3, min.features = 200)

## Number of genes and cells respectively
dim(new) # 26443 14558
dim(old) # 21477  9915

## Summary of cell depth
summary(new$nCount_RNA)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 1512    1913   29693   39173   64806  371445 
summary(old$nCount_RNA)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 500    8601   16980   18374   25674  114632 

## Barcode assignment
newps <- read.table("../Sequencing/Cat/NHA_1_ENRICHED/getproto.txt", header = TRUE)
oldps <- read.table("../../Results/1_GuideAssignment/NHA-1_hnPCR2.txt", header = TRUE)
ps.whitelist <- read.csv("../Whitelists/Protospacer Whitelist.csv")
ps.whitelist <- ps.whitelist[which(ps.whitelist$Celltype == "NHA"),]

    newps$cr_call <- newps$cell %in% colnames(new)
    m <- match(substr(newps$barcode, 2, 1000), substr(ps.whitelist$GuideSequence, 1, 15))
    newps$GuideID <- ps.whitelist$GuideID[m]  
    newps$TargetCat <- ps.whitelist$TargetCat[m]  
    newps$TargetID <- ps.whitelist$TargetID[m]  
    newps$TargetCt <- ps.whitelist$Celltype[m]  
    
    oldps$cr_call <- oldps$cell %in% colnames(new) # note: new, not old
    m <- match(substr(oldps$barcode, 2, 1000), substr(ps.whitelist$GuideSequence, 1, 15))
    oldps$GuideID <- ps.whitelist$GuideID[m]  
    oldps$TargetCat <- ps.whitelist$TargetCat[m]  
    oldps$TargetID <- ps.whitelist$TargetID[m]  
    oldps$TargetCt <- ps.whitelist$Celltype[m]  

    thresh <- 3
    oldpst <- oldps[which(oldps$umi_count >= thresh & !(is.na(oldps$GuideID))),] # 88929 assignments
    newpst <- newps[which(newps$umi_count >= thresh & !(is.na(newps$GuideID))),] # 87864 assignments
    
    
