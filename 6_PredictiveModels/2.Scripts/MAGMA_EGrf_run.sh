#!/bin/bash
#This script runs MAGMA enrichment analyses in EGrf enhancers

##This part runs MAGMA annotation step and Gene level analyses. The MAGMA annotation step is executed using a 1-kb window around each region.
#The default SNP-wise mean model was used for gene level analyses
Rscript RunMAGMA.R -a -l "" -c 3 -d General_Background -w 1 -b /mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/ --gene-analysis --input="General_Background"

#This part run MAGMA Set-level enrichment analysis. 
Rscript RunMAGMA.R -s -b /mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/ -f "EGrf_distance_05RPKM.set"  -i "General_Background" -d "EGrf_Distance_05RPKM"

