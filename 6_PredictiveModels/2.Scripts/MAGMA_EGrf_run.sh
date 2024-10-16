#Final Background w = 1,0
##This part Run MAGMA annotation step and Gene level analyses. The MAGMA annotation step  was executed using a 1-kb window around each region.
#The default SNP-wise mean model was used for gene level analyses
Rscript RunMAGMA.R -a -l "" -c 3 -d General_Background -w 1 -b /mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/ --gene-analysis --input="General_Background"

#This part Run MAGMA Set-level enrichment analyses. 
Rscript RunMAGMA.R -s -b /mnt/Data0/PROJECTS/CROPSeq/Manuscript/GitHub/Scripts/6.Predictive_models/MAGMA/EGrf/ -f "EGrf_distance_05RPKM.set"  -i "General_Background" -d "EGrf_Distance_05RPKM"

