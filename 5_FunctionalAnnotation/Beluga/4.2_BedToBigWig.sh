#!/bin/bash
###### on RNA2

cd /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Beluga/BelugaVariants/

source /home/rna2/PROGRAMS/miniconda3/etc/profile.d/conda.sh
conda activate mosdepth
sizes=http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
for file in BelugaVariants BelugaVariants_Escore
do
  sort -k1,1 -k2,2n $file.bedGraph > $file.sorted.bedGraph 
  bedGraphToBigWig $file.sorted.bedGraph $sizes $file.bw
done