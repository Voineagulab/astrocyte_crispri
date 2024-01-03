#!/bin/bash
#Overlaps our peaks with ENCODE rE2G predicitons
cd /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/
ENHANCERS="Data/PeakData/Tested_Enh.bed"
for BED in /mnt/Data0/PROJECTS/CROPSeq/PublicData/ENCODE_rE2G/*/*.all.bed
do 
  name=`echo $BED | sed -E 's/.*\\/(.*).all.bed/\\1/g'`
  bedtools intersect -wa -wb -a $ENHANCERS -b $BED > "Results/ENCODE_rE2GPredictions/$name.bed"
done 