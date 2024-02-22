#!/bin/bash
#Overlaps our astrocyte peaks with ENCODE rE2G predictions for astrocytes
# all and th referes to all EG pairs and thresholded EG pairs respectively
# NHA astrocytes(Lonza)is the astrocyte_Donor_ENCDO916IIE sample and was used for comparison with the Astrocyte CRISPRi data

#./astrocyte_Donor_ENCDO916IIE:
#ENCFF405VJJ.th.bed
#ENCFF440FMQ.all.bed

#./astrocyte_of_the_cerebellum_Donor_ENCDO227AAA:
#ENCFF203TSB.all.bed
#ENCFF402NLK.th.bed

#./astrocyte_of_the_hippocampus_Donor_ENCDO223AAA:
#ENCFF661LEX.all.bed
#ENCFF836VJL.th.bed

#./astrocyte_of_the_spinal_cord_Donor_ENCDO224AAA:
#ENCFF425PZQ.all.bed
#ENCFF852GRD.th.bed


cd /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/
ENHANCERS="Data/PeakData/Tested_Enh.bed"
for BED in /mnt/Data0/PROJECTS/CROPSeq/PublicData/ENCODE_rE2G/*/*.all.bed
do 
  name=`echo $BED | sed -E 's/.*\\/(.*).all.bed/\\1/g'`
  bedtools intersect -wa -wb -a $ENHANCERS -b $BED > "Results/ENCODE_rE2GPredictions/$name.bed"
done 