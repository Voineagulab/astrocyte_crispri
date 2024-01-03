#!/bin/bash
set -e
source /home/rna2/PROGRAMS/miniconda3/etc/profile.d/conda.sh
conda activate mosdepth
cd  /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Process_Chip

#get mean of BigWigs for data from enformer prediction datasets Astrocytes
for file in /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/AstroctyeDatasets/ENCFF*
do
  out=`echo $file | sed -E 's/.*\\/(.*).bigWig/\\1.tab/g'`
  bigWigAverageOverBed $file /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Tested_Enh.bed $out
  bigWigAverageOverBed $file /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Intergenic_Peaks_chr.bed Intergenic_$out
done

#I used to use this to see if H3K4me3 is enriched in our peaks compared to promoters and intergenic peaks 
# bigWigAverageOverBed /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/AstroctyeDatasets/ENCFF184NZS.bigWig \
# ../../Data/GeneData/Gene_Promoters.bed Promoters_ENCFF184NZS


ENCODE="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/ENCODE_Candidates.bed"
ENCODE_CHR="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/ENCODE_Candidates_chr.bed"
#get mean of bigWigs for ENCODE dataset 
mkdir -p ENCODE
#get depth for data from enformer prediction datasets K562 cells
for file in /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/RF_data/hg38/*.bigWig #TAPseq/
do
  out=`echo $file | sed -E 's/.*\\/(.*).bigWig/\\1.tab/g'`
  bigWigAverageOverBed $file $ENCODE_CHR ENCODE/Encode_$out
done


#mosdepth of DNAse to replace ATACseq.pileup in our data
mosdepth --by $ENCODE_CHR Enhancer_DNAse_Encode /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/DNase/ENCFF205FNC.bam