#!/bin/bash
set -e
source /home/rna2/PROGRAMS/miniconda3/etc/profile.d/conda.sh
conda activate mosdepth
cd  /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Process_Chip

#get coverage of Astrocyte enhancers using BigWigs from ENCODE Astrocyte datasets
for file in /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/AstroctyeDatasets/ENCFF*
do
  out=`echo $file | sed -E 's/.*\\/(.*).bigWig/\\1.tab/g'`
  bigWigAverageOverBed $file /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Tested_Enh.bed $out
  bigWigAverageOverBed $file /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Intergenic_Peaks_chr.bed Intergenic_$out
done

#get coverage of K562 enhancers from the ENCODE benchmarking dataset (Gschwind et al.) using BigWigs from ENCODE K562 datasets 
#These are used for the TAPseq-rf model.
ENCODE="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/ENCODE_Candidates.bed"
ENCODE_CHR="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/ENCODE_Candidates_chr.bed"

mkdir -p ENCODE

for file in /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/RF_data/hg38/*.bigWig #TAPseq/
do
  out=`echo $file | sed -E 's/.*\\/(.*).bigWig/\\1.tab/g'`
  bigWigAverageOverBed $file $ENCODE_CHR ENCODE/Encode_$out
done


#mosdepth of DNAse to replace ATACseq.pileup 
mosdepth --by $ENCODE_CHR Enhancer_DNAse_Encode /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/DNase/ENCFF205FNC.bam