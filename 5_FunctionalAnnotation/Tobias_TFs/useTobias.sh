#!/bin/bash

# This script runs the TOBIAS pipeline for footprinting analysis on ATAC-seq data:  
# 1) ATACorrect: Performs bias correction of ATAC-seq reads.  
# 2) ScoreBigwig: Calculates footprint scores from corrected cut sites.  
# 3) BINDetect: Estimates differentially bound motifs based on scores, sequence, and motif.  
# The script further processes TOBIAS outputs.

#conda install tobias -c bioconda
#conda env create -n tobias python=3.6.0
#pip install tobias
#pip install pyparsing==2.4.7
#pip install certifi
#Original Genome is on RNA1 GENOME=/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/GRCh38_hg38 
#. tobiasHeader.sh
ATACedits=""

#TOBIAS drops CHRM in the previous step but doesn't accept peaks with CHRm anymore ...
if [ ! $# -eq 1  ]
then
    echo "Specify Header file for example see tobiasHeader.sh"
    exit
else 
    . $1
fi

source /home/rna2/PROGRAMS/miniconda3/etc/profile.d/conda.sh
conda activate tobias

#TOBIAS drops CHRM in the first step but doesn't accept peaks with CHRm anymore ...
head $PEAKS
egrep -v "chrM" $PEAKS > $PEAKS.nochrM
NOCHRMPEAKS="$PEAKS.nochrM"


echo $OUTDIR
if [[ !  -d "$OUTDIR"  ]]
then
    mkdir -p $OUTDIR
fi

#VERY SLOW ~5hrs 
echo "Running TOBIAS ATACorrect"
TOBIAS ATACorrect \
--bam $ATAC_BAM  \
--genome $GENOME \
--peaks $PEAKS \
$ATACedits \
--blacklist $BLACKLIST \
--cores 6 \
--outdir $OUTDIR | tee ${OUTDIR}/ATAC.log

ATAC_OUT=`sed -E 's/.*\/(.*).bam/\1_corrected.bw/g' <<< $ATAC_BAM`
echo "Running TOBIAS ScoreBigWig"
#~45 minutes
TOBIAS ScoreBigwig \
--signal ${OUTDIR}/${ATAC_OUT} \
--regions $PEAKS \
--cores 5 \
--output ${OUTDIR}/ATAC_footprints.bw | tee  ${OUTDIR}/Score.log


if [[ !  -d "${OUTDIR}/BINDetect"  ]]
then
    mkdir ${OUTDIR}/BINDetect
fi

echo "Runnning TOBIAS BINDetect"
TOBIAS BINDetect \
--motifs $MOTIFS \
--signals ${OUTDIR}/ATAC_footprints.bw \
--genome $GENOME \
--peaks $NOCHRMPEAKS \
--cores 3 \
--outdir ${OUTDIR}/BINDetect | tee ${OUTDIR}/Detect.log

#Combine all bed files
echo "" >  ${OUTDIR}/BINDetect/combined_all.bed
echo "" >  ${OUTDIR}/BINDetect/combined.bed
echo "" >  ${OUTDIR}/BINDetect/unbound_combined.bed
for dir in `ls -d ${OUTDIR}/BINDetect/*/`
do
    ALL=`echo $dir | sed -E 's/.*\/(.*)\//\1_all.bed/g'`
    cat ${dir}beds/$ALL >> ${OUTDIR}/BINDetect/combined_all.bed
    BOUND=`echo $dir | sed -E 's/.*\/(.*)\//\1_ATAC_footprints_bound.bed/g'`
    cat ${dir}beds/$BOUND >> ${OUTDIR}/BINDetect/combined.bed
    UNBOUND=`echo $dir | sed -E 's/.*\/(.*)\//\1_ATAC_footprints_unbound.bed/g'`
    cat ${dir}beds/$UNBOUND >> ${OUTDIR}/BINDetect/unbound_combined.bed
done

#sort the beds
bedtools sort -i ${OUTDIR}/BINDetect/combined.bed > ${OUTDIR}/BINDetect/combined.sorted.bed 
awk '{ print $1,$2,$3,$4 }' OFS='\t' ${OUTDIR}/BINDetect/combined.sorted.bed > ${OUTDIR}/BINDetect/combined.sorted.bed2
mv ${OUTDIR}/BINDetect/combined.sorted.bed2 ${OUTDIR}/BINDetect/combined.sorted.bed

bedtools sort -i ${OUTDIR}/BINDetect/unbound_combined.bed > ${OUTDIR}/BINDetect/unbound_combined.sorted.bed 
awk '{ print $1,$2,$3,$4 }' OFS='\t' ${OUTDIR}/BINDetect/unbound_combined.sorted.bed > ${OUTDIR}/BINDetect/unbound_combined.sorted.bed2
mv ${OUTDIR}/BINDetect/unbound_combined.sorted.bed2 ${OUTDIR}/BINDetect/unbound_combined.sorted.bed

#Get counts for all peaks (quicker than loading)
bedtools intersect -c -a $PEAKS -b ${OUTDIR}/BINDetect/combined.sorted.bed > ${OUTDIR}/BINDetect/count_overlaps.bed
bedtools intersect -c -a $PEAKS -b ${OUTDIR}/BINDetect/unbound_combined.sorted.bed > ${OUTDIR}/BINDetect/unbound_count_overlaps.bed


############
#TODO
##From here needs to be re-run for Other data sets
bedtools intersect -wa -wb -a  $TESTEDPEAKS -b ${OUTDIR}/BINDetect/combined.sorted.bed > ${OUTDIR}/BINDetect/bound_overlaps.bed
bedtools intersect -wa -wb -a  $TESTEDPEAKS -b ${OUTDIR}/BINDetect/unbound_combined.sorted.bed > ${OUTDIR}/BINDetect/unbound_overlaps.bed

#INTERGENIC=/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Intergenic_Peaks_chr.bed 
bedtools intersect -wa -wb -a  $INTERGENIC -b ${OUTDIR}/BINDetect/combined.sorted.bed > ${OUTDIR}/BINDetect/bound_overlaps_intergenic.bed
bedtools intersect -wa -wb -a  $INTERGENIC -b ${OUTDIR}/BINDetect/unbound_combined.sorted.bed > ${OUTDIR}/BINDetect/unbound_overlaps_intergenic.bed


#Get expressed TFs only 
#If we need this information again re-write in R
#bedtools intersect -c -a $PEAKS -b ${OUTDIR}/BINDetect/exp_bound_overlaps.bed > ${OUTDIR}/BINDetect/exp_bound_overlaps.bed
#bedtools intersect -c -a $PEAKS -b ${OUTDIR}/BINDetect/exp_unbound_overlaps.bed > ${OUTDIR}/BINDetect/exp_unbound_overlaps.bed
egrep -f ${EXP_REGEX} ${OUTDIR}/BINDetect/unbound_overlaps.bed > ${OUTDIR}/BINDetect/exp_unbound_overlaps.bed
egrep -f ${EXP_REGEX} ${OUTDIR}/BINDetect/bound_overlaps.bed > ${OUTDIR}/BINDetect/exp_bound_overlaps.bed

#intersect most highly expressed Gene 5' end with peaks and label as gene promoter, Then take Bound/Unbound from those peaks
# bedtools window -w 500  -b $GENETSS -a $PEAKS > ${OUTDIR}/Peak_TSS_overlaps.bed
# TSSPEAKs="${OUTDIR}/Peak_TSS_overlaps.bed"
# bedtools intersect -wa -wb -a $TSSPEAKs -b ${OUTDIR}/BINDetect/combined.sorted.bed > ${OUTDIR}/BINDetect/gene_bound_overlaps.bed
# bedtools intersect -wa -wb -a $TSSPEAKs -b ${OUTDIR}/BINDetect/unbound_combined.sorted.bed > ${OUTDIR}/BINDetect/gene_unbound_overlaps.bed
# egrep -f ${EXP_REGEX} ${OUTDIR}/BINDetect/gene_bound_overlaps.bed > ${OUTDIR}/BINDetect/exp_gene_bound_overlaps.bed
# 
# bedtools window -w 500  -b $GENETSSALL -a $PEAKS > ${OUTDIR}/Peak_TSS_overlaps_All.bed
# TSSPEAKsALL="${OUTDIR}/Peak_TSS_overlaps_All.bed"
# bedtools intersect -wa -wb -a $TSSPEAKsALL -b ${OUTDIR}/BINDetect/combined.sorted.bed > ${OUTDIR}/BINDetect/gene_bound_overlaps_all.bed
# bedtools intersect -wa -wb -a $TSSPEAKsALL -b ${OUTDIR}/BINDetect/unbound_combined.sorted.bed > ${OUTDIR}/BINDetect/gene_unbound_overlaps_all.bed

