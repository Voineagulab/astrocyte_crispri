#!/bin/bash
#This script extracts FASTA sequences for enhancer regions
#sed "s/^chr//" /mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Variants/DeepLearning/Webtool_Input/Beluga_2000bp.bed > /mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/3_HitEnrichment/Variants/DeepLearning/Webtool_Input/Beluga_2000bp_nochr.bed

OUT="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/"

ENH_BED_2000bp="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Enhancers_2kbwindow.bed"

FASTA="/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/Genome/UCSC/hg38.fa"
#FASTA="/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/Genome/Ensembl/genome.fa"
ENH_BED="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/PeakData/Tested_Enh.bed"
bedtools getfasta  -fi $FASTA -bed $ENH_BED -fo ${OUT}PeaksFasta.tab -name -tab 
bedtools getfasta  -fi $FASTA -bed $ENH_BED_2000bp -fo ${OUT}PeaksFasta_2000bp.tab -name -tab
bedtools getfasta  -fi $FASTA -bed $ENH_BED -fo ${OUT}PeaksFasta.fa -name
bedtools getfasta  -fi $FASTA -bed $ENH_BED_2000bp -fo ${OUT}PeaksFasta_2000bp.fa -name
