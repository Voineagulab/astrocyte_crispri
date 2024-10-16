#!/bin/bash

#using the BAM they said they used for predictions instead
ATAC_BAM="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/DNase/ENCFF205FNC.bam"
GENOME="/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/Genome/UCSC/hg38.fa"
PEAKS="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/DNase/ENCFF185XRG.bed" #This is the hg38 version 
#GENETSS="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/Gene_TSS.bed" # hg19 but not using it

OUTDIR="../../Results/Tobias/ENCODEFootprint"
MOTIFS="/mnt/Data0/PROJECTS/CROPSeq/PublicData/TF_BindingSites/JASPAR2022_homo_sapiens_latest_release_pfms_jaspar.txt"


ATACedits="--read_shift 0 0 --k_flank 6"
EXP_REGEX="../../Results/Tobias/ExpressedTFs/K562_Expressed_TFs_regex.txt"
TESTEDPEAKS="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/ENCODE_rE2G/ENCODE_Candidates_chr.bed"

BLACKLIST="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/Blacklist/hg38-blacklist.v2.bed"
#BLACKLIST="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/Blacklist/hg19-blacklist.v2.bed"

echo Files set
