#!/bin/bash

ATAC_BAM="/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/Bulk_ATAC_RNAseq/GOK8505-GOK8837/GOK8837A1/Mapping/filtered_bam/NHA_ATAC_S3.filtered.bam"
GENOME="/home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/Genome/UCSC/hg38.fa"
PEAKS="/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/Bulk_ATAC_RNAseq/GOK8505-GOK8837/GOK8837A1/Mapping/MACS2/NHA_ATAC_S3.filtered.BAM_peaks.narrowPeak"
OUTDIR="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/Footprint"
#GENETSS="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/GeneData/All_transcript_TSSs.bed"
#GENETSSALL="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/GeneData/All_transcript_TSSs_allGenes.bed"

BLACKLIST="/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/Blacklist/hg38-blacklist.v2.bed"
MOTIFS="/mnt/Data0/PROJECTS/CROPSeq/PublicData/TF_BindingSites/JASPAR2022_homo_sapiens_latest_release_pfms_jaspar.txt"


EXP_REGEX="../Results/Tobias/ExpressedTFs/Expressed_TFs_regex.txt"

TESTEDPEAKS="/mnt/Data0/PROJECTS/CROPSeq/FullLibrary_Selection/Results/Final_List/NHA_Peaks.bed"

