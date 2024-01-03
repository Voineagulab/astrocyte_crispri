#This is for running TTseq analysis on Katana
#hg38
rsync /home/rna2/REFERENCE/HUMAN/GRCh38_hg38/cellranger_refdata-gex-GRCh38-2020-A/genes/genes.gtf z5205171@kdm.restech.unsw.edu.au://srv/scratch/voineagu/PROJECTS/EnhancerPrediction/Data/HG38
rsync /home/rna2/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/Genome/UCSC/hg38.fa z5205171@kdm.restech.unsw.edu.au://srv/scratch/voineagu/PROJECTS/EnhancerPrediction/Data/HG38

#hg19
#rsync /home/rna2/REFERENCE/HUMAN/GRCh37_hg19/GRCh37/Sequence/WholeGenomeFasta/genome.fa z5205171@kdm.restech.unsw.edu.au:/srv/scratch/voineagu/PROJECTS/EnhancerPrediction/Data/hg19/
#rsync  /home/rna2/REFERENCE/HUMAN/GRCh37_hg19/GRCh37/Annotation/Genes/genes.gtf z5205171@kdm.restech.unsw.edu.au:/srv/scratch/voineagu/PROJECTS/EnhancerPrediction/Data/hg19/

#Rsync bam files 

rsync z5205171@kdm.restech.unsw.edu.au:/srv/scratch/voineagu/PROJECTS/EnhancerPrediction/TTseq/Results/star_salmon/RNAseq_R1.markdup.sorted.bam /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/TTseq/
rsync z5205171@kdm.restech.unsw.edu.au:/srv/scratch/voineagu/PROJECTS/EnhancerPrediction/TTseq/Results/star_salmon/TTseq_R1.markdup.sorted.bam /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/TTseq/

#Rsync Feature counts from salmon 
rsync z5205171@kdm.restech.unsw.edu.au:/srv/scratch/voineagu/PROJECTS/EnhancerPrediction/TTseq/Results/star_salmon/salmon.merged.gene_counts.rds /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/TTseq/
rsync z5205171@kdm.restech.unsw.edu.au:/srv/scratch/voineagu/PROJECTS/EnhancerPrediction/TTseq/Results/star_salmon/salmon.merged.gene_counts.tsv /mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Data/K562Data/TTseq/
