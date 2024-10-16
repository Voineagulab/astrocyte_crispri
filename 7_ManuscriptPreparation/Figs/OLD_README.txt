Figure production scripts (in R) are organised by theme, rather than Figure / Supplementary Figure number. Plots are output as .pdf, with corresponding statistics (where relevant) saved as .txt with no delimiter.

File names automatically note from which script they are generated, as listed below. Note that not all figures are produced in these scripts, such as those produced by Cytoscape or BioRender

############################ Figure1
####1B
Script:2_NHACharacterisation_CandidateEnhancerSelection/2f_RNAseq Clustering.R
Output:IV_RNAseqClustering.pdf

####1C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SingleR_annotation.R
Outputs: 
Barplot_SingleR.pdf
UMAP_SingleR.pdf

####1D
Script: LibraryDesign.R
Outputs:
Encode annotation of candidates.pdf
Encode annotation of candidates.txt

#####1E
Script: 7_ManuscriptPreparation/Figs/Manuscript_Final_UpsetPlot.R
Output: IV_Fig1_UpsetPlot_long.pdf

#####1F&G:
Script:5_FunctionalAnnotation/WGCNAenh.R
Outputs:
Fig1.ATAC.large.pdf (F)
Fig1.ModuleSig.pdf (G)
moduleStats_Hits.csv
moduleStats.csv

############################ Figure2
#####2A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Volcano.pdf

#####2B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SanityChecks.R
Outputs: 
K562, Superenhancer, TAD.pdf
K562.txt
Superenhancer.txt
TAD.txt

#####2C
Script: 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output: Nanostring vs scRNAseq.pdf


############################ Figure3

#####3A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs: 
Distance Density.pdf
Distance Density.txt

#####3B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Nearest gene stacked barplot.pdf

#####3D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Outputs: 
Reads in hits versus non-hits.pdf
Reads Reads in hits versus non-hits.txt

#####3F
Script:Script:5_FunctionalAnnotation/WGCNAenh.R
Output:
Fig3.hitsATAC_v2_large.pdf

#####3G
Script:Script:5_FunctionalAnnotation/WGCNAenh.R
Output:
CT_Stage_BarplotData.csv
############################ Figure4
#####4D*

############################ Figure5
#####6D

############################ Figure6
Script:5_FunctionalAnnotation/Beluga/5.ISM_Beluga_CCL2.R
Output:Enh427_ISM.pdf

*1E - Script LibraryDesign - Candidate length histogram.pdf
*1F - Script LibraryDesign - Genes tested per candidate.pdf
2A - Script MainDE - Volcano.pdf
2B - Script SanityCheck - K562, Superenhancer, TAD.pdf
2B - Script SanityCheck - K562.txt
2B - Script SanityCheck - Superenhancer.txt
2B - Script SanityCheck - TAD.txt
2C - Script MainDE - Nanostring vs scRNAseq.pdf
3A - Script MainDE - Distance Density.pdf
3A - Script MainDE - Distance Density.txt
3B - Script MainDE - Nearest gene stacked barplot.pdf
3D - Script tt - eRNA.pdf
3E - Script tt - Reads in hits versus non-hits.pdf
3E - Script tt - Reads in hits versus non-hits.txt
*3F
*3G
4A - Script HumBiol - Functional annotation.pdf
4C - Script TF - binding counts and fractions.pdf
4C - Script TF - Combined.pdf
4C - Script TF - footprinting signal.pdf
4D*
5A - Script HumBiol - Disease annotation.pdf
5C*
6A - Script eqtl - SNPs per enhancer.pdf
6A - Script eqtl - SNPs per enhancer.txt
6B - Script eqtl - eqtl vs crispri (pooled).pdf
6C - Script eqtl - eqtl reproducibility.pdf
6C - Script eqtl - eqtl reproducibility.txt
6D* 
7C - Script SanityCheck EGrf - NHA_EGRF_OverlapBarplot_0.5RPKM.pdf
*SFig1A - Script LibraryDesign - SingleR.pdf
*SFig1B - Script LibraryDesign - NHA ATAC vs Herring snATACseq.pdf
SFig2A - Script LibraryDesign - Candidate annotation upset.pdf
SFig2B - Script LibraryDesign - Candidate annotation barplot.pdf
SFig3A*
SFig3B*
SFig4A - Script qc - nCells per candidate.pdf
SFig4B - Script qc - nCells per guide.pdf
SFig4C - Script qc - MOI.pdf
SFig4D - Script qc - positive control volcano.pdf
SFig4E - Script qc - qqplot.pdf
SFig4F - Script LibraryDesign - dCas9KRAB.pdf
SFig4G - Script MainDE - nEnh Per Gene.pdf
SFig4G - Script MainDE - nGenes Per Enhancer.pdf
SFig4H - Script MainDE - Downregulation Percentage.pdf
SFig4H - Script MainDE - Downregulation percentage.txt
SFig4I - Script MainDE - Downregulation Distribution.pdf
SFig5B - Script tt - Scatterplot RNAseq vs TTseq.pdf
SFig5C - Script tt - eRNA across thresholds.pdf
SFig5D - Script tt - FANTOM5.pdf
SFig7A - Script TF - JASPAR class enrichments.pdf
SFig7B*
SFig8 - Script eqtl - eqtl vs crispri (per study).pdf
SFig10A - Script qc - cells per batch.pdf
SFig10B - Script qc - cell-level qc.pdf
SFig10 - Script qc - scRNAseq qc.txt
SFig11A - Script MainDE - Expression threshold density.pdf
SFig11B - Script MainDE - Expression threshold dropout rate.pdf
SFig11B - Script MainDE - Expression threshold dropout rate.txt
SFig11C - Script MainDE - Expression threshold Seurat vs CPM.pdf
SFig11C - Script MainDE - Expression threshold Seurat vs CPM.txt
SFig12A-C - Script MainDE - Nanostring QC.pdf
SFig12D - Script MainDE - Nanostring QC on genes.pdf