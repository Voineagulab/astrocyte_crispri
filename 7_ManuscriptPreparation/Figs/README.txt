Figure production scripts (in R) are organised by theme, rather than Figure / Supplementary Figure number. Plots are output as .pdf, with corresponding statistics (where relevant) saved as .txt with no delimiter.

File names automatically note from which script they are generated, as listed below. Note that not all figures are produced in these scripts, such as those produced by Cytoscape or BioRender

############################ Figure1
####1B
Script:2_NHACharacterisation_CandidateEnhancerSelection/2f_RNAseq Clustering.R
Output:IV_RNAseqClustering.pdf

####1C
Script:SingleR Annotation.R
Output:
Barplot_SingleR.pdf
UMAP_SingleR.pdf

####1D
Script:LibraryDesign.R
Outputs:
Encode annotation of candidates.pdf
Encode annotation of candidates.txt

#####1E
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_UpsetPlot.R
Output:IV_Fig1_UpsetPlot_long.pdf

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
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Nanostring vs scRNAseq.pdf
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
#####4A 
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_RelevanceToHumanBiology.R
Output:Functional annotation.pdf

#####4C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Outputs:
binding counts and fractions.pdf
Combined.pdf
footprinting signal.pdf

#####4D
Script:5_FunctionalAnnotation/EnhGenePairs/TF_CytoscapeAstronet.R
Output:AstroNet.HeatmapEG_v2.pdf

############################ Figure5
#####5A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_RelevanceToHumanBiology.R 
Output:Disease annotation.pdf

#####5C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:_no_pool.pdf

############################ Figure6
####6A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_eQTL.R
OutputS: 
SNPs per enhancer.pdf
SNPs per enhancer.txt

####6B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_eQTL.R
Output:eqtl vs crispri (pooled).pdf

####6C 
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_eQTL.R
Outputs: 
eqtl reproducibility.pdf
reproducibility.txt

####6D
Script:5_FunctionalAnnotation/Beluga/5.ISM_Beluga_CCL2.R
Output:Enh427_ISM.pdf

############################ Figure7
####7C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SanityChecks_EGrf.R
Output:NHA_EGRF_OverlapBarplot_0.5RPKM.pdf


############################ Figure S2
####2A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:Candidate annotation upset.pdf

####2B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:Candidate annotation barplot.pdf

############################ Figure S4
####4A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:nCells per candidate.pdf

####4B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:nCells per guide.pdf

####4C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:MOI.pdf

####4D 
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:positive control volcano.pdf

####4E
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:qqplot.pdf

####4F
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:dCas9KRAB.pdf

####4G
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:nEnh Per Gene.pdf

####4H
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs:
Downregulation Percentage.pdf
Downregulation Percentage.txt

####4I
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Downregulation Distribution.pdf

############################ Figure S5
####5B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:Scatterplot RNAseq vs TTseq.pdf

####5C
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:eRNA across thresholds.pdf

####5D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:FANTOM5.pdf

############################ Figure S7
####7A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TranscriptionFactors.R
Output:JASPAR class enrichments.pdf

############################ Figure S8
####8A
Script Manuscript_Final_eQTL.R
Output eqtl vs crispri (per study).pdf

############################ Figure S10
####10A
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:cells per batch.pdf

####10B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:cell-level qc.pdf

####10C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:scRNAseq qc.txt

############################ Figure S11
####11A
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Expression threshold density.pdf

####11B
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs
Expression threshold dropout rate.pdf
Output Expression threshold dropout rate.txt

####11C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs: 
Expression threshold Seurat vs CPM.pdf
Output Expression threshold Seurat vs CPM.txt

############################ Figure S12
####12A-C
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Nanostring QC.pdf

####12D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Nanostring QC on genes.pdf