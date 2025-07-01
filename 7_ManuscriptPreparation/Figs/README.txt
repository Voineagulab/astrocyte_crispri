Figure production scripts (in R) are organised by theme, rather than Figure / Supplementary Figure number. Plots are output as .pdf, with corresponding statistics (where relevant) saved as .txt with no delimiter.

File names automatically note from which script they are generated, as listed below. Note that not all figures are produced in these scripts, such as those produced by Cytoscape or BioRender

############################ Figure1
####1B
Script:2_NHACharacterisation_CandidateEnhancerSelection/2f_RNAseq Clustering.R
Output:IV_RNAseqClustering.pdf

####1C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SingleR_annotation.R
Output:
Fig1_Top_Script LibraryDesign_SingleR_Barplot.pdf
Fig1_Bottom_Left_UMAP_CellType.pdf
Fig1_Bottom_Middle_UMAP_DevStage.pdf
Fig1_Bottom_Right_UMAP_CellCycle.pdf

####1D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Outputs:
1D - Script LibraryDesign - Encode annotation of candidates.pdf
1D - Script LibraryDesign - Encode annotation of candidates.txt

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
Output:2A - Script MainDE - Volcano (Revised2).pdf

#####2B
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs:
2B - Script MainDE - nEnh Per Gene (Update for Revision 2).pdf
2B - Script MainDE - nGenes Per Enhancer (Update for Revision 2).pdf

#####2C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SanityChecks.R
Outputs: 
2C - Script SanityCheck - K562, Superenhancer, TAD.pdf
2C - Script SanityCheck - K562.txt
2C - Script SanityCheck - Superenhancer.txt
2C - Script SanityCheck - TAD.txt

#####2D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:2D - Script MainDE - Nanostring vs scRNAseq.pdf

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

#####3E
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output: Reads in hits versus non-hits.pdf

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
####7A
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R.R
Output:7A_BoxplotAstro_AUPRC_1000bootstraps.pdf

####7B
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R.R
Output:7B.Performance_Astrocytes_K562_PR.pdf

####7C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SanityChecks_EGrf.R
Output:NHA_EGRF_OverlapBarplot_0.5RPKM.pdf


############################ Extended Data Figure 1
####1A
Script:2_NHACharacterisation_CandidateEnhancerSelection/2e_dCas9-KRAB.R
Output:dCas9-KRAB.pdf

####1B
Script:2_NHACharacterisation_CandidateEnhancerSelection/2g_PCA_ATAC_NHA_Herring_Nott.R
Output:PCA_NHA_ATAC validation.pdf

############################ Extended Data Figure 2
####2A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:ExtFig2A - Script LibraryDesign - Candidate annotation upset.pdf

####2B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:ExtFig2B - Script LibraryDesign - Candidate annotation barplot.pdf

############################ Extended Data Figure 3
####3A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:ExtFig3A - Script LibraryDesign - Candidate length histogram.pdf

####3B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:ExtFig3B - Script LibraryDesign - Genes tested per candidate.pdf

####3C
Script:7_ManuscriptPreparation/Figs/Manuscript_final_enhancer_coverage_module_significance.R
Output:Nott_StatsCoverage_heatmap.pdf

####3D
Script:7_ManuscriptPreparation/Figs/Manuscript_final_enhancer_coverage_module_significance.R
Output:Nott_Coverage_heatmap.pdf

############################ Extended Data Figure 4
####4A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:ExtFig4A - Script qc - nCells per candidate.pdf

####4B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:ExtFig4b - Script qc - nCells per guide.pdf

####4C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:ExtFig4c - Script qc - MOI.pdf

####4D 
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:ExtFig4d - Script qc - Positive controls volcano.pdf

####4E
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:ExtFig4e - Script qc - qqplot.pdf

####4F
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:ExtFig4F - Script LibraryDesign - dCas9KRAB.pdf

####4G
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs:
ExtFig4G - Script MainDE - Downregulation Percentage.pdf
ExtFig4G - Script MainDE - Downregulation Percentage.txt

####4H
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:ExtFig4H - Script MainDE - Downregulation Distribution.pdf

############################ Extended Data Figure 5
####5B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:ExtFig5B - Script tt - Scatterplot RNAseq vs TTseq.pdf

####5C
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:ExtFig5C - Script tt - eRNA across thresholds.pdf

####5D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:ExtFig5D - Script tt - FANTOM5.pdf

####5E
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:ExtFig5E - Script tt - FANTOM5 vs TTseq.pdf

############################ Extended Data Figure 7
####7A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TranscriptionFactors.R
Output:ExtFig7A - Script TF - JASPAR class enrichments.pdf
####7B
5_FunctionalAnnotation/EnhGenePairs/TF_CytoscapeAstronet.R
Output:AstroNet.HeatmapTF_v2.pdf

############################ Extended Data Figure 8
####8A
7_ManuscriptPreparation/Figs/Script Manuscript_Final_eQTL.R
Output:ExtFig8A - Script eqtl - eqtl vs crispri (per study).pdf

############################Extended Data Figure 9
####9A
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:ExtFig9A.Astrocytes_Variable_Prediction.pdf

####9B
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:Bootstrapping_AstrocytesPlusonly.pdf

####9C
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:ExtFig9C_BoxplotK562_AUPRC_1000bootstraps.pdf

####9D
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:7B.Performance_Astrocytes_K562_PR.pdf

####9E
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:ExtFig9E.K562_Variable_Prediction.pdf

####9F
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:ExtFig9F.Astrocytes_RFimportance.pdf & 8D.K562_RFimportance.pdf

############################ Supplementary Figure 1
####S1A
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:SFig1A - Script qc - cells per batch.pdf

####S1B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:SFig1B - Script qc - cell-level qc.pdf

####S1C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:SFig1C - Script qc - scRNAseq qc.txt

############################ Supplementary Figure 2
####S2A
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:SFig2A - Script MainDE - Expression threshold density

####S2B
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs:
SFig2B - Script MainDE - Expression threshold dropout rate.pdf
SFig2B - Script MainDE - Expression threshold dropout rate.txt

####S2C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs: 
SFig2C - Script MainDE - Expression threshold Seurat vs CPM.pdf
SFig2C - Script MainDE - Expression threshold Seurat vs CPM.txt

############################ Supplementary Figure 3
####S3A-C
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:SFig3A-C - Script MainDE - Nanostring QC.pdf

####S3D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:SFig3D - Script MainDE - Nanostring QC on genes.pdf


############################ Supplementary Figure 4
####S4A
Script: 3_CRISPRi/3f_PowerSimulations.R
Output:Power Simulation - All variables V2.pdf

####S4B
Script: 3_CRISPRi/3f_PowerSimulations.R
Output:Proportion well-powered (wide).pdf

