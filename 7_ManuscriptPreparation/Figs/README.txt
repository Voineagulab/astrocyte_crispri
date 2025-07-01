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
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
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

#####2D
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:nEnh Per Gene.pdf

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
Output:8B.Performance_Astrocytes_K562_PR.pdf

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
Output:Candidate annotation upset.pdf

####2B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:Candidate annotation barplot.pdf

############################ Extended Data Figure 3
####3A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:Sfig3A - Script LibraryDesign - Candidate length histogram.pdf

####3B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_LibraryDesign.R
Output:Sfig3B - Script LibraryDesign - Genes tested per candidate.pdf

####3C
Script:7_ManuscriptPreparation/Figs/Manuscript_final_enhancer_coverage_module_significance.R
Output:Nott_StatsCoverage_heatmap.pdf

####3D
Script:7_ManuscriptPreparation/Figs/Manuscript_final_enhancer_coverage_module_significance.R
Output:Nott_Coverage_heatmap.pdf

############################ Extended Data Figure 4
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
Outputs:
Downregulation Percentage.pdf
Downregulation Percentage.txt

####4H
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Downregulation Distribution.pdf

############################ Extended Data Figure 5
####5B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:5B - Script tt - Scatterplot RNAseq vs TTseq.pdf

####5C
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:5C - Script tt - eRNA across thresholds.pdf

####5D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:5D - Script tt - FANTOM5.pdf

####5E
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TTseq.R
Output:5E - Script tt - FANTOM5 vs TTseq.pdf

############################ Extended Data Figure 7
####7A
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_TranscriptionFactors.R
Output:JASPAR class enrichments.pdf
####7B
5_FunctionalAnnotation/EnhGenePairs/TF_CytoscapeAstronet.R
Output:AstroNet.HeatmapTF_v2.pdf
############################ Extended Data Figure 8
####8A
Script Manuscript_Final_eQTL.R
Output eqtl vs crispri (per study).pdf

############################ Extended Data Figure 9
####9A
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:8C.Astrocytes_Variable_Prediction.pdf

####9B
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:Bootstrapping_AstrocytesPlusonly.pdf

####9C
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:SFig9C_BoxplotK562_AUPRC_1000bootstraps.pdf

####9D
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:8B.Performance_Astrocytes_K562_PR.pdf

####9E
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:8C.Astrocytes_Variable_Prediction.pdf

####9F
Script:6_PredictiveModels/2.Scripts/EvaluatePredictionModels.R
Output:8D.Astrocytes_RFimportance.pdf & 8D.K562_RFimportance.pdf

############################ Supplementary Figure 1
####S1A
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:cells per batch.pdf

####S1B
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:cell-level qc.pdf

####S1C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_SequencingQC.R
Output:scRNAseq qc.txt

############################ Supplementary Figure 2
####S2A
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Expression threshold density.pdf

####S2B
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs
Expression threshold dropout rate.pdf
Output Expression threshold dropout rate.txt

####S2C
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Outputs: 
Expression threshold Seurat vs CPM.pdf
Output Expression threshold Seurat vs CPM.txt

############################ Supplementary Figure 3
####S3A-C
Script:Script 7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Nanostring QC.pdf

####S3D
Script:7_ManuscriptPreparation/Figs/Manuscript_Final_MainDE.R
Output:Nanostring QC on genes.pdf


############################ Supplementary Figure 4
####S4A
Script: 3_CRISPRi/3f_PowerSimulations.R
Output:Power Simulation - All variables V2.pdf

####S4B
Script: 3_CRISPRi/3f_PowerSimulations.R
Output:Proportion well-powered (wide).pdf

