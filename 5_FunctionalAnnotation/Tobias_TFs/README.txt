
#The usage order for running TOBIAS was the following:

TobiasExpressedGenes.R :Processesses TFs from JASPAR and generate files for use in TOBIAS analysis
useTobias.sh tobiasHeader.sh: Runs TOBIAS on in-house ATAC-seq data from NHA cells.
useTobias.sh tobiasHeaderENCODE.sh Runs TOBIAS on in-house ENCODE DNase-seq data from K562 cells (not included in the manuscript).
TobiasAnalysis.R: Analyzes TOBIAS outputs
TobiasFamilies.R: Compares families of TFs
TobiasAggregatePlots.sh: Generates aggregate plots of ATAC-seq signals across all sites using the PlotAggregate function from TOBIAS, comparing chromatin accessibility at enhancer regions between bound and unbound transcription factors (TFs)
Aggregate_Signal.R:Defines a function for aggregating signal footprinting data from TOBIAS and prepares it for visualization in R

TobiasRFunctions.R Defines R Functions used for TOBIAS analyses.
