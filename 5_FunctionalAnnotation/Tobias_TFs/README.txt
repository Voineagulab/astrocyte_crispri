Aggregate_Signal.R Defines a function for aggregating signal footprinting data from TOBIAS and prepares it for visualization in R

TobiasAggregatePlots.sh Generate aggregates footprinting plots comparing bound and unbound transcription factors (TFs) at enhancer regions

TobiasAnalysis.R Analyzes TOBIAS outputs

TobiasExpressedGenes.R Processesses TFs from JASPAR and generate files for use in TOBIAS analysis

TobiasFamilies.R Compares families of TFs

tobiasHeader.sh Defines path directories for files used in TOBIAS  analysis

tobiasHeaderENCODE.sh Defines path directories for files used in TOBIAS  analysis

TobiasRFunctions.R Defines R Functions used for TOBIAS analyses.

useTobias.sh Runs TOBIAS


#The usage order for running TOBIAS was the following:
TobiasExpressedGenes.R 
useTobias.sh tobiasHeader.sh
useTobias.sh tobiasHeaderENCODE.sh
TobiasAnalysis.R
TobiasFamilies.R
TobiasAggregatePlots.sh
Aggregate_Signal.R
